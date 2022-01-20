import argparse
import datetime
import pathlib

import numpy as np
import netCDF4 as nc
import cftime

from common.data import extract_lat_lon, extract_offsets
from common.io import write_tim, get_inputfiles, read_csv, write_pli, BCFileWriter


def extract_qlateral(current_netcdf_filename, idxs):
    with nc.Dataset(current_netcdf_filename) as ncdata:
        # extract qlateral
        q_vals = ncdata['q_lateral'][idxs]
        return q_vals

def create_qlat_bc_file(output_dir: pathlib.Path, data: dict):
    """Create boundary condition file for lateral discharge

    Args:
        output_dir (pathlib.Path): output directory for BoundaryConditions.bc
        data (dict): [description]
    """
    units = [("time", "minutes since 2000-01-01 00:00:00"),
            ("lateral_discharge", "m^3/s")]
    date_index = cftime.date2num(data['row_index'], units[0][1], calendar='julian')
    values = data['qlateral']
    with BCFileWriter(output_dir/"BoundaryConditions.bc") as bcwriter:
        for i, commid in enumerate(data['col_index']):
            v = values[:, i]
            if v.mask.any():
                continue
            ts = zip(date_index, v)
            bcwriter.add_forcing(f"{commid}_0001", "timeseries", units, ts)

def create_qlat_pli_files(output_dir: pathlib.Path, data: dict):
    for i, commid in enumerate(data['col_index']):
        lat = data['lat'][i]
        lon = data['lon'][i]
        _data = {'name': commid, 'values': np.array([[lon, lat]]), 'index': [commid]}
        write_pli(output_dir/f"{commid}.pli", _data)

def create_qlat_tim_files(output_dir: pathlib.Path, data: dict):
    """Create DFlow tim files for q_lateral values

    Args:
        output_dir (pathlib.Path): Output directory
        data (dict):
    """
    date_index = cftime.date2num(data['row_index'], 'seconds since 2000-01-01 00:00:00', calendar='julian')
    values = data['qlateral']
    ts = np.empty((len(date_index), 2), dtype=float)
    ts[:, 0] = date_index
    for i, commid in enumerate(data['col_index']):
        ts[:, 1] = values[:, i]
        write_tim(output_dir/f"{commid}.tim", ts)


def main(args):
    #get the path for the first input file
    files = get_inputfiles(args.input_dir, args.start_time, args.stop_time)
    print("Processing", len(files), "files in", args.input_dir)
    f0 = next(iter(files.values()))

    # get the list of comm ids that need to be read
    print ("Finding comm IDs in CHRT netCDF")
    commdata = read_csv(args.comm_id_path)
    comm_ids = np.asarray(commdata['NWMCommID'])
    boundaryid = np.asarray(commdata["BoundaryID"])

    # mask = nwm streamflow mask for selected commids
    # fidx = array of indices to reorder selected commids to nwm streamflow order
    # (ie comm_ids[fidx] corresponds to data[mask]
    mask, feature_ids, fidx = extract_offsets(f0, comm_ids)
    lat, lon = extract_lat_lon(f0, mask)

    #extract lat long and streamflow for the ids with stored offsets
    qlats = np.ma.masked_array(np.zeros((len(files), len(comm_ids))), fill_value=-9999)
    data = {'lat': lat, 'lon': lon,
            'col_index': comm_ids[fidx],
            'row_index': list(files.keys())}

    print("Reading qlateral...")
    ffiles = len(files)
    for i, f in enumerate(files.values()):
        #extract streamflow for the comm ids with stored offsets
        qlats[i] = extract_qlateral(f, mask)
        print("{}/{} ({:.2%})".format(i, ffiles, i/ffiles).ljust(20), end="\r")

    # Write qlateral data
    data['qlateral'] = qlats
    create_qlat_tim_files(args.output_dir, data)
    create_qlat_pli_files(args.output_dir, data)
    create_qlat_bc_file(args.output_dir, data)


def get_options():
    parser = argparse.ArgumentParser(description='Create nundging files from NWM streamflow output.')
    parser.add_argument('--comm_ids', dest='comm_id_path', required=True,
                    help='The file list the comm ids to be extracted. Two columns map the comm id to an identifier string (used for boundary condition output).')
    parser.add_argument('--start', dest='start_time', required=True,
                    help='The date and time to begin making timeslice files for YYYY-mm-dd_HH:MM:SS')
    parser.add_argument('--stop', dest='stop_time', required=True,
                    help='The date and time to stop making timeslice files for YYYY-mm-dd_HH:MM:SS')
    parser.add_argument('--input', dest='input_dir', required=True, type=pathlib.Path,
                    help='The directory that stores NWM streamflow input files')
    parser.add_argument('--output', dest='output_dir', required=True, type=pathlib.Path,
                    help='The directory to write timeslice files to')

    args = parser.parse_args()
    args.start_time = datetime.datetime.strptime(args.start_time, '%Y-%m-%d_%H:%M:%S')
    args.stop_time = datetime.datetime.strptime(args.stop_time, '%Y-%m-%d_%H:%M:%S')
    return args


# Run main when this file is run
if __name__ == "__main__":
    args = get_options()
    main(args)
