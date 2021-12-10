# This script is intended to extract streamflow data from NWM v2.1 results
# to be used to provide streamflow forcing for a subset domain
#
# usage

import datetime
import pathlib
import numpy as np
import netCDF4 as nc
import cftime
import argparse
import operator

from common.io import BCFileWriter, read_csv, write_tim


def binary_isin(elements, test_elements, assume_sorted=False, return_indices=False):
    """Test if values of elements are present in test_elements.
    Returns a boolean array: True if element is present in test_elements, False otherwise.

    If return_indices=True, return an array of indexes that transforms the elements of
    test_elements to same order of elements (for the elements that are present in test_elements).
    ie, for every True in the returned mask, also return the index such that test_elements[indices] == elements[mask]
    """
    elements = np.asarray(elements)
    test_elements = np.asarray(test_elements)

    if assume_sorted:
        idxs = np.searchsorted(test_elements, elements)
    else:
        asorted = test_elements.argsort()
        idxs = np.searchsorted(test_elements, elements, sorter=asorted)

    valid = idxs != len(test_elements)
    test_selector = idxs[valid]
    if not assume_sorted:
        test_selector = asorted[test_selector]

    mask = np.zeros(elements.shape, dtype=bool)
    mask[valid] = test_elements[test_selector] == elements[valid]

    #indices are the array of indexes that transform the elements in
    # test_elements to match the order of elements
    if return_indices:
        return mask, test_selector[mask[valid]]
    else:
        return mask



# this function extracts the offsest values for the requested ids
def extract_offsets(current_netcdf_filename, selected_ids):
    """Return a mask on feature_id for the selected_ids

    Args:
        current_netcdf_filename (str): NetCDF file to read
        selected_ids (list_like): ids to select.

    Returns:
        (ndarray): Indexes *in order* of selected comm_ids.
    """
    with nc.Dataset(current_netcdf_filename) as ncdata:
        feature_id_index = ncdata['feature_id'][:]
        print("Size of feature id data is ", len(feature_id_index))
        fmask, fidx = binary_isin(feature_id_index, selected_ids, return_indices=True)
        fvals = feature_id_index[fmask]
        return fmask, fvals, fidx

# extract the lat long and streamflow for each comm ids
def extract_lat_lon(current_netcdf_filename, idxs):

    with nc.Dataset(current_netcdf_filename) as ncdata:
        # extract the streamflow
        lat_vals = ncdata['latitude'][idxs]
        lon_vals = ncdata['longitude'][idxs]
        return (lat_vals, lon_vals)

# extract the lat long and streamflow for each comm ids
def extract_streamflow(current_netcdf_filename, idxs):

    with nc.Dataset(current_netcdf_filename) as ncdata:
        # extract the streamflow
        stream_flow_vals = ncdata['streamflow'][idxs]
        return stream_flow_vals

def extract_qlateral(current_netcdf_filename, idxs):
    with nc.Dataset(current_netcdf_filename) as ncdata:
        # extract qlateral
        q_vals = ncdata['q_lateral'][idxs]
        return q_vals

# parse the comm ids from the first column of the com-id text file
def read_comm_ids(comm_id_file: str):
    id_list = []
    with open(comm_id_file) as comm_ids:
        for line in comm_ids:
            part = line.split()
            if len(part) == 2:
                num = int(part[0])
                id_list.append((num, part[1]))
            else:
                raise RuntimeError("Invalid comm id mapping detected")
    return id_list

def create_boundary_files(output_dir: pathlib.Path, data: dict):
    """Create DFlow boundary files in bc text file format

    Args:
        output_dir (pathlib.Path): Output directory
        data (dict):

    Raises:
        FileNotFoundError: [description]

    Returns:
        [type]: [description]
    """
    units = [("time", "seconds since 2000-01-01 00:00:00"),
             ("dischargebnd", "m^3/s")]
    date_index = cftime.date2num(data['row_index'], units[0][1], calendar='julian')
    values = data['streamflow']
    with BCFileWriter(output_dir) as bcwriter:
        for i, commid in enumerate(data['col_index']):
            v = values[:, i]
            if v.mask.any():
                continue
            ts = zip(date_index, v)
            print("Writing station:", commid, end='\r')
            bcwriter.add_forcing(commid, "timeseries", units, ts)


def create_qlat_tim_files(output_dir: pathlib.Path, data: dict):
    """Create DFlow tim files for q_lateral values

    Args:
        output_dir (pathlib.Path): Output directory
        data (dict):
    """
    date_index = cftime.date2num(data['row_index'], 'seconds since 2000-01-01 00:00:00', calendar='julian')
    values = data['qlateral']
    for i, commid in enumerate(data['col_index']):
        ts = zip(date_index, values[:, i])
        write_tim(output_dir/f"{commid}.tim", ts)


def get_inputfiles(input_path, start, end):
    """Generate list of files between start and end timestamps.

    Args:
        input_path (pathlib.Path): Directory to process
        start (datetime.datetime): Start date
        end (datetime.datetime): End date

    Raises:
        FileNotFoundError: Raised when a missing file is encountered

    Returns:
        (list): List of files matching date pattern
    """
    files = {}
    fdate = start
    one_hour = datetime.timedelta(hours=1)
    chrt = next(input_path.glob("*.CHRTOUT_DOMAIN1"))
    while fdate <= end:
        cand = chrt.with_name(fdate.strftime('%Y%m%d%H%M.CHRTOUT_DOMAIN1'))
        if cand.exists():
            files[fdate] = cand
        else:
            raise FileNotFoundError(cand)
        fdate += one_hour
    return files

# this is the main function of the script
# usage:
# python streamflow2nundging --comm_ids COMM_IDS --start_time YYYY-mm-dd_HH:MM:SS --stop_time YYYY-mm-dd_HH:MM:SS
#        [--input INPUT_DIR] [--output OUTPUT_DIR]
#
# if INPUT_DIR is not set ./input is used
# if OUTPUT_DIR is not set ./output is used

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
    parser.add_argument('--qlat', dest='qlat', action='store_true',
                    help='Extract qlateral values along with streamflow.')
    #parser.add_argument('--route-link', dest='route_link', action='store', type=pathlib.Path,
    #                help='The routlink file for the NWM run which provided input data')
    args = parser.parse_args()
    args.start_time = datetime.datetime.strptime(args.start_time, '%Y-%m-%d_%H:%M:%S')
    args.stop_time = datetime.datetime.strptime(args.stop_time, '%Y-%m-%d_%H:%M:%S')
    return args


def main(args):
    #get the path for the first input file
    files = get_inputfiles(args.input_dir, args.start_time, args.stop_time)
    print("Processing", len(files), "files in", args.input_dir)
    f0 = next(iter(files.values()))

    # get the list of comm ids that need to be read
    print ("Finding comm IDs in CHRT netCDF")
    commdata = read_csv(args.comm_id_path)
    comm_ids = np.asarray(commdata['NWMCommID'], dtype=int)
    boundaryid = np.asarray(commdata["BoundaryID"])

    # mask = nwm streamflow mask for selected commids
    # fidx = array of indices to reorder selected commids to nwm streamflow order
    # (ie comm_ids[fidx] corresponds to data[mask]
    mask, feature_ids, fidx = extract_offsets(f0, comm_ids)
    lat, lon = extract_lat_lon(f0, mask)

    #extract lat long and streamflow for the ids with stored offsets
    streamflow = np.ma.masked_array(np.zeros((len(files), len(comm_ids))), fill_value=-9999)
    if args.qlat:
        qlats = np.ma.masked_all_like(streamflow)
    data = {'lat': lat, 'lon': lon,
            'col_index': boundaryid[fidx],
            'row_index': list(files.keys())}

    print("Reading streamflow...")
    ffiles = len(files)
    for i, f in enumerate(files.values()):
        #extract streamflow for the comm ids with stored offsets
        streamflow[i] = extract_streamflow(f, mask)
        if args.qlat:
            qlats[i] = extract_qlateral(f, mask)
        print("{}/{} ({:.2%})".format(i, ffiles, i/ffiles).ljust(20), end="\r")

    # Reorder columns of streamflow to match col_index
    #streamflow = streamflow[:, selector_idx]

    # Multiply by flow direction
    if "FlowDir" in commdata:
        print("Detected flow direction...")
        fd = np.asarray(commdata["FlowDir"], dtype=int)[fidx]
        streamflow *= fd

    # Set streamflow data
    data["streamflow"] = streamflow
    create_boundary_files(args.output_dir, data)

    # Write qlateral data
    if args.qlat:
        data['qlateral'] = qlats
        create_qlat_tim_files(args.output_dir, data)

# Run main when this file is run
if __name__ == "__main__":
    args = get_options()
    main(args)
