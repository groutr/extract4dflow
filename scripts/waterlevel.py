import argparse
import pathlib
import time

import netCDF4
#import h5netcdf.legacyapi as netCDF4
import numpy as np

from common.io import BCFileWriter, read_csv
from common.geometry import kd_nearest_neighbor


def invalid_mask(var, chunksize=2**18):
    """Generate a mask of a masked array for the columns that only have completely valid observations.

    An optional chunksize can be set to avoid reading the entire variable mask array into memory.

    Note that the masked array mask is inverted, ie a True value indicates a missing value.
    The mask returned from this function, a True indicates a column with no missing values.

    Args:
        var (netCDF4.Variable): Variable to consider
        chunksize (int, optional): Number of columns to consider at a time. Defaults to 2**18.

    Returns:
        np.ndarray: Mask of valid columns
    """
    rv = np.empty(var.shape[1], dtype=bool)
    buf = np.empty((var.shape[0], chunksize), dtype='bool')
    for i in range(0, var.shape[1], chunksize):
        print(i)
        rv_view = rv[i:i+chunksize]
        buf_view = buf[:, :rv_view.shape[0]]
        #np.equal(var[:, i:i+chunksize], var._FillValue, out=buf_view)
        cols = np.ma.getmaskarray(var[:, i:i+chunksize])
        np.any(buf_view, axis=0, out=rv_view)
        np.logical_not(rv_view, out=rv_view)
    return rv


def main(args):
    print("Reading Boundary CSV")
    csv_data = read_csv(args.boundary_csv)

    with netCDF4.Dataset(args.fort63, mode='r') as ds:
        zeta = ds.variables['zeta']
        print("Masking invalid stations")
        mask = invalid_mask(zeta)
        valid_stations = mask.nonzero()[0]
        print("Masking coordinates")
        #adlons = ds.variables['x'][:][mask]
        #adlats = ds.variables['y'][:][mask]
        adlons = np.ma.getdata(ds.variables['x'][mask])
        adlats = np.ma.getdata(ds.variables['y'][mask])

        adpts = np.column_stack([adlons, adlats])
        csvpts = np.column_stack([csv_data['long'], csv_data['lat']])
        print("Querying nearest points")
        _, CN = kd_nearest_neighbor(adpts, csvpts)
        stations = valid_stations[CN]

        time_col = ds.variables['time']
        ref_time = time_col.units.rstrip('UTC').rstrip()
        units = [('time', ref_time),
                 ('waterlevelbnd', 'm')]
        out_buf = np.empty((len(time_col), 2), dtype='float64')
        out_buf[:, 0] = time_col[:]
        if args.output.is_dir():
            args.output = args.output/f"waterlevel.bc"

        with BCFileWriter(args.output) as bc_out:
            print("Writing BC output", bc_out.filename)
            for name, station in zip(csv_data['NWMCommID'], stations):
                print(f"Station {name} ({station})".ljust(50), end="\r")
                out_buf[:, 1] = np.ma.getdata(zeta[:, station])
                bc_out.add_forcing(name, 'timeseries', units, out_buf)
        print()

def get_options():
    parser = argparse.ArgumentParser()

    parser.add_argument('fort63', type=pathlib.Path, help="Adcirc fort.63.nc path")
    parser.add_argument('boundary_csv', type=pathlib.Path, help="Path to boundary csv file")
    parser.add_argument("-o", "--output", type=pathlib.Path, default=pathlib.Path('.'), help="Path to bc output directory. Default is current directory")

    return parser.parse_args()

if __name__ == "__main__":
    args = get_options()
    main(args)
