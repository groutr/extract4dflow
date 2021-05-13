import argparse
import pathlib
import time

import netCDF4
#import h5netcdf.legacyapi as netCDF4
import numpy as np
from scipy.spatial import cKDTree


class BCFileWriter:
    functions = {'timeseries'}

    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        self._filehandle = open(self.filename, 'w')
        return self

    def __exit__(self, type, value, traceback):
        self._filehandle.close()

    def add_forcing(self, name, function, units, data):
        """Add forcing
        Args:
            name (str): Name
            function (str): One of BCFileWriter.functions
            units (list[tuples]): A list of tuples mapping column name to column units.
                The ordering should match the ordering of the data columns.
            data (Iterable of lists): Number of columns in data and len(units) must match.
                Data will be iterated thru row by row
        Returns:
            None
        """
        if function not in self.functions:
            raise ValueError("Invalid function")

        fh = self._filehandle
        fh.write("[forcing]\n")
        fh.write(f"Name = {name}\n")
        fh.write(f"Function = {function}\n")
        fh.write(f"Time-interpolation = linear\n")
        for i, (col, unit) in enumerate(units):
            fh.write(f"Quantity = {col}\n")
            fh.write(f"Unit = {unit}\n")

        if isinstance(data, np.ndarray):
            np.savetxt(fh, data, fmt='%f', delimiter=' ')
        else:
            for row in data:
                fh.write(" ".join(map(str, row)))
                fh.write("\n")

        fh.write("\n")



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


def read_pli(path):
    index = []
    with open(path) as fp:
        name = next(fp).strip()  # Name row
        shape = tuple(map(int, next(fp).split()))  # shape
        values = np.empty(shape, dtype='float64')
        for i, L in enumerate(fp):
            L = L.split()
            if L:
                values[i] = L[:2]
                index.append(L[2])
    return {'name': name, 'index': index, 'values': values}

def read_waterlevel(fort63, pli, bc_output):
    print("Reading PLI")
    pli_data = read_pli(pli)

    with netCDF4.Dataset(fort63, mode='r') as ds:
        zeta = ds.variables['zeta']
        print("Masking invalid stations")
        mask = invalid_mask(zeta)
        valid_stations = mask.nonzero()[0]
        print("Masking coordinates")
        #adlons = ds.variables['x'][:][mask]
        #adlats = ds.variables['y'][:][mask]
        adlons = np.ma.getdata(ds.variables['x'][mask])
        adlats = np.ma.getdata(ds.variables['y'][mask])
        print("Querying nearest points")
        tree = cKDTree(np.column_stack([adlons, adlats]))
        CN = tree.query(pli_data['values'])
        stations = valid_stations[CN[1]]

        time_col = ds.variables['time']
        ref_time = time_col.units.rstrip('UTC').rstrip()
        units = [('time', ref_time),
                 ('waterlevelbnd', 'm')]
        out_buf = np.empty((len(time_col), 2), dtype='float64')
        out_buf[:, 0] = time_col[:]
        if bc_output.is_dir():
            bc_output = bc_output/f"{pli_data['name']}.bc"

        with BCFileWriter(bc_output) as bc_out:
            print("Writing BC output", bc_out.filename)
            for name, station in zip(pli_data['index'], stations):
                print(f"Station {name} ({station})".ljust(50), end="\r")
                #out_buf[:, 1] = zeta[:, station]
                out_buf[:, 1] = np.ma.getdata(zeta[:, station])
                bc_out.add_forcing(name, 'timeseries', units, out_buf)
        print()

def get_options():
    parser = argparse.ArgumentParser()

    parser.add_argument('fort63', type=pathlib.Path, help="Adcirc fort.63.nc path")
    parser.add_argument('pli', type=pathlib.Path, help="Path to PLI boundary file")
    parser.add_argument("-o", "--output", type=pathlib.Path, default=pathlib.Path('.'), help="Path to bc output directory. Default is current directory")

    return parser.parse_args()

if __name__ == "__main__":
    args = get_options()
    read_waterlevel(args.fort63, args.pli, args.output)


