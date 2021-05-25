from itertools import islice
import math
import numpy as np

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


def read_pli(path, step=1):
    index = []
    with open(path) as fp:
        name = next(fp).strip()  # Name row
        shape = list(map(int, next(fp).split()))  # shape
        assert shape[1] == 2
        if step > 1:
            # adjust shape by dividing by step
            shape[0] = math.ceil(shape[0] / step)
        values = np.empty(shape, dtype='float64')
        for i, L in enumerate(islice(fp, 0, None, step)):
            L = L.split()
            if L:
                values[i] = L[:2]
                index.append(L[2])
    return {'name': name, 'index': index, 'values': values}


def read_polygon(path):
    with open(path) as fp:
        # Skip all lines that begin with '*' character at beginning of file
        L = next(fp)
        while L.lstrip().startswith('*'):
            L = next(fp)

        shape = list(map(int, next(fp).split()))
        assert shape[1] == 2
        pts = np.empty(shape, dtype='float64')
        for i, L in enumerate(fp):
            L = L.split()
            if L:
                pts[i] = L
    return pts


def write_pli(path, pli):
    """Write a PLI file located at <path>

    Args:
        path (pathlib.Path): File path
        pli (dict): PLI contents
    """
    with open(path, 'w') as out:
        out.write(pli['name'])
        out.write("\n")
        out.write(" ".join(map(str, pli['values'].shape)))
        out.write("\n")
        for (lon, lat), name in zip(pli['values'], pli['index']):
            out.write(f"{lon} {lat}  {name}\n")

def read_csv(path):
    return np.recfromcsv(path, encoding=None)

def read_ext(path):
    """Read ext file and return list of forcing sections"""
    sections = []
    with open(path, 'r') as fin:
        s = []
        for L in map(str.strip, fin):
            if not L:
                continue

            if L == "[boundary]" and s:
                yield Block(s)
                s = [L]
            else:
                s.append(L)
    return sections

class Block:
    def __init__(self, block_lines):
        self.data = {}

        if block_lines[0][0] == '[' and block_lines[0][-1] == ']':
            self._type = block_lines[0]
        else:
            raise ValueError("Invalid block")

        for line in block_lines[1:]:
            key, value = line.split('=')
            self.data[key.strip()] = value.strip()

    def __str__(self):
        lines = [self._type]
        for key, value in self.data.items():
            lines.append(f"{key}={value}")
        return "\n".join(lines)

    @property
    def type(self):
        return self._type[1:-1]


