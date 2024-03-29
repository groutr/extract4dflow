from itertools import islice
import math
import numpy as np
import string
import csv
import operator
import datetime
from textwrap import dedent
from io import StringIO
from tlz import take, drop
from tlz.curried import get

class Fort53Parser:
    def __init__(self, fort15, fort53):
        self.fort15 = fort15
        self.fort53 = fort53

        self.freq = None
        self.nodes = None

    def _parse_freqs(self):
        if self.freq is not None:
            return

        freq = []
        with open(self.fort53) as fh:
            nfreq = int(next(fh))
            for i, line in enumerate(fh):
                if i == nfreq:
                    break
                f = line.rsplit(None, 1)
                freq.append(f[1])
        self.freq = tuple(freq)

    def _parse_grd(self):
        if self.nodes is not None:
            return

        with open(self.fort15) as fh:
            next(fh)
            NE, NP = map(int, next(fh).split())
            self.nodes = np.recfromtxt(fh, dtype='float64', max_rows=NP, names=('jn', 'lon', 'lat', 'depth'))

    def node_coords(self):
        self._parse_grd()
        return np.asarray([self.nodes['lon'], self.nodes['lat']]).T

    def read_freqs(self, nodes=None, freqs=None):
        """Read tidal frequencies for each node

        Nodes is assumed to be sorted.
        """
        # Find the indexes for each freq
        self._parse_grd()
        self._parse_freqs()

        if freqs is None:
            freqs = self.freq

        if nodes is None:
            nodes = self.nodes['jn'].astype(int)

        idxget = get(list(map(self.freq.index, freqs)))

        with open(self.fort53) as fh:
            nfreq = int(next(fh))
            # skip the next nfreq lines
            fh = drop(nfreq, fh)
            NP = int(next(fh))
            assert NP == len(self.nodes)

            assert len(self.freq) == nfreq
            node = int(next(fh))
            for i, n in enumerate(nodes):
                if node - 1 == n:
                    yield (node, block)
                    continue

                if node > n:
                    raise RuntimeError

                if node < n:
                    # How many lines to skip to node
                    fh = drop((n-node) * (nfreq + 1) - 1, fh)
                    node = int(next(fh))
                if node == n:
                    if i == 0 or nodes[i-1] != node:
                        data = tuple(take(nfreq, fh))
                        block = np.loadtxt(idxget(data), dtype=float)
                    yield (node, block)
                    node = int(next(fh))


class BCFileWriter:
    functions = {'timeseries', 'astronomic'}

    forcing_template = (
        "[forcing]\n"
        "Name = $name\n"
        "Function = $function\n"
        "Time-interpolation = linear"
    )

    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        self._filehandle = open(self.filename, 'w')
        return self

    def __exit__(self, type, value, traceback):
        self._filehandle.close()

    def add_forcing(self, name, function, units, data, vectors=None):
        """Add forcing
        Args:
            name (str): Name
            function (str): One of BCFileWriter.functions
            units (list[tuples]): A list of tuples mapping column name to column units.
                The ordering should match the ordering of the data columns.
            data (Iterable of lists): Number of columns in data and len(units) must match.
                Data will be iterated thru row by row
            vectors (list[str]):
        Returns:
            None
        """
        T = [string.Template(self.forcing_template).substitute(name=name, function=function)]
        if vectors is not None:
            for v in vectors:
                T.append(f"Vector = {v}")

        for i, (col, unit) in enumerate(units):
            T.append(f"Quantity = {col}")
            T.append(f"Unit = {unit}")
        
        arr = StringIO()
        if isinstance(data, np.ndarray):
            np.savetxt(arr, data, fmt='%f', delimiter=' ')
        else:
            for row in data:
                arr.write(" ".join(map(str, row)))
                arr.write("\n")
        
        T.append(arr.getvalue())
        self._write_str(T)

    def _write_str(self, lines):
        for L in lines:
            self._filehandle.write(L)
            self._filehandle.write("\n")
        self._filehandle.write("\n")


def write_tim(path, data):
    with open(path, 'w') as fh:
        np.savetxt(fh, data, fmt='%f', delimiter=' ')

def read_pli(path, step=1):
    with open(path) as fp:
        name = next(fp).strip()  # Name row
        shape = list(map(int, next(fp).split()))  # shape
        assert shape[1] == 2
        index, values = _read_xyn(fp, step=step)
        values = np.asarray(values, dtype='float64')
    return {'name': name, 'index': index, 'values': values}


def read_xyn(path, step=1):
    with open(path) as fp:
        index, values = _read_xyn(fp, step=step)
        values = np.asarray(values, dtype='float64')
    return {'name': None, 'index': index, 'values': values}


def _read_xyn(fp, step=1):
    index = []
    values = []
    for L in islice(fp, 0, None, step):
        L = L.strip().split(maxsplit=2)
        if L:
            values.append(L[:2])
            index.append(L[2])
    return index, values


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
            PLI dictionary is keyed by 'name', 'values', and 'index'
    """
    with open(path, 'w') as out:
        out.write(pli['name'])
        out.write("\n")
        out.write(" ".join(map(str, pli['values'].shape)))
        out.write("\n")
        for (lon, lat), name in zip(pli['values'], pli['index']):
            out.write(f"{lon} {lat}  {name}\n")


def write_xyn(path, xyn):
    with open(path, 'w') as out:
        for (lon, lat), name in zip(xyn['values'], xyn['index']):
            out.write(f"{lon} {lat}  {name}\n")


def read_csv(path, cols=None):
    # We could use pandas here, but we need more of a
    # justification than just reading a csv.
    with open(path) as fin:
        csvr = csv.reader(fin)
        columns = next(csvr)

        if not cols:
            cols = columns
            col_slicer = operator.itemgetter(slice(0, None))
            row_append = [].append
        else:
            col_slicer = operator.itemgetter(*(columns.index(i) for i in cols))
            row_append = [].append

        for row in csvr:
            row_append(col_slicer(row))

        return dict(zip(cols, zip(*row_append.__self__)))


def read_ext(path):
    """Read ext file and return list of forcing sections"""
    sections = []
    with open(path, 'r') as fin:
        s = []
        for L in map(str.strip, fin):
            if not L:
                continue
            
            if L.startswith('[') and L.endswith(']') and s:
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


def write_ext_v2(path, data, mode='w'):
    """Write ext format v2

    Arguments:
        path: Path of file to write
        data: List of dictionaries having the following keys:
            id: id of node
            name: name of node
            x: Longitude coordinate
            y: Latitude coordinate
    
    """
    header = "[General]\nfileVersion = 2.00\nfileType = extForce\n"
                    
    template = string.Template(dedent("""
            [lateral]
            id = $id
            name = $name
            type = discharge
            LocationType = all
            numCoordinates = 1
            xCoordinates = $x
            yCoordinates = $y
            discharge = BoundaryConditions.bc
            """))

    with open(path, mode=mode) as ext_out:
        ext_out.write(header)
        for block in data:
            ext_out.write(template.substitute(block))
            ext_out.write("\n\n")

def write_ext(path, data, mode='w'):
    """Write legacy ext format

    Arguments:
        path: Path of file to write
        data: list of dictionaries having the following keys:
            filename: filename
    """

    template = string.Template(dedent("""
    QUANTITY=discharge_salinity_temperature_sorsin
    FILENAME=$filename
    FILETYPE=9
    METHOD=1
    OPERAND=O
    AREA=1"""))
    with open(path, mode=mode) as ext_out:
        ext_out.write("\n")
        for d in data:
            ext_out.write(template.substitute(d))
            ext_out.write("\n\n")


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
        files[fdate] = chrt.with_name(fdate.strftime('%Y%m%d%H%M.CHRTOUT_DOMAIN1'))
        fdate += one_hour
    return files


    


