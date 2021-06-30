import argparse
import pathlib
from tlz import cons

from common.io import Fort53Parser, BCFileWriter, read_pli
from common.geometry import kd_nearest_neighbor


def read_freq_map(path):
    freq_map = {}
    with open(path) as fin:
        for L in fin:
            k, v = L.split()
            freq_map[k] = v
    return freq_map


def main(args):
    print("Reading input files...", end='')
    src = Fort53Parser(args.fort15, args.fort53)
    src._parse_freqs()
    src._parse_grd()
    pli = read_pli(args.pli)
    print("done")

    units = [('astronomic component', '-'),
             ('waterlevelbnd amplitude', 'm'),
             ('waterlevelbnd phase', 'deg')]

    with BCFileWriter(args.output) as bc_out:
        print("Finding nearest stations...", end='')
        _, idx = kd_nearest_neighbor(src.node_coords(), pli['values'])
        idx_s = idx.argsort()
        nodes = src.nodes[idx[idx_s]]['jn'].astype(int)
        print("done")
        if args.frequencies is None:
            freqs = {k: k for k in src.freq}
        else:
            freqs = read_freq_map(args.frequencies)
        print("Reading frequencies...", len(nodes), len(freqs))
        pli_names = pli['index']
        for i, (n, F) in enumerate(src.read_freqs(nodes=nodes, freqs=list(freqs.keys()))):
            data = map(cons, freqs.values(), F)
            name = pli_names[idx_s[i]]
            bc_out.add_forcing(name, "astronomic", units, data)
            print("Station", name, end='\r')
        print()


def get_options():
    parser = argparse.ArgumentParser()

    parser.add_argument('fort15', type=pathlib.Path, help="Adcirc fort15 grid file")
    parser.add_argument('fort53', type=pathlib.Path, help="Adcirc fort.63.nc path")
    parser.add_argument('pli', type=pathlib.Path, help="Path to PLI boundary file")
    parser.add_argument("-o", "--output", type=pathlib.Path, default=pathlib.Path('.'), help="Path to bc output directory. Default is current directory")
    parser.add_argument("-f", "--frequencies", type=pathlib.Path, default=None, help="Tide frequencies to select")

    return parser.parse_args()

if __name__ == "__main__":
    args = get_options()
    main(args)
