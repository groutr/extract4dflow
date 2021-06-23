import argparse
import pathlib
import re
from tlz import cons

from common.io import Fort53Parser, BCFileWriter, read_pli
from common.geometry import kd_nearest_neighbor

dflow_re = re.compile(r'\W')

def main(args):
    src = Fort53Parser(args.fort15, args.fort53)
    src._parse_freqs()
    src._parse_grd()
    pli = read_pli(args.pli)

    units = [('astronomic component', None),
             ('waterlevelbnd amplitude', 'm'),
             ('waterlevelbnd phase', 'deg')]
    
    with BCFileWriter(args.output) as bc_out:
        _, idx = kd_nearest_neighbor(src.node_coords(), pli['values'])
        idx_s = idx.argsort()
        nodes = src.nodes[idx[idx_s]]['jn']
        if args.frequencies is None:
            freqs = src.freq
        else:
            freqs = args.frequencies.read_text().splitlines()
        outfreqs = tuple(dflow_re.sub("", x) for x in freqs)
        for i, (n, F) in enumerate(src.read_freqs(nodes=nodes, freqs=freqs)):
            data = map(cons, outfreqs, F)
            bc_out.add_forcing(pli['index'][idx_s[i]], "astronomic", units, data)


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