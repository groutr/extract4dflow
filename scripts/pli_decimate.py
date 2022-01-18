import argparse
import pathlib
import numpy as np
from itertools import islice, compress
from common.io import read_pli, write_pli, read_polygon
from common.geometry import clip_point_to_roi


def rename_points(pts):
    """Rename pts so they are in ascending order

    Assumes that index names are in the form of *_<number>.
    This will completely replace number left padded to original width.

    Modifies index in-place.

    Args:
        pts (dict): Contents of PLI boundary file
    """
    index = pts['index']
    tmp = index[0].rsplit('_', 1)
    assert len(tmp) == 2
    pad_width = len(tmp[1])
    for i, name in enumerate(index):
        i_str = str(i+1).zfill(pad_width)
        index[i] = f"{name.rsplit('_', 1)[0]}_{i_str}"
    return pts


def main(args):
    pli = read_pli(args.pli)

    if args.polygon:
        P = read_polygon(args.polygon)
        idxs = clip_point_to_roi(P, pli['values'])
        pli['values'] = pli['values'][idxs]
        pli['index'] = list(compress(pli['index'], idxs))

    new_index, new_values = zip(*islice(zip(pli['index'], pli['values']), 0, None, args.n))
    pli['index'] = list(new_index)
    pli['values'] = np.array(new_values)
    write_pli(args.output, pli)


def get_options():
    parser = argparse.ArgumentParser()

    parser.add_argument('pli', type=pathlib.Path, help="Path to PLI boundary file")
    parser.add_argument('-n', type=int, default=1, help="Decimation factor")
    parser.add_argument('-o', '--output', type=pathlib.Path, help="Output file path")
    parser.add_argument("-c", "--clip", dest="polygon", default=None, type=pathlib.Path, help="Use a polygon to select sites for output")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_options()
    main(args.pli, args.output, args.n, polygon=args.clip)
