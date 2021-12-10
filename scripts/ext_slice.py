import argparse
import pathlib
import numpy as np
import csv
import itertools

from common.geometry import clip_point_to_roi
from common.io import read_polygon, read_csv, read_ext


def removesuffix(s, suffix):
    if suffix and s.endswith(suffix):
        return s[:-len(suffix)]
    else:
        return s

def create_ext_subdomain(src_ext, dst_ext, mask):
    with open(dst_ext, mode='w') as fout:
        for block in read_ext(src_ext):
            # check if location is masked
            if removesuffix(block.data['locationfile'], '.pli') not in mask:
                fout.write(str(block))
                fout.write("\n\n")

def create_csv_subdomain(src_csv, dst_csv, mask):
    """Read a source csv and write masked output to dst_csv.

    Reading the source from disk is necessary because NumPy will lowercase
    all column names
    """
    with open(dst_csv, 'w', newline='') as fout:
        csvw = csv.writer(fout)
        with open(src_csv, 'r') as fin:
            csvr = csv.reader(fin)

            csvw.writerow(next(csvr))

            for row in itertools.compress(csvr, mask):
                csvw.writerow(row)


# user defined options
def get_options():
    parser = argparse.ArgumentParser(description='Create ext DFlow subdomain file based on user specified Polygon.')
    parser.add_argument('boundary_csv', type=pathlib.Path,
                    help='CSV file mapping boundary ids to geospatial position.')
    parser.add_argument('polygon', type=pathlib.Path,
                    help='The path of the polygon file defining the region of interest.')
    parser.add_argument('--ext', dest='ext', type=pathlib.Path,
                    help='The path of the original DFlow continetal mesh boundary ID ext file to extract boundary segments within polygon')
    parser.add_argument('--streamlines', type=pathlib.Path,
                    help='CSV file georeferencing boundary ids for lateral dicharges. This will be clipped to region of interest.')
    parser.add_argument('-o', '--output', dest='output_dir', default=pathlib.Path('.'), type=pathlib.Path,
                    help='The directory to write DFlow subdomain ext file to')
    args = parser.parse_args()

    # Validate that output_dir is an existing directory
    if not args.output_dir.exists():
        raise FileNotFoundError(args.output_dir)
    elif not args.output_dir.is_dir():
        raise NotADirectoryError(args.output_dir)

    return args

def main(args):

    # read Boundary ID dataframe containing geospatial information
    bnd_id_data = read_csv(args.boundary_csv, cols=["long", "lat", "boundaryid"])
    # extract polygon coordinate info from user defined file
    polygon_coords = read_polygon(args.polygon)

    # Extract boundary point coordinates
    bnd_pts = np.column_stack([bnd_id_data["long"], bnd_id_data["lat"]]).astype(float)
    poly_mask = np.array(clip_point_to_roi(polygon_coords, bnd_pts), dtype='bool')

    if args.ext:
        ext_name = f"{args.ext.stem}_slice_{args.polygon.stem}{args.ext.suffix}"
        exclude_blocks = set(itertools.compress(bnd_id_data["boundaryid"], ~poly_mask))
        create_ext_subdomain(args.ext, args.output_dir.joinpath(ext_name), exclude_blocks)

    if args.streamlines:
        stl_name = f"{args.streamlines.stem}_slice_{args.polygon.stem}{args.streamlines.suffix}"
        _data = read_csv(args.streamlines, usecols=['lat', 'lon'])
        _pts = np.column_stack([_data['lon'], _data['lat']]).astype(float)
        stl_mask = np.array(clip_point_to_roi(polygon_coords, _pts), dtype=bool)
        create_csv_subdomain(args.streamlines, args.output_dir.joinpath(stl_name), stl_mask)

    csv_name = f"{args.boundary_csv.stem}_slice_{args.polygon.stem}{args.boundary_csv.suffix}"
    create_csv_subdomain(args.boundary_csv, args.output_dir.joinpath(csv_name), poly_mask)


##### User command line example  ########
#./python3.7 ../../DFlow_polygon_slice_ext.py --bnd_ids_csv="./InflowBoundaryConditions.csv" --polygon_file="./Irma_Enclosure.pol" --flowfm_bnd_ext="./FlowFM_bnd.ext" --output="./"


# Run main when this file is run
if __name__ == "__main__":
    args = get_options()
    main(args)

