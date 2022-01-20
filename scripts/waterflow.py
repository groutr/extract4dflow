import argparse
import json
import pathlib
import itertools
import random
import string
import datetime

import pli_decimate
import ext_slice
import waterlevel
import streamflow
import qlateral
from common.io import read_ext

FILE_KEYS = (
        "global_pli",
        "region",
        "fort63",
        "boundary_csv",
        "global_boundary",
        "streamlines"
)

DIR_KEYS = ("output_directory", "streamflow_input")


def update_boundary(src_ext, wl, sf, ld, pli):
    fname = random.choices(string.ascii_letters, k=8)
    tmp = src_ext.parent.joinpath(''.join(fname))
    with open(tmp, 'w') as fout:
        for block in read_ext(src_ext):
            if 'quantity' in block.data:
                if block.data['quantity'] == "waterlevelbnd":
                    block.data["locationfile"] = pli
                    block.data["forcingfile"] = wl
                elif block.data['quantity'] == "dischargebnd":
                    block.data["forcingfile"] = sf
            if block._type == "[lateral]":
                block.data["discharge"] = ld
            fout.write(str(block))
            fout.write("\n\n")

    # Overwrite src_ext
    tmp.rename(src_ext)



def get_options():
    parser = argparse.ArgumentParser()

    parser.add_argument('config', type=pathlib.Path, help='Configuration JSON')

    args = parser.parse_args()
    return args

def validate_directory(path):
    """Ensure that path is a directory"""
    if not path.is_dir():
        raise NotADirectoryError(path)

def validate_file(path):
    """Ensure that path to a file exists"""
    if not path.is_file():
        raise FileNotFoundError(path)

def main(args):
    with open(args.config, 'r') as fh:
        config = json.load(fh)

    #parse dates
    start_time = datetime.datetime.strptime(config["start_time"], '%Y-%m-%d_%H:%M')
    stop_time = datetime.datetime.strptime(config["stop_time"], '%Y-%m-%d_%H:%M')

    for p in FILE_KEYS:
        config[p] = pathlib.Path(config[p])
        validate_file(config[p])

    for p in DIR_KEYS:
        config[p] = pathlib.Path(config[p])
        validate_directory(config[p])

    region_stem = config["region"].stem
    out_dir = config["output_directory"]

    print("PLI Decimate...")
    pli_decimate_args = argparse.Namespace()
    pli_decimate_args.n = 1
    pli_decimate_args.pli = config["global_pli"]
    pli_decimate_args.polygon = config["region"]
    pli_out = out_dir.joinpath(f"{pli_decimate_args.pli.stem}_slice_{region_stem}.pli")
    pli_decimate_args.output = pli_out
    print(pli_decimate_args)
    pli_decimate.main(pli_decimate_args)
    assert pli_out.exists()
    print(pli_out)

    print("EXT Slice...")
    ext_slice_args = argparse.Namespace()
    ext_slice_args.boundary_csv = config["boundary_csv"]
    ext_slice_args.polygon = config["region"]
    ext_slice_args.ext = config["global_boundary"]
    ext_slice_args.streamlines = config["streamlines"]
    ext_slice_args.output_dir = out_dir
    ext_slice.main(ext_slice_args)
    [ext_out] = list(out_dir.glob(f"*_slice_{region_stem}.ext"))
    [boundary_csv] = list(out_dir.glob(f"{config['boundary_csv'].stem}_slice_{region_stem}.csv"))
    print(ext_out)

    print("Waterlevel extraction...")
    waterlevel_args = argparse.Namespace()
    waterlevel_args.fort63 = config["fort63"]
    waterlevel_args.boundary_csv = boundary_csv
    wl_out = out_dir.joinpath(f"waterlevel_slice_{region_stem}.bc")
    waterlevel_args.output = wl_out
    waterlevel.main(waterlevel_args)
    assert wl_out.exists()
    print(wl_out)

    print("Streamflow extraction...")
    streamflow_args = argparse.Namespace()
    streamflow_args.comm_id_path = boundary_csv
    streamflow_args.start_time = start_time
    streamflow_args.stop_time = stop_time
    streamflow_args.input_dir = config["streamflow_input"]
    sf_out = out_dir.joinpath(f"streamflow_slice_{region_stem}.bc")
    streamflow_args.output_dir = sf_out
    print(streamflow_args)
    streamflow.main(streamflow_args)
    assert sf_out.exists()
    print(sf_out)

    print("Lateral flow extraction...")
    qlateral_args = argparse.Namespace()
    qlateral_args.comm_id_path = boundary_csv
    qlateral_args.start_time = start_time
    qlateral_args.stop_time = stop_time
    qlateral_args.input_dir = config["streamflow_input"]
    ql_out = out_dir.joinpath(f"qlateral_slice_{region_stem}.bc")
    qlateral_args.output_dir = ql_out
    print(qlateral_args)
    qlateral.main(qlateral_args)
    assert ql_out.exists()
    print(ql_out)

    print("Updating boundary", ext_out)
    update_boundary(ext_out, wl_out.name, sf_out.name, ql_out.name, pli_out.name)


if __name__ == "__main__":
    args = get_options()
    main(args)


