"""
Compute tangent and normal velocity

Extracts U/V wind velocities and projects to tangent and normal vectors

Usage:
python velocity.py fort64.nc boundary.pli


Returns:
    _type_: _description_
"""

import argparse
import pathlib
import numpy as np
import xarray as xr

from common.io import read_pli, BCFileWriter
from common.geometry import kd_nearest_neighbor

def tangent_normal_vector(pts):
    diff = np.diff(pts, axis=0)
    dists = (diff ** 2).sum(axis=1) ** 0.5
    utan = diff/dists[:, np.newaxis]
    unorm = np.column_stack((utan[:,1], -utan[:,0]))
    return utan, unorm

def project_unit(unit, U, V):
    return unit[:,0] * U + unit[:,1] * V

def get_options():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('fort64', type=pathlib.Path, help="Adcirc fort.64.nc path")
    parser.add_argument('pli', type=pathlib.Path, help="Path to boundary csv file")
    parser.add_argument("-o", "--output", type=pathlib.Path, default=pathlib.Path('.'), help="Path to bc output directory. Default is current directory")
    return parser.parse_args()

def main(args):
    pli_values = read_pli(args.pli)
    pts = pli_values['values']

    utan, unorm = tangent_normal_vector(pts)
    utan = np.append(utan, np.atleast_2d(utan[-1]), axis=0)
    unorm = np.append(unorm, np.atleast_2d(unorm[-1]), axis=0)

    with xr.open_dataset(args.fort64) as DS:
        vertices = np.column_stack((DS.x.values, DS.y.values))
        _, nodes = kd_nearest_neighbor(vertices, pts)
        U = DS['u-vel']
        V = DS['v-vel']

        ref_time = DS.variables['time'].base_date.rstrip('UTC').rstrip()
        tunits = [('time', ref_time),
                ('tangentialvelocitybnd', 'm/s')]
        nunits = [('time', ref_time),
                ('normalvelocitybnd', 'm/s')]

        out_buf = np.empty((len(DS.variables['time']), 2), dtype='float64')
        out_buf[:, 0] = DS.variables['time'].values
        with BCFileWriter(args.output/"TangentVelocity.bc") as tan_out:
            with BCFileWriter(args.output/"NormalVelocity.bc") as nor_out:
                for name, n in zip(pli_values['index'], nodes):
                    print(f"Station {name}".ljust(50), end="\r")
                    uv = U[:, n].values
                    vv = V[:, n].values
                    out_buf[:, 1] = project_unit(utan, uv, vv)
                    tan_out.add_forcing(name, 'timeseries', tunits, out_buf)

                    out_buf[:, 1] = project_unit(unorm, uv, vv)
                    nor_out.add_forcing(name, 'timeseries', nunits, out_buf)
                print()


if __name__ == "__main__":
    args = get_options()
    main(args)
