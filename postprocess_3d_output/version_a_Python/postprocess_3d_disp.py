#!/usr/bin/env python3
"""

Postprocess 3D displacement output to NetCDF.

Input: a parameter file with 8 lines in fixed order:
1) PREFIX
2) Step_id
3) N_CPU_surface
4) N_CPU_Z
5) n_surface_nodes_per_cpu
6) n_z_nodes_per_cpu
7) output_dir
8) output_prefix

Output: NetCDF file with variables
- Vr_cumu (r, loc): cumulative radial velocity and etc
- dr_incr (r, loc): incremental displacement for this step and etc

also output an interpolated grid NetCDF file with variables
- Vr_cumu (r, lat, lon): cumulative radial velocity and etc
- dr_incr (r, lat, lon): incremental displacement for this step and etc
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Tuple

import numpy as np
import xarray as xr
import pygmt


def read_params(path: str) -> Tuple[str, int, int, int, int, int, str, str]:
    with open(path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]

    if len(lines) < 8:
        raise ValueError(f"Expected at least 8 non-empty lines in {path}, got {len(lines)}")

    prefix = lines[0]
    step_id = int(lines[1])
    n_cpu_surface = int(lines[2])
    n_cpu_z = int(lines[3])
    n_surface_nodes_per_cpu = int(lines[4])
    n_z_nodes_per_cpu = int(lines[5])
    output_dir = lines[6]
    output_prefix = lines[7]

    return (
        prefix,
        step_id,
        n_cpu_surface,
        n_cpu_z,
        n_surface_nodes_per_cpu,
        n_z_nodes_per_cpu,
        output_dir,
        output_prefix,
    )


# def infer_and_convert_deg(theta: np.ndarray, phi: np.ndarray, units: str) -> Tuple[np.ndarray, np.ndarray]:
#     if units == "deg":
#         return theta, phi
#     if units == "rad":
#         return np.degrees(theta), np.degrees(phi)

#     # auto: if ranges look like radians, convert to degrees
#     theta_max = float(np.nanmax(np.abs(theta)))
#     phi_max = float(np.nanmax(np.abs(phi)))
#     if theta_max <= np.pi * 1.1 and phi_max <= (2.0 * np.pi) * 1.1:
#         return np.degrees(theta), np.degrees(phi)

#     return theta, phi


def read_r_coords(prefix: str, n_cpu_z: int, n_z_nodes_per_cpu: int) -> np.ndarray:
    r_all = []
    for cpuid in range(n_cpu_z):
        path = f"{prefix}.horiz_ave.{cpuid}.1"
        data = np.loadtxt(path, skiprows=1)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        r = data[:, 0]
        if r.size != n_z_nodes_per_cpu:
            raise ValueError(
                f"Unexpected r size in {path}: {r.size} (expected {n_z_nodes_per_cpu})"
            )
        r_all.append(r)
    return np.concatenate(r_all, axis=0)


def read_surface_coords(
    prefix: str,
    n_cpu_surface: int,
    n_cpu_z: int,
    n_surface_nodes_per_cpu: int,
) -> Tuple[np.ndarray, np.ndarray]:
    '''
    read surface coordinates (theta, phi) from all surface cpus, and convert to degrees and lat/lon.
    '''
    theta_all = []
    phi_all = []
    for k in range(1, n_cpu_surface + 1):  # cpu hori_id, cpu_id_xy, is from 1 to n_cpu_surface+1, surface cpuid is cpu_id_xy * n_cpu_z - 1
        cpuid = n_cpu_z * k - 1
        path = f"{prefix}.coord_s.{cpuid}"
        data = np.loadtxt(path, skiprows=1)  ##### one header line !
        if data.ndim == 1:  
            data = data.reshape(1, -1)  # 1 row, x column
        if data.shape[1] < 2:
            raise ValueError(f"Invalid coord data in {path}: need 2 columns")
        if data.shape[0] != n_surface_nodes_per_cpu:
            raise ValueError(
                f"Unexpected surface node count in {path}: {data.shape[0]} "
                f"(expected {n_surface_nodes_per_cpu})"
            )
        theta, phi = data[:, 0], data[:, 1]
        theta, phi = np.degrees(theta), np.degrees(phi)

        # convert theta to lat
        theta = 90.0 - theta
        theta_all.append(theta)
        phi_all.append(phi)

    return np.concatenate(theta_all, axis=0), np.concatenate(phi_all, axis=0)


def read_disp_block(
    path: str,
    n_surface_nodes_per_cpu: int,
    n_z_nodes_per_cpu: int,
) -> Tuple[np.ndarray, float, float]:
    '''
    read data for one block (from one cpu)
    assume header line being "step num_nodes time dt", and data lines being "vtheta vphi vr dtheta dphi dr ...".
    '''
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline().strip()

    # header: step, num_nodes, time, dt
    parts = header.split()
    if len(parts) < 4:
        raise ValueError(f"Unexpected header format in {path}: {header}")

    time = float(parts[2])
    dt = float(parts[3])

    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 6:
        raise ValueError(f"Expected 6 columns in {path}, got {data.shape[1]}")

    expected_rows = n_surface_nodes_per_cpu * n_z_nodes_per_cpu
    if data.shape[0] != expected_rows:
        raise ValueError(
            f"Unexpected row count in {path}: {data.shape[0]} (expected {expected_rows})"
        )

    data = data[:, :6]
    # reshape: (surface_node, z_node, component), z changes fastest
    data = data.reshape((n_surface_nodes_per_cpu, n_z_nodes_per_cpu, 6))
    return data, time, dt


def build_regular_grid(step_deg: float) -> Tuple[np.ndarray, np.ndarray]:
    theta_grid = np.arange(0.0, 180.0 + 0.5 * step_deg, step_deg)
    phi_grid = np.arange(0.0, 360.0 + 0.5 * step_deg, step_deg)

    return theta_grid, phi_grid


def main() -> int:
    parser = argparse.ArgumentParser(description="Postprocess 3D displacement to NetCDF")
    parser.add_argument("param_file", help="Path to parameter file")

    args = parser.parse_args()

    (
        prefix,
        step_id,
        n_cpu_surface,
        n_cpu_z,
        n_surface_nodes_per_cpu,
        n_z_nodes_per_cpu,
        output_dir,
        output_prefix,
    ) = read_params(args.param_file)

    n_loc_total = n_cpu_surface * n_surface_nodes_per_cpu
    n_r_total = n_cpu_z * n_z_nodes_per_cpu

    r = read_r_coords(prefix, n_cpu_z, n_z_nodes_per_cpu)
    if r.size != n_r_total:
        raise ValueError(f"Unexpected r size: {r.size} (expected {n_r_total})")

    theta, phi = read_surface_coords(
        prefix, n_cpu_surface, n_cpu_z, n_surface_nodes_per_cpu
    )
    if theta.size != n_loc_total or phi.size != n_loc_total:
        raise ValueError(
            f"Unexpected loc size: theta {theta.size}, phi {phi.size} (expected {n_loc_total})"
        )

    # Preallocate output arrays
    vr_cumu = np.zeros((n_r_total, n_loc_total), dtype=float)
    vtheta_cumu = np.zeros((n_r_total, n_loc_total), dtype=float)
    vphi_cumu = np.zeros((n_r_total, n_loc_total), dtype=float)
    dr_incr = np.zeros((n_r_total, n_loc_total), dtype=float)
    dtheta_incr = np.zeros((n_r_total, n_loc_total), dtype=float)
    dphi_incr = np.zeros((n_r_total, n_loc_total), dtype=float)

    time_val = None
    dt_val = None

    for k in range(1, n_cpu_surface + 1):
        loc_start = (k - 1) * n_surface_nodes_per_cpu
        loc_end = k * n_surface_nodes_per_cpu

        for z in range(n_cpu_z):
            cpuid = n_cpu_z * (k - 1) + z
            path = f"{prefix}.nodal_disp.{cpuid}.{step_id}"
            data, time, dt = read_disp_block(path, n_surface_nodes_per_cpu, n_z_nodes_per_cpu)

            if time_val is None:
                time_val = time
                dt_val = dt
            else:
                if time != time_val or dt != dt_val:
                    print(
                        f"Warning: time/dt mismatch in {path}: time={time}, dt={dt}",
                        file=sys.stderr,
                    )

            r_start = z * n_z_nodes_per_cpu
            r_end = (z + 1) * n_z_nodes_per_cpu

            # Convert to (z_node, surface_node)
            vtheta = data[:, :, 0].T
            vphi = data[:, :, 1].T
            vr = data[:, :, 2].T
            dtheta = data[:, :, 3].T
            dphi = data[:, :, 4].T
            dr = data[:, :, 5].T

            vtheta_cumu[r_start:r_end, loc_start:loc_end] = vtheta
            vphi_cumu[r_start:r_end, loc_start:loc_end] = vphi
            vr_cumu[r_start:r_end, loc_start:loc_end] = vr
            dtheta_incr[r_start:r_end, loc_start:loc_end] = dtheta
            dphi_incr[r_start:r_end, loc_start:loc_end] = dphi
            dr_incr[r_start:r_end, loc_start:loc_end] = dr

    if time_val is None or dt_val is None:
        raise RuntimeError("No displacement files were read")

    ds = xr.Dataset(
        data_vars={
            "Vr_cumu": (("r", "loc"), vr_cumu),
            "Vtheta_cumu": (("r", "loc"), vtheta_cumu),
            "Vphi_cumu": (("r", "loc"), vphi_cumu),
            "dr_incr": (("r", "loc"), dr_incr),
            "dtheta_incr": (("r", "loc"), dtheta_incr),
            "dphi_incr": (("r", "loc"), dphi_incr),
        },
        coords={
            "r": ("r", r),
            "lat": ("loc", theta),
            "phi": ("loc", phi),
        },
        attrs={
            "time": float(time_val),
            "dt": float(dt_val),
        },
    )

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{output_prefix}.{step_id}.nc")
    ds.to_netcdf(out_path)
    print(f"Wrote: {out_path}")

    #####################################################################
    # Interpolate to a regular NxN degree grid (nearest neighbor).
    #####################################################################
    

    

    def grid_by_nearneighbor(value_2d, lon1=0, lon2=360, lat1=-90, lat2=90, grid_step=1):
        '''
        using pygmt.nearneighbor to grid the 2D (r, loc) array to (r, theta_grid, phi_grid) array.
        by looping throught each r
        '''
        nr = value_2d.shape[0]
        output = np.zeros((nr, int((lat2-lat1)/grid_step)+1, int((lon2-lon1)/grid_step)+1), dtype=float)
        for i in range(nr):
            print("Gridding r index {}/{}...".format(i+1, nr))
            # create a temporary dataframe for pygmt
            df = np.column_stack((phi, theta, value_2d[i, :])) # lon, lat, data
            # grid the data using pygmt.nearneighbor
            grid = pygmt.nearneighbor(
                data=df,
                region=[lon1, lon2, lat1, lat2],
                spacing=grid_step,
                search_radius="200k",
            )
            # output.append(grid.values)
            output[i, :, :] = grid.values
        lon = grid.lon.values
        lat = grid.lat.values
        return output, lon, lat
    
    lon1 = 200
    lon2 = 330
    lat1 = 0
    lat2 = 90
    grid_step = 1  # the resolution for interpolated regular grid

    vr_grid, longi, lati = grid_by_nearneighbor(vr_cumu, lon1, lon2, lat1, lat2, grid_step)
    # vtheta_grid,_, _ = grid_by_nearneighbor(vtheta_cumu, lon1, lon2, lat1, lat2, grid_step)
    # vphi_grid,_, _ = grid_by_nearneighbor(vphi_cumu, lon1, lon2, lat1, lat2, grid_step)
    dr_grid,_, _ = grid_by_nearneighbor(dr_incr, lon1, lon2, lat1, lat2, grid_step)
    dtheta_grid,_, _ = grid_by_nearneighbor(dtheta_incr, lon1, lon2, lat1, lat2, grid_step)
    dphi_grid,_, _ = grid_by_nearneighbor(dphi_incr, lon1, lon2, lat1, lat2, grid_step)

    ds_grid = xr.Dataset(
        data_vars={
            "Vr_cumu": (("r", "lat", "lon"), vr_grid),
            "dr_incr": (("r", "lat", "lon"), dr_grid),
            # "Vtheta_cumu": (("r", "lat", "lon"), vtheta_grid),
            # "Vphi_cumu": (("r", "lat", "lon"), vphi_grid),
            "dtheta_incr": (("r", "lat", "lon"), dtheta_grid),
            "dphi_incr": (("r", "lat", "lon"), dphi_grid),
        },
        coords={
            "r": ("r", r),
            "lat": ("lat", lati),
            "lon": ("lon", longi),
        },
        attrs={
            "time": float(time_val),
            "dt": float(dt_val),
            "interp_method": "nearest from pygmt",
            "grid_step_deg": grid_step,
        },
    )

    out_path_grid = os.path.join(output_dir, f"{output_prefix}.{step_id}.grid_lon{lon1}-{lon2}_lat{lat1}-{lat2}_res{grid_step:.1f}.nc")
    ds_grid.to_netcdf(out_path_grid)
    print(f"Wrote: {out_path_grid}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
