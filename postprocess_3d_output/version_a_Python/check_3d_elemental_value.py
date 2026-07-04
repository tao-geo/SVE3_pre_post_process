'''
check mat group debug file
'''

import numpy as np
from typing import Tuple
import xarray as xr
import sys
import os

# CONFIG = {
#     prefix : ,
#     step_id : ,
#     n_cpu_surface : ,
#     n_cpu_z : ,
#     nxy : , # node in x or y direction per cpu
#     n_z_nodes_per_cpu : ,
#     output_dir : ,
#     output_prefix : ,
# }

prefix = '/glade/derecho/scratch/taoyuan/Proj_RSL_3Dvisc/cases_nonlinear/case_a1/case_a1'
step_id = None
n_cpu_surface = 192
n_cpu_z = 1
nxy = 33 # node in x or y direction per cpu
n_z_nodes_per_cpu = 65
output_dir = '/glade/work/taoyuan/Proj_RSL_3Dvisc/CASES_nonNewtonian/'
output_prefix = 'case_a1.matgroup.'

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

def read_element_r_coords(prefix: str, n_cpu_z: int, n_z_nodes_per_cpu: int) -> np.ndarray:
    '''
    get the r coord of element center
    '''
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
        r_all.append((r[:-1] + r[1:]) / 2.0)
    return np.concatenate(r_all, axis=0)

def read_surface_coords(
    prefix: str,
    n_cpu_surface: int,
    n_cpu_z: int,
    nxy: int,
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
        if data.shape[0] != (nxy**2):
            raise ValueError(
                f"Unexpected surface node count in {path}: {data.shape[0]} "
                f"(expected {(nxy**2)})"
            )
        theta, phi = data[:, 0], data[:, 1]
        theta, phi = np.degrees(theta), np.degrees(phi)

        # convert theta to lat
        theta = 90.0 - theta
        theta_all.append(theta)
        phi_all.append(phi)

    return np.concatenate(theta_all, axis=0), np.concatenate(phi_all, axis=0)

def read_element_surface_coords(
    prefix: str,
    n_cpu_surface: int,
    n_cpu_z: int,
    nxy: int,
) -> Tuple[np.ndarray, np.ndarray]:
    '''
    read element's surface coordinates (theta, phi) from all surface cpus, and convert to degrees and lat/lon.
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
        if data.shape[0] != (nxy**2):
            raise ValueError(
                f"Unexpected surface node count in {path}: {data.shape[0]} "
                f"(expected {(nxy**2)})"
            )
        # first, convert to 2-D array
        data = data.reshape((nxy, nxy, 2))  
        data_center = (data[:-1, :-1, :] + data[1:, :-1, :] + data[:-1, 1:, :] + data[1:, 1:, :]) / 4.0
        theta = data_center[:, :, 0].flatten()
        phi = data_center[:, :, 1].flatten()

        # theta, phi = data[:, 0], data[:, 1]
        theta, phi = np.degrees(theta), np.degrees(phi)

        # convert theta to lat
        theta = 90.0 - theta
        theta_all.append(theta)
        phi_all.append(phi)

    return np.concatenate(theta_all, axis=0), np.concatenate(phi_all, axis=0)

def read_disp_block(
    path: str,
    nxy: int,
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

    expected_rows = (nxy**2) * n_z_nodes_per_cpu
    if data.shape[0] != expected_rows:
        raise ValueError(
            f"Unexpected row count in {path}: {data.shape[0]} (expected {expected_rows})"
        )

    data = data[:, :6]
    # reshape: (surface_node, z_node, component), z changes fastest
    data = data.reshape(((nxy**2), n_z_nodes_per_cpu, 6))
    return data, time, dt

def read_matgroup_block(
    path: str,
    nxy: int,
    n_z_nodes_per_cpu: int,
):
    '''
    read data for one block (from one cpu)
    assume header line being "step num_nodes time dt", and data lines being "vtheta vphi vr dtheta dphi dr ...".
    '''
    # with open(path, "r", encoding="utf-8") as f:
    #     header = f.readline().strip()

    # # header: step, num_nodes, time, dt
    # parts = header.split()
    # if len(parts) < 4:
    #     raise ValueError(f"Unexpected header format in {path}: {header}")

    # time = float(parts[2])
    # dt = float(parts[3])

    data = np.loadtxt(path)
    # if data.ndim == 1:
    #     data = data.reshape(1, -1)
    # if data.shape[1] != 1:
    #     raise ValueError(f"Expected 1 columns in {path}, got {data.shape[1]}")

    expected_rows = ((nxy-1)**2) * (n_z_nodes_per_cpu-1)
    if data.shape[0] != expected_rows:
        raise ValueError(
            f"Unexpected row count in {path}: {data.shape[0]} (expected {expected_rows})"
        )

    # data = data[:, :6]
    # reshape: (surface_node, z_node, component), z changes fastest
    data = data.reshape((((nxy-1)**2), (n_z_nodes_per_cpu-1)))
    return data #, time, dt

def save_matgroup() -> int:

    n_loc_total = n_cpu_surface * ((nxy-1)**2)
    n_r_total = n_cpu_z * (n_z_nodes_per_cpu - 1)

    r = read_element_r_coords(prefix, n_cpu_z, n_z_nodes_per_cpu)
    if r.size != n_r_total:
        raise ValueError(f"Unexpected r size: {r.size} (expected {n_r_total})")

    theta, phi = read_element_surface_coords(
        prefix, n_cpu_surface, n_cpu_z, nxy
    )
    if theta.size != n_loc_total or phi.size != n_loc_total:
        raise ValueError(
            f"Unexpected loc size: theta {theta.size}, phi {phi.size} (expected {n_loc_total})"
        )

    # Preallocate output arrays
    matgroup = np.zeros((n_r_total, n_loc_total), dtype=float)


    for k in range(1, n_cpu_surface + 1):
        loc_start = (k - 1) * ((nxy-1)**2)
        loc_end = k * ((nxy-1)**2)

        for z in range(n_cpu_z):
            cpuid = n_cpu_z * (k - 1) + z
            path = f"{prefix}.mat_group.{cpuid}"
            # data, time, dt = read_disp_block(path, (nxy**2), n_z_nodes_per_cpu)
            data = read_matgroup_block(path, nxy, n_z_nodes_per_cpu)

            r_start = z * n_z_nodes_per_cpu
            r_end = (z + 1) * n_z_nodes_per_cpu

            # Convert to (z_node, surface_node)
            # vtheta = data[:, :, 0].T
            # vphi = data[:, :, 1].T
            # vr = data[:, :, 2].T
            # dtheta = data[:, :, 3].T
            # dphi = data[:, :, 4].T
            # dr = data[:, :, 5].T

            matgroup[r_start:r_end, loc_start:loc_end] = data[:,:].T

    ds = xr.Dataset(
        data_vars={
            "matgroup": (("r", "loc"), matgroup),
        },
        coords={
            "r": ("r", r),
            "lat": ("loc", theta),
            "phi": ("loc", phi),
        },
    )

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{output_prefix}.nc")
    ds.to_netcdf(out_path)
    print(f"Wrote: {out_path}")

if __name__ == '__main__':
    save_matgroup()