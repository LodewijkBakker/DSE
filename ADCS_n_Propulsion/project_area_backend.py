import argparse
import os
import cProfile

import numpy as np
# import stl
# from stl import mesh

from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely import centroid # correct?

from scipy.spatial.transform import Rotation as R
import multiprocessing as mp


def _area2(a, b, c):
    """Calculate twice the area of a 2D triangle

    Parameters
    ----------
    a, b, c : tuple
        vertex tuples containing two floats (x, y)

    Returns
    -------
    float
        twice the triangle area
    References
    ----------
    [1] Computational Geometry in C by Joseph O'Rourke
    """
    # --- constants ----------------------------------------------------------------+

    x, y = [0, 1]
    return (
            (b[x] - a[x]) * (c[y] - a[y]) -
            (c[x] - a[x]) * (b[y] - a[y]))


def _munge_triangle(v1, v2, v3, min_area=0):
    """MUNGE TRIANGLE

    For every new triangle that gets added, it  needs to pass some basic
    quality controls and if need be, fixed.

    1) Triangles with zero area are automatically rejected.

    2) Triangles must have three unique vertices.

    3) All triangles will be returned with a CCW winding regardless of the
    winding upon input.

    4) Final output is a list of three edges, consisting of vertex master list
    indices.

    Inputs
    ------
    v1, v2, v3 : tuples

    Returns
    -------
    t : a list consisting of three edge tuples
    """

    # first check if the triangle has three unique vertices
    if v1 == v2 or v2 == v3 or v3 == v1:
        return []

    a2 = _area2(v1, v2, v3)

    # if the triangle has cw winding, reverse the winding order, this perhaps can be removed
    if a2 < 0:
        t_verts = [v1, v3, v2]

    # the winding order is correct, leave vertices alone
    if a2 > 0:
        t_verts = [v1, v2, v3]

    # minimum triangle area check
    if abs(a2 * 0.5) <= min_area:
        return []

    return t_verts


def _project_triangle_to_plane(tri_3d, rot_matrix):
    """PROJECT TRIANGLE TO ARBITRARY PLANE
    Parameters
    ----------
    tri_3d : list

    Returns
    -------
    list
    """
    return [(rot_matrix.dot(v)[1], rot_matrix.dot(v)[2]) for v in tri_3d]


def _read_stl(stl_model):
    """
    Reads both ASCII and Binary STL Files
    """

    for i in stl_model.points:
        triangle = i.reshape(3, 3).astype(np.float16)
        yield triangle


def _load_stl(stl_model, rot_angles):
    """LOAD STL
    Parameters
    ----------
    stl_model : numpy mesh
    pvec : str
    Returns
    -------
    list of shapely polygon objects
    """

    # assemble stl
    triangles = []
    # The following matrix is used to rotate the triangle as to represent a rotated object
    rot_matrix = R.from_euler('zyx', rot_angles, degrees=True).as_matrix()

    for i, t in enumerate(_read_stl(stl_model)):
        v1, v2, v3 = _project_triangle_to_plane(t, rot_matrix)
        verts = _munge_triangle(v1, v2, v3)
        if len(verts) == 3:
            triangles.append(Polygon(verts))

    return triangles


def _calculate_projected_area(rot_angles, stl_model):
    """CALCULATE PROJECTED AREA
    Parameters
    ----------
    stl_model : numpy mesh
    pvec : str
    Returns
    -------
    float
    """
    # bulk load all triangles from stl
    triangles = _load_stl(stl_model, rot_angles)

    # merge triangles into single polygon
    sub_poly = unary_union(triangles)

    grav_center = [1, 1, 1]
    geo_center = centroid(sub_poly)
    arm = _calculate_arm(rot_angles, grav_center, geo_center)

    return [sub_poly.area, rot_angles]


def _calculate_arm(rot_angles, grav_center, geo_center):
    rot_matrix = R.from_euler('zyx', rot_angles, degrees=True).as_matrix()
    grav_center_projected = [(rot_matrix.dot(grav_center)[1], rot_matrix.dot(grav_center)[2])]
    arm = ((grav_center_projected[0] - geo_center[0])**2 + (grav_center_projected[1] - geo_center[1])**2)**0.5
    return arm

# --- MAIN ---------------------------------------------------------------------+


# def calculate_area_multi_processing_start(stl_path, rot_angles):
#     #stl_path = "satellite2.stl" # "testfileformat.stl"
#     stl_model_main = mesh.Mesh.from_file(stl_path)
#
#     ng = 1
#     list_rot_angles = []
#     for i in rot_angles:
#         list_rot_angles.append([i, stl_model_main])
#
#     pool = mp.Pool()
#     pool = mp.Pool(processes=6)
#     res = pool.starmap(_calculate_projected_area, list_rot_angles)
#     pool.close()
#     #print(res)
#     return res
#     # 24.2, 0 degree angles for fokker g1

# get centroid, distance from 0.0