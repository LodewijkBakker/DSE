import argparse
import os
import cProfile

import numpy as np
# import stl
# from stl import mesh

from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely import centroid # correct?
from shapely import Point

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
    """PROJECT TRIANGLE TO ARBITRARY PLANE -> AFTER ROTATION MATRIX APPLICATION
     PROJECT ON YZ PLANE
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


def _calculate_projected_area(rot_angle, stl_model, list_grav_centers):
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
    triangles = _load_stl(stl_model, rot_angle)

    # merge triangles into single polygon
    sub_poly = unary_union(triangles)  # runtime warning but should still store it correctly?

    geo_center = centroid(sub_poly)

    area_arm_list = []
    for grav_center in list_grav_centers:
        arm = _calculate_arm(rot_angle, grav_center, geo_center)
        area_arm_list.append(arm)

    return [sub_poly.area, area_arm_list]


def _calculate_arm(rot_angle, grav_center, geo_center):
    rot_matrix = R.from_euler('zyx', rot_angle, degrees=True).as_matrix()
    grav_center_projected = rot_matrix.dot(grav_center)
    #print(grav_center_projected)
    arm_projected = [0, grav_center_projected[1] - geo_center.x, grav_center_projected[2] - geo_center.y]
    # geo_center.x is not actually x and neither is y. Reason is because it is projected on
    # yz axis and shapely thinks this is xy axis
    #print(arm_projected)
    rot_matrix_transpose = rot_matrix.transpose()  # test this works
    arm_body_frame = rot_matrix_transpose.dot(arm_projected)  # matrix multiplication

    #print(grav_center_projected)
    #arm = ((grav_center_projected[0] - geo_center.x)**2 + (grav_center_projected[1] - geo_center.y)**2)**0.5
    # arm = []
    return arm_body_frame


def unit_tests():
    geo_center_test = Point(1, 0)
    print(_calculate_arm([0, 0, 0], [1, 1, 0], geo_center_test))

    geo_center_test = Point(0, 0)
    assert(np.all([0, 0, 0] == _calculate_arm([0, 0, 0], [0, 0, 0], geo_center_test)))
    assert(np.all([0, 1, 0] == _calculate_arm([0, 0, 0], [1, 1, 0], geo_center_test)))
    # for distance of center of pressure to cg z axis doesn't matter (should be x axis)
    assert(np.all([0, 1, 0] == _calculate_arm([0, 0, 0], [0, 1, 0], geo_center_test)))
    assert(np.all([0, 1, 1] == _calculate_arm([0, 0, 0], [0, 1, 1], geo_center_test)))

    # rotated 360 should still be the same
    assert(np.all(np.isclose([0, 0, 0],_calculate_arm([0, 0, 360], [0, 0, 0], geo_center_test))))
    assert(np.all(np.isclose([0, 1, 0], _calculate_arm([0, 360, 0], [1, 1, 0], geo_center_test))))
    # for distance of center of pressure to cg z axis doesn't matter (should be x axis)
    assert(np.all(np.isclose([0, 1, 0], _calculate_arm([360, 0, 0], [0, 1, 0], geo_center_test))))
    assert(np.all(np.isclose([0, 1, 1], _calculate_arm([360, 0, 0], [0, 1, 1], geo_center_test))))


if __name__ == "__main__":
    unit_tests()

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