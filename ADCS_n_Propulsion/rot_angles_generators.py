
import numpy as np
from scipy.spatial.transform import Rotation as R


def in_val_octant(point_coord, octant_sel):
    # Z is a coordinate
    # Octants selected = array[0,1,2,3,4,5,6,7] where 1 2 3 are filled with booleans
    # octant 1 is +++ then clockwise where z is constant + and then when rotated in full
    # do the same but with z = -

    # Octant 1 X+Y+Z+
    if octant_sel[0] and point_coord[0] >= 0 and point_coord[1] >= 0 and point_coord[2] >= 0:
        return True
    # Octant 2 X-Y+Z+
    elif octant_sel[1] and point_coord[0] <= 0 <= point_coord[1] and point_coord[2] >= 0:
        return True
    # Octant 3 X-Y-Z+
    elif octant_sel[2] and point_coord[0] <= 0 and point_coord[1] <= 0 and point_coord[2] >= 0:
        return True
    # Octant 4 X+Y-Z+
    elif octant_sel[3] and point_coord[0] >= 0 >= point_coord[1] and point_coord[2] >= 0:
        return True
    # Octant 5 X+Y+Z-
    elif octant_sel[4] and point_coord[0] >= 0 and point_coord[1] >= 0 and point_coord[2] <= 0:
        return True
    # Octant 6 X-Y+Z-
    elif octant_sel[5] and point_coord[0] <= 0 <= point_coord[1] and point_coord[2] <= 0:
        return True
    # Octant 7 X-Y-Z-
    elif octant_sel[6] and point_coord[0] <= 0 and point_coord[1] <= 0 and point_coord[2] <= 0:
        return True
    # Octant 8 X+Y-Z-
    elif octant_sel[7] and point_coord[0] >= 0 >= point_coord[1] and point_coord[2] <= 0:
        return True
    else:
        return False

def pole_accurate_angles(ng, octant_sel):
    # it is difficult to do accurately number of point so this is about precise
    # shitty prime numbers
    # You cannot return a value when using a callback function unfortunately
    # so it has to be updated this way

    rot_angles = []
    angle_step_size = 360 / ng
    # needed to satisfy ng for the octants selected
    for angle_theta_sphere in np.arange(0, 360, angle_step_size):
        polar_rotation = R.from_euler('zyx', [angle_theta_sphere, 0, 0], degrees=True)
        # taitbryan
        if in_val_octant(polar_rotation.apply([1, 0, 0]), octant_sel):
            rot_angles.append([angle_theta_sphere, 0, 0])
            #rot_angles.append([0, 0, 30])

    return rot_angles


def sphere_fibonacci_grid_points(ng: int):
    # *****************************************************************************80
    #
    # SPHERE_FIBONACCI_GRID_POINTS: Fibonacci spiral gridpoints on a sphere.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    15 May 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Richard Swinbank, James Purser,
    #    Fibonacci grids: A novel approach to global modelling,
    #    Quarterly Journal of the Royal Meteorological Society,
    #    Volume 132, Number 619, July 2006 Part B, pages 1769-1793.
    #
    #  Parameters:
    #
    #    Input, integer NG, the number of points.
    #
    #    Output, real XG(3,N), the grid points.
    #

    phi = (1.0 + np.sqrt(5.0)) / 2.0

    theta = np.zeros(ng)
    sphi = np.zeros(ng)
    cphi = np.zeros(ng)

    # creating points
    for i in range(0, ng):
        i2 = 2 * i - (ng - 1)
        theta[i] = 2.0 * np.pi * float(i2) / phi
        sphi[i] = float(i2) / float(ng)
        cphi[i] = np.sqrt(float(ng + i2) * float(ng - i2)) / float(ng)

    return theta, sphi, cphi


def sphere_fib_octant(ng, octant_sel):

    n_select_octant = 0
    for i in octant_sel:
        if i:
            n_select_octant += 1

    if n_select_octant != 0:
        n_points_fibo = int(8/n_select_octant * ng)

    theta, sphi, cphi = sphere_fibonacci_grid_points(n_points_fibo)

    #rot_angles = np.zeros((n_points_fibo, 3))
    rot_angles = []
    # transforming points in rotation and discarding if not within range
    for i in range(0, n_points_fibo):
        # Cross product of 1,0,0 and point on sphere = vector u
        # Norm of vector u
        # This vector * invcos(vector_u_dot x_axis_vector) (basically vector * angle between
        # vectors and as vector_u and x_axis vector are going to be norms no length must be
        # taken into account)

        # Get euler angles from this using scipy
        # Write it to xg

        u1 = cphi[i] * np.sin(theta[i])
        u2 = cphi[i] * np.cos(theta[i])
        u3 = sphi[i]
        z = [u1, u2, u3]

        if in_val_octant(z, octant_sel):  # rot angles still initialised in total!!!
            # rotation vector should be in radians,
            fib_point_angle = np.arccos(np.dot([u1, u2, u3], [1, 0, 0]))
            # z_rot_vector = np.cross([1, 0, 0], [u1, u2, u3]) / np.linalg.norm(
            #     z) * fib_point_angle

            q_vector = np.cross([1, 0, 0], [u1, u2, u3])
            q_rot = np.sqrt(np.linalg.norm([1, 0, 0]) * np.linalg.norm([u1, u2, u3])) \
                    + np.dot([1, 0, 0], [u1, u2, u3])
            quat_total = [q_vector[0], q_vector[1], q_vector[2], q_rot]
            # print(quat_total) # this would show that the above will not work if
            # specific indixes are not taken "(
            rotation_from_quat = R.from_quat(quat_total)

            # rotation_from_rot_vector = R.from_rotvec(z_rot_vector)
            # euler_rotation_z = rotation_from_rot_vector.as_euler('zyx', degrees=True)
            euler_rotation_z = rotation_from_quat.as_euler('zyx', degrees=True)
            rot_angles.append(euler_rotation_z)
            #[euler_rotation_z[0], euler_rotation_z[1], euler_rotation_z[2]]
            # rot_angles[i, 0] = euler_rotation_z[0]
            # rot_angles[i, 1] = euler_rotation_z[1]
            # rot_angles[i, 2] = euler_rotation_z[2]

    return np.array(rot_angles)

def sphere_fib_user_angle(ng):
    # less efficient
    #
    yaw_range = [-1, 1]
    pitch_range = [-1, 1]
    roll_range = [-31, 31]

    rot_angles = []
    for i in np.linspace(yaw_range[0], yaw_range[1], 10):
        for j in np.linspace(pitch_range[0], pitch_range[1], 10):
            for k in np.linspace(roll_range[0], roll_range[1], 30):
                rot_angles.append([i, j, k])

    return rot_angles

