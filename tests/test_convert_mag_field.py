import unittest
import numpy as np
import scipy
import matplotlib.pyplot as plt
from DSE.ADCS_n_Propulsion.convert_magnetic_field import avg_torque_calc, sizing_dipole, res_torques_calc, \
    integrate_torques, angular_momentum_realism_creator, get_sizing_from_angular_momentum, mag_field_creator

class TestConvertMagField(unittest.TestCase):
    def test_avg_torque_calc(self):
        np.testing.assert_almost_equal(avg_torque_calc(np.array([0, 0, 0]), np.array([0, 0, 0]), 0, 1000), np.array([0, 0, 0]))
        np.testing.assert_almost_equal(avg_torque_calc(np.array([1, 1, 1]), np.array([0.5, 0, 0.2]), 200, 1000), np.array([1, 1, 1]))
        np.testing.assert_almost_equal(avg_torque_calc(np.array([1, 1, 1]), np.array([2, 0, 2]), 1000, 1000), np.array([2, 1, 2]))

    def test_sizing_dipole(self):
        t_orbits_detumbling = 40*5423
        m = 0.5  # tesla
        v_r = 300 / 180 * np.pi  # rad /s
        torque_required = v_r / (t_orbits_detumbling * m)

        np.testing.assert_almost_equal(sizing_dipole(np.array([[0, 0, 0], [2, 2, 2]]), np.array([1, 1, 1]), np.array([1, 1, 1]),
                      t_orbits_detumbling), np.array([1, 1, 1]))
        np.testing.assert_almost_equal(sizing_dipole(np.array([[0, 0, 0], [1, 1, 1]]), np.array([0, 0, 0]), np.array([1, 1, 1]),
                                           t_orbits_detumbling), np.array([torque_required, torque_required,
                                                                           torque_required]))

    def test_res_torques_calc(self):
        np.testing.assert_almost_equal(res_torques_calc(np.array([4,5,6]), np.array([0.5, 10, 0]), np.array([1,2,3])),
                                       np.array([-3.5, 0, -18]))

    def test_integrate_torques(self):
        print(integrate_torques(np.array([[1, 2, 3], [4, 5, 6]]).T))
        np.testing.assert_almost_equal(integrate_torques(np.array([[1, 2, 3], [4, 5, 6]])), np.array([[2.5, 3.5,
                                                                                                       4.5]]))

    def test_angular_momentum_realism_creator(self):
        x = np.array([[0, 1, 2, 0, -1, -3, -2, 0, 1], [0,0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0]]).T
        res = np.array([[0, 1, 2, 0, 0, 0, 1, 3, 4], [0,0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0]]).T
        np.testing.assert_almost_equal(angular_momentum_realism_creator(x), res)

    def test_get_sizing_from_angular_momentum(self):
        x = np.array([[1,76,3], [3,6,8], [8,3,6]])
        np.testing.assert_almost_equal(get_sizing_from_angular_momentum(x), np.array([16, 152, 16]))

    def test_angular_momentum_calc(self):
        # x = np.array([[1, 3, 3], [3, 6, 8], [8, 3, 6]])
        # plt.plot(x)
        # plt.plot(scipy.integrate.cumtrapz(x, axis=0)*2)
        # plt.show()
        pass



