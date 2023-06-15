import unittest
import numpy as np
import scipy
import matplotlib.pyplot as plt
from DSE.ADCS_n_Propulsion.convert_magnetic_field import avg_torque_calc, sizing_minimum_dipole, res_torques_calc, \
    integrate_torques, angular_momentum_realism_creator, get_sizing_from_angular_momentum, \
    sizing_angular_momentum_calc, sizing_cmg, sizing_magnetorquer, angular_momentum_calc, optimum_sizer

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

        np.testing.assert_almost_equal(sizing_minimum_dipole(np.array([[0, 0, 0], [2, 2, 2]]), np.array([1, 1, 1]), np.array([1, 1, 1]),
                      t_orbits_detumbling), np.array([1, 1, 1]))
        np.testing.assert_almost_equal(sizing_minimum_dipole(np.array([[0, 0, 0], [1, 1, 1]]), np.array([0, 0, 0]), np.array([1, 1, 1]),
                                           t_orbits_detumbling), np.array([torque_required, torque_required,
                                                                           torque_required]))

    def test_res_torques_calc(self):
        np.testing.assert_almost_equal(res_torques_calc(np.array([4,5,6]), np.array([0.5, 10, 0]), np.array([1,2,3])),
                                       np.array([-3.5, 0, -18]))

    def test_integrate_torques(self):
        np.testing.assert_almost_equal(integrate_torques(np.array([[1, 2, 3], [4, 5, 6]])), np.array([[2.5, 3.5,4.5]]))

    def test_angular_momentum_realism_creator(self):
        x = np.array([[0, 1, 2, 0, -1, -3, -2, 0, 1], [0,0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0]]).T
        res = np.array([[0, 1, 2, 0, 0, 0, 1, 3, 4], [0,0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0,0]]).T
        np.testing.assert_almost_equal(angular_momentum_realism_creator(x), res)

    def test_get_sizing_from_angular_momentum(self):
        x = np.array([[1,76,3], [3,6,8], [8,3,6]])
        np.testing.assert_almost_equal(get_sizing_from_angular_momentum(x), np.array([8, 76, 8]))

    def test_angular_momentum_calc(self):
        a = np.array([[1, 1, 4, 1, 0, 1, 0, 8, 4, 6], [1, 1, 4, 1, 0, 1, 0, 8, 4, 6], [1, 1, 4, 1, 0, 1, 0, 8, 4, 6]]).T
        d = np.ones((10, 3))
        m = np.array([[1, 1, 2, 1, 4, -1, 6, 0, 8, 0], [1, 1, 2, 1, 4, -1, 6, 0, 8, 0], [1, 1, 2, 1, 4, -1, 6, 0, 8, 0]]).T
        np.testing.assert_almost_equal(angular_momentum_calc(m, a, d), [4,4,4])

    def test_sizing_angular_momentum_calc(self):
        max_ang_mom = [10, 20, 30]
        self.assertEqual(sizing_angular_momentum_calc(max_ang_mom), 125/16)

    def test_sizing_cmg(self):
        sam, radius, max_ang_mom = 0.005, 0.1, [1,2,1]
        np.testing.assert_almost_equal(sizing_cmg(max_ang_mom, r_wheel=radius, sizing_angular_momentum=sam),
                                       [0.0127324, 0.0114459])
        np.testing.assert_almost_equal(sizing_cmg(max_ang_mom, r_wheel=radius), [1.352817, 1.2161307])

    def test_sizing_magnetorquer(self):
        np.testing.assert_almost_equal(sizing_magnetorquer(np.array([1,1,1])), [0.045, 0.189])
        self.assertRaises(AssertionError, sizing_magnetorquer, np.array([0,0,0]))

    def test_optimum_sizer(self):
        a = np.ones((50, 3))*2
        d = np.ones((50, 3))*1
        m = np.ones((50, 3))*1
        np.testing.assert_almost_equal(optimum_sizer(m, a, d), [0.5, 0.5, 0.5])
