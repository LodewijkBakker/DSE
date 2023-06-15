import unittest
import numpy as np
import DSE.Thermal.thermal_dyn_sim as tds
from DSE.Thermal.materials import *
from DSE.Thermal.env import Env


class TestThermalSimulation(unittest.TestCase):
    def test_thermal_node(self):
        TN = tds.ThermalNode(0, 'node', {'area': 1, 'radiation_area': 2, 'contact_area': 3, 'mass': 4,
                                         'heat_generated': 0}, aluminium(), bare_Al())
        self.assertEqual(TN.node_id, 0)
        self.assertEqual(TN.name, 'node')
        self.assertEqual(TN.area, 1)
        self.assertEqual(TN.radiation_area, 2)
        self.assertEqual(TN.contact_area, 3)
        self.assertEqual(TN.mass, 4)
        self.assertEqual(TN.heat_generated, 0)
        self.assertEqual(TN.thermal_capacitance, 897)
        self.assertEqual(TN.emissivity, 0.15)

    def test_thermal_model(self):
        TM = tds.ThermalModel([], [], Env(), 200, n_orbits=10, sun_shield=False, unit='C', ESATAN=False)
        self.assertEqual(TM.nodes, [])
        self.assertEqual(TM.connections, [])
        self.assertEqual(TM.init_temp, 200)
        self.assertEqual(TM.n_orbits, 10)
        self.assertEqual(TM.sun_shield, False)
        self.assertEqual(TM.unit, 'C')
        self.assertEqual(TM.ESATAN, False)
        self.assertEqual(TM.solution, None)

    def test_thermal_solution(self):
        TN1 = tds.ThermalNode(0, 'node1', {'area': 0.5, 'radiation_area': 1, 'contact_area': [0, 0.1], 'mass': 1,
                                           'heat_generated': 10}, aluminium(), white_paint())
        TN2 = tds.ThermalNode(1, 'node2', {'area': 1, 'radiation_area': 1, 'contact_area': [0.1, 0], 'mass': 1,
                                           'heat_generated': 10}, aluminium(), bare_Al())
        nodes = [TN1, TN2]
        TC = np.array([[1], [0]])
        Q = [np.ones(Env().t_orbit), np.ones(Env().t_orbit)]
        TM = tds.ThermalModel(nodes, TC, Env(), [250, 260], n_orbits=1, sun_shield=False, unit='K', ESATAN=True,
                              Q_ESATAN=Q)
        TM.solve()
        self.assertAlmostEqual(TM.Tdot_arr[0][0], 0.05154, delta=1e-4)


