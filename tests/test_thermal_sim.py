import unittest
import numpy as np
import DSE.Thermal.thermal_dyn_sim as tds
from DSE.Thermal.materials import Material
from DSE.Thermal.env import Env


class TestThermalSimulation(unittest.TestCase):
    def test_thermal_node(self):
        TN = tds.ThermalNode(0, 'node', {'area': 1, 'radiation_area': 2, 'contact_area': 3, 'mass': 4, 'heat_generated':0}, Material().aluminium())
        self.assertEqual(TN.node_id, 0)
        self.assertEqual(TN.name, 'node')
        self.assertEqual(TN.area, 1)
        self.assertEqual(TN.radiation_area, 2)
        self.assertEqual(TN.contact_area, 3)
        self.assertEqual(TN.mass, 4)
        self.assertEqual(TN.heat_generated, 0)

    def test_thermal_model(self):
        TM = tds.ThermalModel([], [], Env(), 1000, 200)
        self.assertEqual(TM.nodes, [])
        self.assertEqual(TM.connections, [])
        self.assertEqual(TM.t_sim, 1000)
        self.assertEqual(TM.init_temp, 200)
        self.assertEqual(TM.solution, None)

    def test_thermal_solution(self):
        TN1 = tds.ThermalNode(0, 'node1', {'area': 1, 'radiation_area': 1, 'contact_area': 2, 'mass': 4, 'heat_generated': 10}, Material().aluminium())
        TN2 = tds.ThermalNode(1, 'node2', {'area': 2, 'radiation_area': 2, 'contact_area': 2, 'mass': 4, 'heat_generated': 10}, Material().aluminium())
        nodes = [TN1, TN2]
        TC = np.array([[1],
                       [0]])
        TM = tds.ThermalModel(nodes, TC, Env(), 1000, [200, 200])
        TM.solve()
        self.assertAlmostEqual(TM.solution.y[0][-1], 280, delta=1)

