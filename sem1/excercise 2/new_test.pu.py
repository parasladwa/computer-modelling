import unittest
import numpy as np
from particle3D import Particle3D

class TestParticle3D(unittest.TestCase):

    def test_init(self):
        p = Particle3D("A", 1.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

        self.assertEqual(p.label, "A")
        self.assertAlmostEqual(p.mass, 1.0, places=6)
        np.testing.assert_array_equal(p.position, np.array([0.0, 0.0, 0.0]))
        np.testing.assert_array_equal(p.velocity, np.array([0.0, 0.0, 0.0]))

    def test_str(self):
        p = Particle3D("A", 1.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
        expected_str = "A 0.0 0.0 0.0"

        self.assertEqual(str(p), expected_str)

    def test_kinetic_energy(self):
        p = Particle3D("A", 1.0, [0.0, 0.0, 0.0], [1.0, 2.0, 3.0])

        ke = p.kinetic_energy()

        self.assertAlmostEqual(ke, 14.5, places=6)

    def test_momentum(self):
        p = Particle3D("A", 2.0, [1.0, 2.0, 3.0], [1.0, 1.0, 1.0])

        momentum = p.momentum()

        np.testing.assert_array_equal(momentum, np.array([2.0, 2.0, 2.0]))

    def test_update_position_1st(self):
        p = Particle3D("A", 1.0, [1.0, 2.0, 3.0], [1.0, 1.0, 1.0])
        dt = 0.5

        p.update_position_1st(dt)

        np.testing.assert_array_equal(p.position, np.array([1.5, 2.5, 3.5]))

    def test_update_velocity(self):
        p = Particle3D("A", 2.0, [1.0, 2.0, 3.0], [1.0, 1.0, 1.0])
        dt = 0.5
        force = np.array([0.5, 0.5, 0.5])

        p.update_velocity(dt, force)

        np.testing.assert_array_equal(p.velocity, np.array([1.25, 1.25, 1.25]))

    def test_read_line(self):
        line = "B 2.0 1.0 2.0 3.0 0.5 0.5 0.5"
        p = Particle3D.read_line(line)

        self.assertEqual(p.label, "B")
        self.assertAlmostEqual(p.mass, 2.0, places=6)
        np.testing.assert_array_equal(p.position, np.array([1.0, 2.0, 3.0]))
        np.testing.assert_array_equal(p.velocity, np.array([0.5, 0.5, 0.5]))

    def test_total_kinetic_energy(self):
        particles = [
            Particle3D("A", 1.0, [0.0, 0.0, 0.0], [1.0, 2.0, 3.0]),
            Particle3D("B", 2.0, [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]),
            Particle3D("C", 3.0, [2.0, 2.0, 2.0], [2.0, 2.0, 2.0])
        ]

        total_ke = Particle3D.total_kinetic_energy(particles)

        self.assertAlmostEqual(total_ke, 58.0, places=6)

    def test_com_velocity(self):
        particles = [
            Particle3D("A", 1.0, [1.0, 2.0, 3.0], [1.0, 1.0, 1.0]),
            Particle3D("B", 2.0, [2.0, 2.0, 2.0], [2.0, 2.0, 2.0])
        ]

        com_vel = Particle3D.com_velocity(particles)

        np.testing.assert_array_equal(com_vel, np.array([1.33333333, 1.33333333, 1.33333333]))

if __name__ == '__main__':
    unittest.main()
