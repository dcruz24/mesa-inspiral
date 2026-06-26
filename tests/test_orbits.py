import unittest
from types import SimpleNamespace

import numpy as np

from src.config import constants as const
from src.physics import orbits


class OrbitIntegratorTests(unittest.TestCase):
    def test_initialize_two_body_state_from_circular_scalars(self):
        sim = SimpleNamespace(a=10.0, v=2.0)
        orbits.initialize_orbit_state(sim)
        self.assertEqual((sim.x, sim.y, sim.vx, sim.vy), (10.0, 0.0, 0.0, 2.0))

    def test_circular_step_preserves_existing_euler_update(self):
        sim = SimpleNamespace(m1=1.0, a=10.0, th=0.0, Om=1.0, da_dt=-0.5, th_dt=2.0)
        accretor = SimpleNamespace(dm_dt=0.1, dOm_dt=0.2)
        orbits.step_circular_euler(sim, accretor, 4.0)
        self.assertAlmostEqual(sim.m1, 1.4)
        self.assertAlmostEqual(sim.a, 8.0)
        self.assertAlmostEqual(sim.th, 8.0)
        self.assertAlmostEqual(sim.Om, 1.8)

    def test_two_body_step_updates_relative_state(self):
        sim = SimpleNamespace(
            m1=1.4 * const.Msun_cgs,
            m2=10.0 * const.Msun_cgs,
            a=10.0 * const.Rsun_cgs,
            v=np.sqrt(const.G_cgs * 11.4 * const.Msun_cgs / (10.0 * const.Rsun_cgs)),
            Fd=0.0,
            Om=0.0,
        )
        accretor = SimpleNamespace(dm_dt=0.0, dOm_dt=0.0)
        orbits.initialize_orbit_state(sim)
        old_x = sim.x
        orbits.step_two_body_euler(sim, accretor, 1.0)
        self.assertTrue(np.isfinite(sim.a))
        self.assertTrue(np.isfinite(sim.v))
        self.assertLess(sim.vx, 0.0)
        self.assertAlmostEqual(sim.x, old_x)


if __name__ == '__main__':
    unittest.main()
