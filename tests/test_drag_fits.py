import unittest

import numpy as np

from src.config import constants as const
from src.physics import basic, drag


class DragFitTests(unittest.TestCase):
    def test_cad_fits_are_positive_on_reference_grid(self):
        q = np.linspace(0.05, 0.5, 8)
        mach = np.linspace(1.1, 5.0, 8)
        Q, M = np.meshgrid(q, mach)
        for ca, cd in [drag.Cad43(Q, M), drag.Cad53(Q, M)]:
            self.assertTrue(np.all(np.isfinite(ca)))
            self.assertTrue(np.all(np.isfinite(cd)))
            self.assertTrue(np.all(ca > 0.0))
            self.assertTrue(np.all(cd > 0.0))

    def test_macleod_2015_fit_mode_returns_dimensionless_curves(self):
        eps = np.linspace(0.3, 3.0, 20)
        mdot, fd, ldot = drag.getNS2015Drag(1.4 * const.Msun_cgs, 1e-7, 1e7, 1.0, eps, fit=True)
        self.assertEqual(mdot.shape, eps.shape)
        self.assertEqual(fd.shape, eps.shape)
        self.assertEqual(ldot.shape, eps.shape)
        self.assertTrue(np.all(mdot > 0.0))
        self.assertTrue(np.all(fd > 0.0))

    def test_hoyle_lyttleton_drag_matches_mdot_times_velocity(self):
        mass = 1.4 * const.Msun_cgs
        rho = 1e-7
        velocity = 1e7
        self.assertAlmostEqual(
            basic.getFdHL(mass, rho, velocity),
            basic.getMdotHL(mass, rho, velocity) * velocity,
        )


if __name__ == '__main__':
    unittest.main()
