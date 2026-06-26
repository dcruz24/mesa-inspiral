import os
import runpy
import unittest
from pathlib import Path


class FitCheckScriptTests(unittest.TestCase):
    def run_script_and_check_outputs(self, script_name, expected_outputs):
        try:
            import matplotlib  # noqa: F401
            import numpy  # noqa: F401
        except ImportError as exc:
            self.skipTest(f'Missing plotting dependency: {exc}')

        repo_root = Path(__file__).resolve().parents[1]
        script = repo_root / 'tests' / 'fit_checks' / script_name
        plots = repo_root / 'tests' / 'fit_checks' / 'plots'
        env_backend = os.environ.get('MPLBACKEND')
        os.environ['MPLBACKEND'] = 'Agg'
        try:
            runpy.run_path(str(script), run_name='__main__')
        finally:
            if env_backend is None:
                os.environ.pop('MPLBACKEND', None)
            else:
                os.environ['MPLBACKEND'] = env_backend

        for output in expected_outputs:
            path = plots / output
            self.assertTrue(path.exists(), f'Missing generated plot: {path}')
            self.assertGreater(path.stat().st_size, 0, f'Generated plot is empty: {path}')

    def test_macleod_fit_curve_plots(self):
        self.run_script_and_check_outputs(
            'macleod_fit_curves.py',
            ['WindTunnelFit.png', 'WindTunnel2.png'],
        )

    def test_macleod_fit_surface_plots(self):
        self.run_script_and_check_outputs(
            'macleod_fit_surfaces.py',
            ['A1_logCa43.png', 'A2_logCd43.png', 'A3_logCa53.png', 'A4_logCd53.png'],
        )

    def test_neutron_star_eos_plot(self):
        self.run_script_and_check_outputs(
            'neutron_star_eos_curves.py',
            ['SLY4_eos_curves.png'],
        )


if __name__ == '__main__':
    unittest.main()
