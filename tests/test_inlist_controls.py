import shutil
import tempfile
import unittest
from pathlib import Path

from src.io import inlist as ih
from src.mesa import inlist_controls


class InlistControlTests(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        workspace_root = Path(__file__).resolve().parents[1]
        src_inlist = workspace_root / "inlist_project"
        self.test_inlist = Path(self.tmpdir.name) / "inlist_project"
        shutil.copy2(src_inlist, self.test_inlist)

    def tearDown(self):
        self.tmpdir.cleanup()

    def _line_value(self, variable_name):
        lines = ih.getTxtArrays(str(self.test_inlist))
        idx = ih.getTxtIndex(lines, variable_name)
        return lines[idx].split()[-1]

    def test_set_drag_energy_injection_updates_controls_uniform(self):
        inlist_controls.setDragEnergyInjection(
            str(self.test_inlist),
            True,
            energy_erg=2.5e43,
            r_inner_cm=2.0e12,
            r_outer_cm=6.0e12,
            duration_sec=1.25e4,
            kernel="Uniform",
            gaussian_sigma_fraction=0.35,
        )

        self.assertEqual(self._line_value("use_other_energy"), ".true.")
        self.assertEqual(float(self._line_value("x_ctrl(3)")), 2.5e43)
        self.assertEqual(float(self._line_value("x_ctrl(4)")), 2.0e12)
        self.assertEqual(float(self._line_value("x_ctrl(5)")), 6.0e12)
        self.assertEqual(float(self._line_value("x_ctrl(6)")), 1.25e4)
        self.assertEqual(float(self._line_value("x_ctrl(7)")), 0.35)
        self.assertEqual(self._line_value("x_integer_ctrl(1)"), "0")

    def test_set_drag_energy_injection_gaussian_reorders_radii(self):
        inlist_controls.setDragEnergyInjection(
            str(self.test_inlist),
            True,
            energy_erg=-3.0e40,
            r_inner_cm=9.0e11,
            r_outer_cm=3.0e11,
            duration_sec=500.0,
            kernel="Gauss",
            gaussian_sigma_fraction=0.2,
        )

        self.assertEqual(float(self._line_value("x_ctrl(3)")), 3.0e40)
        self.assertEqual(float(self._line_value("x_ctrl(4)")), 3.0e11)
        self.assertEqual(float(self._line_value("x_ctrl(5)")), 9.0e11)
        self.assertEqual(self._line_value("x_integer_ctrl(1)"), "1")

    def test_set_drag_energy_injection_disable_turns_hook_off(self):
        inlist_controls.setDragEnergyInjection(
            str(self.test_inlist),
            False,
            energy_erg=1.0,
            r_inner_cm=1.0,
            r_outer_cm=2.0,
            duration_sec=1.0,
        )

        self.assertEqual(self._line_value("use_other_energy"), ".false.")
        self.assertEqual(self._line_value("x_ctrl(3)"), "0d0")

    def test_set_drag_energy_injection_rejects_unknown_kernel(self):
        with self.assertRaises(ValueError):
            inlist_controls.setDragEnergyInjection(
                str(self.test_inlist),
                True,
                energy_erg=1.0,
                r_inner_cm=1.0,
                r_outer_cm=2.0,
                duration_sec=1.0,
                kernel="Triangle",
            )


if __name__ == "__main__":
    unittest.main()
