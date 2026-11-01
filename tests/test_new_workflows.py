"""新写回归检查：新版案例、压力定义、周期显示和全资源身份。"""
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from gear_contact import analyze_pair
from gear_transmission.case_workflow import prepare_case_request
from gear_transmission.plotting import contact_pass_data, plot_result
from gear_transmission.protocol import prepare_request, read_json
from gear_transmission.upstream import snapshot_root
from gear_transmission.workflow import identity, replay_run, run_request
from gear_ltca_snapshot.examples.case_contact_plotting import repeat_cycle_for_display
from gear_ltca_snapshot.examples.case_native_contact_plotting import project_contact_pass_samples

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'gear-skills.project.json').is_file())


class NewWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('GEAR_SKILLS_TEST_TMP'))
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def request(self, extended=False):
        request = read_json(ROOT/'skills/cylindrical-gear-analysis/assets/single_case.json')
        request['case'].update(phases=4, strips=5)
        if extended:
            request['case']['pair']['pinion']['helix_angle_deg'] = 15
            request['case']['pair']['wheel']['helix_angle_deg'] = -15
            request['case'].update(elastic_mode='iso_love', local_contact=False,
                iso_settings={'nominal_pinion_torque_n_m': 80}, phase_offset=.13,
                extended_contact_settings={'face_extension_factor': .14, 'roll_extension_factor': .4})
        return request

    def result(self, extended=False):
        prepared = prepare_request(self.request(extended))
        return analyze_pair(prepared['pair'], **prepared['kwargs'])

    def test_complete_resources_and_tamper_detection(self):
        self.assertGreater(len(identity()['runtime_resources']), 40)
        with patch('gear_transmission.upstream.snapshot_root', return_value=self.directory):
            with self.assertRaises(FileNotFoundError):
                identity()

    def test_full_case_executes_and_replays_strength_and_arrays(self):
        request = read_json(snapshot_root()/'examples/helical_18_77_case.json')
        request['analysis'].update(phases=4, strips=5)
        prepared = prepare_case_request(request)
        self.assertAlmostEqual(prepared['kwargs']['torque_n_m'],
            request['loading']['nominal_pinion_torque_n_m']*request['loading']['ltca_torque_multiplier'])
        run = run_request(request, self.directory/'case')
        self.assertTrue(run['analysis_completed'])
        self.assertFalse(run['validation']['same_inputs_as_reference'])
        self.assertFalse(run['validation']['strength_comparison_passed'])
        self.assertGreater(run['strength']['pinion']['root_stress_pa'], 0)
        replay = replay_run(self.directory/'case', self.directory/'replay')
        self.assertTrue(replay['passed'])
        self.assertTrue(replay['strength_equal'])

    def test_endpoint_pass_samples_are_original_values_and_flags(self):
        result = self.result(extended=True)
        data = contact_pass_data(result)
        np.testing.assert_array_equal(data['samples'], project_contact_pass_samples(result))
        self.assertEqual(data['definition'], '独立端点压力')
        self.assertTrue(np.isnan(data['pressure_pa_x_y']).any())
        self.assertTrue(np.any(data['samples'][:, 3:5]))
        self.assertFalse(np.any(data['ordinary_samples'][:, 3:5]))

    def test_new_generic_figures_and_unmodified_cycle(self):
        result = self.result()
        before = [r.pressure_pa.copy() for r in result.rotations]
        metrics = plot_result(result, self.directory)
        self.assertFalse(metrics['reference_used'])
        self.assertIn('中心', metrics['pressure_definition'])
        for name in ('01_case_geometry', '03_ltca_cycle', '04_contact_pattern',
                     '05_contact_coordinates', '07_project_contact_pass'):
            self.assertTrue((self.directory/(name+'.png')).is_file())
        for original, rotation in zip(before, result.rotations):
            np.testing.assert_array_equal(original, rotation.pressure_pa)
        cycle = np.loadtxt(self.directory/'03_ltca_cycle.csv', delimiter=',', skiprows=1, usecols=range(1, 8))
        displayed = np.loadtxt(self.directory/'03_ltca_cycle_display.csv', delimiter=',', skiprows=1, usecols=range(1, 8))
        np.testing.assert_array_equal(displayed, repeat_cycle_for_display(cycle, 3))

    def test_missing_strength_factors_fail_before_output(self):
        request = read_json(snapshot_root()/'examples/helical_18_77_case.json')
        del request['strength']['mesh_factors']['dynamic_factor']
        with self.assertRaises(TypeError):
            run_request(request, self.directory/'invalid')
        self.assertFalse((self.directory/'invalid').exists())


if __name__ == '__main__':
    unittest.main()
