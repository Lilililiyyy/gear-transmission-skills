"""压力指标从请求贯穿报告、完整复算及只读重绘的接入验收。"""
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from gear_transmission.protocol import prepare_request, read_json
from gear_transmission.replot import replot_run
from gear_transmission.workflow import artifact_hashes, replay_run, run_request


ROOT = next(p for p in Path(__file__).resolve().parents if (p/'gear-skills.project.json').is_file())


class DesignMetricWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('GEAR_SKILLS_TEST_TMP'))
        self.addCleanup(self.temp.cleanup)
        self.output = Path(self.temp.name)

    def request(self, extended=False):
        request = read_json(ROOT/'skills/cylindrical-gear-analysis/assets/crown_study.json')
        request['case'].update(phases=4, strips=5)
        request['study']['designs'] = request['study']['designs'][:2]
        request['study']['operating_points'] = request['study']['operating_points'][:1]
        if extended:
            request['case']['pair']['pinion']['helix_angle_deg'] = 15
            request['case']['pair']['wheel']['helix_angle_deg'] = -15
            request['case'].update(elastic_mode='iso_love', local_contact=False, phase_offset=.13,
                iso_settings={'nominal_pinion_torque_n_m': 80},
                extended_contact_settings={'face_extension_factor': .14, 'roll_extension_factor': .4})
            request['study']['pressure_metric'] = 'endpoint'
        return request

    def test_default_and_invalid_metric_are_explicit_before_output(self):
        request = self.request()
        prepared = prepare_request(request)
        self.assertEqual(prepared['pressure_metric'], 'center')
        self.assertIn('study.pressure_metric', prepared['defaults'])
        for metric in ('endpoint', 'unknown', [], True):
            request['study']['pressure_metric'] = metric
            with self.assertRaises(ValueError):
                run_request(request, self.output/'invalid')
            self.assertFalse((self.output/'invalid').exists())

    def test_endpoint_report_and_replay_use_original_endpoint_array(self):
        run = run_request(self.request(extended=True), self.output/'run')
        self.assertEqual(run['pressure_metric'], 'endpoint')
        study = read_json(self.output/'run/study.json')
        self.assertEqual(study['pressure_metric'], 'endpoint')
        for case, row in zip(run['cases'], study['rows']):
            with np.load(self.output/'run'/case['directory']/'arrays.npz', allow_pickle=False) as arrays:
                expected = float(np.max(arrays['endpoint_pressure_pa']))
                self.assertAlmostEqual(case['selected_pressure_pa']/expected, 1.)
                self.assertAlmostEqual(row['max_pressure_mpa']*1e6/expected, 1.)
                self.assertAlmostEqual(case['summary']['max_pressure_pa']/float(np.max(arrays['pressure_pa'])), 1.)
        report = (self.output/'run/report.md').read_text(encoding='utf-8')
        self.assertIn('本次压力指标：独立端点压力', report)
        replay = replay_run(self.output/'run', self.output/'replay')
        self.assertTrue(replay['passed'])
        self.assertTrue(replay['pressure_selection_equal'])
        self.assertTrue(replay['study_equal'])

    def test_endpoint_replot_preserves_selected_metric_without_solving(self):
        run_request(self.request(extended=True), self.output/'run')
        before = artifact_hashes(self.output/'run')
        with patch('gear_transmission.workflow.analyze_pair', side_effect=AssertionError('重绘不能求解')):
            result = replot_run(self.output/'run', self.output/'replot')
        self.assertFalse(result['solver_executed'])
        self.assertEqual(result['pressure_metric'], 'endpoint')
        self.assertEqual(artifact_hashes(self.output/'run'), before)
        self.assertTrue((self.output/'replot/plots/study/robust_design_plot_data.csv').is_file())


if __name__ == '__main__':
    unittest.main()
