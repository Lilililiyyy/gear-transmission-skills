"""新写适配验收：验证单位、原数组、真实参数变化、失败及独立复算，不认证物理精度。"""
import csv
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from gear_contact import analyze_config, export_result
from gear_transmission.protocol import prepare_request, read_json
from gear_transmission.workflow import identity, replay_run, run_request

ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'gear-skills.project.json').is_file())
ASSETS = ROOT / 'skills/cylindrical-gear-analysis/assets'


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('GEAR_SKILLS_TEST_TMP'))
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def request(self, study=False):
        request = read_json(ASSETS / ('crown_study.json' if study else 'single_case.json'))
        request['case'].update(phases=4, strips=5)
        return request

    def test_explicit_units_and_conflicts(self):
        request = self.request()
        prepared = prepare_request(request)
        self.assertAlmostEqual(prepared['pair'].pinion.normal_module_m, 0.0025)
        self.assertAlmostEqual(prepared['kwargs']['modification_pinion'].lead_crown_m, 8e-6)
        self.assertAlmostEqual(prepared['kwargs']['material_pinion'].youngs_modulus_pa, 210e9)
        request['case']['pair']['pinion']['normal_module_m'] = 0.0025
        with self.assertRaises(ValueError):
            prepare_request(request)

    def test_reject_invalid_input_before_output(self):
        for value in (-1, True, float('inf')):
            request = self.request()
            request['case']['torque_n_m'] = value
            output = self.directory / 'invalid'
            with self.assertRaises(ValueError):
                run_request(request, output)
            self.assertFalse(output.exists())
        request = self.request()
        request['case']['material_pinion']['unknown_factor'] = 1
        with self.assertRaises(TypeError):
            prepare_request(request)

    def test_reject_duplicate_json_keys(self):
        path = self.directory / 'duplicate.json'
        path.write_text('{"schema":"first","schema":"second"}', encoding='utf-8')
        with self.assertRaises(ValueError):
            read_json(path)

    def test_single_matches_direct_kernel_and_replays(self):
        request = self.request()
        run = run_request(request, self.directory / 'run')
        direct = analyze_config(prepare_request(request)['normalized']['case'])
        export_result(direct, self.directory / 'direct')
        with np.load(self.directory / 'run/cases/case_001/arrays.npz', allow_pickle=False) as a, \
             np.load(self.directory / 'direct/arrays.npz', allow_pickle=False) as b:
            self.assertEqual(set(a.files), set(b.files))
            for key in a.files:
                np.testing.assert_array_equal(a[key], b[key], err_msg=key)
        self.assertEqual(run['status'], 'completed')
        self.assertFalse(run['validation']['physical_validation'])
        replay = replay_run(self.directory / 'run', self.directory / 'replay')
        self.assertTrue(replay['passed'])
        self.assertGreater(len(replay['cases'][0]['arrays']), 15)

    def test_comparison_uses_each_candidate_and_operating_point(self):
        request = self.request(study=True)
        request['study']['designs'] = request['study']['designs'][:2]
        request['study']['operating_points'] = request['study']['operating_points'][:2]
        run = run_request(request, self.directory / 'run')
        self.assertEqual(len(run['cases']), 4)
        self.assertEqual([r['summary']['torque_n_m'] for r in run['cases']], [20, 80, 20, 80])
        a, b = run['cases'][1]['summary'], run['cases'][3]['summary']
        self.assertNotEqual(a['max_pressure_pa'], b['max_pressure_pa'])
        replay = replay_run(self.directory / 'run', self.directory / 'replay')
        self.assertTrue(replay['passed'])
        self.assertTrue(replay['study_equal'])

    def test_no_silent_study_override(self):
        request = self.request(study=True)
        request['case']['torque_n_m'] = 80
        with self.assertRaises(ValueError):
            prepare_request(request)
        request = self.request(study=True)
        request['study']['designs'][1]['name'] = request['study']['designs'][0]['name']
        with self.assertRaises(ValueError):
            prepare_request(request)

    def test_iso_rejects_duplicate_local_deformation(self):
        request = self.request()
        request['case']['elastic_mode'] = 'iso'
        with self.assertRaises(ValueError):
            prepare_request(request)

    def test_iso_explicit_settings_run_and_replay(self):
        request = self.request()
        request['case'].update(elastic_mode='iso', local_contact=False,
                               iso_settings={'nominal_pinion_torque_n_m': 80})
        run = run_request(request, self.directory / 'iso')
        self.assertEqual(run['cases'][0]['summary']['elastic_mode'], 'iso')
        self.assertTrue(replay_run(self.directory / 'iso', self.directory / 'iso_replay')['passed'])

    def test_zero_load_internal_iso_keeps_unloaded_curvature_out_of_pressure(self):
        """零载荷内啮合 ISO 请求不应因未承载点的奇异曲率而失败。"""
        request = self.request()
        for gear in request['case']['pair'].values():
            gear.update(normal_module_mm=2, face_width_mm=20)
        request['case']['pair']['pinion']['teeth'] = 20
        request['case']['pair']['wheel'].update(teeth=60, internal=True)
        request['case'].update(torque_n_m=0, elastic_mode='iso', local_contact=False,
                               iso_settings={'nominal_pinion_torque_n_m': 0})
        run = run_request(request, self.directory / 'zero_iso')
        self.assertEqual(run['cases'][0]['summary']['max_pressure_pa'], 0)
        with np.load(self.directory / 'zero_iso/cases/case_001/arrays.npz', allow_pickle=False) as arrays:
            np.testing.assert_array_equal(arrays['forces_n'], 0)
        replay = replay_run(self.directory / 'zero_iso', self.directory / 'zero_iso_replay')
        self.assertTrue(replay['passed'])
        self.assertTrue(all(row['maximum_absolute_difference'] == 0 for row in replay['cases'][0]['arrays']))

    def test_extended_contact_csv_preserves_flags_and_endpoint_layout(self):
        """延伸接触的 CSV 与原数组保持逐行一致，并能独立复算。"""
        request = self.request()
        request['case']['pair']['pinion']['helix_angle_deg'] = 15
        request['case']['pair']['wheel']['helix_angle_deg'] = -15
        request['case'].update(elastic_mode='iso_love', local_contact=False, phase_offset=0.13,
                               iso_settings={'nominal_pinion_torque_n_m': 80},
                               extended_contact_settings={'face_extension_factor': 0.14, 'roll_extension_factor': 0.4})
        run_request(request, self.directory / 'extended')
        case = self.directory / 'extended/cases/case_001'
        with (case / 'contact_cells.csv').open(encoding='utf-8-sig', newline='') as stream:
            rows = list(csv.DictReader(stream))
        with np.load(case / 'arrays.npz', allow_pickle=False) as arrays:
            self.assertEqual(len(rows), len(arrays['forces_n']))
            for name in ('eligible', 'extended_pinion', 'extended_wheel'):
                np.testing.assert_array_equal([row[name] == 'True' for row in rows], arrays[name])
            np.testing.assert_array_equal([float(row['additional_gap_m']) for row in rows], arrays['additional_gap_m'])
            np.testing.assert_array_equal(np.diff(arrays['line_endpoint_offsets']), np.diff(arrays['line_cell_offsets']) + 1)
            self.assertTrue(np.any(arrays['extended_pinion'] | arrays['extended_wheel']))
            self.assertTrue(np.all(arrays['forces_n'][~arrays['eligible']] == 0))
        self.assertTrue((case / 'endpoint_cells.csv').is_file())
        self.assertTrue(replay_run(self.directory / 'extended', self.directory / 'extended_replay')['passed'])

    def test_wrong_material_and_boolean_geometry_rejected(self):
        request = self.request()
        request['case']['material_pinion'] = 'steel'
        with self.assertRaises(ValueError):
            prepare_request(request)
        request = self.request()
        request['case']['pair']['wheel']['internal'] = 0
        with self.assertRaises(ValueError):
            prepare_request(request)

    def test_internal_helical_and_nondefault_material(self):
        request = self.request()
        request['case']['pair']['pinion']['helix_angle_deg'] = 20
        request['case']['pair']['wheel'].update(teeth=84, internal=True, helix_angle_deg=20)
        request['case']['material_pinion']['youngs_modulus_gpa'] = 200
        run = run_request(request, self.directory / 'internal')
        self.assertTrue(run['cases'][0]['summary']['all_converged'])
        request['case']['pair']['wheel']['helix_angle_deg'] = -20
        with self.assertRaises(ValueError):
            prepare_request(request)

    def test_zero_torque_and_output_preservation(self):
        request = self.request()
        request['case']['torque_n_m'] = 0
        run = run_request(request, self.directory / 'zero')
        self.assertEqual(run['cases'][0]['summary']['max_pressure_pa'], 0)
        with self.assertRaises(FileExistsError):
            run_request(request, self.directory / 'zero')
        self.assertEqual(read_json(self.directory / 'zero/run.json')['status'], 'completed')

    def test_tampered_result_is_not_replayed(self):
        run_request(self.request(), self.directory / 'run')
        target = self.directory / 'run/cases/case_001/summary.json'
        target.write_text('{}', encoding='utf-8')
        with self.assertRaises(ValueError):
            replay_run(self.directory / 'run', self.directory / 'replay')
        self.assertFalse((self.directory / 'replay').exists())

    def test_failure_keeps_diagnostic_status(self):
        with patch('gear_transmission.workflow.analyze_pair', side_effect=ArithmeticError('测试求解失败')):
            with self.assertRaises(ArithmeticError):
                run_request(self.request(), self.directory / 'failed')
        result = read_json(self.directory / 'failed/run.json')
        self.assertEqual(result['status'], 'failed')
        self.assertFalse(result['analysis_completed'])
        self.assertFalse((self.directory / 'failed/report.md').exists())

    def test_installed_kernel_identity(self):
        self.assertGreater(len(identity()['kernel_files']), 10)


if __name__ == '__main__':
    unittest.main()
