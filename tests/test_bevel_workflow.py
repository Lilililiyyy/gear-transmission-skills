"""新写行为验收：真实直／弧齿强度、来源绑定、错误输入及独立复算。"""
import json
from pathlib import Path
import tempfile
import unittest

from bevel_gear.workflow import load_bevel_case
from gear_transmission.bevel_workflow import identity, prepare_case, read_json, run_case, replay_run


ROOT = next(p for p in Path(__file__).resolve().parents if (p/'gear-skills.project.json').is_file())
ASSETS = ROOT/'skills/bevel-gear-design/assets'


class BevelWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.output = Path(self.temp.name)

    def test_straight_actual_strength_and_replay(self):
        case = load_bevel_case(ASSETS/'straight_16_32_case.json')
        result = run_case(case, self.output/'run')
        self.assertEqual(result['status'], 'completed')
        data = read_json(self.output/'run/calculation/workflow.json')
        self.assertEqual(data['geometry']['inputs']['wheel_teeth'], 32)
        self.assertIsNotNone(data['strength_assessment'])
        self.assertTrue(replay_run(self.output/'run', self.output/'replay')['passed'])

    def test_spiral_actual_strength_matches_bound_reference(self):
        case = load_bevel_case(ASSETS/'spiral_17_47_case.json')
        reference = read_json(ASSETS/'spiral_17_47_reference.json')
        run_case(case, self.output/'run', reference=reference)
        result = read_json(self.output/'run/calculation/workflow.json')
        self.assertTrue(result['masta_benchmark']['same_case_inputs'])
        self.assertTrue(result['masta_benchmark']['comparison_passed'])
        self.assertTrue(replay_run(self.output/'run', self.output/'replay')['passed'])

    def test_hypoid_and_inconsistent_manufacturing_rejected(self):
        case = load_bevel_case(ASSETS/'spiral_17_47_case.json')
        case['geometry']['offset_m'] = .001
        with self.assertRaisesRegex(ValueError, '零偏置'):
            prepare_case(case, 'geometry')
        case = load_bevel_case(ASSETS/'spiral_17_47_case.json')
        case['geometry']['wheel_outer_pitch_diameter_m'] *= 1.1
        with self.assertRaises(ValueError):
            prepare_case(case, 'strength')

    def test_existing_output_and_modified_source_data_rejected(self):
        case = load_bevel_case(ASSETS/'straight_16_32_case.json')
        run_case(case, self.output/'run')
        with self.assertRaises(FileExistsError):
            run_case(case, self.output/'run')
        path = self.output/'run/calculation/workflow.json'
        value = json.loads(path.read_text(encoding='utf-8'))
        value['forces'] = {}
        path.write_text(json.dumps(value), encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'SHA256'):
            replay_run(self.output/'run', self.output/'replay')

    def test_identity_is_pinned(self):
        self.assertEqual(identity()['source_commit'], '05da3ac166a750a6326907af5bd0dba8aa6ab034')


if __name__ == '__main__':
    unittest.main()
