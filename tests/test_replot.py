"""绘图故障回归：标准比较入口、只读重绘及损坏结果拒绝。"""
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from gear_contact.contact_display import contact_pass_data as source_contact_pass_data
from gear_transmission.plotting import contact_pass_data
from gear_transmission.protocol import read_json
from gear_transmission.replot import replot_run
from gear_transmission.workflow import artifact_hashes, run_request

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'gear-skills.project.json').is_file())


class ReplotTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('GEAR_SKILLS_TEST_TMP'))
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def request(self):
        value = read_json(ROOT/'skills/cylindrical-gear-analysis/assets/crown_study.json')
        value['case'].update(phases=4, strips=5)
        value['study']['designs'] = value['study']['designs'][:2]
        value['study']['operating_points'] = value['study']['operating_points'][:1]
        return value

    def test_standard_study_outputs_shared_style_comparison(self):
        self.assertIs(contact_pass_data, source_contact_pass_data)
        out = self.directory/'run'
        run_request(self.request(), out, plots=True)
        figure = out/'plots/point_001/07_contact_pass_comparison.png'
        self.assertTrue(figure.is_file())
        record = read_json(figure.parent/'contact_pass_comparison.json')
        self.assertEqual(len(record['records']), 2)
        self.assertFalse(record['raw_results_modified'])

    def test_replot_never_solves_and_preserves_all_original_artifacts(self):
        out = self.directory/'run'
        run_request(self.request(), out)
        before = artifact_hashes(out)
        with patch('gear_transmission.workflow.analyze_pair', side_effect=AssertionError('重绘不应求解')):
            record = replot_run(out, self.directory/'plots')
        self.assertEqual(record['status'], 'completed')
        self.assertFalse(record['solver_executed'])
        self.assertFalse(record['numerical_replay_performed'])
        self.assertEqual(before, artifact_hashes(out))
        self.assertTrue((self.directory/'plots/plots/point_001/07_contact_pass_comparison.png').is_file())

    def test_modified_source_is_rejected_before_creating_output(self):
        out = self.directory/'run'
        run_request(self.request(), out)
        (out/'cases/case_001/summary.json').write_text('{}', encoding='utf-8')
        with self.assertRaises(ValueError):
            replot_run(out, self.directory/'plots')
        self.assertFalse((self.directory/'plots').exists())


if __name__ == '__main__':
    unittest.main()
