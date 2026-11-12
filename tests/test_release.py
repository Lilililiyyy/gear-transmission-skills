"""版本策略回归：单一来源、双语内核表、Skill无独立编号及历史保留。"""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'gear-skills.project.json').is_file())
SPEC = importlib.util.spec_from_file_location('library_release', ROOT/'scripts/release.py')
RELEASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RELEASE)


class ReleaseTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('GEAR_SKILLS_TEST_TMP'))
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        files = ['pyproject.toml', 'gear-skills.project.json', '.codex-plugin/plugin.json', 'README.md', 'README_EN.md',
                 'src/gear_transmission/kernel_snapshot.json', 'src/gear_transmission/bevel_snapshot.json',
                 'skills/cylindrical-gear-analysis/SKILL.md', 'skills/bevel-gear-design/SKILL.md']
        for name in files:
            (self.root/name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT/name, self.root/name)

    def test_only_authoritative_version_changes_all_generated_metadata(self):
        path = self.root/'pyproject.toml'
        previous = RELEASE.project_version(self.root)
        path.write_text(path.read_text(encoding='utf-8').replace(f'version = "{previous}"', 'version = "9.8.7"', 1), encoding='utf-8')
        before = (self.root/'src/gear_transmission/kernel_snapshot.json').read_bytes()
        with self.assertRaises(ValueError):
            RELEASE.synchronize(self.root)
        RELEASE.synchronize(self.root, write=True)
        RELEASE.synchronize(self.root)
        for name in ('gear-skills.project.json', '.codex-plugin/plugin.json'):
            self.assertEqual(json.loads((self.root/name).read_text(encoding='utf-8'))['version'], '9.8.7')
        for name in ('README.md', 'README_EN.md'):
            text = (self.root/name).read_text(encoding='utf-8')
            self.assertNotIn('library-version:', text)
            self.assertIn('img.shields.io/github/v/release/', text)
            self.assertIn('v0.4.1', text)
        self.assertEqual(before, (self.root/'src/gear_transmission/kernel_snapshot.json').read_bytes())

    def test_new_kernel_pin_updates_both_tables_without_changing_library_version(self):
        path = self.root/'src/gear_transmission/bevel_snapshot.json'
        manifest = json.loads(path.read_text(encoding='utf-8'))
        manifest.update(source_version='99.1.2', source_commit='1'*40)
        path.write_text(json.dumps(manifest), encoding='utf-8')
        before = RELEASE.project_version(self.root)
        with self.assertRaises(ValueError):
            RELEASE.synchronize(self.root)
        RELEASE.synchronize(self.root, write=True)
        self.assertEqual(RELEASE.project_version(self.root), before)
        for name in ('README.md', 'README_EN.md'):
            text = (self.root/name).read_text(encoding='utf-8')
            self.assertIn('99.1.2', text)
            self.assertNotIn('/commit/'+'1'*40, text)
        self.assertIn('/commit/'+'1'*40, RELEASE.release_notes(self.root, before))

    def test_current_readme_layout_is_preserved_and_sync_is_idempotent(self):
        """徽章与四列表格可直接校验；同步不恢复已删除的版本段或状态列。"""
        before = {name: (self.root/name).read_bytes() for name in ('README.md', 'README_EN.md')}
        RELEASE.synchronize(self.root)
        self.assertEqual(RELEASE.synchronize(self.root, write=True)['updated'], [])
        for name, original in before.items():
            self.assertEqual((self.root/name).read_bytes(), original)
            rows = [line for line in original.decode('utf-8').splitlines()
                    if line.startswith('|') and '/SKILL.md' in line]
            self.assertEqual(len(rows), len(RELEASE.KERNELS))
            self.assertTrue(all(len(line.strip().strip('|').split('|')) == 4 for line in rows))

    def test_stale_kernel_version_is_rejected_without_writing(self):
        """只读检查发现过期内核版本后不改文件；显式同步仅更新版本单元格。"""
        path = self.root/'README.md'
        records = RELEASE.kernel_records(self.root)
        current = records['bevel-gear-design']['version']
        path.write_text(path.read_text(encoding='utf-8').replace('| '+current+' |', '| 99.9.9 |'), encoding='utf-8')
        before = {p: p.read_bytes() for p in self.root.rglob('*') if p.is_file()}
        with self.assertRaisesRegex(ValueError, '发布元数据不同步'):
            RELEASE.synchronize(self.root)
        self.assertEqual(before, {p: p.read_bytes() for p in self.root.rglob('*') if p.is_file()})
        self.assertEqual(RELEASE.synchronize(self.root, write=True)['updated'], ['README.md'])
        RELEASE.synchronize(self.root)

    def test_unexpected_table_column_is_rejected(self):
        """表格仍严格检查结构，避免在错误列中写入内核版本。"""
        path = self.root/'README_EN.md'
        lines = path.read_text(encoding='utf-8').splitlines()
        lines = [line+' unexpected |' if line.startswith('|') and '/SKILL.md' in line else line for line in lines]
        path.write_text('\n'.join(lines)+'\n', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, '四列'):
            RELEASE.synchronize(self.root)

    def test_independent_skill_version_is_rejected(self):
        path = self.root/'skills/bevel-gear-design/SKILL.md'
        text = path.read_text(encoding='utf-8').replace('license: AGPL-3.0-only', 'license: AGPL-3.0-only\nmetadata:\n  version: "99.1.2"')
        path.write_text(text, encoding='utf-8')
        with self.assertRaises(ValueError):
            RELEASE.synchronize(self.root)

    def test_actual_installation_must_match_authoritative_version(self):
        with patch.object(RELEASE, 'installed_version', return_value='0.0.0'):
            with self.assertRaises(ValueError):
                RELEASE.synchronize(self.root, check_installed=True)

    def test_release_notes_include_both_languages_and_actual_pins_only(self):
        notes = RELEASE.release_notes(self.root, RELEASE.project_version(self.root))
        self.assertIn('## 中文', notes)
        self.assertIn('## English', notes)
        self.assertNotIn('v0.4.1', notes)
        for record in RELEASE.kernel_records(self.root).values():
            self.assertIn(record['version'], notes)
            self.assertIn(record['commit'], notes)


if __name__ == '__main__':
    unittest.main()
