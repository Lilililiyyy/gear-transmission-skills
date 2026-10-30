"""新写安装验收：验证完整资源、更新边界、只读检查及真实 Git 快进。"""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'gear-skills.project.json').is_file())
SPEC = importlib.util.spec_from_file_location('gear_skill_installer', ROOT / 'scripts/update_codex_skills.py')
INSTALLER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(INSTALLER)


class InstallerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('GEAR_SKILLS_TEST_TMP'))
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name).resolve()
        self.root = self.directory / 'project'
        self.source = self.root / 'skills/test-gear'
        (self.source / 'assets').mkdir(parents=True)
        (self.source / 'SKILL.md').write_text('---\nname: test-gear\ndescription: 齿轮计算\n---\n', encoding='utf-8')
        (self.source / 'assets/sample.bin').write_bytes(bytes(range(256)))
        (self.root / 'pyproject.toml').write_text('[project]\nname="fixture"\n', encoding='utf-8')
        (self.root / 'gear-skills.project.json').write_text(
            json.dumps({'skills': ['skills/test-gear']}), encoding='utf-8')
        self.destination = self.directory / 'installed'
        self.runtime = self.directory / 'runtime'

    def invoke(self, *args):
        with patch.object(INSTALLER, 'project_root', return_value=self.root), \
             patch.object(INSTALLER, 'install_runtime', return_value=Path('fixture-python')) as install, \
             patch.object(INSTALLER, 'check_runtime', return_value={'python': 'fixture-python'}):
            result = INSTALLER.main(['--dest', str(self.destination), '--runtime', str(self.runtime), *args])
        return result, install.call_count

    def test_complete_install_and_read_only_check(self):
        self.assertEqual(self.invoke()[0], 0)
        target = self.destination / 'test-gear'
        self.assertEqual((target / 'assets/sample.bin').read_bytes(), bytes(range(256)))
        before = {str(p): (p.read_bytes(), p.stat().st_mtime_ns) for p in self.destination.rglob('*') if p.is_file()}
        self.assertEqual(self.invoke('--check'), (0, 0))
        after = {str(p): (p.read_bytes(), p.stat().st_mtime_ns) for p in self.destination.rglob('*') if p.is_file()}
        self.assertEqual(before, after)

    def test_update_removes_obsolete_owned_resource(self):
        self.invoke()
        (self.source / 'assets/sample.bin').unlink()
        (self.source / 'assets/new.json').write_text('{}', encoding='utf-8')
        with self.assertRaises(ValueError):
            self.invoke('--check')
        self.invoke()
        self.assertFalse((self.destination / 'test-gear/assets/sample.bin').exists())
        self.assertTrue((self.destination / 'test-gear/assets/new.json').is_file())
        self.invoke('--check')

    def test_local_edits_are_preserved(self):
        self.invoke()
        target = self.destination / 'test-gear/SKILL.md'
        target.write_text('本地修改', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, '本地修改'):
            self.invoke()
        self.assertEqual(target.read_text(encoding='utf-8'), '本地修改')

    def test_unrelated_skill_is_untouched_and_collision_is_rejected(self):
        unrelated = self.destination / 'other-skill/SKILL.md'
        unrelated.parent.mkdir(parents=True)
        unrelated.write_text('其他技能', encoding='utf-8')
        self.invoke()
        self.assertEqual(unrelated.read_text(encoding='utf-8'), '其他技能')
        other_dest = self.directory / 'collision'
        (other_dest / 'test-gear').mkdir(parents=True)
        (other_dest / 'test-gear/SKILL.md').write_text('第三方技能', encoding='utf-8')
        with self.assertRaises(ValueError):
            INSTALLER.preflight(INSTALLER.sources(self.root), other_dest, {})

    def test_runtime_failure_does_not_register_skill(self):
        with patch.object(INSTALLER, 'project_root', return_value=self.root), \
             patch.object(INSTALLER, 'install_runtime', side_effect=ValueError('依赖安装失败')):
            with self.assertRaises(ValueError):
                INSTALLER.main(['--dest', str(self.destination)])
        self.assertFalse(self.destination.exists())

    def test_missing_install_and_extra_resource_fail_check(self):
        with self.assertRaises(ValueError):
            self.invoke('--check')
        self.invoke()
        (self.destination / 'test-gear/extra.txt').write_text('额外文件', encoding='utf-8')
        with self.assertRaises(ValueError):
            self.invoke('--check')

    def test_source_directory_cannot_be_install_destination(self):
        with patch.object(INSTALLER, 'project_root', return_value=self.root):
            with self.assertRaises(ValueError):
                INSTALLER.main(['--dest', str(self.root / 'skills')])

    def test_source_registration_cannot_escape_project(self):
        (self.root / 'gear-skills.project.json').write_text(
            json.dumps({'skills': ['../test-gear']}), encoding='utf-8')
        with self.assertRaises(ValueError):
            INSTALLER.sources(self.root)

    def test_pull_fast_forwards_and_rejects_local_changes(self):
        def git(*args):
            return subprocess.run(['git', *map(str, args)], check=True, capture_output=True, text=True)
        git('init', self.root)
        git('-C', self.root, 'config', 'user.name', 'Installer Test')
        git('-C', self.root, 'config', 'user.email', 'installer@example.invalid')
        git('-C', self.root, 'config', 'core.autocrlf', 'false')
        git('-C', self.root, 'add', '.')
        git('-C', self.root, 'commit', '-m', 'fixture')
        remote = self.directory / 'remote.git'
        clone = self.directory / 'clone'
        git('clone', '--bare', self.root, remote)
        git('clone', remote, clone)
        git('-C', self.root, 'remote', 'add', 'origin', remote)
        (self.source / 'assets/new.json').write_text('{}', encoding='utf-8')
        git('-C', self.root, 'add', '.')
        git('-C', self.root, 'commit', '-m', 'update fixture')
        git('-C', self.root, 'push', 'origin', 'HEAD')
        INSTALLER.pull(clone)
        self.assertTrue((clone / 'skills/test-gear/assets/new.json').is_file())
        (clone / 'skills/test-gear/SKILL.md').write_text('本地修改', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, '未提交修改'):
            INSTALLER.pull(clone)


if __name__ == '__main__':
    unittest.main()
