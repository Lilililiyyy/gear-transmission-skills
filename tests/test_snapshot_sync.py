"""来源提交字节验收：检出换行转换不应进入快照，历史CRLF证据不得规范化。"""
import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = next(p for p in Path(__file__).resolve().parents if (p/'gear-skills.project.json').is_file())
SPEC = importlib.util.spec_from_file_location('sync_ltca_snapshot', ROOT/'scripts/sync_ltca_snapshot.py')
SYNC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SYNC)


class SnapshotSyncTests(unittest.TestCase):
    def test_committed_bytes_ignore_checkout_conversion_and_preserve_raw_evidence(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('GEAR_SKILLS_TEST_TMP')) as directory:
            repo = Path(directory)
            def git(*args):
                return subprocess.check_output(['git', '-C', str(repo), *args], stderr=subprocess.DEVNULL)
            git('init')
            git('config', 'core.autocrlf', 'false')
            (repo/'.gitattributes').write_bytes(b'*.py text eol=lf\n*.json -text\n')
            (repo/'kernel.py').write_bytes(b'value = 1\n')
            original_evidence = b'{"index": [0, 1]}\r\n'
            (repo/'evidence.json').write_bytes(original_evidence)
            git('add', '.')
            git('-c', 'user.name=Snapshot Test', '-c', 'user.email=snapshot@example.invalid',
                'commit', '-m', '固定字节夹具')
            commit = git('rev-parse', 'HEAD').decode().strip()
            (repo/'kernel.py').write_bytes(b'value = 1\r\n')
            # 旧版Git需刷新索引才能确认仅换行转换；入库内容仍为原LF字节。
            git('add', 'kernel.py')
            self.assertEqual(git('status', '--porcelain'), b'')
            result = SYNC.committed_files(repo, commit, ['kernel.py', 'evidence.json'])
            self.assertEqual(result['kernel.py'], b'value = 1\n')
            self.assertEqual(result['evidence.json'], original_evidence)


if __name__ == '__main__':
    unittest.main()
