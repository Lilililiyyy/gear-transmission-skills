"""新写安装工具：配置独立环境，同步完整技能，并核对文件及内核身份。"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


RECEIPT = '.gear-transmission-skills-install.json'
IGNORED = {'__pycache__', '.pytest_cache', '.DS_Store'}


def project_root(start):
    """通过项目标记定位根目录，不依赖当前工作目录或固定父目录层数。"""
    for candidate in (start, *start.parents):
        if (candidate / 'gear-skills.project.json').is_file() and (candidate / 'pyproject.toml').is_file():
            return candidate.resolve()
    raise ValueError('找不到完整项目，请保留 gear-skills.project.json 和 pyproject.toml')


def file_hashes(directory):
    """比较完整资源，排除 Python 运行缓存；链接不作为可复制的资源接受。"""
    result = {}
    if not directory.is_dir():
        raise ValueError(f'技能目录不存在：{directory}')
    for item in sorted(directory.rglob('*')):
        relative = item.relative_to(directory)
        if any(part in IGNORED for part in relative.parts) or item.suffix in {'.pyc', '.pyo'}:
            continue
        if item.is_symlink() or item.resolve() != directory.resolve() / relative:
            raise ValueError(f'技能内部包含链接，请先改为完整资源文件：{item}')
        if item.is_file():
            result[relative.as_posix()] = hashlib.sha256(item.read_bytes()).hexdigest()
    return result


def sources(root):
    marker = json.loads((root / 'gear-skills.project.json').read_text(encoding='utf-8'))
    result = {}
    for relative in marker['skills']:
        path = root / relative
        name = path.name
        if not re.fullmatch(r'[a-z0-9][a-z0-9-]{0,63}', name) or path.resolve().parent != root / 'skills':
            raise ValueError(f'非法技能登记路径：{relative}')
        if name in result or not (path / 'SKILL.md').is_file():
            raise ValueError(f'技能重复或缺少 SKILL.md：{relative}')
        result[name] = (path, file_hashes(path))
    if not result:
        raise ValueError('项目未登记技能')
    return result


def command(args, root):
    print('执行：' + subprocess.list2cmdline([str(arg) for arg in args]), flush=True)
    subprocess.run([str(arg) for arg in args], cwd=root, check=True,
                   env={**os.environ, 'PYTHONUTF8': '1', 'PIP_DISABLE_PIP_VERSION_CHECK': '1'})


def pull(root):
    """只允许干净 Git 检出的快进更新，不覆盖本地修改。"""
    status = subprocess.run(['git', '-C', str(root), 'status', '--porcelain'],
                            capture_output=True, text=True, check=True)
    if status.stdout.strip():
        raise ValueError('工作区有未提交修改；请先处理修改，再执行 --pull')
    command(['git', '-C', root, 'pull', '--ff-only'], root)


def runtime_python(runtime):
    return runtime / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')


def install_runtime(root, runtime):
    """始终使用目标虚拟环境；包含绘图依赖，不改动全局 Python 环境。"""
    python = runtime_python(runtime)
    if not python.is_file():
        if runtime.exists():
            raise ValueError(f'目标目录存在但不是可用虚拟环境：{runtime}')
        command([sys.executable, '-m', 'venv', runtime], root)
    if not (runtime / 'pyvenv.cfg').is_file():
        raise ValueError(f'请指定独立虚拟环境：{runtime}')
    command([python, '-m', 'pip', 'install', 'setuptools>=77', 'wheel'], root)
    command([python, '-m', 'pip', 'install', '--no-build-isolation', '-c', root / 'constraints.txt',
             '-e', f'{root}[plots]'], root)
    return python


def check_runtime(root, python):
    """只读检查实际导入路径和内核 SHA256；不执行求解或写入 pyc。"""
    code = '''import json, sys
from importlib.metadata import version
from pathlib import Path
import gear_transmission, gear_contact
from gear_transmission.workflow import identity
root = Path(sys.argv[1]).resolve()
assert Path(gear_transmission.__file__).resolve().parent == root / 'src/gear_transmission', '工作流来自其他项目'
assert Path(gear_contact.__file__).resolve().parent == root / 'packages/cylindrical-gear-ltca/gear_contact', '内核来自其他项目'
print(json.dumps({'python': sys.executable, 'identity': identity(), 'matplotlib': version('matplotlib')}, ensure_ascii=False))
'''
    result = subprocess.run([str(python), '-B', '-c', code, str(root)], cwd=root,
                            capture_output=True, text=True, encoding='utf-8',
                            env={**os.environ, 'PYTHONUTF8': '1', 'PYTHONDONTWRITEBYTECODE': '1'})
    if result.returncode:
        raise ValueError('运行环境校验失败：' + result.stderr.strip())
    data = json.loads(result.stdout)
    return {'python': data['python'], 'versions': data['identity']['versions'],
            'snapshot_sha256': data['identity']['snapshot_sha256'], 'matplotlib': data['matplotlib']}


def load_receipt(destination):
    path = destination / RECEIPT
    if not path.is_file():
        return {}
    data = json.loads(path.read_text(encoding='utf-8'))
    if data.get('project') != 'gear-transmission-skills' or data.get('schema') != 1:
        raise ValueError(f'无法识别已有安装记录：{path}')
    return data


def preflight(source, destination, previous):
    """仅更新本库拥有且未被本地修改的目录；同名第三方目录不会被覆盖。"""
    for name, (path, expected) in source.items():
        target = destination / name
        if not target.exists() and not target.is_symlink():
            continue
        if target.resolve() == path.resolve():
            continue
        if target.is_symlink() or target.resolve() != destination / name:
            raise ValueError(f'目标为链接，拒绝替换：{target}')
        actual = file_hashes(target)
        if actual != expected and actual != previous.get('skills', {}).get(name):
            raise ValueError(f'同名目录不是本库安装内容，或含本地修改：{target}')


def safe_child(path, parent):
    """移动或递归清理前验证绝对路径，限定为明确目标目录的直接子项。"""
    if path.is_symlink() or path.resolve().parent != parent.resolve():
        raise ValueError(f'目录边界校验失败：{path}')


def sync_skills(source, destination, receipt):
    destination.mkdir(parents=True, exist_ok=True)
    for name, (path, expected) in source.items():
        target = destination / name
        if target.exists() and file_hashes(target) == expected:
            continue
        stage = Path(tempfile.mkdtemp(prefix='.gear-skills-stage-', dir=destination))
        backup = stage / 'previous'
        staged_skill = stage / name
        try:
            shutil.copytree(path, staged_skill,
                            ignore=shutil.ignore_patterns(*IGNORED, '*.pyc', '*.pyo'))
            if file_hashes(staged_skill) != expected:
                raise ValueError(f'暂存资源校验失败：{name}')
            safe_child(target, destination)
            safe_child(stage, destination)
            if target.exists():
                target.rename(backup)
            try:
                staged_skill.rename(target)
            except BaseException:
                if backup.exists():
                    backup.rename(target)
                raise
        finally:
            safe_child(stage, destination)
            shutil.rmtree(stage)
    receipt_path = destination / RECEIPT
    temporary = destination / (RECEIPT + '.tmp')
    temporary.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    temporary.replace(receipt_path)


def verify(root, source, destination, runtime, previous):
    if previous.get('root') != str(root) or previous.get('runtime') != str(runtime):
        raise ValueError('安装记录缺失或项目／环境位置不同，请先运行安装命令')
    for name, (_, expected) in source.items():
        if previous.get('skills', {}).get(name) != expected or file_hashes(destination / name) != expected:
            raise ValueError(f'完整技能目录与当前检出不一致：{name}，请重新同步')
    return check_runtime(root, runtime_python(runtime))


def main(argv=None):
    parser = argparse.ArgumentParser(description='安装／更新本库完整 Codex 技能与 Python 环境')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--pull', action='store_true', help='先执行 git pull --ff-only，再安装')
    mode.add_argument('--check', action='store_true', help='只读校验技能资源与运行环境')
    parser.add_argument('--dest', type=Path, default=Path(os.environ.get('CODEX_SKILLS_DIR', '~/.agents/skills')),
                        help='技能安装目录，默认 ~/.agents/skills')
    parser.add_argument('--runtime', type=Path, help='虚拟环境目录，默认项目 .venv')
    args = parser.parse_args(argv)
    root = project_root(Path(__file__).resolve().parent)
    destination = args.dest.expanduser().resolve()
    runtime = (args.runtime or root / '.venv').expanduser().resolve()
    if args.pull:
        pull(root)
    source = sources(root)
    if destination == root / 'skills' or any(destination == p.resolve() or p.resolve() in destination.parents
                                            for p, _ in source.values()):
        raise ValueError('不能将安装目标设为源技能目录或其内部')
    previous = load_receipt(destination)
    if not args.check:
        preflight(source, destination, previous)
        python = install_runtime(root, runtime)
        check_runtime(root, python)
        receipt = {'schema': 1, 'project': 'gear-transmission-skills', 'root': str(root),
                   'runtime': str(runtime), 'python': str(python),
                   'skills': {name: hashes for name, (_, hashes) in source.items()}}
        sync_skills(source, destination, receipt)
        previous = load_receipt(destination)
    checked = verify(root, source, destination, runtime, previous)
    print(json.dumps({'passed': True, 'mode': 'check' if args.check else 'install',
                      'destination': str(destination), 'skills': list(source), **checked},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f'安装／校验失败：{error}', file=sys.stderr)
        raise SystemExit(1)
