"""新写快照工具：从干净的公开锥齿轮工程复制原字节并登记文件及函数来源。"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def root_path():
    return next(p for p in Path(__file__).resolve().parents if (p / 'gear-skills.project.json').is_file())


def sync(source):
    root = root_path()
    source = Path(source).resolve()
    def git(*args):
        return subprocess.check_output(['git', '-C', str(source), *args], text=True, encoding='utf-8').strip()
    if git('status', '--porcelain'):
        raise ValueError('来源工程含未提交修改，不能登记为干净提交快照')
    commit = git('rev-parse', 'HEAD')
    version = next(ast.literal_eval(node.value) for node in ast.parse((source/'bevel_gear/__init__.py').read_text(encoding='utf-8-sig')).body
                   if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == '__version__' for t in node.targets))
    destination = root / 'packages/bevel-gear-design'
    selected = sorted(p for p in source.rglob('*') if p.is_file() and
        (p.relative_to(source).parts[0] == 'bevel_gear' and p.suffix in ('.py', '.html') or
         p.relative_to(source).parts[0] == 'docs' and p.suffix == '.md' or
         p.parent == source and p.name in ('LICENSE', 'pyproject.toml')) and
        '__pycache__' not in p.parts)
    files, functions = {}, []
    for path in selected:
        relative = path.relative_to(source).as_posix()
        data = path.read_bytes()
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        files[relative] = hashlib.sha256(data).hexdigest()
        if path.suffix == '.py':
            tree = ast.parse(data.decode('utf-8-sig'))
            for node in ast.walk(tree):
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    functions.append(dict(source_file=relative, source_sha256=files[relative],
                        source_function=node.name, start_line=node.lineno, end_line=node.end_lineno,
                        snapshot_file='packages/bevel-gear-design/' + relative,
                        implementation='原字节保留的公开重建计算内核；未在适配层改写公式'))
    assets = root / 'skills/bevel-gear-design/assets'
    assets.mkdir(parents=True, exist_ok=True)
    for name in ('straight_16_32_case.json', 'spiral_17_47_case.json',
                 'spiral_17_47_pinion_manufacturing.json', 'spiral_17_47_wheel_manufacturing.json',
                 'spiral_17_47_reference.json'):
        shutil.copyfile(source / 'examples' / name, assets / name)
    identity = dict(schema=1, source_repository='https://github.com/Lilililiyyy/Bevel-Gear-Design',
        source_commit=commit, source_clean=True, source_version=version,
        snapshot_path='packages/bevel-gear-design', files=files,
        implementation='公开重建计算内核快照；Skill工作流为新写适配，非厂商原始源码')
    for path, value in ((root/'src/gear_transmission/bevel_snapshot.json', identity),
                        (root/'docs/bevel-source-map.json', dict(**identity, functions=functions))):
        path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    subprocess.run([sys.executable, str(root/'scripts/release.py'), '--sync'], check=True)
    print(json.dumps(dict(commit=commit, files=len(files), functions=len(functions)), ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='同步公开锥齿轮计算内核快照')
    parser.add_argument('--source', required=True)
    sync(parser.parse_args().source)
