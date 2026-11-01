"""维护工具：同步独立项目的计算、案例、后处理及对标输入，保留原始字节。"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def selected_files(source):
    """从已跟踪文件选取运行资源；不把历史结果、缓存和 Git 目录打包。"""
    tracked = subprocess.check_output(['git', '-C', str(source), 'ls-files', '-z']).decode().split('\0')
    return sorted(name for name in tracked if name and (
        name in {'LICENSE', 'pyproject.toml', 'constraints.txt'}
        or name.startswith('gear_contact/') and name.endswith('.py')
        or name.startswith('docs/') and name.endswith('.md')
        or name.startswith('examples/masta_reference/')
        or name.startswith('examples/') and len(Path(name).parts) == 2 and Path(name).suffix in {'.py', '.json'}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True, help='独立 LTCA 项目的干净检出目录')
    parser.add_argument('--check', action='store_true', help='只校验完整资源清单及字节，不写入')
    args = parser.parse_args()
    root = next(p for p in Path(__file__).resolve().parents if (p/'gear-skills.project.json').is_file())
    source = args.source.resolve()
    git = lambda *a: subprocess.check_output(['git', '-C', str(source), *a]).decode().strip()
    if git('status', '--porcelain'):
        raise ValueError('来源项目含未提交修改，请先固定来源版本')
    files = selected_files(source)
    target = root/'packages/cylindrical-gear-ltca'
    manifest_path = root/'src/gear_transmission/kernel_snapshot.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    hashes = {name: hashlib.sha256((source/name).read_bytes()).hexdigest() for name in files}
    mapping_path = root/'docs/source-map.json'
    mapping = json.loads(mapping_path.read_text(encoding='utf-8'))
    example_symbols = {
        'examples/run_cylindrical_workflow.py': ['prepare_case', 'calculate_case', 'run'],
        'examples/case_contact_plotting.py': ['repeat_cycle_for_display', '_project_line_samples', 'plot_case_contact'],
        'examples/case_native_contact_plotting.py': ['project_contact_pass_samples', 'plot_native_contact_comparison'],
        'examples/case_geometry_plotting.py': ['plot_case_geometry'],
        'examples/case_sensitivity.py': ['run_case_sensitivity'],
    }
    records = {(r['source_file'], r['symbol']): r for r in mapping['symbols']}
    for name, symbols in example_symbols.items():
        for symbol in symbols:
            records.setdefault((name, symbol), {'source_project': 'cylindrical-gear-ltca',
                'source_file': name, 'integrated_file': 'packages/cylindrical-gear-ltca/'+name,
                'symbol': symbol, 'implementation': '原字节保留的现有案例／后处理实现'})
    for (name, symbol), record in records.items():
        node = next(n for n in ast.parse((source/name).read_text(encoding='utf-8')).body
                    if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.name == symbol)
        record.update(lines=[node.lineno, node.end_lineno], sha256=hashes[name])
    mapping.update(symbols=list(records.values()), source_base_commit=git('rev-parse', 'HEAD'))
    if args.check:
        assert manifest['source_base_commit'] == git('rev-parse', 'HEAD'), '来源提交不一致'
        assert manifest['files'] == hashes, '登记资源不完整或来源哈希变化'
        assert all((target/name).read_bytes() == (source/name).read_bytes() for name in files), '快照字节不一致'
        assert json.loads(mapping_path.read_text(encoding='utf-8')) == mapping, '函数来源映射未更新'
    else:
        for name in set(manifest['files'])-set(files):
            stale = (target/name).resolve()
            if not stale.is_relative_to(target.resolve()):
                raise ValueError('旧快照文件路径超出集成目录：'+name)
            stale.unlink()
        for name in files:
            (target/name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source/name, target/name)
        manifest.update(source_base_commit=git('rev-parse', 'HEAD'), files=hashes,
                        runtime_scope='完整计算模块、案例编排、绘图后处理、案例输入及公开对标资源；不含生成结果和缓存')
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        mapping_path.write_text(json.dumps(mapping, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'checked': args.check, 'source_commit': git('rev-parse', 'HEAD'), 'files': len(files)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
