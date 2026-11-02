"""发布元数据工具：只从根pyproject读取库版本，从快照读取内核版本。"""
import argparse
from importlib.metadata import version as installed_version
import json
from pathlib import Path
import re


KERNELS = {
    'cylindrical-gear-analysis': ('kernel_snapshot.json', 'kernel_version', 'source_base_commit'),
    'bevel-gear-design': ('bevel_snapshot.json', 'source_version', 'source_commit'),
}


def project_version(root):
    """兼容Python3.10，仅读取[project]内唯一的X.Y.Z字面值，完整TOML由构建工具校验。"""
    text = (root/'pyproject.toml').read_text(encoding='utf-8')
    section = re.search(r'(?ms)^\[project\]\s*\n(.*?)(?=^\[|\Z)', text)
    values = re.findall(r'''(?m)^version\s*=\s*["'](\d+\.\d+\.\d+)["']\s*$''', section[1] if section else '')
    if len(values) != 1:
        raise ValueError('根pyproject.toml的[project]须声明唯一的X.Y.Z版本')
    return values[0]


def kernel_records(root):
    marker = json.loads((root/'gear-skills.project.json').read_text(encoding='utf-8'))
    skills = [Path(path).name for path in marker['skills']]
    if set(skills) != set(KERNELS):
        raise ValueError('技能登记与发布工具的内核映射不一致，请登记新增技能的快照来源')
    result = {}
    for name in skills:
        file, version_key, commit_key = KERNELS[name]
        snapshot = json.loads((root/'src/gear_transmission'/file).read_text(encoding='utf-8'))
        result[name] = {'version': snapshot[version_key], 'commit': snapshot[commit_key],
                        'repository': snapshot['source_repository'].rstrip('/')}
    return result


def block(text, name, replacement):
    pattern = rf'(?s)<!-- {name}:start -->\n.*?\n<!-- {name}:end -->'
    if len(re.findall(pattern, text)) != 1:
        raise ValueError('文档必须保留唯一的自动维护区块：'+name)
    return re.sub(pattern, lambda match: f'<!-- {name}:start -->\n{replacement}\n<!-- {name}:end -->', text)


def document(root, name, version, records):
    text = (root/name).read_text(encoding='utf-8')
    introduction = (f'当前库版本为 **{version}**。所有 Skill 随库统一发布；所用内核版本见下表。'
                    if name == 'README.md' else
                    f'Library version: **{version}**. All skills ship together; integrated kernel versions are listed below.')
    text = block(text, 'library-version', introduction)
    table_match = re.search(r'(?s)<!-- skills:start -->\n(.*?)\n<!-- skills:end -->', text)
    if table_match is None:
        raise ValueError(name+'缺少技能表格自动维护区块')
    lines, found = [], set()
    for line in table_match[1].splitlines():
        match = re.search(r'skills/([^/]+)/SKILL\.md', line)
        if match:
            skill = match[1]
            if skill not in records or skill in found:
                raise ValueError(name+'含未登记或重复技能：'+skill)
            found.add(skill)
            columns = [part.strip() for part in line.strip().strip('|').split('|')]
            if len(columns) != 5:
                raise ValueError(name+'技能表格应有五列：技能、能力、内核版本、状态、独立项目')
            record = records[skill]
            columns[2] = f"{record['version']} ([{record['commit'][:7]}]({record['repository']}/commit/{record['commit']}))"
            line = '| '+' | '.join(columns)+' |'
        lines.append(line)
    if found != set(records):
        raise ValueError(name+'技能表格与登记清单不完整对应')
    return block(text, 'skills', '\n'.join(lines))


def check_skills(root):
    for name in KERNELS:
        content = (root/'skills'/name/'SKILL.md').read_text(encoding='utf-8')
        front = content.split('---', 2)[1]
        if re.search(r'(?m)^\s*version\s*:', front):
            raise ValueError(name+'仍有独立Skill版本，请删除；技能随整库发布')


def synchronize(root, *, write=False, check_installed=False):
    root = Path(root).resolve()
    version = project_version(root)
    check_skills(root)
    records = kernel_records(root)
    outputs = {}
    for name in ('gear-skills.project.json', '.codex-plugin/plugin.json'):
        data = json.loads((root/name).read_text(encoding='utf-8'))
        data['version'] = version
        outputs[name] = json.dumps(data, ensure_ascii=False, indent=2)+'\n'
    for name in ('README.md', 'README_EN.md'):
        outputs[name] = document(root, name, version, records)
    different = [name for name, value in outputs.items() if (root/name).read_text(encoding='utf-8') != value]
    if different and not write:
        raise ValueError('发布元数据不同步：'+', '.join(different)+'；请运行 python scripts/release.py --sync')
    for name in different:
        (root/name).write_text(outputs[name], encoding='utf-8', newline='\n')
    if check_installed:
        from gear_transmission import __version__
        if installed_version('gear-transmission-skills') != version or __version__ != version:
            raise ValueError('实际安装版本与根pyproject不一致，请重新安装本项目')
    return {'library_version': version, 'kernels': records, 'updated': different}


def release_notes(root, version):
    """从两份README当前版本条目生成Release正文，历史记录不改写。"""
    sections = []
    for name, language in (('README.md', '中文'), ('README_EN.md', 'English')):
        text = (root/name).read_text(encoding='utf-8')
        pattern = rf'(?ms)^- \*\*[^\n]*v{re.escape(version)}\*\*\n(.*?)(?=^- \*\*|\Z)'
        match = re.search(pattern, text)
        if not match:
            raise ValueError(name+'缺少当前版本更新记录')
        sections += [f'## {language}', '', '\n'.join(line[2:] if line.startswith('  ') else line
                     for line in match[1].strip('\n').splitlines()), '']
    sections += ['## 集成内核 / Integrated kernels', '']
    for skill, record in kernel_records(root).items():
        sections.append(f"- `{skill}`: {record['version']} ([{record['commit'][:7]}]({record['repository']}/commit/{record['commit']}))")
    return '\n'.join(sections)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--sync', action='store_true', help='同步自动维护的库版本和README内核表格')
    action.add_argument('--check', action='store_true', help='只读校验，默认操作')
    action.add_argument('--version', action='store_true', help='输出根pyproject版本')
    action.add_argument('--notes', action='store_true', help='输出当前版本中英文Release正文')
    parser.add_argument('--check-installed', action='store_true', help='同时校验实际安装和导入版本')
    args = parser.parse_args()
    root = next(p for p in Path(__file__).resolve().parents if (p/'gear-skills.project.json').is_file())
    if args.version:
        print(project_version(root))
    elif args.notes:
        print(release_notes(root, project_version(root)), end='')
    else:
        print(json.dumps(synchronize(root, write=args.sync, check_installed=args.check_installed),
                         ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
