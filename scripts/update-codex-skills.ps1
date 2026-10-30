# 新写 PowerShell 入口；参数原样传给跨平台安装程序。
$ErrorActionPreference = 'Stop'
$gearSkillsPython = if ($env:GEAR_SKILLS_PYTHON) { $env:GEAR_SKILLS_PYTHON } else { 'python' }
& $gearSkillsPython -B (Join-Path $PSScriptRoot 'update_codex_skills.py') @args
exit $LASTEXITCODE
