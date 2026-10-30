#!/usr/bin/env bash
# 新写 Bash 入口，安装逻辑与 Windows 共用同一 Python 实现。
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${GEAR_SKILLS_PYTHON:-}" ]]; then
  PYTHON="$GEAR_SKILLS_PYTHON"
elif [[ "${OSTYPE:-}" == msys* || "${OSTYPE:-}" == cygwin* ]]; then
  # Windows 的 python3 可能只是商店占位程序，优先使用已安装的 python。
  PYTHON=python
elif command -v python3 >/dev/null 2>&1; then
  PYTHON=python3
else
  PYTHON=python
fi
exec "$PYTHON" -B "$SCRIPT_DIR/update_codex_skills.py" "$@"
