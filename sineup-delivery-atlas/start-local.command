#!/bin/zsh

set -e

PROJECT_DIR="${0:A:h}"
CODEX_RUNTIME="$HOME/.cache/codex-runtimes/codex-primary-runtime/dependencies"

if [[ -x "$CODEX_RUNTIME/bin/fallback/pnpm" ]]; then
  export PATH="$CODEX_RUNTIME/node/bin:$CODEX_RUNTIME/bin/override:$CODEX_RUNTIME/bin/fallback:$PATH"
fi

cd "$PROJECT_DIR"
exec pnpm dev
