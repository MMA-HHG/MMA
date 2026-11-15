#!/usr/bin/env bash
# Recipe: Claude Code. Contract: docker/agents/README.md.
# Installation method: Anthropic's native installer, https://code.claude.com/docs/en/setup
#
# Usage (through mma-agents): claude.sh [VERSION]
#   (none)          the current release, which then updates itself in the background
#   stable, latest  that release channel, which then also governs the updates
#   X.Y.Z           that version, kept: background updates are turned off
set -euo pipefail

request="${1:-}"
export CLAUDE_CONFIG_DIR="${CLAUDE_CONFIG_DIR:-${MMA_AGENTS_STATE:-${HOME}/.mma-agents}/claude}"
settings="${CLAUDE_CONFIG_DIR}/settings.json"
command="${HOME}/.local/bin/claude"
request_file="${HOME}/.local/state/mma-agents/claude.request"

say() { echo "[mma-agents:claude] $*"; }
fail() { echo "[mma-agents:claude:error] $*" >&2; exit 1; }

case "${request}" in
    ""|stable|latest) ;;
    *) [[ "${request}" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] \
           || fail "Version must be X.Y.Z, 'stable' or 'latest', not '${request}'." ;;
esac
wanted="${request:-current}"
pinned=0
if [[ "${request}" =~ ^[0-9] ]]; then
    pinned=1
fi

installed_version() {
    [[ -x "${command}" ]] || return 0
    "${command}" --version 2>/dev/null | awk 'NR == 1 { print $1 }' || true
}

# Background updates are off for a pinned version and otherwise as Claude Code intends.
# The switch is DISABLE_AUTOUPDATER in the env block of Claude Code's own settings.json.
configure_updates() {
    python3 - "${settings}" "${pinned}" <<'PY'
import json
import os
import sys

path, pinned = sys.argv[1], sys.argv[2] == "1"
data = {}
if os.path.exists(path):
    with open(path, encoding="utf-8") as handle:
        text = handle.read().strip()
    try:
        data = json.loads(text) if text else {}
    except json.JSONDecodeError as exc:
        sys.exit(f"[mma-agents:claude:error] {path} is not valid JSON ({exc}); fix it and retry.")
elif not pinned:
    sys.exit(0)
env = data.get("env") or {}
if pinned:
    env["DISABLE_AUTOUPDATER"] = "1"
else:
    env.pop("DISABLE_AUTOUPDATER", None)
if env:
    data["env"] = env
else:
    data.pop("env", None)
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path + ".tmp", "w", encoding="utf-8") as handle:
    json.dump(data, handle, indent=2)
    handle.write("\n")
os.replace(path + ".tmp", path)
PY
}

# Idempotent and offline when the request is already satisfied.
current="$(installed_version)"
if [[ -n "${current}" && "$(cat "${request_file}" 2>/dev/null || true)" == "${wanted}" ]] \
        && { (( ! pinned )) || [[ "${current}" == "${request}" ]]; }; then
    configure_updates
    say "Claude Code ${current} is already installed (${wanted})."
    exit 0
fi

command -v curl >/dev/null || fail "curl is missing from the image; it belongs in the Dockerfile's package list."
mkdir -p "${CLAUDE_CONFIG_DIR}" "${HOME}/.local/bin" "$(dirname "${request_file}")"

say "Installing Claude Code (${wanted}) with the native installer."
installer=(bash)
if [[ -n "${request}" ]]; then
    installer+=(-s "${request}")
fi
curl -fsSL https://claude.ai/install.sh | "${installer[@]}" \
    || fail "The native installer failed; see its messages above."

current="$(installed_version)"
[[ -n "${current}" ]] || fail "The installer finished, but ${command} reports no version."
(( ! pinned )) || [[ "${current}" == "${request}" ]] \
    || fail "Requested ${request}, but ${command} reports ${current}."

configure_updates
echo "${wanted}" > "${request_file}"
if (( pinned )); then
    say "Installed Claude Code ${current}; background updates are off."
else
    say "Installed Claude Code ${current} (${wanted}); it updates itself as Claude Code intends."
fi
