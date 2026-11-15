# AI coding agents in the Docker image

The MMA Docker image can contain AI coding agents, such as Codex CLI and Claude Code, which then work in the mounted repository together with MMA's AI support. Agents are optional: without the options below, the image is the one described in [README.md, Reference Docker installation](../../README.md#reference-docker-installation-and-running-the-code). This file describes how to use the agents and how to add or maintain a recipe for one.

## Using an agent

Build the image once, as in the README, and choose the agents when you start a container. On top of the usual `docker run` command, add the agents' volume and your selection:

```bash
docker build -t mma .
docker run --name CONTAINER_NAME --hostname CONTAINER_NAME -v .:/MMA -w /MMA -p 8888:8888 \
    -v mma-agents:/home/active_electron/.mma-agents -e AI_AGENTS=claude,codex -it mma
```

The first start downloads the selected agents, so it needs internet access; `-e AI_AGENTS=claude` or `-e AI_AGENTS=codex` installs only one of them. The volume `mma-agents` keeps logins and settings for later containers; replace the home path if you build with a different `USERNAME`. `docker start -ai CONTAINER_NAME` restarts the same container with everything in place, without downloading again.

In the container, start an agent in `/MMA` from the main shell or a JupyterLab terminal, and sign in:

- Claude Code: run `claude` and follow the sign-in prompt.
- Codex CLI: run `codex login --device-auth`, then `codex`. Device-code login avoids a browser callback into the container. It must be enabled for your OpenAI account or workspace; otherwise use another method from [OpenAI's authentication documentation](https://developers.openai.com/codex/auth/).

Shells opened with `docker exec` lack the variables that the entrypoint sets. Neither the image nor the repository contains credentials, and there is no sign-in during an image build.

To have agents ready in the image instead, so that containers start without downloading, select them at build time and run the image with the same command:

```bash
docker build --build-arg AI_AGENTS=claude,codex -t mma-ai .
```

The same options apply when you build from `docker/Dockerfile_Windows`.

`AI_AGENTS` takes one agent or a comma-separated list. Without a version, an agent's current release is installed; subsequent update behaviour is agent-specific (see below). An entry may also name a version, such as `codex@X.Y.Z` or `claude@X.Y.Z`. The recipe then verifies that version and suppresses automatic update checks or updates as appropriate for that agent. For Claude Code, `claude@stable` and `claude@latest` choose a release channel instead. Inside a container, `mma-agents list` shows the supported agents, and `mma-agents install claude@X.Y.Z` installs an agent or changes its version. An unknown name or a failed requested installation stops the build or container start with an error. To start a container without agents, omit `AI_AGENTS` or set it to `none`.

### Git

`/MMA` is your checkout on the host, so changes and commits appear there immediately. Set your Git identity once in the container; it is stored in the repository's local Git configuration:

```
git config user.name "Your Name"
git config user.email "you@example.org"
```

The container has no Git credentials; push from the host.

### Autonomy and isolation

The image leaves each agent's approval and sandbox settings unchanged. Its normal rules determine which edits and commands can run without another prompt. Choose the collaboration mode inside the agent. In Codex, inspect `/status` and use `/permissions`; in Claude Code, for example, use `Shift+Tab` or `claude --permission-mode acceptEdits`, while `claude --dangerously-skip-permissions` removes all prompts. Before you choose, consider what the container isolates and what it does not:

- Isolated: files on the host outside the repository, and the host's processes.
- Not isolated: the repository, including its Git history, which an agent can change or delete; network access, which is unrestricted; and the agents' logins in `~/.mma-agents`, which every process in the container can read.
- The container user has `sudo` without a password, so an agent can act as root in the container. On Linux hosts, files it creates as root in `/MMA` then belong to root on the host.

Do not mount host secrets, such as `~/.ssh` or cloud credentials, into the container.

## How the installation works

- `docker/mma-agents` is installed in the image as the command `mma-agents`. It reads a list such as `claude,other@1.2.0`, checks every name against the recipes present here, runs a shared setup once and then each recipe. An unknown name stops with a message pointing to this file.
- The Dockerfile runs `mma-agents install "$AI_AGENTS"` at the end of the build, as the container user. The entrypoint runs the same command at container start when the variable `AI_AGENTS` is set, before compiling MMA, so that a failed request stops the start quickly. Without `AI_AGENTS`, it re-applies the image's own build-time selection, offline when nothing changed; this keeps a pinned version's settings even when a reused volume hides the image's copy, and a failure there only produces a warning. `mma-agents list` shows the supported agents.
- The shared setup creates the state directory `~/.mma-agents` and lets Git work in the mounted repository (`safe.directory`). Agent-neutral system packages, such as `curl`, are in the Dockerfile's package list.
- The set of supported agents is simply the set of `*.sh` files here: the file name is the agent's name.

## Recipe contract

A recipe `docker/agents/<name>.sh` must:

1. **Be named** with lowercase letters, digits and single hyphens; the name is what users write in `AI_AGENTS`.
2. **Accept an optional version** as its only argument. Without one, install the vendor's current release. With one, verify that exact version and suppress unsolicited updates through a documented setting or launcher option. Document the agent's update behaviour, including how to request a newer version. Leave user configuration intact where possible.
3. **Work without root.** It runs as the container user, at build time or at container start, possibly where `sudo` is unavailable. Install under `$HOME` only; never call `sudo` or a system package manager. It must not rely on its own file location.
4. **Be idempotent and offline when there is nothing to do.** If the request is already satisfied, exit successfully without network access; the entrypoint runs recipes on every container start while `AI_AGENTS` is set.
5. **Put its command on `PATH` for every shell**, in `~/.local/bin`, which the image adds to `PATH`; prefer the command exactly as the vendor's installer leaves it. If the agent needs an environment variable, for example to find its state, add one `ENV` line to the agents block of both Dockerfiles; it is inert while the agent is absent. Never use shell start-up files, so that `docker exec` and non-interactive shells behave the same.
6. **Keep the agent's login and settings** in `$MMA_AGENTS_STATE/<name>`, through the agent's own setting or environment variable. One optional volume then keeps every agent's login. Document any user resources, such as personal skills, that use another location.
7. **Stay neutral.** No credentials or API keys, no permission or autonomy settings, no MCP servers: the user decides these at run time.
8. **Use the vendor's documented installation method**, check the result with the agent's `--version`, prefix its messages with `[mma-agents:<name>]` and exit non-zero on failure.

## Adding an agent

1. Copy `claude.sh` to `<name>.sh` and adapt it to the agent's documented installer.
2. Test it in a running container without rebuilding, from the mounted repository:

       MMA_AGENTS_RECIPES=/MMA/docker/agents mma-agents install <name>

   Run it twice: the second run must report that nothing needed doing, without network access.
3. Build with `--build-arg AI_AGENTS=<name>` and start a container with `-e AI_AGENTS=<name>` from the plain image; both must give a working `<name> --version`.
4. A recipe only installs the agent. Whether the agent finds `AGENTS.md` and the skills by itself, or needs a signpost, is decided in the AI-support handbook, `docs/ai/README.md`, section 10.

## Claude Code (`claude.sh`)

- Installed with Anthropic's native installer. `claude` installs the current release, which then updates itself in the background, as Claude Code intends. `claude@stable` or `claude@latest` installs that release channel, which then also governs the updates. `claude@X.Y.Z` installs that version and keeps it: the recipe sets `DISABLE_AUTOUPDATER` to `1` in the `env` block of Claude Code's own `settings.json`, and removes it again when a later install requests no version.
- `~/.local/bin/claude` is the installer's own link, so `claude update` and `claude doctor` behave as documented; `claude doctor` also shows whether auto-updates are on.
- The image sets `CLAUDE_CONFIG_DIR` to `~/.mma-agents/claude`, with one `ENV` line in the agents block that is inert while Claude Code is absent. The login, settings and personal skills therefore live there, not in `~/.claude`.
- Claude Code finds the project's instructions and skills in the mounted repository through `CLAUDE.md` and `.claude/skills/`. The image adds nothing to the repository's `.claude/` folder.


## Codex CLI (`codex.sh`)

- Installed without root from the Linux standalone binaries in
  [OpenAI's releases](https://github.com/openai/codex/releases), as documented in
  [OpenAI's repository](https://github.com/openai/codex#quickstart). Linux x86-64
  and ARM64 are supported, including containers built with the Windows Dockerfile.
  No Node.js or npm dependency is added to the scientific image.
- `codex` and `codex@latest` install the latest stable release on first use.
  Repeating a satisfied request stays offline and keeps the installed version.
  To select a newer release, run `mma-agents install codex@X.Y.Z` with its real
  version. The recipe does not run an update service.
- `codex@X.Y.Z` verifies the downloaded executable's version before replacing a
  working installation. A small launcher adds the documented
  `check_for_update_on_startup=false` CLI override for pinned installs, without
  editing `config.toml`. A later unversioned installation removes this override.
  Pinning controls the installed CLI release, not the hosted model or service.
- The image sets `CODEX_HOME` to the `codex` subdirectory of the state volume.
  File-based credentials, settings and sessions live there. If you select an OS
  keyring or ephemeral credential storage, that volume alone does not persist the
  login. Codex's normal personal-skill location is `~/.agents/skills`; mount it
  separately if you want personal skills to survive container replacement.
- The executable and pin marker stay under the user's `~/.local/share/mma-agents`
  directory, outside the credentials volume. Reusing a login volume therefore
  does not hide the Codex version policy. The launcher is in `~/.local/bin`.
- Codex uses the existing root `AGENTS.md`, module-local instructions and
  `.agents/skills/` directly. Run `/skills` to inspect discovery. No duplicated
  Codex-specific project instructions or permission defaults are generated.
- Start from `/MMA` in the main shell or a JupyterLab terminal, so commands inherit
  the MMA environment prepared by the entrypoint. User settings remain under the
  user's control. See [Codex authentication](https://developers.openai.com/codex/auth/)
  and [configuration](https://developers.openai.com/codex/config-reference/).

### Installation checks

Run `python3 tools/test_ai_agents.py` for offline tests of installation selection,
version handling, failure propagation and launcher behaviour. These tests use
small fake release archives; they do not test a real Codex download or login.
Validate both Docker build paths and authentication in an actual Docker environment
before calling the integration tested end to end.
