# AI coding agents in the Docker image

The MMA Docker image can contain AI coding agents: Claude Code, Codex CLI and Antigravity CLI. They work in the mounted repository together with MMA's AI support. Agents are optional: MMA builds, runs and is tested without them, and without the options below, the image is the one described in [README.md, Reference Docker installation](../../README.md#reference-docker-installation-and-running-the-code). This file describes how agents work with MMA, how to use them in the container and how to add or maintain a recipe for one.

## How agents work with MMA

- **Entry point.** Agents start from [AGENTS.md](../../AGENTS.md), which routes them to MMA's workflows and skills. It and the skills in `.agents/skills/` are where the agents installed here look for them, Claude Code through the generated `CLAUDE.md`; section 10 of the handbook records which discovery has been verified.
- **Autonomy is your choice.** You decide inside the agent how freely it edits files and runs commands; see [Autonomy and isolation](#autonomy-and-isolation). Unless you authorise them in advance, agents ask before the consequential actions listed in `AGENTS.md`, such as pushing or deleting simulation data.
- **AI-assisted work is traced.** Every AI-assisted commit carries `Assisted-by:` and `Numerical-impact:` trailers, and pull requests state the degree of assistance, as defined in [CONTRIBUTING.md](../../CONTRIBUTING.md).
- **Open to extension.** Anyone can add skills, workflows, records or agent recipes, on their own responsibility. The conventions are in the AI-support handbook, [docs/ai/README.md](../../docs/ai/README.md).

## Using an agent

Build the image once, as in the README, and choose the agents when you start a container. On top of the usual `docker run` command, add the agents' volume and your selection:

```bash
docker build -t mma .
docker run --name CONTAINER_NAME --hostname CONTAINER_NAME -v .:/MMA -w /MMA -p 8888:8888 \
    -v mma-agents:/home/active_electron/.mma-agents -e AI_AGENTS=claude,codex -it mma
```

The first start downloads the selected agents, so it needs internet access; `-e AI_AGENTS=claude`, `-e AI_AGENTS=codex` or `-e AI_AGENTS=antigravity` installs only one of them. The volume `mma-agents` keeps logins and settings for later containers; replace the home path if you build with a different `USERNAME`. `docker start -ai CONTAINER_NAME` restarts the same container with everything in place, without downloading again.

In the container, start an agent in `/MMA` from the main shell or a JupyterLab terminal, and sign in:

- Claude Code: run `claude` and follow the sign-in prompt.
- Codex CLI: run `codex login --device-auth`, then `codex`. Device-code login avoids a browser callback into the container. It must be enabled for your OpenAI account or workspace; otherwise use another method from [OpenAI's authentication documentation](https://developers.openai.com/codex/auth/).
- Antigravity CLI: run `agy` and choose a sign-in method. The container has no keyring for the login; see [Antigravity CLI](#antigravity-cli-antigravitysh) below.

Shells opened with `docker exec` lack the variables that the entrypoint sets. Neither the image nor the repository contains credentials, and there is no sign-in during an image build.

To have agents ready in the image instead, so that containers start without downloading, select them at build time and run the image with the same command:

```bash
docker build --build-arg AI_AGENTS=claude,codex -t mma-ai .
```

The same options apply when you build from `docker/Dockerfile_Windows`.

`AI_AGENTS` takes one agent or a comma-separated list. Without a version, an agent's current release is installed; subsequent update behaviour is agent-specific (see below). An entry may also name a version, such as `codex@X.Y.Z` or `claude@X.Y.Z`. The recipe then verifies that version and suppresses automatic update checks or updates as appropriate for that agent. For Claude Code, `claude@stable` and `claude@latest` choose a release channel instead. For Antigravity CLI, whose installer offers only the current release, `antigravity@X.Y.Z` keeps a version that is installed or current. Inside a container, `mma-agents list` shows the supported agents, and `mma-agents install claude@X.Y.Z` installs an agent or changes its version. An unknown name or a failed requested installation stops the build or container start with an error. To start a container without agents, omit `AI_AGENTS` or set it to `none`.

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
6. **Keep the agent's login and settings** in `$MMA_AGENTS_STATE/<name>`, through the agent's own setting or environment variable; only where the agent offers neither, through a link from its fixed location, documented in the agent's section. One optional volume then keeps every agent's login. Document any user resources, such as personal skills, that use another location.
7. **Stay neutral.** No credentials or API keys, no permission or autonomy settings, no MCP servers: the user decides these at run time.
8. **Use the vendor's documented installation method**, check the result with the agent's `--version`, prefix its messages with `[mma-agents:<name>]` and exit non-zero on failure.

## Adding an agent

1. Copy `claude.sh` to `<name>.sh` and adapt it to the agent's documented installer.
2. Test it in a running container without rebuilding, from the mounted repository:

       MMA_AGENTS_RECIPES=/MMA/docker/agents mma-agents install <name>

   Run it twice: the second run must report that nothing needed doing, without network access.
3. Add offline tests to `ai/tests/test_ai_agents.py` that replace the vendor download with a fake, as the tests of the existing recipes do, and run `python3 ai/tests/test_ai_agents.py`.
4. Build with `--build-arg AI_AGENTS=<name>` and start a container with `-e AI_AGENTS=<name>` from the plain image; both must give a working `--version` of the agent's command.
5. A recipe only installs the agent. Whether the agent finds `AGENTS.md` and the skills by itself, or needs a signpost, is decided in the AI-support handbook, `docs/ai/README.md`, section 10.

## Claude Code (`claude.sh`)

- Installed with Anthropic's native installer. `claude` installs the current release, which then updates itself in the background, as Claude Code intends. `claude@stable` or `claude@latest` installs that release channel, which then also governs the updates. `claude@X.Y.Z` installs that version and keeps it: the recipe sets `DISABLE_AUTOUPDATER` to `1` in the `env` block of Claude Code's own `settings.json`, and removes it again when a later install requests no version.
- `~/.local/bin/claude` is the installer's own link, so `claude update` and `claude doctor` behave as documented; `claude doctor` also shows whether auto-updates are on.
- The image sets `CLAUDE_CONFIG_DIR` to `~/.mma-agents/claude`, with one `ENV` line in the agents block that is inert while Claude Code is absent. The login, settings and personal skills therefore live there, not in `~/.claude`.
- Claude Code finds the project's instructions and skills in the mounted repository through `CLAUDE.md` and `.claude/skills/`. The image adds nothing to the repository's `.claude/` folder.


## Codex CLI (`codex.sh`)

- Installed without root from the complete Linux `codex-package` archives in
  [OpenAI's releases](https://github.com/openai/codex/releases), as documented in
  [OpenAI's repository](https://github.com/openai/codex#quickstart). The whole
  package is kept, so the helper executables next to `codex` are available; the
  recipe checks them before it switches to a new package. Linux x86-64 and ARM64
  are supported, including containers built with the Windows Dockerfile.
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

## Antigravity CLI (`antigravity.sh`)

- Installed with Google's installer, `curl -fsSL https://antigravity.google/cli/install.sh | bash`, as on [Antigravity's download page](https://antigravity.google/download/). The command is `agy`, which the installer leaves in `~/.local/bin`. The installer also adds `~/.local/bin` to `PATH` in `~/.bashrc` and `~/.profile`; the image has it on `PATH` already, so these lines change nothing, and the recipe does not rely on them. Google's authentication guide mentions installer options `--skip-path` and `--skip-aliases` against this, but the installer script rejects them (checked with version 1.5.18), so the recipe does not pass them. As with the documented `curl | bash`, the installer gets no terminal input. The executable is large, about 270 MB, which matters when it is built into an image.
- `antigravity` and `antigravity@latest` install the current release, which then updates itself in the background, as the Antigravity CLI intends. Repeating a satisfied request stays offline.
- `antigravity@X.Y.Z` keeps that version. Google's installer offers only the current release, so X.Y.Z must already be installed or be the current release; otherwise the recipe stops, names the current release and leaves the installation unchanged. The kept executable moves to `~/.local/share/mma-agents/antigravity/`, and a small launcher in `~/.local/bin/agy` sets the documented `AGY_CLI_DISABLE_AUTO_UPDATE=true`. A later unversioned request restores the installer's layout, offline, and background updates resume. A version older than the installed one is refused without a download, since the CLI only updates forwards.
- An image built with `antigravity@X.Y.Z` keeps that version in every container started from it: the entrypoint re-applies the image's selection at every start, so an unversioned request inside such a container lasts only until the next start. To let the CLI update itself there, start the container with `-e AI_AGENTS=antigravity`.
- Before the installer runs, an existing installation waits in `~/.local/share/mma-agents/antigravity.previous`. It is put back if the installation fails or is interrupted, for example with `Ctrl+C`, and at the latest by the next run if the container was stopped meanwhile.
- The Antigravity CLI keeps its settings, logs and conversations in `~/.gemini/antigravity-cli` and offers no setting or variable to move them. The recipe therefore makes that path a link to `~/.mma-agents/antigravity`, so that the state volume keeps them, as item 6 of the recipe contract allows in this case. State found at the old path moves into the volume, unless the volume already holds state; then it is set aside as `~/.gemini/antigravity-cli.before-mma-agents.<date>`. A link pointing elsewhere is redirected to the volume, with a message. Other files in `~/.gemini`, such as personal skills in `~/.gemini/config/skills`, are not kept; mount them separately if you need them.
- Sign-in: the CLI keeps its login in the operating system's keyring (Secret Service over D-Bus on Linux), which the container lacks. Whether it then falls back to a file in its state directory, which the volume would keep, has not been verified. Google documents two alternatives: a Gemini API key in `GEMINI_API_KEY`, used when `settings.json` sets `"modelProvider": "gemini"`, and Application Default Credentials with `AGY_ADC_AUTH=true`; see [Antigravity's authentication guide](https://antigravity.google/docs/install/). The recipe sets neither.
- The Antigravity CLI reads `AGENTS.md` and the skills in `.agents/skills/` directly, from its working directory and the folders above it up to the repository root (see [its discovery rules](https://antigravity.google/docs/customization-discovery/)). Start it in `/MMA` for work spanning modules; started in a module directory, it also reads that module's `AGENTS.md`.
- Tested on 2026-11-27 with Antigravity CLI 1.5.18 in a running MMA container, with a separate home directory: installation, offline repeat, keeping and releasing a version, the refusal of unavailable versions, and `Ctrl+C` during the download. Image builds with `--build-arg AI_AGENTS=antigravity` and sign-in have not been tested.

## Installation checks

Run `python3 ai/tests/test_ai_agents.py` for offline tests of the Codex CLI and Antigravity CLI recipes, the `mma-agents` dispatcher and the entrypoint: installation selection, version handling, failure propagation, launcher behaviour and, for Antigravity CLI, the state link. The tests replace vendor downloads with small fakes, such as release archives and installer scripts; they test neither real downloads nor sign-in. Validate both Docker build paths and sign-in in an actual Docker environment before calling an integration tested end to end. `pytest ai -v` runs these tests as well, together with the structural checks of the AI support; MMA's default test run does not include them.
