#!/usr/bin/env python3
"""Offline installer checks; vendor downloads and the root identity are stubbed."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class AgentInstallationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.home = self.work / "user home"
        self.home.mkdir()
        self.stubs = self.work / "stubs"
        self.stubs.mkdir()
        self.env = dict(os.environ)
        for key in ("CODEX_HOME", "MMA_AGENTS_STATE", "AI_AGENTS", "FAIL_CURL", "CORRUPT_ARCHIVE"):
            self.env.pop(key, None)
        self.env.update(HOME=str(self.home),
                        PATH=f"{self.stubs}:{self.home / '.local/bin'}:{os.environ['PATH']}",
                        MMA_AGENTS_RECIPES=str(ROOT / "docker/agents"),
                        GIT_CONFIG_GLOBAL=str(self.home / ".gitconfig"),
                        GIT_CONFIG_NOSYSTEM="1", TEST_ARCH="x86_64",
                        TEST_SYSTEM="Linux", RELEASE_VERSION="1.2.3",
                        FETCH_LOG=str(self.work / "fetches"),
                        ARGS_LOG=str(self.work / "arguments.json"))
        self.stub("uname", '#!/bin/sh\ncase "$1" in\n-s) echo "${TEST_SYSTEM:-Linux}";;\n-m) echo "${TEST_ARCH:-x86_64}";;\nesac\n')
        self.stub("id", '#!/bin/sh\nif [ "$1" = "-u" ]; then echo 1000; else /usr/bin/id "$@"; fi\n')
        self.stub("curl", f"#!{sys.executable}\n" + '''import io, os, pathlib, sys, tarfile
args = sys.argv[1:]
url = next(a for a in args if a.endswith('.tar.gz'))
with open(os.environ['FETCH_LOG'], 'a') as log:
    log.write(url + '\\n')
if os.environ.get('FAIL_CURL'):
    sys.exit(22)
output = pathlib.Path(args[args.index('-o') + 1])
if os.environ.get('CORRUPT_ARCHIVE'):
    output.write_bytes(b'not a gzip archive')
    sys.exit(0)
version = os.environ.get('RELEASE_VERSION', '1.2.3')
program = ('#!' + sys.executable + '\\nimport json, os, pathlib, sys\\n'
           + 'if "--version" in sys.argv:\\n    print("codex-cli ' + version + '")\\n'
           + 'elif "--help" in sys.argv:\\n'
           + '    assert not ' + repr(bool(os.environ.get('NO_NO_DAEMON'))) + '\\n'
           + '    assert (pathlib.Path(__file__).resolve().parent / "resources/bundle.txt").read_text() == "resource"\\n'
           + 'else:\\n    with open(os.environ["ARGS_LOG"], "w") as out:\\n'
           + '        json.dump([sys.argv[1:], os.environ.get("CODEX_HOME")], out)\\n')
layout = os.environ.get('TEST_LAYOUT', 'bin').strip('/')
prefix = layout + '/' if layout else ''
files = {'codex': (program.encode(), 0o755),
         'resources/bundle.txt': (b'resource', 0o644)}
for name in ('codex-code-mode-host', 'codex-exec-server', 'codex-app-server'):
    if name == 'codex-code-mode-host' and os.environ.get('OMIT_HOST'):
        continue
    broken = (name == 'codex-code-mode-host' and os.environ.get('BROKEN_HOST')) or (name == 'codex-exec-server' and os.environ.get('BROKEN_EXEC'))
    body = b'#!/bin/sh\\nexit 126\\n' if broken else b'#!/bin/sh\\n[ "$1" = "--help" ]\\n'
    mode = 0o644 if name == 'codex-code-mode-host' and os.environ.get('NONEXEC_HOST') else 0o755
    files[name] = (body, mode)
if os.environ.get('OMIT_CODEX'):
    del files['codex']
if os.environ.get('DUPLICATE_CODEX'):
    files['other/codex'] = files['codex']
with tarfile.open(output, 'w:gz') as archive:
    for name, (payload, mode) in files.items():
        entry = tarfile.TarInfo(prefix + name)
        entry.mode, entry.size = mode, len(payload)
        archive.addfile(entry, io.BytesIO(payload))
''')
        self.recipe = ROOT / "docker/agents/codex.sh"
        self.dispatcher = ROOT / "docker/mma-agents"
        self.launcher = self.home / ".local/bin/codex"
        self.binary = self.home / ".local/share/mma-agents/codex/codex"

    def stub(self, name, body):
        path = self.stubs / name
        path.write_text(body)
        path.chmod(0o755)

    def run_script(self, script, *args, expected=0, **overrides):
        result = subprocess.run(["bash", str(script), *args],
                                env=dict(self.env, **overrides), text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if expected == 0:
            self.assertEqual(result.returncode, 0, result.stdout)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def launch(self, *args, **overrides):
        result = subprocess.run([str(self.launcher), *args],
                                env=dict(self.env, **overrides), text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertEqual(result.returncode, 0, result.stdout)
        return result

    def test_current_install_and_offline_repeat(self):
        self.run_script(self.recipe)
        self.assertEqual(self.launch("--version").stdout.strip(), "codex-cli 1.2.3")
        self.run_script(self.recipe, FAIL_CURL="1")
        self.run_script(self.recipe, "latest", FAIL_CURL="1")
        fetches = Path(self.env["FETCH_LOG"]).read_text().splitlines()
        self.assertEqual(len(fetches), 1)
        self.assertIn("/latest/download/codex-package-x86_64-unknown-linux-musl.tar.gz", fetches[0])

    def test_pin_unpin_argument_forwarding_and_existing_config(self):
        config = self.home / ".mma-agents/codex/config.toml"
        config.parent.mkdir(parents=True)
        config.write_text('model = "chosen-model"\napproval_policy = "on-request"\n')
        original = config.read_bytes()
        self.run_script(self.recipe, "1.2.3")
        self.assertIn("/download/rust-v1.2.3/", Path(self.env["FETCH_LOG"]).read_text())
        self.run_script(self.recipe, "1.2.3", FAIL_CURL="1")
        args = ("a prompt with spaces", "--", "$literal")
        external_state = str(self.work / "persistent state")
        self.launch(*args, CODEX_HOME=external_state)
        observed, state = json.loads(Path(self.env["ARGS_LOG"]).read_text())
        self.assertEqual(observed, ["--no-daemon", "-c", "check_for_update_on_startup=false", *args])
        self.assertEqual(state, external_state)
        self.run_script(self.recipe, RELEASE_VERSION="1.2.4")
        self.launch(*args)
        observed, state = json.loads(Path(self.env["ARGS_LOG"]).read_text())
        self.assertEqual(observed, ["--no-daemon", *args])
        self.assertEqual(state, str(config.parent))
        self.assertEqual(config.read_bytes(), original)

    def test_failed_upgrade_preserves_working_install(self):
        config = self.home / ".mma-agents/codex/config.toml"
        config.parent.mkdir(parents=True)
        config.write_text('model = "chosen-model"\n')
        self.run_script(self.recipe, "1.2.3")
        originals = {p: p.read_bytes() for p in (self.binary, self.launcher,
                     self.binary.parent / "request", config)}
        original_target = self.binary.resolve()
        original_packages = set((self.binary.parent / "packages").iterdir())
        for overrides in ({"FAIL_CURL": "1"}, {"CORRUPT_ARCHIVE": "1"},
                          {"RELEASE_VERSION": "9.9.9"}, {"RELEASE_VERSION": "invalid"},
                          {"OMIT_HOST": "1"}, {"BROKEN_HOST": "1"},
                          {"NONEXEC_HOST": "1"}, {"BROKEN_EXEC": "1"},
                          {"NO_NO_DAEMON": "1"}, {"OMIT_CODEX": "1"},
                          {"DUPLICATE_CODEX": "1"}):
            with self.subTest(overrides=overrides):
                self.run_script(self.recipe, "2.0.0", expected=1,
                                **dict({"RELEASE_VERSION": "2.0.0"}, **overrides))
                for path, contents in originals.items():
                    self.assertEqual(path.read_bytes(), contents)
                self.assertEqual(list(self.binary.parent.glob(".install.*")), [])
                self.assertEqual(self.binary.resolve(), original_target)
                self.assertEqual(set((self.binary.parent / "packages").iterdir()), original_packages)
        self.assertEqual(self.launch("--version").stdout.strip(), "codex-cli 1.2.3")

    def test_supported_architectures(self):
        for arch, target in (("x86_64", "x86_64"), ("amd64", "x86_64"),
                             ("aarch64", "aarch64"), ("arm64", "aarch64")):
            with self.subTest(arch=arch):
                if self.binary.exists():
                    self.binary.unlink()
                self.run_script(self.recipe, TEST_ARCH=arch)
                url = Path(self.env["FETCH_LOG"]).read_text().splitlines()[-1]
                self.assertTrue(url.endswith(f"codex-package-{target}-unknown-linux-musl.tar.gz"))

    def test_package_layout_and_resources_survive_installation(self):
        for layout in ("", "bin", "codex-package/bin"):
            with self.subTest(layout=layout):
                if self.binary.exists():
                    self.binary.unlink()
                self.run_script(self.recipe, TEST_LAYOUT=layout)
                directory = self.binary.resolve().parent
                for helper in ("codex-code-mode-host", "codex-exec-server", "codex-app-server"):
                    self.assertTrue(os.access(directory / helper, os.X_OK))
                self.assertEqual((directory / "resources/bundle.txt").read_text(), "resource")
                self.launch("--help")

    def test_legacy_binary_only_install_is_repaired_at_same_version(self):
        self.run_script(self.recipe)
        original = self.binary.read_bytes()
        self.binary.unlink()
        self.binary.write_bytes(original)
        self.binary.chmod(0o755)
        # This still passes --version, reproducing the old false health check.
        self.assertEqual(self.launch("--version").stdout.strip(), "codex-cli 1.2.3")
        self.run_script(self.recipe)
        self.assertTrue(self.binary.is_symlink())
        self.assertTrue((self.binary.resolve().parent / "codex-code-mode-host").exists())
        url = Path(self.env["FETCH_LOG"]).read_text().splitlines()[-1]
        self.assertIn("/download/rust-v1.2.3/", url)
        self.run_script(self.recipe, FAIL_CURL="1")

    def test_missing_or_broken_installed_host_triggers_repair(self):
        self.run_script(self.recipe)
        for damage in ("missing", "nonexecutable", "broken"):
            with self.subTest(damage=damage):
                helper = self.binary.resolve().parent / "codex-code-mode-host"
                if damage == "missing":
                    helper.unlink()
                elif damage == "nonexecutable":
                    helper.chmod(0o644)
                else:
                    helper.write_text("#!/bin/sh\nexit 126\n")
                old_package = self.binary.resolve().parent
                self.run_script(self.recipe)
                self.assertNotEqual(self.binary.resolve().parent, old_package)
                self.run_script(self.recipe, FAIL_CURL="1")

    def test_failed_repair_preserves_legacy_executable(self):
        self.binary.parent.mkdir(parents=True)
        self.binary.write_text("#!/bin/sh\necho 'codex-cli 1.2.3'\n")
        self.binary.chmod(0o755)
        (self.binary.parent / "request").write_text("current\n")
        original = self.binary.read_bytes()
        self.run_script(self.recipe, expected=1, FAIL_CURL="1")
        self.assertFalse(self.binary.is_symlink())
        self.assertEqual(self.binary.read_bytes(), original)

    def test_successful_upgrade_keeps_previous_package_for_running_sessions(self):
        self.run_script(self.recipe, "1.2.3")
        previous = self.binary.resolve()
        self.run_script(self.recipe, "1.2.4", RELEASE_VERSION="1.2.4")
        self.assertNotEqual(self.binary.resolve(), previous)
        self.assertTrue(previous.exists())
        self.assertTrue((previous.parent / "codex-code-mode-host").exists())
        self.assertEqual(self.launch("--version").stdout.strip(), "codex-cli 1.2.4")

    def test_unsupported_platform_and_invalid_versions_do_not_download(self):
        for overrides, request in (({"TEST_ARCH": "riscv64"}, ""),
                                   ({"TEST_SYSTEM": "Darwin"}, ""),
                                   ({}, "stable"), ({}, "1.2"), ({}, "../bad")):
            with self.subTest(overrides=overrides, request=request):
                self.run_script(self.recipe, request, expected=1, **overrides)
        self.assertFalse(Path(self.env["FETCH_LOG"]).exists())

    def test_dispatcher_preflight_defaults_and_listing(self):
        self.run_script(self.dispatcher, "install", "")
        self.run_script(self.dispatcher, "install", "none")
        self.assertFalse((self.home / ".mma-agents").exists())
        self.assertEqual(self.run_script(self.dispatcher, "list").stdout.split(),
                         ["claude", "codex"])
        for selection in ("codex,unknown", "codex,codex", "codex@", "../codex"):
            with self.subTest(selection=selection):
                self.run_script(self.dispatcher, "install", selection, expected=1)
        self.assertFalse(Path(self.env["FETCH_LOG"]).exists())
        self.assertFalse((self.home / ".mma-agents").exists())

    def test_dispatcher_runs_crlf_recipe(self):
        recipes = self.work / "recipes"
        recipes.mkdir()
        (recipes / "sample.sh").write_bytes(b'#!/bin/bash\r\nset -euo pipefail\r\nprintf "%s" "$1" > "$HOME/recipe-result"\r\n')
        self.run_script(self.dispatcher, "install", "sample@1.2.3",
                        MMA_AGENTS_RECIPES=str(recipes))
        self.assertEqual((self.home / "recipe-result").read_text(), "1.2.3")

    def test_entrypoint_propagates_requested_install_failure(self):
        self.stub("mpirun", "#!/bin/sh\necho node\necho node\n")
        self.stub("mma-agents", "#!/bin/sh\nexit 1\n")
        marker = self.work / "command-ran"
        result = self.run_script(ROOT / "docker/entrypoint.sh", "touch", str(marker),
                                 expected=1, MMA_AUTO_BUILD="0", AI_AGENTS="codex")
        self.assertIn("Requested AI agent installation failed", result.stdout)
        self.assertFalse(marker.exists())


if __name__ == "__main__":
    unittest.main()
