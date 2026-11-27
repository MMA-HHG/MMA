"""The structural checks of MMA's AI support, run as a pytest test.

Like the AI support itself, these tests are optional: "pytest ai -v" runs them together
with the installer tests, while MMA's default test run covers the core modules only.
CI runs ai/check_ai_layer.py directly. Conventions: docs/ai/README.md, section 11.
"""

import sys
from pathlib import Path

import pytest

pytest.importorskip("yaml", reason="the checks of the AI support need PyYAML")

AI_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(AI_DIR))
import check_ai_layer  # noqa: E402  (after the PyYAML check, which it needs)


def test_structural_checks_pass(capsys):
    status = check_ai_layer.main(["--root", str(AI_DIR.parent)])
    report = capsys.readouterr().out
    assert status == 0, report
