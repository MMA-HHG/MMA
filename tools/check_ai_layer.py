#!/usr/bin/env python3
"""Structural checks for MMA's AI support (docs/ai/README.md, section 11).

Usage:
    python3 tools/check_ai_layer.py [--with-readme] [--root PATH]

Checks S1-S8 are listed in
.agents/skills/mma-documentation-consistency/references/verification-summary.md.
Exit status 1 if any error is found; warnings do not fail the run.
With --with-readme, README.md is checked too (links, anchors and skill references) and its
findings are reported as warnings.

Needs Python 3.8 or later and PyYAML.
"""

import argparse
import datetime
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # keep the checkout free of __pycache__
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ai_signposts as signposts  # noqa: E402  (also reports a missing PyYAML)

NAME_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
ENTITY_RE = re.compile(r"^mma-[a-z0-9]+(?:-[a-z0-9]+)*$")
STATUSES = {"maintained", "community", "experimental", "generated"}
TRACKS = {"operation", "applications", "core", "documentation"}
WORKFLOW_KEYS = ["description", "mma-status", "mma-track", "mma-owner"]
WORKFLOW_HEADINGS = ["Purpose", "Entry conditions", "Route",
                     "Decision points and transitions", "Exit criteria"]
AGENTS_HEADINGS = ["Collaboration and autonomy", "Environment", "Repository map", "Routing",
                   "Repository rules", "Findings", "Extending the AI support"]
RECORD_KEYS = ["machine", "hpc-entry", "record-status", "last-verified", "mma-commit", "verified-by"]
RECORD_STATUSES = {"unverified", "verified"}
PLACEHOLDER_OWNER = "@REPLACE-WITH-OWNER"
AGENTS_MAX_LINES = 150
WORKFLOW_MAX_LINES = 80
SKILL_MAX_LINES = 500
DESCRIPTION_MAX = 1024
RECORD_STALE_DAYS = 365
HANDBOOK = "docs/ai/README.md"

CODE_SPAN = re.compile(r"`([^`\n]+)`")
LINK = re.compile(r"!?\[[^\]\n]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
FENCE = re.compile(r"^\s*(```|~~~)")
HEADING = re.compile(r"^(#{1,6})\s+(.+?)\s*#*\s*$")
NOT_A_PATH = set("<>*$~{}@=|…")
# A skill reference: "the skill `mma-io`" or "skills `mma-a`, `mma-b` and `mma-c`".
SKILL_REFERENCE = re.compile(r"\bskills?\s+((?:`[^`\n]+`(?:\s*(?:,|and|or)\s*)?)+)", re.IGNORECASE)


class Report:
    def __init__(self):
        self.errors, self.warnings = [], []

    def add(self, as_warning, check, where, message):
        line = f"[{check}] {where}: {message}"
        (self.warnings if as_warning else self.errors).append(line)

    def error(self, check, where, message):
        self.add(False, check, where, message)

    def warn(self, check, where, message):
        self.add(True, check, where, message)


def read(path):
    return path.read_text(encoding="utf-8")


def rel(root, path):
    return path.relative_to(root).as_posix()


def outside_fences(text):
    """Yield (line number, line) for lines that are not inside fenced code blocks."""
    in_fence = False
    for number, line in enumerate(text.splitlines(), 1):
        if FENCE.match(line):
            in_fence = not in_fence
            continue
        if not in_fence:
            yield number, line


def headings(text, level=None):
    found = []
    for _, line in outside_fences(text):
        match = HEADING.match(line)
        if match and (level is None or len(match.group(1)) == level):
            found.append(match.group(2).strip())
    return found


def slugify(heading):
    """GitHub-style anchor of a heading."""
    text = heading.strip().lower()
    text = re.sub(r"[^\w\- ]", "", text)
    return text.replace(" ", "-")


_anchor_cache = {}


def anchors(path):
    if path not in _anchor_cache:
        seen, result = {}, set()
        for heading in headings(read(path)):
            slug = slugify(heading)
            count = seen.get(slug, 0)
            result.add(slug if count == 0 else f"{slug}-{count}")
            seen[slug] = count + 1
        _anchor_cache[path] = result
    return _anchor_cache[path]


def parse_frontmatter(root, path, report, check):
    try:
        frontmatter, body = signposts.split_frontmatter(read(path))
    except Exception as exc:  # YAML errors come in several types
        report.error(check, rel(root, path), f"invalid front matter: {exc}")
        return None, ""
    if frontmatter is None:
        report.error(check, rel(root, path), "missing YAML front matter")
    return frontmatter, body


def check_tracks(root, path, value, report, check):
    tracks = [t.strip() for t in str(value or "").split(",") if t.strip()]
    if not tracks:
        report.error(check, rel(root, path), "mma-track is empty")
    for track in tracks:
        if track not in TRACKS:
            report.error(check, rel(root, path),
                         f"unknown track '{track}'; use one of {', '.join(sorted(TRACKS))}")


# S1 -----------------------------------------------------------------------------------
def check_skills(root, report):
    names, maintained, handwritten = set(), set(), set()
    base = root / ".agents" / "skills"
    if not base.is_dir():
        report.error("S1", ".agents/skills", "directory missing")
        return names, maintained, handwritten
    for folder in sorted(p for p in base.iterdir() if p.is_dir()):
        skill_md = folder / "SKILL.md"
        where = rel(root, skill_md)
        if not skill_md.is_file():
            report.error("S1", rel(root, folder), "skill folder without SKILL.md")
            continue
        frontmatter, body = parse_frontmatter(root, skill_md, report, "S1")
        if frontmatter is None:
            continue
        names.add(folder.name)
        name, description = frontmatter.get("name"), frontmatter.get("description")
        metadata = frontmatter.get("metadata") or {}
        if name != folder.name:
            report.error("S1", where, f"name '{name}' differs from the folder name '{folder.name}'")
        if not isinstance(name, str) or not NAME_RE.match(name) or len(name) > 64:
            report.error("S1", where, "name must use lowercase letters, digits and single hyphens, "
                                      "at most 64 characters")
        elif not name.startswith("mma-"):
            report.error("S1", where, "repository skills must use the reserved prefix 'mma-'")
        if not isinstance(description, str) or not description.strip():
            report.error("S1", where, "description is missing")
        elif len(description) > DESCRIPTION_MAX:
            report.error("S1", where, f"description has {len(description)} characters; "
                                      f"the limit is {DESCRIPTION_MAX}")
        if not isinstance(metadata, dict):
            report.error("S1", where, "metadata must be a mapping")
            continue
        for key, value in metadata.items():
            if not isinstance(value, str):
                report.error("S1", where, f"metadata value of '{key}' must be a string")
        status = metadata.get("mma-status")
        if status not in STATUSES:
            report.error("S1", where, f"metadata.mma-status must be one of {', '.join(sorted(STATUSES))}")
        if "mma-track" in metadata:
            check_tracks(root, skill_md, metadata.get("mma-track"), report, "S1")
        if status == "maintained":
            maintained.add(folder.name)
        if status != "generated":
            handwritten.add(folder.name)
        if len(body.splitlines()) > SKILL_MAX_LINES:
            report.warn("S1", where, f"body is longer than {SKILL_MAX_LINES} lines; "
                                     "move material to references/")
    return names, maintained, handwritten


def check_name_clashes(root, report, handwritten, workflow_names):
    """A hand-written skill must not take the name of a workflow's signpost skill."""
    for workflow in sorted(workflow_names):
        name = f"mma-{workflow}"
        if name in handwritten:
            report.error("S1", f".agents/skills/{name}/SKILL.md",
                         f"hand-written skill has the name of the signpost skill of workflow "
                         f"'{workflow}'; rename one of them")


# S3 -----------------------------------------------------------------------------------
def check_agents(root, report, maintained_skills):
    path = root / "AGENTS.md"
    if not path.is_file():
        report.error("S3", "AGENTS.md", "file missing")
        return ""
    text = read(path)
    if len(text.splitlines()) > AGENTS_MAX_LINES:
        report.warn("S3", "AGENTS.md", f"longer than {AGENTS_MAX_LINES} lines; "
                                       "move detail into workflows, skills or the handbook")
    present = set(headings(text, level=2))
    for heading in AGENTS_HEADINGS:
        if heading not in present:
            report.error("S3", "AGENTS.md", f"required section '## {heading}' missing")
    for skill in sorted(maintained_skills):
        if skill not in text:
            report.error("S3", "AGENTS.md", f"maintained skill '{skill}' is not routed")
    return text


# S2 -----------------------------------------------------------------------------------
def check_workflows(root, report, agents_text):
    names = set()
    base = root / "docs" / "ai" / "workflows"
    if not base.is_dir():
        report.error("S2", "docs/ai/workflows", "directory missing")
        return names
    for path in sorted(base.glob("*.md")):
        where = rel(root, path)
        names.add(path.stem)
        if not NAME_RE.match(path.stem):
            report.error("S2", where, "file name must use lowercase letters, digits and single hyphens")
        frontmatter, body = parse_frontmatter(root, path, report, "S2")
        if frontmatter is None:
            continue
        for key in WORKFLOW_KEYS:
            if not frontmatter.get(key):
                report.error("S2", where, f"front matter key '{key}' missing")
        if frontmatter.get("mma-status") and frontmatter["mma-status"] not in STATUSES:
            report.error("S2", where, f"mma-status must be one of {', '.join(sorted(STATUSES))}")
        if frontmatter.get("mma-track"):
            check_tracks(root, path, frontmatter["mma-track"], report, "S2")
        if len(str(frontmatter.get("description", ""))) > DESCRIPTION_MAX:
            report.error("S2", where, f"description longer than {DESCRIPTION_MAX} characters")
        present = set(headings(body, level=2))
        for heading in WORKFLOW_HEADINGS:
            if heading not in present:
                report.error("S2", where, f"required section '## {heading}' missing")
        if len(read(path).splitlines()) > WORKFLOW_MAX_LINES:
            report.warn("S2", where, f"longer than {WORKFLOW_MAX_LINES} lines; a workflow is a route, "
                                     "not a manual")
        if frontmatter.get("mma-status") == "maintained" and \
                f"docs/ai/workflows/{path.name}" not in agents_text:
            report.error("S2", where, "maintained workflow is not routed from AGENTS.md")
    return names


# S7 -----------------------------------------------------------------------------------
def check_codeowners(root, report):
    path = root / ".github" / "CODEOWNERS"
    if not path.is_file():
        report.error("S7", ".github/CODEOWNERS", "file missing")
        return
    text = read(path)
    if signposts.CORE_BEGIN not in text or signposts.CORE_END not in text:
        report.error("S7", ".github/CODEOWNERS", "core block markers missing")
        return
    paths = signposts.core_paths(root)
    if not paths:
        report.error("S7", ".github/CODEOWNERS", "core block is empty")
    for entry in paths:
        if not (root / entry.strip("/")).exists():
            report.error("S7", ".github/CODEOWNERS", f"core path '{entry}' does not exist")
    if PLACEHOLDER_OWNER in text:
        report.warn("S7", ".github/CODEOWNERS", f"replace the placeholder owner {PLACEHOLDER_OWNER}")


# S8 -----------------------------------------------------------------------------------
def check_records(root, report):
    base = root / "docs" / "ai" / "records"
    if not base.is_dir():
        return
    today = datetime.date.today()
    for path in sorted(base.rglob("*.md")):
        if path.name == "README.md":
            continue
        where = rel(root, path)
        frontmatter, _ = parse_frontmatter(root, path, report, "S8")
        if frontmatter is None:
            continue
        for key in RECORD_KEYS:
            if key not in frontmatter:
                report.error("S8", where, f"front matter key '{key}' missing")
        if frontmatter.get("record-status") not in RECORD_STATUSES:
            report.error("S8", where, "record-status must be 'unverified' or 'verified'")
        try:
            verified = datetime.date.fromisoformat(str(frontmatter.get("last-verified")))
            if (today - verified).days > RECORD_STALE_DAYS:
                report.warn("S8", where, f"last verified more than {RECORD_STALE_DAYS} days ago")
        except ValueError:
            report.error("S8", where, "last-verified must be a date in the form YYYY-MM-DD")


# S4, S5 -------------------------------------------------------------------------------
def planned_names(root):
    """Entity names reserved in the handbook's section 'Planned entities'."""
    path = root / HANDBOOK
    names, inside = set(), False
    if not path.is_file():
        return names
    for _, line in outside_fences(read(path)):
        match = HEADING.match(line)
        if match and len(match.group(1)) == 2:
            inside = re.sub(r"^\d+\.\s*", "", match.group(2)).strip() == "Planned entities"
            continue
        if inside:
            names.update(s for s in CODE_SPAN.findall(line) if ENTITY_RE.match(s))
    return names


def files_to_scan(root, with_readme):
    files = [root / "AGENTS.md", root / "CONTRIBUTING.md", root / "docs" / "deviations-from-paper.md",
             root / ".github" / "pull_request_template.md", root / "docker" / "agents" / "README.md"]
    files += sorted((root / "docs" / "ai").rglob("*.md"))
    files += sorted((root / ".agents" / "skills").rglob("*.md"))
    if with_readme:
        files.append(root / "README.md")
    local = root / "docs" / "ai" / "local"
    return [f for f in files if f.is_file() and local not in f.parents]


def looks_like_path(span):
    """A relative repository path; absolute paths, such as /MMA in the container, are skipped."""
    return ("/" in span and not span.startswith("/") and not any(c.isspace() for c in span)
            and not any(c in NOT_A_PATH for c in span)
            and "://" not in span and not span.startswith(("-", "//")))


def check_text(root, report, files, known_names, planned):
    readme = root / "README.md"
    for path in files:
        where_file = rel(root, path)
        as_warning = path == readme
        in_planned = False
        for number, line in outside_fences(read(path)):
            where = f"{where_file}:{number}"
            match = HEADING.match(line)
            if match and len(match.group(1)) == 2:
                in_planned = (where_file == HANDBOOK and
                              re.sub(r"^\d+\.\s*", "", match.group(2)).strip() == "Planned entities")
            for target in LINK.findall(line):
                check_link(root, report, path, where, target, as_warning)
            for group in SKILL_REFERENCE.findall(line):
                for span in CODE_SPAN.findall(group):
                    if ENTITY_RE.match(span) and span not in known_names | planned:
                        report.add(as_warning, "S5", where, f"unknown skill '{span}'")
            for span in (s.strip() for s in CODE_SPAN.findall(line)):
                if in_planned or as_warning or not looks_like_path(span):
                    continue
                candidates = [path.parent / span, root / span.lstrip("/")]
                if not any(c.exists() for c in candidates):
                    report.add(as_warning, "S4", where, f"path '{span}' not found")


def check_link(root, report, path, where, target, as_warning):
    if target.startswith(("http://", "https://", "mailto:", "ftp://")) or "<" in target or "\\" in target:
        return
    file_part, _, fragment = target.partition("#")
    destination = path if not file_part else (path.parent / file_part).resolve()
    if file_part and not destination.exists():
        report.add(as_warning, "S4", where, f"link target '{target}' not found")
        return
    if fragment and destination.is_file() and destination.suffix.lower() == ".md":
        if fragment.lower() not in anchors(destination):
            report.add(as_warning, "S4", where, f"anchor '#{fragment}' not found in "
                                                f"{rel(root, destination)}")


# S6 -----------------------------------------------------------------------------------
def check_signposts(root, report):
    for problem in signposts.check(root):
        report.error("S6", "signposts", problem)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Structural checks for MMA's AI support.")
    parser.add_argument("--with-readme", action="store_true",
                        help="also check README.md; its findings are reported as warnings")
    parser.add_argument("--root", type=Path, default=signposts.repo_root(),
                        help="repository root (default: the parent of tools/)")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    report = Report()

    skill_names, maintained, handwritten = check_skills(root, report)
    agents_text = check_agents(root, report, maintained)
    workflow_names = check_workflows(root, report, agents_text)
    check_name_clashes(root, report, handwritten, workflow_names)
    check_codeowners(root, report)
    check_records(root, report)
    known = skill_names | {f"mma-{w}" for w in workflow_names}
    check_text(root, report, files_to_scan(root, args.with_readme), known, planned_names(root))
    check_signposts(root, report)

    for line in report.errors:
        print(f"ERROR   {line}")
    for line in report.warnings:
        print(f"WARNING {line}")
    print(f"{len(report.errors)} error(s), {len(report.warnings)} warning(s).")
    return 1 if report.errors else 0


if __name__ == "__main__":
    sys.exit(main())
