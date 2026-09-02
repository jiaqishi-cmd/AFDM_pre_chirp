"""Install project research skills into the local Codex skills directory.

The repository keeps communication research skills in a shared-reference
layout:

    skills/communication-research-toolkit/<skill>/SKILL.md
    skills/communication-research-toolkit/references/*.md

For Codex global installation, this script creates self-contained skill
folders under ~/.codex/skills, copying the shared references into each skill.
"""

import argparse
import os
import re
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOOLKIT = ROOT / "skills" / "communication-research-toolkit"
REFERENCES = TOOLKIT / "references"
NAME_PATTERN = re.compile(r"^name:\s*([A-Za-z0-9_-]+)\s*$", re.MULTILINE)


def default_codex_skills_dir():
    codex_home = os.environ.get("CODEX_HOME")
    if codex_home:
        return Path(codex_home) / "skills"
    return Path.home() / ".codex" / "skills"


def discover_skills():
    skills = []
    for skill_md in sorted(TOOLKIT.glob("*/SKILL.md")):
        text = skill_md.read_text(encoding="utf-8")
        match = NAME_PATTERN.search(text)
        if not match:
            raise ValueError(f"Missing frontmatter name in {skill_md}")
        skills.append((match.group(1), skill_md.parent))
    return skills


def install_skill(name, source_dir, dest_root, force):
    dest_dir = dest_root / name
    if dest_dir.exists():
        if not force:
            return f"skip existing {dest_dir}"
        shutil.rmtree(dest_dir)

    shutil.copytree(source_dir, dest_dir)

    dest_skill_md = dest_dir / "SKILL.md"
    text = dest_skill_md.read_text(encoding="utf-8")
    text = text.replace("../references/", "references/")
    with dest_skill_md.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)

    shutil.copytree(REFERENCES, dest_dir / "references")
    return f"installed {dest_dir}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dest", type=Path, default=default_codex_skills_dir())
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if not TOOLKIT.exists():
        raise FileNotFoundError(f"Toolkit not found: {TOOLKIT}")
    if not REFERENCES.exists():
        raise FileNotFoundError(f"References not found: {REFERENCES}")

    skills = discover_skills()
    print(f"Destination: {args.dest}")
    for name, source_dir in skills:
        if args.dry_run:
            print(f"would install {name} from {source_dir}")
        else:
            args.dest.mkdir(parents=True, exist_ok=True)
            print(install_skill(name, source_dir, args.dest, args.force))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
