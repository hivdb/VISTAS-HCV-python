#!/usr/bin/env python3
"""Reject staged files that contain absolute filesystem paths."""

from __future__ import annotations

import re
import subprocess
import sys


UNIX_ABSOLUTE_PATH = re.compile(
    r"(?<![:\w])/"
    r"(?:Applications|Library|System|Users|Volumes|bin|etc|home|lib|opt|private|root|sbin|tmp|usr|var)"
    r"(?:/[^\s'\"`<>),;]*)?"
)
WINDOWS_ABSOLUTE_PATH = re.compile(
    r"(?<![\w])(?:[A-Za-z]:[\\/][^\s'\"`<>),;]*|\\\\[^\\/\s]+[\\/][^\\/\s]+[\\/][^\s'\"`<>),;]*)"
)


def git(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", *args],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def staged_files() -> list[str]:
    result = git("diff", "--cached", "--name-only", "--diff-filter=ACMR")
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr, end="")
        sys.exit(result.returncode)
    return [line for line in result.stdout.splitlines() if line]


def staged_text(path: str) -> str | None:
    result = subprocess.run(
        ["git", "show", f":{path}"],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    if result.returncode != 0:
        return None
    try:
        return result.stdout.decode("utf-8")
    except UnicodeDecodeError:
        return None


def findings_for(path: str, text: str) -> list[tuple[int, str]]:
    findings = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        if line_number == 1 and line.startswith("#!"):
            continue
        for pattern in (UNIX_ABSOLUTE_PATH, WINDOWS_ABSOLUTE_PATH):
            for match in pattern.finditer(line):
                findings.append((line_number, match.group(0)))
    return findings


def main() -> int:
    all_findings: list[tuple[str, int, str]] = []
    for path in staged_files():
        text = staged_text(path)
        if text is None:
            continue
        all_findings.extend((path, line_number, value) for line_number, value in findings_for(path, text))

    if not all_findings:
        return 0

    print("Absolute filesystem paths are not allowed in committed files:", file=sys.stderr)
    for path, line_number, value in all_findings:
        print(f"  {path}:{line_number}: {value}", file=sys.stderr)
    print("\nUse relative paths or configurable paths instead.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
