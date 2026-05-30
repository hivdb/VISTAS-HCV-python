#!/usr/bin/env python3
"""Run run_hcv_vistas.py for each FASTA file directly under a folder."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
CONVERTER = SCRIPT_DIR / "run_hcv_vistas.py"
FASTA_EXTENSIONS = {".fa", ".fas", ".fasta", ".fna"}


def find_fasta_files(folder: Path) -> list[Path]:
    return sorted(
        path
        for path in folder.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    )


def ref_id_from_filename(path: Path) -> str:
    return path.stem.split("_", 1)[0]


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Scan a folder for FASTA files and create one Step 2 CSV per file. "
            "Extra options after the folder are passed to run_hcv_vistas.py."
        )
    )
    parser.add_argument("folder", type=Path, help="Folder containing FASTA files")
    return parser


def main() -> int:
    parser = get_parser()
    args, converter_args = parser.parse_known_args()

    if not args.folder.is_dir():
        parser.exit(1, f"error: folder not found: {args.folder}\n")
    if not CONVERTER.is_file():
        parser.exit(1, f"error: converter not found: {CONVERTER}\n")

    fasta_files = find_fasta_files(args.folder)
    if not fasta_files:
        parser.exit(1, f"error: no FASTA files found directly under {args.folder}\n")

    failures = 0
    for fasta_path in fasta_files:
        output_dir = args.folder / ref_id_from_filename(fasta_path)
        output_dir.mkdir(exist_ok=True)
        output_path = output_dir / fasta_path.with_suffix(".csv").name
        moved_fasta_path = output_dir / fasta_path.name
        existing_csv_path = fasta_path.with_suffix(".csv")
        if existing_csv_path.is_file():
            existing_csv_path.replace(output_path)
            fasta_path.replace(moved_fasta_path)
            print(f"Moved existing CSV {existing_csv_path} -> {output_path}", file=sys.stderr, flush=True)
            print(f"Moved FASTA {fasta_path} -> {moved_fasta_path}", file=sys.stderr, flush=True)
            continue

        command = [
            sys.executable,
            str(CONVERTER),
            str(fasta_path),
            "-o",
            str(output_path),
            *converter_args,
        ]
        print(f"Processing {fasta_path} -> {output_path}", file=sys.stderr, flush=True)
        completed = subprocess.run(command)
        if completed.returncode != 0:
            failures += 1
            print(f"Failed: {fasta_path}", file=sys.stderr, flush=True)
        else:
            fasta_path.replace(moved_fasta_path)
            print(f"Moved FASTA {fasta_path} -> {moved_fasta_path}", file=sys.stderr, flush=True)

    if failures:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
