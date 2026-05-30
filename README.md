# VISTAS-HCV-python

## Batch FASTA conversion

Run Step 2 CSV conversion for every FASTA file directly under a folder:

```sh
python scan_hcv_vistas_folder.py path/to/folder
```

For each `*.fa`, `*.fas`, `*.fasta`, or `*.fna` file, the script writes a CSV under a `RefID` subfolder and moves the FASTA there too. The `RefID` is the first part of the FASTA filename before `_`, for example `ABC123_pair.fasta` becomes `ABC123/ABC123_pair.fasta` and `ABC123/ABC123_pair.csv`. If `ABC123_pair.csv` already exists next to the FASTA file, the script moves both files into `ABC123/` and skips rerunning conversion.

Options supported by `run_hcv_vistas.py` can be appended after the folder:

```sh
python scan_hcv_vistas_folder.py path/to/folder --country USA --year 2024
```

## Pre-commit Checks

This repository uses `.githooks/pre-commit` to block commits that contain absolute filesystem paths in staged text files.

Enable it for a fresh clone with:

```sh
git config core.hooksPath .githooks
```
