# VISTAS-HCV-python

## Batch FASTA conversion

Run Step 2 CSV conversion for every FASTA file directly under a folder:

```sh
python scan_hcv_vistas_folder.py path/to/folder
```

For each `*.fa`, `*.fas`, `*.fasta`, or `*.fna` file, the script writes a CSV next to it with the same base name, for example `sample.fasta` becomes `sample.csv`.

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
