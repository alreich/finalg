# Changelog

## 0.1.0 — Repackaging

- Restructured the repository into a `src/` layout (`src/finalg/`) so the
  package is installable via `pip install -e .` and behaves correctly once
  built as a wheel.
- Added `pyproject.toml` (setuptools backend), making `finalg` a proper
  installable Python package for the first time.
- Fixed a broken import in `finalg/__init__.py`
  (`from permutation import Perm` → `from finalg.permutation import Perm`)
  that only "worked" when running from an unbuilt source checkout.
- Moved the example-algebra JSON files into `src/finalg/data/algebras/` and
  switched `Examples`/`finalg.__init__` to locate them via
  `importlib.resources` instead of a path computed relative to `__file__`,
  so example algebras are found correctly from an installed wheel, not just
  a source checkout.
- Excluded the Tesseract example files (`tesseract.json`,
  `Tesseract.group.SCRATCHWORK`) and all `F*.json` field-order example
  files (F3/F7/F11/F19/F23 ACN/SQR algebras) from the packaged data set, to
  keep the package small.
- Renamed `test/` to `tests/`.
- Moved `experimental.py` and `work_in_progress.py` out of the installable
  package into a top-level `scratch/` directory.
- Removed `trash/`, `misc/`, and `anaconda_projects/` from the repository;
  added `anaconda_projects/` to `.gitignore`.
- Brought over the remaining Sphinx documentation from `abstract_algebra`:
  `conf.py`, `index.rst`, `api.rst`, `README.rst`, `license.rst`,
  `Makefile`, `make.bat`, `requirements.in`/`requirements.txt`,
  `00_definitions.rst`, `30_vectors.rst`, `40_reg_rep.rst`,
  `50_matrices.rst`, `55_cayley_dickson.rst`, `90_resources.rst`, and
  `_static/`.
- Added `.readthedocs.yaml` at the repo root so Read the Docs can build the
  Sphinx site.
- Added a GitHub Actions CI workflow (`.github/workflows/ci.yml`) that
  installs the package and runs `pytest` across Python 3.10–3.12.
