# Publishing to PyPI

## Prerequisites
- A PyPI account
- A PyPI API token (create at https://pypi.org/manage/account/token/)

## Steps

1. Build the package:
   ```bash
   pip install build
   python -m build
   ```
   This creates `dist/biocompiler-0.9.0-py3-none-any.whl` and `dist/biocompiler-0.9.0.tar.gz`

2. Check the package:
   ```bash
   pip install twine
   twine check dist/*
   ```

3. Upload to TestPyPI (optional, recommended for first time):
   ```bash
   twine upload --repository testpypi dist/*
   ```
   Verify at https://test.pypi.org/project/biocompiler/

4. Upload to PyPI:
   ```bash
   twine upload dist/*
   ```

5. Verify installation:
   ```bash
   pip install biocompiler
   python -c "import biocompiler; print(biocompiler.__version__)"
   ```

## Version Bumping

Before publishing a new version:
1. Update `version` in `pyproject.toml`
2. Update `__version__` and `SAFETY_VERSION` in `src/biocompiler/__init__.py`
3. Update the version in `README.md`
4. Commit and tag: `git tag v0.9.3`
5. Build and upload
