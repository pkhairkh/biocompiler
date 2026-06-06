# Task 2.3 — BLAST+ Integration for Biosecurity Screening

## Agent: main

## Summary

Implemented BLAST+ integration for biosecurity screening, extending the existing exact/fuzzy substring matching with homology-based search against hazardous sequence databases.

## Files Created

1. **`src/biocompiler/blast_integration.py`** — Main module (~530 lines)
   - `BlastScanner` class: Wrapper for BLAST+ operations (blastn/blastp)
   - `build_hazard_db()`: Build BLAST databases from hazardous sequences
   - `screen_blast()`: Convenience function for BLAST screening
   - `check_biosecurity_blast()`: Combined exact/fuzzy + BLAST screening
   - `is_blast_available()` / `find_blast_bin()`: BLAST+ detection
   - Risk level computation from identity/coverage thresholds
   - Hit classification into hazard categories
   - Graceful fallback when BLAST+ is not installed

2. **`tests/test_blast_integration.py`** — Comprehensive test suite (77 tests)
   - All tests use mocks (no BLAST+ installation required)
   - Covers: binary detection, XML parsing, risk levels, DB building,
     report merging, API integration, edge cases

## Files Modified

3. **`pyproject.toml`** — Added `blast = []` optional dependency group
4. **`src/biocompiler/api.py`** — Added blast_screening, blast_db_path, blast_evalue to ProteinInput; added blast fields to OptimizeResponse; updated /optimize endpoint
5. **`src/biocompiler/__init__.py`** — Added BLAST+ imports with try/except guards

## Test Results

- 77/77 new BLAST integration tests PASSED
- 62/62 existing biosecurity tests PASSED (no regressions)

## Key Design Decisions

- BLAST+ is optional: When not installed, screening falls back to exact/fuzzy only with a warning
- BLAST hits are merged with exact matching results using highest-risk-wins
- Environment variables (BIOCOMPILER_BLAST_PATH, BIOCOMPILER_BLAST_DB_PATH) for custom paths
- Import guards follow the same pattern as Z3/ViennaRNA/MHCflurry
