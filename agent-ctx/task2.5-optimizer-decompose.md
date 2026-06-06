# Task 2.5 — Decompose optimization.py into focused modules

## Summary
Successfully decomposed the monolithic 10,986-line optimization.py into a modular subpackage (biocompiler.optimizer) with 5 focused modules, maintaining full backward compatibility.

## Files Created
- `src/biocompiler/optimizer/__init__.py` — Re-exports all public symbols
- `src/biocompiler/optimizer/__init__.pyi` — Type stub for package
- `src/biocompiler/optimizer/cai.py` — CAI computation (~440 lines)
- `src/biocompiler/optimizer/cai.pyi` — Type stub
- `src/biocompiler/optimizer/constraints.py` — Constraint handling (~460 lines)
- `src/biocompiler/optimizer/greedy.py` — Greedy optimizer (~3,300 lines)
- `src/biocompiler/optimizer/greedy.pyi` — Type stub
- `src/biocompiler/optimizer/pipeline.py` — BioOptimizer + high-level API (~7,100 lines)
- `src/biocompiler/optimizer/pipeline.pyi` — Type stub
- `src/biocompiler/optimizer/utils.py` — Data classes and utilities (~170 lines)
- `src/biocompiler/optimizer/py.typed` — PEP 561 marker

## Files Modified
- `src/biocompiler/optimization.py` — Replaced with thin shim + deprecation warning
- `src/biocompiler/__init__.py` — Changed `.optimization` → `.optimizer` import

## Test Results
- 283/283 optimization-related tests PASSED
- 383/383 broader tests PASSED (1 pre-existing failure)
- Zero new failures introduced

## Backward Compatibility
- `from biocompiler.optimization import BioOptimizer` still works (deprecation warning)
- `from biocompiler.optimizer import BioOptimizer` works (new canonical import)
- `from biocompiler import BioOptimizer` works (top-level, no warning)
