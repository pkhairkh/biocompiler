"""
Backward-compatibility shim — the canonical module is now ``escherichia.py``.

All other organism files use the genus name only (human.py, mouse.py,
yeast.py, bacillus.py, etc.).  This shim re-exports everything so that
existing ``from biocompiler.organisms.e_coli import ...`` statements
continue to work.
"""

from .escherichia import *  # noqa: F401,F403
from .escherichia import __all__  # noqa: F401
