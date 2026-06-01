"""Bootstrap runner for comprehensive_benchmark from the src package."""
import importlib
import sys
from pathlib import Path

# Ensure src/ package is found first (before root biocompiler/)
src_dir = str(Path(__file__).resolve().parent.parent / "src")
sys.path.insert(0, src_dir)

# Remove the root biocompiler from cached modules so we get the src version
if 'biocompiler' in sys.modules:
    # Remove all biocompiler submodules
    to_remove = [k for k in sys.modules if k.startswith('biocompiler')]
    for k in to_remove:
        del sys.modules[k]

# Now import fresh from the src package
from biocompiler.comprehensive_benchmark import main

if __name__ == "__main__":
    main()
