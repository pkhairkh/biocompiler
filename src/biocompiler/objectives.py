# Backward compatibility shim - use biocompiler.optimizer.objectives instead
from biocompiler.optimizer.objectives import *  # noqa: F401,F403
import warnings
warnings.warn(
    "Import from biocompiler.objectives is deprecated. "
    "Use biocompiler.optimizer.objectives instead.",
    DeprecationWarning,
    stacklevel=2,
)
