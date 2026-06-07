# Backward compatibility shim - use biocompiler.optimizer.incremental instead
from biocompiler.optimizer.incremental import *  # noqa: F401,F403
import warnings
warnings.warn(
    "Import from biocompiler.incremental is deprecated. "
    "Use biocompiler.optimizer.incremental instead.",
    DeprecationWarning,
    stacklevel=2,
)
