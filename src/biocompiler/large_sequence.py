# Backward compatibility shim - use biocompiler.optimizer.large_sequence instead
from biocompiler.optimizer.large_sequence import *  # noqa: F401,F403
import warnings
warnings.warn(
    "Import from biocompiler.large_sequence is deprecated. "
    "Use biocompiler.optimizer.large_sequence instead.",
    DeprecationWarning,
    stacklevel=2,
)
