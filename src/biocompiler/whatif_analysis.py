# Backward compatibility shim - use biocompiler.optimizer.whatif_analysis instead
from biocompiler.optimizer.whatif_analysis import *  # noqa: F401,F403
import warnings
warnings.warn(
    "Import from biocompiler.whatif_analysis is deprecated. "
    "Use biocompiler.optimizer.whatif_analysis instead.",
    DeprecationWarning,
    stacklevel=2,
)
