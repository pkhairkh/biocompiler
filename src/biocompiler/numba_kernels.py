# Backward compatibility shim - use biocompiler.optimizer.numba_kernels instead
from biocompiler.optimizer.numba_kernels import *  # noqa: F401,F403
import warnings
warnings.warn(
    "Import from biocompiler.numba_kernels is deprecated. "
    "Use biocompiler.optimizer.numba_kernels instead.",
    DeprecationWarning,
    stacklevel=2,
)
