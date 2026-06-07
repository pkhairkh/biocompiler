# Backward compatibility shim - use biocompiler.optimizer.multigene instead
from biocompiler.optimizer.multigene import *  # noqa: F401,F403
from biocompiler.optimizer.multigene import (  # noqa: F401
    _assemble_construct,
    _infer_construct_type,
    _protein_to_2a_dna,
    _generate_genbank_multigene,
)
import warnings
warnings.warn(
    "Import from biocompiler.multigene is deprecated. "
    "Use biocompiler.optimizer.multigene instead.",
    DeprecationWarning,
    stacklevel=2,
)
