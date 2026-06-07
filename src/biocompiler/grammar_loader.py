# Backward compatibility shim - use biocompiler.optimizer.grammar_loader instead
from biocompiler.optimizer.grammar_loader import *  # noqa: F401,F403
from biocompiler.optimizer.grammar_loader import (  # noqa: F401
    _DEFAULT_CAI_THRESHOLD,
    _DEFAULT_CELLULAR_CONTEXT,
    _DEFAULT_CRYPTIC_SPLICE_THRESHOLD,
    _DEFAULT_ENZYMES,
    _DEFAULT_EXON_BOUNDARIES,
    _DEFAULT_GC_HI,
    _DEFAULT_GC_LO,
    _DEFAULT_ORGANISM,
    _DEFAULT_UNCERTAIN_LO,
)
# Re-export 'yaml' attribute for monkeypatching in tests
from biocompiler.optimizer import grammar_loader as _gl
yaml = _gl.yaml
_check_yaml_available = _gl._check_yaml_available
import warnings
warnings.warn(
    "Import from biocompiler.grammar_loader is deprecated. "
    "Use biocompiler.optimizer.grammar_loader instead.",
    DeprecationWarning,
    stacklevel=2,
)
