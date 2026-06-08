# Backward compatibility shim - use biocompiler.optimizer.* instead
# Import the module to make all names accessible
import biocompiler.optimizer as _opt

# Map to the correct submodule
import importlib as _il
_mod_map = {
    'optimization_helpers': 'constraint_helpers',
    'optimization_gc': 'gc_adjustment',
    'optimization_cpg': 'cpg_disruption',
    'optimization_splice': 'splice_elimination',
}
_submod_name = _mod_map['optimization_gc']
_submod = _il.import_module(f'biocompiler.optimizer.{_submod_name}')

# Copy all public names
for _name in dir(_submod):
    if not _name.startswith('__'):
        globals()[_name] = getattr(_submod, _name)

import warnings
warnings.warn(
    "Import from biocompiler.optimization_gc is deprecated. "
    f"Use biocompiler.optimizer.{_submod_name} instead.",
    DeprecationWarning,
    stacklevel=2,
)
