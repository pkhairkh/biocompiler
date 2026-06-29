"""
BioCompiler IR — Multi-Level Intermediate Representation
=======================================================

Five IR levels mirroring the Central Dogma:

    IR-L0  GenomicDNA      →  raw DNA + region annotations
    IR-L1  PreMRNA         →  transcribed RNA (T → U)
    IR-L2  MatureMRNA      →  spliced: 5'UTR + CDS + 3'UTR
    IR-L3  Polypeptide     →  translated amino acid sequence
    IR-L4  FoldedProtein   →  3D structure (oracle — Phase 2)

Each level is a dataclass defined in :mod:`biocompiler.ir.types`.
Lowering passes (L0→L1, L1→L2, L2→L3, L3→L4) are pure functions in
:mod:`biocompiler.ir.passes`.  Semantic invariants for each level are
checked by :mod:`biocompiler.ir.invariants`.

Typical usage::

    from biocompiler.ir import (
        IR_L0_GenomicDNA, IRLevel, compile_gene,
    )

    gene = IR_L0_GenomicDNA(
        sequence="ATGGCTAAGTAA",
        regions=[],
        organism="e_coli",
        gene_name="test",
    )
    protein = compile_gene(gene, IRLevel.L3)
    print(protein.sequence)  # → "MAK*"
"""

from .types import (
    IRLevel,
    GeneRegion,
    IR_L0_GenomicDNA,
    IR_L1_PreMRNA,
    IR_L2_MatureMRNA,
    IR_L3_Polypeptide,
    IR_L4_FoldedProtein,
    IRError,
)
from .passes import (
    transcribe,
    splice,
    translate,
    fold,
    compile_gene,
    IRObject,
)
from .invariants import (
    check_l0_invariants,
    check_l1_invariants,
    check_l2_invariants,
    check_l3_invariants,
    check_l4_invariants,
    check_invariants,
)

__all__ = [
    # IR level enum
    "IRLevel",
    # IR dataclasses
    "GeneRegion",
    "IR_L0_GenomicDNA",
    "IR_L1_PreMRNA",
    "IR_L2_MatureMRNA",
    "IR_L3_Polypeptide",
    "IR_L4_FoldedProtein",
    # Error type
    "IRError",
    # Lowering passes
    "transcribe",
    "splice",
    "translate",
    "fold",
    "compile_gene",
    "IRObject",
    # Invariant checkers
    "check_l0_invariants",
    "check_l1_invariants",
    "check_l2_invariants",
    "check_l3_invariants",
    "check_l4_invariants",
    "check_invariants",
]

__version__ = "0.1.0"
