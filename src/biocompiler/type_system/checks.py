"""
BioCompiler Type System — Predicate Check Functions (re-export shim)
====================================================================
Historically this module was a ~3,200-LOC monolith implementing every
``check_*`` predicate. As of the W8-b refactor it is a thin re-export
layer: the actual implementations live in the predicate-category
submodules:

  - sequence_checks    : shared helpers + DNA/codon-level checks
  - mrna_checks        : mRNA structure / stability / folding / cryptic ORF
  - splice_checks      : cryptic splice + GT dinucleotide family
  - cpg_checks         : CpG island detection
  - restriction_checks : restriction-site detection
  - mirna_checks       : miRNA binding-site detection (+ private helpers)
  - immuno_checks      : m6A, PolyA, nucleoside-mod guidance, RQC trigger
  - alu_checks         : Alu repeat detection (W5-a fix preserved)

All historical imports continue to work unchanged, e.g.::

    from biocompiler.type_system.checks import check_no_stop_codons
    from biocompiler.type_system.checks import _count_dinucs_fast
    from biocompiler.type_system.checks import *            # noqa: F401,F403
    from biocompiler.type_system import check_no_stop_codons

This is a pure code-move refactor — no function signatures or logic were
modified. Each submodule declares ``__all__`` (including its underscore
helpers) so that the star imports below transitively re-export every
historical public and private name.
"""

# Star-import from each submodule. The submodules' __all__ explicitly
# include their underscore-prefixed helpers, so this re-exports both
# public and private names — preserving the historical public surface.
from .sequence_checks import *      # noqa: F401,F403
from .mrna_checks import *          # noqa: F401,F403
from .splice_checks import *        # noqa: F401,F403
from .cpg_checks import *           # noqa: F401,F403
from .restriction_checks import *   # noqa: F401,F403
from .mirna_checks import *         # noqa: F401,F403
from .immuno_checks import *        # noqa: F401,F403
from .alu_checks import *           # noqa: F401,F403


__all__ = [
    # ── sequence_checks ────────────────────────────────────────────────
    "_count_dinucs_fast",
    "_is_prokaryotic_organism",
    "_compute_max_gt_count",
    "_translate_dna_to_aa",
    "_resolve_species_cai",
    "_compute_codon_ramp_score",
    "find_cross_codon_gt",
    "find_cross_codon_cg",
    "find_cross_codon_restriction",
    "check_no_stop_codons",
    "check_no_cryptic_promoter",
    "check_valid_coding_seq",
    "check_conservation_score",
    "check_codon_optimality",
    "check_no_blast_matches",
    "check_primer_compatibility",
    # ── mrna_checks ────────────────────────────────────────────────────
    "check_no_unexpected_tm_domain",
    "check_mrna_secondary_structure",
    "check_co_translational_folding",
    "check_mrna_stability",
    "check_no_cryptic_orf",
    # ── splice_checks ──────────────────────────────────────────────────
    "check_no_cryptic_splice",
    "check_no_gt_dinucleotide",
    "check_no_avoidable_gt",
    "check_no_gt_dinucleotide_soft",
    # ── cpg_checks ─────────────────────────────────────────────────────
    "check_no_cpg_island",
    # ── restriction_checks ─────────────────────────────────────────────
    "check_no_restriction_site",
    # ── mirna_checks ───────────────────────────────────────────────────
    "_rna_revcomp_to_dna",
    "_mirna_context_score",
    "check_no_mirna_binding_site",
    # ── immuno_checks ──────────────────────────────────────────────────
    "check_no_m6a_site",
    "check_no_polya_signal",
    "check_nucleoside_modification_guidance",
    "check_no_rqc_trigger",
    # ── alu_checks ─────────────────────────────────────────────────────
    "check_no_alu_repeat",
]
