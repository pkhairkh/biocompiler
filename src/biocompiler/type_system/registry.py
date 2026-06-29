"""
BioCompiler Type System — Predicate Registry
=============================================
Named dispatch for certificate generation and verification.
"""

import inspect
from typing import Dict, List

from biocompiler.shared.types import TypeCheckResult
from biocompiler.shared.exceptions import UnknownPredicateError
from .predicates import (
    evaluate_co_translational_folding,
    evaluate_codon_adapted,
    evaluate_codon_optimality,
    evaluate_conservation_score,
    evaluate_gc_in_range,
    evaluate_in_frame,
    evaluate_mrna_secondary_structure,
    evaluate_mrna_stability,
    evaluate_no_alu_repeat,
    evaluate_no_blast_matches,
    evaluate_no_cpg_island,
    evaluate_no_cryptic_orf,
    evaluate_no_cryptic_promoter,
    evaluate_no_cryptic_splice,
    evaluate_no_m6a_site,
    evaluate_no_mirna_binding_site,
    evaluate_no_gt_dinucleotide,
    evaluate_no_instability_motif,
    evaluate_no_polya_signal,
    evaluate_no_restriction_site,
    evaluate_no_rqc_trigger,
    evaluate_no_stop_codons,
    evaluate_no_unexpected_tm_domain,
    evaluate_nucleoside_modification_guidance,
    evaluate_primer_compatibility,
    evaluate_splice_correct,
    evaluate_valid_coding_seq,
)


class PredicateRegistry:
    """Registry of named type predicates with evaluate() and verify() dispatch.

    The registry provides a single entry point for certificate generation
    and verification, mapping predicate names to their evaluate functions.
    It supports both evaluation (with default parameters) and verification
    (re-running a predicate with specific parameters from a certificate).

    Diagnostic predicates
    ---------------------
    A small number of predicate names are emitted by the optimization
    pipeline but are NOT part of the 43-predicate canonical contract
    (``PREDICATE_NAMES`` in ``type_system/codon_tables.py`` and the
    ``test_predicate_registry_has_43`` regression test). The notable case
    is ``MRNAStability`` (emitted by ``check_mrna_stability`` on both the
    fast and slow paths via ``evaluate_extended_predicates``), which is a
    composite diagnostic score computed differently from the 43 canonical
    predicates.

    Such predicates are registered via :meth:`register_diagnostic`. They
    are dispatchable through :meth:`evaluate` / :meth:`verify` /
    :meth:`__contains__` (so the certificate verifier can re-evaluate
    them) but are NOT listed by :meth:`names` (so the 43-predicate
    contract is preserved). :meth:`diagnostic_names` returns them
    separately for callers that need the full set.
    """

    def __init__(self) -> None:
        self._predicates: Dict[str, callable] = {}
        self._diagnostic_predicates: Dict[str, callable] = {}
        self._verify_params: Dict[str, Dict[str, str]] = {}

    def register(
        self,
        name: str,
        fn: callable,
        verify_param_map: Dict[str, str] | None = None,
    ) -> None:
        """Register a predicate evaluation function.

        Args:
            name: Predicate name (e.g., 'NoCrypticSplice').
            fn: Callable that returns TypeCheckResult.
            verify_param_map: Optional mapping from certificate param names
                to function kwarg names.
        """
        self._predicates[name] = fn
        if verify_param_map:
            self._verify_params[name] = verify_param_map

    def register_diagnostic(
        self,
        name: str,
        fn: callable,
        verify_param_map: Dict[str, str] | None = None,
    ) -> None:
        """Register a *diagnostic* predicate that is dispatchable but not
        counted in the canonical 43-predicate contract.

        Diagnostic predicates participate in ``evaluate`` / ``verify`` /
        ``__contains__`` exactly like canonical predicates, but are
        excluded from :meth:`names` (use :meth:`diagnostic_names` to list
        them). This lets the certificate verifier re-evaluate predicates
        such as ``MRNAStability`` without breaking the
        ``test_predicate_registry_has_43`` contract.

        Args:
            name: Predicate name (e.g., 'MRNAStability').
            fn: Callable that returns TypeCheckResult.
            verify_param_map: Optional mapping from certificate param names
                to function kwarg names.
        """
        self._diagnostic_predicates[name] = fn
        if verify_param_map:
            self._verify_params[name] = verify_param_map

    def _ensure_protein_predicates(self) -> None:
        """Lazily register protein-level predicates on first access."""
        if not _protein_predicates_registered:
            _register_protein_predicates()

    def names(self) -> List[str]:
        """Return sorted list of *canonical* registered predicate names.

        Excludes diagnostic predicates (registered via
        :meth:`register_diagnostic`) so that the 43-predicate contract
        (``test_predicate_registry_has_43``) is preserved. Use
        :meth:`diagnostic_names` to retrieve the diagnostic subset.
        """
        self._ensure_protein_predicates()
        return sorted(self._predicates.keys())

    def diagnostic_names(self) -> List[str]:
        """Return sorted list of diagnostic predicate names.

        These predicates are dispatchable via :meth:`evaluate` /
        :meth:`verify` / :meth:`__contains__` but not listed in
        :meth:`names` (so the canonical 43-predicate contract is
        preserved).
        """
        return sorted(self._diagnostic_predicates.keys())

    def _resolve_fn(self, name: str):
        """Return the evaluation function for ``name`` (canonical or
        diagnostic), or ``None`` if not registered."""
        self._ensure_protein_predicates()
        if name in self._predicates:
            return self._predicates[name]
        if name in self._diagnostic_predicates:
            return self._diagnostic_predicates[name]
        return None

    def evaluate(self, name: str, **kwargs) -> TypeCheckResult:
        """Evaluate a named predicate.

        Args:
            name: Predicate name.
            **kwargs: Arguments to pass to the predicate function.

        Returns:
            TypeCheckResult from the predicate.

        Raises:
            UnknownPredicateError: if name is not registered.
        """
        self._ensure_protein_predicates()
        fn = self._resolve_fn(name)
        if fn is None:
            raise UnknownPredicateError(name)
        return fn(**kwargs)

    def verify(self, name: str, **kwargs) -> TypeCheckResult:
        """Verify a predicate — same as evaluate but with certificate params.

        The verify method maps certificate parameter names to the function's
        expected kwargs before calling evaluate.

        Args:
            name: Predicate name.
            **kwargs: Certificate parameters, possibly with different names
                than the evaluate function expects.

        Returns:
            TypeCheckResult from re-evaluation.

        Raises:
            UnknownPredicateError: if name is not registered.
        """
        self._ensure_protein_predicates()
        fn = self._resolve_fn(name)
        if fn is None:
            raise UnknownPredicateError(name)

        # Map certificate-style kwargs to evaluate-style kwargs
        param_map = self._verify_params.get(name, {})
        mapped_kwargs = {}
        for cert_key, fn_key in param_map.items():
            if cert_key in kwargs:
                mapped_kwargs[fn_key] = kwargs[cert_key]

        # Pass through any kwargs that match the function's signature directly
        sig = inspect.signature(fn)
        for key, val in kwargs.items():
            if key in sig.parameters:
                mapped_kwargs[key] = val

        return fn(**mapped_kwargs)

    def __contains__(self, name: str) -> bool:
        self._ensure_protein_predicates()
        return (
            name in self._predicates
            or name in self._diagnostic_predicates
        )


# ────────────────────────────────────────────────────────────
# Global registry instance with all predicates registered
# ────────────────────────────────────────────────────────────
registry = PredicateRegistry()

registry.register(
    "NoCrypticSplice",
    evaluate_no_cryptic_splice,
    verify_param_map={
        "seq": "seq",
        "known_exon_boundaries": "boundaries",
        "organism": "organism",
        "cryptic_splice_threshold": "cryptic_threshold",
    },
)

registry.register(
    "SpliceCorrect",
    evaluate_splice_correct,
    verify_param_map={
        "seq": "seq",
        "known_exon_boundaries": "boundaries",
        "cellular_context": "cellular_context",
    },
)

registry.register(
    "GCInRange",
    evaluate_gc_in_range,
    verify_param_map={
        "seq": "seq",
        "gc_lo": "gc_lo",
        "gc_hi": "gc_hi",
    },
)

# SlidingGC predicate (local/sliding-window GC constraint)
from biocompiler.sequence.sliding_gc import evaluate_sliding_gc as _evaluate_sliding_gc

registry.register(
    "SlidingGC",
    _evaluate_sliding_gc,
    verify_param_map={
        "seq": "seq",
        "window_size": "window_size",
        "gc_min": "gc_min",
        "gc_max": "gc_max",
    },
)

registry.register(
    "CodonAdapted",
    evaluate_codon_adapted,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "threshold": "threshold",
    },
)

registry.register(
    "NoRestrictionSite",
    evaluate_no_restriction_site,
    verify_param_map={
        "seq": "seq",
        "enzyme_set": "enzymes",
    },
)

registry.register(
    "InFrame",
    evaluate_in_frame,
    verify_param_map={
        "seq": "seq",
        "exon_boundaries": "boundaries",
    },
)

registry.register(
    "NoInstabilityMotif",
    evaluate_no_instability_motif,
    verify_param_map={
        "seq": "seq",
    },
)

registry.register(
    "NoCpGIsland",
    evaluate_no_cpg_island,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
    },
)

registry.register(
    "NoStopCodons",
    evaluate_no_stop_codons,
    verify_param_map={
        "seq": "seq",
    },
)

registry.register(
    "NoGTDinucleotide",
    evaluate_no_gt_dinucleotide,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
    },
)

registry.register(
    "ValidCodingSeq",
    evaluate_valid_coding_seq,
    verify_param_map={
        "seq": "seq",
    },
)

registry.register(
    "ConservationScore",
    evaluate_conservation_score,
    verify_param_map={
        "seq": "seq",
        "protein": "protein",
        "min_score": "min_score",
    },
)

registry.register(
    "CodonOptimality",
    evaluate_codon_optimality,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "threshold": "threshold",
    },
)

registry.register(
    "NoCrypticPromoter",
    evaluate_no_cryptic_promoter,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "threshold": "threshold",
    },
)

registry.register(
    "CoTranslationalFolding",
    evaluate_co_translational_folding,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "domain_boundaries": "domain_boundaries",
        "min_pause_cai": "min_pause_cai",
    },
)

registry.register(
    "NoUnexpectedTMDomain",
    evaluate_no_unexpected_tm_domain,
    verify_param_map={
        "seq": "seq",
        "is_cytosolic": "is_cytosolic",
        "threshold": "threshold",
    },
)

registry.register(
    "mRNASecondaryStructure",
    evaluate_mrna_secondary_structure,
    verify_param_map={
        "seq": "seq",
        "window_start": "window_start",
        "window_end": "window_end",
        "dg_threshold": "dg_threshold",
    },
)

registry.register(
    "NoBlastMatches",
    evaluate_no_blast_matches,
    verify_param_map={
        "seq": "seq",
        "reference_sequences": "reference_sequences",
        "k": "k",
    },
)

registry.register(
    "PrimerCompatibility",
    evaluate_primer_compatibility,
    verify_param_map={
        "seq": "seq",
        "region_start": "region_start",
        "region_end": "region_end",
        "min_tm": "min_tm",
        "max_tm": "max_tm",
    },
)

# ── Extended Diagnostic Predicates ──

registry.register(
    "NoCrypticORF",
    evaluate_no_cryptic_orf,
    verify_param_map={
        "seq": "seq",
        "min_orf_length": "min_orf_length",
        "organism": "organism",
    },
)

registry.register(
    "NoRQCTrigger",
    evaluate_no_rqc_trigger,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "poly_a_min_length": "poly_a_min_length",
    },
)

registry.register(
    "NoAluRepeat",
    evaluate_no_alu_repeat,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "min_match_score": "min_match_score",
    },
)

registry.register(
    "NoMiRNABindingSite",
    evaluate_no_mirna_binding_site,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "min_seed_match": "min_seed_match",
        "tissue": "tissue",
    },
)

registry.register(
    "NoM6ASite",
    evaluate_no_m6a_site,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "scan_mode": "scan_mode",
    },
)

registry.register(
    "NoPolyASignal",
    evaluate_no_polya_signal,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "scan_cds_only": "scan_cds_only",
    },
)

registry.register(
    "NucleosideModificationGuidance",
    evaluate_nucleoside_modification_guidance,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "modification_type": "modification_type",
    },
)


# ────────────────────────────────────────────────────────────
# Diagnostic predicates
#
# These predicate names are emitted by the optimization pipeline
# (``check_mrna_stability`` on both the fast and slow paths via
# ``evaluate_extended_predicates``) but are NOT part of the 43-predicate
# canonical contract (``PREDICATE_NAMES`` in ``type_system/codon_tables.py``
# and the ``test_predicate_registry_has_43`` regression test). They are
# registered via ``register_diagnostic`` so that the certificate verifier
# (``verify_certificate`` in ``provenance/certificate.py``) can re-evaluate
# them, while ``registry.names()`` continues to return exactly the 43
# canonical predicate names.
# ────────────────────────────────────────────────────────────

registry.register_diagnostic(
    "MRNAStability",
    evaluate_mrna_stability,
    verify_param_map={
        "seq": "seq",
        "organism": "organism",
        "threshold": "threshold",
    },
)


# ────────────────────────────────────────────────────────────
# Protein-level predicates (stability, solubility,
# immunogenicity, structure)
#
# Deferred registration to avoid circular imports — these
# modules import from type_system themselves, so we register
# them lazily after the package is fully loaded.
# ────────────────────────────────────────────────────────────

_protein_predicates_registered = False


def _register_protein_predicates() -> None:
    """Register protein-level predicates (deferred to break circular imports)."""
    global _protein_predicates_registered
    if _protein_predicates_registered:
        return
    _protein_predicates_registered = True

    from biocompiler.type_system.stability_predicates import (
        evaluate_stable_folding,
        evaluate_no_destabilizing_mutation,
        evaluate_disulfide_bond_integrity,
        evaluate_hydrophobic_core_quality,
    )
    from biocompiler.type_system.solubility_predicates import (
        evaluate_soluble_expression,
        evaluate_no_aggregation_prone_region,
        evaluate_charge_composition,
        evaluate_no_long_hydrophobic_stretch,
    )
    from biocompiler.immunogenicity.predicates import (
        evaluate_low_immunogenicity,
        evaluate_no_strong_t_cell_epitope,
        evaluate_no_dominant_b_cell_epitope,
        evaluate_population_coverage_safe,
    )
    from ..structure.predicates import (
        evaluate_structure_confidence,
        evaluate_no_misfolding_risk,
        evaluate_correct_fold_topology,
        evaluate_no_unexpected_interaction,
    )

    # Structure predicates
    registry.register(
        "StructureConfidence",
        evaluate_structure_confidence,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "pdb_string": "pdb_string",
        },
    )
    registry.register(
        "NoMisfoldingRisk",
        evaluate_no_misfolding_risk,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "pdb_string": "pdb_string",
        },
    )
    registry.register(
        "CorrectFoldTopology",
        evaluate_correct_fold_topology,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "pdb_string": "pdb_string",
        },
    )
    registry.register(
        "NoUnexpectedInteraction",
        evaluate_no_unexpected_interaction,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "pdb_string": "pdb_string",
        },
    )

    # Stability predicates
    registry.register(
        "StableFolding",
        evaluate_stable_folding,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "stability_threshold": "stability_threshold",
            "pdb_string": "pdb_string",
        },
    )
    registry.register(
        "NoDestabilizingMutation",
        evaluate_no_destabilizing_mutation,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "max_ddg": "max_ddg",
            "pdb_string": "pdb_string",
        },
    )
    registry.register(
        "DisulfideBondIntegrity",
        evaluate_disulfide_bond_integrity,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "pdb_string": "pdb_string",
        },
    )
    registry.register(
        "HydrophobicCoreQuality",
        evaluate_hydrophobic_core_quality,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "pdb_string": "pdb_string",
        },
    )

    # Solubility predicates
    registry.register(
        "SolubleExpression",
        evaluate_soluble_expression,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "min_solubility_score": "min_solubility_score",
            "pdb_string": "pdb_string",
        },
    )
    registry.register(
        "NoAggregationProneRegion",
        evaluate_no_aggregation_prone_region,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "pdb_string": "pdb_string",
        },
    )
    registry.register(
        "ChargeComposition",
        evaluate_charge_composition,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "pdb_string": "pdb_string",
        },
    )
    registry.register(
        "NoLongHydrophobicStretch",
        evaluate_no_long_hydrophobic_stretch,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "pdb_string": "pdb_string",
        },
    )

    # Immunogenicity predicates
    registry.register(
        "LowImmunogenicity",
        evaluate_low_immunogenicity,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "max_immunogenicity_score": "max_immunogenicity_score",
        },
    )
    registry.register(
        "NoStrongTCellEpitope",
        evaluate_no_strong_t_cell_epitope,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "mhc_alleles": "mhc_alleles",
        },
    )
    registry.register(
        "NoDominantBCellEpitope",
        evaluate_no_dominant_b_cell_epitope,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
        },
    )
    registry.register(
        "PopulationCoverageSafe",
        evaluate_population_coverage_safe,
        verify_param_map={
            "sequence": "sequence",
            "protein": "protein",
            "organism": "organism",
            "mhc_alleles": "mhc_alleles",
        },
    )


# Protein-level predicates are registered lazily via
# PredicateRegistry._ensure_protein_predicates() to break
# circular imports with stability_predicates, solubility_predicates,
# immuno_predicates, and structure.predicates (which all import from
# type_system).
