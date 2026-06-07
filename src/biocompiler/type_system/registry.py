"""
BioCompiler Type System — Predicate Registry
=============================================
Named dispatch for certificate generation and verification.
"""

import inspect
from typing import Dict, List

from ..types import TypeCheckResult
from ..exceptions import UnknownPredicateError
from .predicates import (
    evaluate_co_translational_folding,
    evaluate_codon_adapted,
    evaluate_codon_optimality,
    evaluate_conservation_score,
    evaluate_gc_in_range,
    evaluate_in_frame,
    evaluate_mrna_secondary_structure,
    evaluate_no_cpg_island,
    evaluate_no_cryptic_promoter,
    evaluate_no_cryptic_splice,
    evaluate_no_gt_dinucleotide,
    evaluate_no_instability_motif,
    evaluate_no_restriction_site,
    evaluate_no_stop_codons,
    evaluate_no_unexpected_tm_domain,
    evaluate_splice_correct,
    evaluate_valid_coding_seq,
)


class PredicateRegistry:
    """Registry of named type predicates with evaluate() and verify() dispatch.

    The registry provides a single entry point for certificate generation
    and verification, mapping predicate names to their evaluate functions.
    It supports both evaluation (with default parameters) and verification
    (re-running a predicate with specific parameters from a certificate).
    """

    def __init__(self) -> None:
        self._predicates: Dict[str, callable] = {}
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

    def _ensure_protein_predicates(self) -> None:
        """Lazily register protein-level predicates on first access."""
        if not _protein_predicates_registered:
            _register_protein_predicates()

    def names(self) -> List[str]:
        """Return sorted list of registered predicate names."""
        self._ensure_protein_predicates()
        return sorted(self._predicates.keys())

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
        if name not in self._predicates:
            raise UnknownPredicateError(name)
        return self._predicates[name](**kwargs)

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
        if name not in self._predicates:
            raise UnknownPredicateError(name)

        # Map certificate-style kwargs to evaluate-style kwargs
        param_map = self._verify_params.get(name, {})
        mapped_kwargs = {}
        for cert_key, fn_key in param_map.items():
            if cert_key in kwargs:
                mapped_kwargs[fn_key] = kwargs[cert_key]

        # Pass through any kwargs that match the function's signature directly
        fn = self._predicates[name]
        sig = inspect.signature(fn)
        for key, val in kwargs.items():
            if key in sig.parameters:
                mapped_kwargs[key] = val

        return fn(**mapped_kwargs)

    def __contains__(self, name: str) -> bool:
        return name in self._predicates


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
from ..sliding_gc import evaluate_sliding_gc as _evaluate_sliding_gc

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

    from ..stability_predicates import (
        evaluate_stable_folding,
        evaluate_no_destabilizing_mutation,
        evaluate_disulfide_bond_integrity,
        evaluate_hydrophobic_core_quality,
    )
    from ..solubility_predicates import (
        evaluate_soluble_expression,
        evaluate_no_aggregation_prone_region,
        evaluate_charge_composition,
        evaluate_no_long_hydrophobic_stretch,
    )
    from ..immuno_predicates import (
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
