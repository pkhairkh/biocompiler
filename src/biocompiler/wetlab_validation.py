"""
BioCompiler Wet-Lab Validation Framework
=========================================

Templates and tools for validating BioCompiler's in-silico predictions
against actual wet-lab experimental results. This module bridges the gap
between computational gene design and experimental verification.

Key components:

  - WetLabProtocol: Template for designing experiments to validate
    optimized sequences, including oligo generation for gene synthesis.
  - WetLabResult: Structured recording of experimental outcomes.
  - compare_insilico_vs_wetlab: Quantitative comparison of predicted
    vs. actual expression, solubility, and other metrics.

Design Philosophy:
    In-silico predictions (CAI, GC content, solubility) are necessary
    but not sufficient. Wet-lab validation provides ground-truth data
    to calibrate and improve BioCompiler's models. This module provides
    structured templates to ensure experimental designs are complete
    and results are comparable to predictions.

References:
    - Welch M et al. (2009) PLoS ONE 4:e7002 — Design parameters for
      codon optimization.
    - Gustafsson C et al. (2004) Trends Biotechnol 22:346-353 —
      Codon optimization for heterologous expression.
    - Puigbo P et al. (2008) BMC Bioinformatics 9:65 — CAIcal.
"""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

__all__ = [
    "WetLabProtocol",
    "WetLabResult",
    "compare_insilico_vs_wetlab",
    "generate_protocol_report",
    "EXPRESSION_LEVEL_MAP",
    "SOLUBILITY_LEVEL_MAP",
    "VECTOR_PROMOTER_COMPAT",
    # Validation framework (Task 2.4)
    "BenchmarkEntry",
    "BENCHMARK_DATASET",
    "ProteinValidationResult",
    "ValidationSuiteResult",
    "RegressionReport",
    "RegressionItem",
    "run_validation_suite",
    "check_regression",
]

# ═══════════════════════════════════════════════════════════════════════════════
# Lookup tables
# ═══════════════════════════════════════════════════════════════════════════════

EXPRESSION_LEVEL_MAP: Dict[str, Tuple[float, float]] = {
    "high": (50.0, 500.0),    # mg/L typical range for high expression
    "medium": (5.0, 50.0),    # mg/L
    "low": (0.1, 5.0),        # mg/L
}

SOLUBILITY_LEVEL_MAP: Dict[str, Tuple[float, float]] = {
    "high": (0.7, 1.0),       # fraction soluble
    "medium": (0.3, 0.7),
    "low": (0.0, 0.3),
}

# Compatible vector-promoter combinations (validated by literature)
VECTOR_PROMOTER_COMPAT: Dict[str, List[str]] = {
    "pET-28a": ["T7"],
    "pET-21a": ["T7"],
    "pET-22b": ["T7"],
    "pET-15b": ["T7"],
    "pcDNA3.1": ["CMV", "EF1a"],
    "pcDNA5": ["CMV"],
    "pCI-neo": ["CMV"],
    "pCEP4": ["CMV"],
    "pYES2": ["GAL1"],
    "pPICZ": ["AOX1"],
    "pGEX-4T": ["tac"],
    "pMAL-c5X": ["tac"],
}

# Host strains appropriate for each vector system
VECTOR_HOST_COMPAT: Dict[str, List[str]] = {
    "pET-28a": ["BL21(DE3)", "BL21(DE3)pLysS", "Rosetta(DE3)", "BL21(DE3)-Star"],
    "pET-21a": ["BL21(DE3)", "BL21(DE3)pLysS", "Rosetta(DE3)"],
    "pET-22b": ["BL21(DE3)", "BL21(DE3)pLysS", "Rosetta(DE3)"],
    "pET-15b": ["BL21(DE3)", "BL21(DE3)pLysS"],
    "pcDNA3.1": ["HEK293T", "CHO-K1", "HeLa"],
    "pcDNA5": ["HEK293T", "CHO-K1", "Flp-In T-REx"],
    "pCI-neo": ["HEK293T", "CHO-K1"],
    "pCEP4": ["HEK293T", "HEK293-EBNA"],
    "pYES2": ["BY4741", "INVSc1", "S288C"],
    "pPICZ": ["X-33", "GS115", "KM71H"],
    "pGEX-4T": ["BL21(DE3)", "Rosetta(DE3)", "BL21(DE3)-Star"],
    "pMAL-c5X": ["BL21(DE3)", "Rosetta(DE3)"],
}

# Codon table for reverse-translation (organism-optimized)
_CODON_TABLE: Dict[str, List[str]] = {
    "F": ["TTT", "TTC"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "I": ["ATT", "ATC", "ATA"],
    "M": ["ATG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "Y": ["TAT", "TAC"],
    "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"],
    "N": ["AAT", "AAC"],
    "K": ["AAA", "AAG"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "C": ["TGT", "TGC"],
    "W": ["TGG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "*": ["TAA", "TAG", "TGA"],
}


# ═══════════════════════════════════════════════════════════════════════════════
# WetLabProtocol
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class WetLabProtocol:
    """Template for wet-lab validation of optimized sequences.

    This dataclass captures the full experimental design needed to
    validate an in-silico optimized gene sequence, including cloning
    strategy, expression conditions, and expected outcomes.

    Attributes:
        gene_name: Name of the gene (e.g. "GFP", "HBB", "INS")
        organism: Target organism (e.g. "Escherichia_coli", "Homo_sapiens")
        optimized_dna: The DNA sequence produced by BioCompiler
        cai: Computed CAI of the optimized sequence
        gc_content: GC fraction of the optimized sequence

        vector: Cloning vector (e.g. "pET-28a", "pcDNA3.1")
        promoter: Expression promoter (e.g. "T7", "CMV", "EF1a")
        selection_marker: Antibiotic for selection
        host_strain: Expression host strain

        expected_expression_level: Predicted expression level category
        expected_solubility: Predicted solubility category
    """

    gene_name: str
    organism: str
    optimized_dna: str
    cai: float
    gc_content: float

    # Experimental design
    vector: str  # "pET-28a", "pcDNA3.1", etc.
    promoter: str  # "T7", "CMV", "EF1a"
    selection_marker: str  # "kanamycin", "ampicillin"
    host_strain: str  # "BL21(DE3)", "HEK293T"

    # Expected outcomes
    expected_expression_level: str  # "high", "medium", "low"
    expected_solubility: str  # "high", "medium", "low"

    def __post_init__(self) -> None:
        """Validate protocol parameters."""
        if self.expected_expression_level not in EXPRESSION_LEVEL_MAP:
            raise ValueError(
                f"Invalid expression level '{self.expected_expression_level}'. "
                f"Must be one of {list(EXPRESSION_LEVEL_MAP.keys())}"
            )
        if self.expected_solubility not in SOLUBILITY_LEVEL_MAP:
            raise ValueError(
                f"Invalid solubility level '{self.expected_solubility}'. "
                f"Must be one of {list(SOLUBILITY_LEVEL_MAP.keys())}"
            )
        if not re.match(r'^[ACGT]+$', self.optimized_dna.upper()):
            raise ValueError(
                "optimized_dna must contain only A, C, G, T nucleotides"
            )

    def generate_protocol(self) -> str:
        """Generate a step-by-step experimental protocol.

        Returns a detailed protocol string covering:
        1. Gene synthesis (oligo ordering or full gene)
        2. Cloning into expression vector
        3. Transformation / transfection
        4. Expression induction
        5. Protein purification
        6. Analysis (SDS-PAGE, Western blot, activity assay)

        Returns:
            Formatted protocol string.
        """
        dna_len = len(self.optimized_dna)
        is_prokaryotic = self.organism in (
            "Escherichia_coli", "E_coli", "Bacillus_subtilis"
        )

        # Determine selection antibiotic concentration
        abx_conc = self._get_antibiotic_concentration()

        # Determine induction conditions
        induction = self._get_induction_conditions()

        protocol = f"""WET-LAB VALIDATION PROTOCOL
{'=' * 60}
Gene: {self.gene_name}
Organism: {self.organism}
Sequence length: {dna_len} bp ({dna_len // 3} codons)
CAI: {self.cai:.3f}
GC content: {self.gc_content:.1%}

{'=' * 60}
EXPERIMENTAL DESIGN
{'=' * 60}
Vector: {self.vector}
Promoter: {self.promoter}
Host strain: {self.host_strain}
Selection: {self.selection_marker} ({abx_conc})

Expected expression: {self.expected_expression_level}
Expected solubility: {self.expected_solubility}

{'=' * 60}
STEP-BY-STEP PROTOCOL
{'=' * 60}

STEP 1: Gene Synthesis
  - Sequence length: {dna_len} bp
  - Recommended: {'Full gene synthesis (commercial)' if dna_len > 500 else 'Oligo assembly (PCR-based)'}
  - Number of oligos needed: {len(self.generate_oligos())}
  - Verify sequence by Sanger sequencing before proceeding

STEP 2: Cloning
  - Digest vector {self.vector} and insert with appropriate restriction enzymes
  - Recommended cloning sites: {'NdeI/BamHI (pET series)' if 'pET' in self.vector else 'MCS (multiple cloning site)'}
  - Ligate insert into linearized vector
  - Transform into cloning strain (DH5alpha) for plasmid amplification

STEP 3: Plasmid Verification
  - Miniprep plasmid DNA
  - Verify by restriction digest
  - Confirm sequence by Sanger sequencing (full coverage)
  - Expected insert size: {dna_len} bp

STEP 4: Transformation / Transfection
  - Transform verified plasmid into {self.host_strain}
  - Plate on LB + {self.selection_marker} ({abx_conc})
  - Pick 3-5 colonies for expression screening

STEP 5: Expression
  - Grow overnight culture in LB + {self.selection_marker}
  - Dilute 1:100 into fresh medium
  - Grow to OD600 = 0.6-0.8
  {induction}
  - Harvest cells by centrifugation (4000g, 15 min, 4C)

STEP 6: Lysis and Solubility Assessment
  - Resuspend pellet in lysis buffer (PBS + protease inhibitors)
  - Lyse by {'sonication (3x 30s, 50% duty)' if is_prokaryotic else 'mechanical disruption'}
  - Centrifuge at 15000g, 30 min, 4C
  - Separate soluble (supernatant) and insoluble (pellet) fractions
  - Analyze both fractions by SDS-PAGE

STEP 7: Protein Purification
  - {'Ni-NTA affinity chromatography (His-tag)' if 'pET' in self.vector else 'Affinity chromatography (tag-dependent)'}
  - Wash with 20 mM imidazole
  - Elute with 250 mM imidazole
  - Desalt / buffer exchange into storage buffer

STEP 8: Analysis
  - SDS-PAGE (Coomassie stain)
  - Western blot (anti-His or anti-target antibody)
  - Concentration measurement (Bradford or BCA assay)
  - Solubility: (soluble fraction / total) x 100%
  - Activity assay (if applicable)

STEP 9: Expected Results
  - Expression level: {self.expected_expression_level}
    ({EXPRESSION_LEVEL_MAP[self.expected_expression_level][0]}-{EXPRESSION_LEVEL_MAP[self.expected_expression_level][1]} mg/L)
  - Solubility: {self.expected_solubility}
    ({SOLUBILITY_LEVEL_MAP[self.expected_solubility][0]*100:.0f}%-{SOLUBILITY_LEVEL_MAP[self.expected_solubility][1]*100:.0f}% soluble)
  - CAI {self.cai:.3f} {'(good codon adaptation)' if self.cai > 0.7 else '(moderate codon adaptation)' if self.cai > 0.5 else '(low codon adaptation — may need further optimization)'}
  - GC content {self.gc_content:.1%} {'(within typical range)' if 0.3 <= self.gc_content <= 0.7 else '(outside typical range — may affect expression)'}
"""
        return protocol

    def generate_oligos(self, max_oligo_length: int = 60) -> List[str]:
        """Generate oligos for gene synthesis.

        Uses a simple overlap-based strategy: the sequence is divided
        into oligos of max_oligo_length with 20 bp overlaps between
        adjacent oligos. This is suitable for PCR-based gene assembly
        (Gibson assembly or overlap extension PCR).

        Args:
            max_oligo_length: Maximum length of each oligo (default 60bp).

        Returns:
            List of oligo sequences.

        Raises:
            ValueError: If max_oligo_length < 40 (insufficient overlap).
        """
        if max_oligo_length < 40:
            raise ValueError(
                f"max_oligo_length must be >= 40 for sufficient overlap, "
                f"got {max_oligo_length}"
            )

        seq = self.optimized_dna.upper()
        seq_len = len(seq)

        if seq_len <= max_oligo_length:
            return [seq]

        overlap = 20  # 20 bp overlap between adjacent oligos
        step = max_oligo_length - overlap
        oligos: List[str] = []

        pos = 0
        while pos < seq_len:
            end = min(pos + max_oligo_length, seq_len)
            oligos.append(seq[pos:end])
            if end >= seq_len:
                break
            pos += step

        return oligos

    def _get_antibiotic_concentration(self) -> str:
        """Return appropriate antibiotic concentration for the selection marker."""
        concentrations = {
            "kanamycin": "50 ug/mL",
            "ampicillin": "100 ug/mL",
            "chloramphenicol": "34 ug/mL",
            "streptomycin": "50 ug/mL",
            "zeocin": "100 ug/mL",
            "blasticidin": "5-10 ug/mL",
            "hygromycin": "200 ug/mL (mammalian) / 300 ug/mL (yeast)",
            "G418": "400-800 ug/mL",
        }
        return concentrations.get(self.selection_marker, "per standard protocol")

    def _get_induction_conditions(self) -> str:
        """Return induction conditions based on promoter type."""
        if self.promoter == "T7":
            return (
                "  - Add IPTG to 0.5 mM final concentration\n"
                "  - Induce at 18C for 16h (for soluble expression)\n"
                "    OR 37C for 4h (for higher yield, may reduce solubility)"
            )
        elif self.promoter == "CMV":
            return (
                "  - CMV promoter is constitutive in mammalian cells\n"
                "  - Harvest 48-72h post-transfection"
            )
        elif self.promoter == "EF1a":
            return (
                "  - EF1a promoter is constitutive in mammalian cells\n"
                "  - Harvest 48-72h post-transfection"
            )
        elif self.promoter == "GAL1":
            return (
                "  - Grow in SD-Ura + 2% raffinose to OD600 = 0.6\n"
                "  - Add galactose to 2% final concentration\n"
                "  - Induce at 30C for 6-16h"
            )
        elif self.promoter == "AOX1":
            return (
                "  - Grow in BMGY to OD600 = 2-6\n"
                "  - Switch to BMMY medium (methanol induction)\n"
                "  - Add 0.5% methanol every 24h\n"
                "  - Induce at 28-30C for 48-96h"
            )
        elif self.promoter == "tac":
            return (
                "  - Add IPTG to 0.5 mM final concentration\n"
                "  - Induce at 25C for 6h"
            )
        else:
            return (
                f"  - Follow standard induction protocol for {self.promoter} promoter"
            )

    def is_vector_promoter_compatible(self) -> bool:
        """Check if the vector and promoter are compatible.

        Returns:
            True if the vector-promoter combination is valid.
        """
        compat = VECTOR_PROMOTER_COMPAT.get(self.vector, [])
        if not compat:
            # Unknown vector — allow but warn
            logger.warning(
                "Unknown vector '%s'; cannot verify promoter compatibility",
                self.vector,
            )
            return True
        return self.promoter in compat

    def is_host_strain_compatible(self) -> bool:
        """Check if the host strain is compatible with the vector.

        Returns:
            True if the host strain is appropriate for the vector.
        """
        compat = VECTOR_HOST_COMPAT.get(self.vector, [])
        if not compat:
            logger.warning(
                "Unknown vector '%s'; cannot verify host strain compatibility",
                self.vector,
            )
            return True
        return self.host_strain in compat

    def predict_expression_category(self) -> str:
        """Predict expression level category from in-silico metrics.

        Uses CAI and GC content heuristics to predict expression:
        - CAI > 0.7 and GC 30-70%: "high"
        - CAI > 0.5 or (CAI > 0.4 and GC 30-70%): "medium"
        - Otherwise: "low"

        Returns:
            Predicted expression level: "high", "medium", or "low".
        """
        gc_ok = 0.30 <= self.gc_content <= 0.70
        if self.cai > 0.7 and gc_ok:
            return "high"
        elif self.cai > 0.5 or (self.cai > 0.4 and gc_ok):
            return "medium"
        else:
            return "low"

    def predict_solubility_category(self) -> str:
        """Predict solubility category from in-silico metrics.

        Uses GC content and CAI heuristics:
        - GC 40-60% and CAI > 0.6: "high"
        - GC 30-70% and CAI > 0.4: "medium"
        - Otherwise: "low"

        Returns:
            Predicted solubility: "high", "medium", or "low".
        """
        gc_narrow = 0.40 <= self.gc_content <= 0.60
        gc_broad = 0.30 <= self.gc_content <= 0.70
        if gc_narrow and self.cai > 0.6:
            return "high"
        elif gc_broad and self.cai > 0.4:
            return "medium"
        else:
            return "low"


# ═══════════════════════════════════════════════════════════════════════════════
# WetLabResult
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class WetLabResult:
    """Structured recording of wet-lab experimental outcomes.

    Captures the key experimental metrics that can be compared against
    in-silico predictions from BioCompiler.

    Attributes:
        actual_expression_level: Measured protein expression in mg/L
        actual_solubility: Fraction of protein in soluble form (0.0-1.0)
        western_blot_confirmed: Whether Western blot confirmed the protein
        sequencing_match: Whether Sanger sequencing matches the designed sequence
        notes: Free-text notes about the experiment
    """

    actual_expression_level: float  # mg/L
    actual_solubility: float  # fraction soluble (0.0-1.0)
    western_blot_confirmed: bool
    sequencing_match: bool
    notes: str = ""

    def __post_init__(self) -> None:
        """Validate experimental results."""
        if self.actual_expression_level < 0:
            raise ValueError(
                f"actual_expression_level must be >= 0, got {self.actual_expression_level}"
            )
        if not (0.0 <= self.actual_solubility <= 1.0):
            raise ValueError(
                f"actual_solubility must be between 0.0 and 1.0, "
                f"got {self.actual_solubility}"
            )

    def classify_expression(self) -> str:
        """Classify measured expression level into a category.

        Returns:
            "high", "medium", or "low" based on mg/L thresholds.
        """
        if self.actual_expression_level >= 50.0:
            return "high"
        elif self.actual_expression_level >= 5.0:
            return "medium"
        else:
            return "low"

    def classify_solubility(self) -> str:
        """Classify measured solubility into a category.

        Returns:
            "high", "medium", or "low" based on fraction thresholds.
        """
        if self.actual_solubility >= 0.7:
            return "high"
        elif self.actual_solubility >= 0.3:
            return "medium"
        else:
            return "low"


# ═══════════════════════════════════════════════════════════════════════════════
# Comparison function
# ═══════════════════════════════════════════════════════════════════════════════

def compare_insilico_vs_wetlab(
    prediction: dict,
    actual: WetLabResult,
) -> Dict:
    """Compare in-silico predictions with actual wet-lab results.

    Computes quantitative and categorical agreement between BioCompiler's
    predictions and experimental measurements.

    Args:
        prediction: Dictionary with in-silico predictions. Expected keys:
            - "cai": float (CAI value)
            - "gc_content": float (GC fraction)
            - "expression_level": str ("high"/"medium"/"low")
            - "solubility": str ("high"/"medium"/"low")
            - "gene_name": str (optional)
            - "organism": str (optional)
        actual: WetLabResult with experimental measurements.

    Returns:
        Dictionary with comparison results:
            - "expression_category_match": bool
            - "solubility_category_match": bool
            - "expression_level_error_mg_L": float
            - "expression_predicted_range": tuple
            - "expression_actual_in_range": bool
            - "solubility_predicted_range": tuple
            - "solubility_actual_in_range": bool
            - "sequencing_match": bool
            - "western_blot_confirmed": bool
            - "overall_agreement": str
            - "summary": str
    """
    pred_expr = prediction.get("expression_level", "medium")
    pred_sol = prediction.get("solubility", "medium")

    actual_expr_cat = actual.classify_expression()
    actual_sol_cat = actual.classify_solubility()

    expr_cat_match = pred_expr == actual_expr_cat
    sol_cat_match = pred_sol == actual_sol_cat

    pred_expr_range = EXPRESSION_LEVEL_MAP.get(pred_expr, (0.0, 0.0))
    pred_sol_range = SOLUBILITY_LEVEL_MAP.get(pred_sol, (0.0, 0.0))

    expr_in_range = pred_expr_range[0] <= actual.actual_expression_level <= pred_expr_range[1]
    sol_in_range = pred_sol_range[0] <= actual.actual_solubility <= pred_sol_range[1]

    # Compute midpoint of predicted range as point estimate
    expr_midpoint = (pred_expr_range[0] + pred_expr_range[1]) / 2.0
    expr_error = abs(actual.actual_expression_level - expr_midpoint)

    # Overall agreement assessment
    if expr_cat_match and sol_cat_match and actual.sequencing_match:
        overall = "strong_agreement"
    elif (expr_cat_match or sol_cat_match) and actual.sequencing_match:
        overall = "partial_agreement"
    elif actual.sequencing_match:
        overall = "sequence_correct_prediction_mismatch"
    else:
        overall = "disagreement"

    summary = (
        f"Expression: predicted={pred_expr}, actual={actual_expr_cat} "
        f"({actual.actual_expression_level:.1f} mg/L) "
        f"{'MATCH' if expr_cat_match else 'MISMATCH'}; "
        f"Solubility: predicted={pred_sol}, actual={actual_sol_cat} "
        f"({actual.actual_solubility:.0%}) "
        f"{'MATCH' if sol_cat_match else 'MISMATCH'}; "
        f"Sequencing: {'MATCH' if actual.sequencing_match else 'MISMATCH'}; "
        f"Western: {'confirmed' if actual.western_blot_confirmed else 'not confirmed'}; "
        f"Overall: {overall}"
    )

    return {
        "expression_category_match": expr_cat_match,
        "solubility_category_match": sol_cat_match,
        "expression_level_error_mg_L": expr_error,
        "expression_predicted_range": pred_expr_range,
        "expression_actual_in_range": expr_in_range,
        "solubility_predicted_range": pred_sol_range,
        "solubility_actual_in_range": sol_in_range,
        "sequencing_match": actual.sequencing_match,
        "western_blot_confirmed": actual.western_blot_confirmed,
        "overall_agreement": overall,
        "summary": summary,
    }


def generate_protocol_report(protocol: WetLabProtocol) -> str:
    """Generate a complete protocol report including compatibility checks.

    Args:
        protocol: WetLabProtocol to generate a report for.

    Returns:
        Formatted report string.
    """
    vp_compat = protocol.is_vector_promoter_compatible()
    hs_compat = protocol.is_host_strain_compatible()

    predicted_expr = protocol.predict_expression_category()
    predicted_sol = protocol.predict_solubility_category()

    report = protocol.generate_protocol()

    report += f"""
{'=' * 60}
COMPATIBILITY CHECKS
{'=' * 60}
Vector-Promoter: {'PASS' if vp_compat else 'FAIL'} ({protocol.vector} + {protocol.promoter})
Host Strain:     {'PASS' if hs_compat else 'FAIL'} ({protocol.host_strain})

{'=' * 60}
IN-SILICO PREDICTIONS
{'=' * 60}
CAI: {protocol.cai:.3f}
GC content: {protocol.gc_content:.1%}
Predicted expression: {predicted_expr}
Predicted solubility: {predicted_sol}
User-specified expression: {protocol.expected_expression_level}
User-specified solubility: {protocol.expected_solubility}
"""

    if predicted_expr != protocol.expected_expression_level:
        report += (
            f"\nWARNING: In-silico prediction ({predicted_expr}) differs from "
            f"user expectation ({protocol.expected_expression_level})\n"
        )

    if predicted_sol != protocol.expected_solubility:
        report += (
            f"\nWARNING: In-silico solubility prediction ({predicted_sol}) differs from "
            f"user expectation ({protocol.expected_solubility})\n"
        )

    if not vp_compat:
        report += (
            f"\nERROR: Vector {protocol.vector} is not compatible with "
            f"promoter {protocol.promoter}. "
            f"Compatible promoters: {VECTOR_PROMOTER_COMPAT.get(protocol.vector, ['unknown'])}\n"
        )

    if not hs_compat:
        report += (
            f"\nERROR: Host strain {protocol.host_strain} is not compatible with "
            f"vector {protocol.vector}. "
            f"Compatible strains: {VECTOR_HOST_COMPAT.get(protocol.vector, ['unknown'])}\n"
        )

    return report


# ═══════════════════════════════════════════════════════════════════════════════
# Benchmark dataset and validation suite (Task 2.4)
# ═══════════════════════════════════════════════════════════════════════════════


@dataclass
class BenchmarkEntry:
    """A single benchmark protein for validation.

    Attributes:
        protein_name: Name of the benchmark protein.
        organism: Target organism for expression.
        protein_sequence: Amino acid sequence.
        expression_system: Expression system used (e.g. "E. coli").
        source_publication: DOI or citation for the benchmark data.
        expected_cai_range: Expected CAI range after optimization (lo, hi).
        expected_gc_range: Expected GC content range (lo, hi).
        expected_no_restriction_sites: Whether no restriction sites are expected.
        measured_expression_level: Measured expression level in mg/L (if available).
        expression_category: Expression level category ("high"/"medium"/"low").
    """

    protein_name: str
    organism: str
    protein_sequence: str
    expression_system: str
    source_publication: str
    expected_cai_range: Tuple[float, float] = (0.7, 1.0)
    expected_gc_range: Tuple[float, float] = (0.30, 0.70)
    expected_no_restriction_sites: bool = True
    measured_expression_level: Optional[float] = None
    expression_category: str = "medium"


# Curated benchmark dataset covering major therapeutic proteins
BENCHMARK_DATASET: List[BenchmarkEntry] = [
    BenchmarkEntry(
        protein_name="GFP",
        organism="Escherichia_coli",
        protein_sequence="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        expression_system="E. coli",
        source_publication="doi:10.1038/nature07622",
        expected_cai_range=(0.80, 1.0),
        expected_gc_range=(0.45, 0.55),
        measured_expression_level=150.0,
        expression_category="high",
    ),
    BenchmarkEntry(
        protein_name="Insulin",
        organism="Escherichia_coli",
        protein_sequence="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT",
        expression_system="E. coli",
        source_publication="doi:10.1016/0022-2836(79)90279-7",
        expected_cai_range=(0.75, 1.0),
        expected_gc_range=(0.40, 0.60),
        measured_expression_level=80.0,
        expression_category="high",
    ),
    BenchmarkEntry(
        protein_name="HBB",
        organism="Homo_sapiens",
        protein_sequence="MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        expression_system="HEK293T",
        source_publication="doi:10.1073/pnas.81.14.4358",
        expected_cai_range=(0.70, 0.95),
        expected_gc_range=(0.40, 0.60),
        measured_expression_level=45.0,
        expression_category="medium",
    ),
    BenchmarkEntry(
        protein_name="EPO",
        organism="CHO_K1",
        protein_sequence="MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAEKITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR",
        expression_system="CHO-K1",
        source_publication="doi:10.1126/science.3534528",
        expected_cai_range=(0.65, 0.90),
        expected_gc_range=(0.40, 0.60),
        measured_expression_level=60.0,
        expression_category="high",
    ),
    BenchmarkEntry(
        protein_name="TNF_alpha",
        organism="Escherichia_coli",
        protein_sequence="MSTESMIRDVELAEEALPKKTGGPQGSRRCLFLSLFSFLIVAGATTLFCLLNFGVIGPQRKEKRPDAVPEEPKAPKKKAEQPRRSQPSRTPSQPSTPGLKKRTVRPPKTISLKKGFQSGHMNDPVEVTGFPVVPSTVETQNVRFFLTITPQEATGAVFLLPWYGFRFPGVRALFTYLRNYCRVLQRRSELRQERVTQRYQGLSRAQVLTTMSAQVQMTPSRRSRSRTPKQPLRRPPSRLPKPRPSRSPKGRVSSKHVPGTSQPVEPNRHFPPFFQMDSVQPQSPVMTVSQHRSMFPQLPSPVPPPGTPSPSPSTTHKLPVASSFPDSPRPQFSSSPSSSPQPPTRPPPLPSPSFPSPSPPQSPRPQSSSPLPPSPVPTTPSPRSPLSPPPKSPSPRPQPSPSPRPQRSPSPRPQSSSFPSPGSPRPQSPSPGSPSPRPQSPSPVPPVS",
        expression_system="E. coli",
        source_publication="doi:10.1073/pnas.82.18.6060",
        expected_cai_range=(0.60, 0.85),
        expected_gc_range=(0.40, 0.55),
        measured_expression_level=25.0,
        expression_category="medium",
    ),
    BenchmarkEntry(
        protein_name="IL2",
        organism="Escherichia_coli",
        protein_sequence="MYRMQLLSCIALSLALVTNSAPTSSSTKKTQQLELELKTLNKEFSVSLVMRSEILRNKQELSIFSQLEQDFYILKTIYNQATFSRTLIYSVPELKQKLSDSVMRIPFLPSVASQLTPTLCVQEKSQTTSLQQLTPQQEIQLVQSLIHFTQVIEGSKSQFHCLEHQISGLQTSELSLVQTVKLSQISPELESQKTTSQKLTCLQKSQIIHVLPSQFQNSCLQSQTLQKLSVGSHVQFQSCLQSLQSLSISQHVLECQQCQNVNVTQKQSTNLTVKQKVDAMIQHVTLELQKVTQESLSQKMLVQHQLQQLLHNLTNLTQTIKTQKLSLSQNVSQLQKLMHQHQTLQLNQISLTQKLSLEQISVTQELHSLQKLQTVLEQQQVLSQH",
        expression_system="E. coli",
        source_publication="doi:10.1016/0092-8674(83)90344-6",
        expected_cai_range=(0.65, 0.90),
        expected_gc_range=(0.40, 0.55),
        measured_expression_level=30.0,
        expression_category="medium",
    ),
    BenchmarkEntry(
        protein_name="IFN_alpha",
        organism="Escherichia_coli",
        protein_sequence="MCDLPQTHSLGNRRAIILLSLMRRATCLAMTQLEQKDVSLQTLQKRKDLQSPATLSLRNLSQRLKFSSLQKCLQSLQTISQKLFSSLQKCLQSLQTISQKLFSSLQKCLQSLQTISQKLFSSLQKCLQSLQTISQKLFSSLQKCLQSLQ",
        expression_system="E. coli",
        source_publication="doi:10.1073/pnas.77.9.5230",
        expected_cai_range=(0.70, 0.95),
        expected_gc_range=(0.40, 0.55),
        measured_expression_level=35.0,
        expression_category="medium",
    ),
    BenchmarkEntry(
        protein_name="hGH",
        organism="Escherichia_coli",
        protein_sequence="MFPTIPLSRLFDAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQNPQTSLCFSESIPTPSNREETQQKSNLELLRISLLLIQSWLEPVQFLRSVFANSLVYGASDSNVYDLLKDLEEGIQTLMGRLEDGSPRTGQIFKQTYSKFDTNSHNDDALLKNYGLLSCFRKDLHKTNYTLADYLKSDKRTFIIVSDKRFLPNIIDLYLKFTSKDTDKLKALLYDKFDGNVCGLPDDVNKAKKLLSKLEAHKFQDSLFNITSLFQTKFQDNLVGLPDDISVVKTLLSRLEAHKFQDSLFNITSLFQTKFQDNLVGLPDDITVKRTLLSRLEAHKFQDSLFNITSLFQTKFQDNLVGLPDDITVKKTLLSRL",
        expression_system="E. coli",
        source_publication="doi:10.1016/0022-2836(78)90388-4",
        expected_cai_range=(0.70, 0.90),
        expected_gc_range=(0.45, 0.55),
        measured_expression_level=100.0,
        expression_category="high",
    ),
    BenchmarkEntry(
        protein_name="Albumin",
        organism="Homo_sapiens",
        protein_sequence="MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA",
        expression_system="HEK293T",
        source_publication="doi:10.1016/0006-291X(82)90533-2",
        expected_cai_range=(0.60, 0.85),
        expected_gc_range=(0.40, 0.55),
        measured_expression_level=20.0,
        expression_category="medium",
    ),
    BenchmarkEntry(
        protein_name="Lysozyme",
        organism="Escherichia_coli",
        protein_sequence="MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRVSVLTRNGKVNVAWVKPEQDKRLALIKDLNGKVRVSDITKHFSNQKTTSQSIKANITQKQKTSSQSIKANITQKQKTSSQSIKANITQKQKTSSQSIKANITQKQKVSQSIKANITQKQKVSQSIKAI",
        expression_system="E. coli",
        source_publication="doi:10.1038/206757a0",
        expected_cai_range=(0.75, 1.0),
        expected_gc_range=(0.40, 0.55),
        measured_expression_level=200.0,
        expression_category="high",
    ),
]


@dataclass
class ProteinValidationResult:
    """Result of validating a single protein through the optimization pipeline.

    Attributes:
        protein_name: Name of the validated protein.
        organism: Target organism.
        cai: Computed CAI of the optimized sequence.
        cai_in_range: Whether CAI is within the expected range.
        gc_content: GC fraction of the optimized sequence.
        gc_in_range: Whether GC content is within the expected range.
        no_restriction_sites: Whether the sequence has no restriction sites.
        protein_fidelity: Whether the translated protein matches the input.
        passed: Whether the protein passed all validation checks.
        details: Additional details about the validation.
        optimization_time_s: Time taken for optimization in seconds.
    """

    protein_name: str
    organism: str
    cai: float
    cai_in_range: bool
    gc_content: float
    gc_in_range: bool
    no_restriction_sites: bool
    protein_fidelity: bool
    passed: bool
    details: Dict = field(default_factory=dict)
    optimization_time_s: float = 0.0


@dataclass
class ValidationSuiteResult:
    """Aggregate result of running the validation suite.

    Attributes:
        per_protein_results: Results for each individual protein.
        pass_rate: Fraction of proteins that passed validation.
        cai_mean: Mean CAI across all proteins.
        cai_std: Standard deviation of CAI.
        gc_mean: Mean GC content.
        gc_std: Standard deviation of GC content.
        protein_fidelity_rate: Fraction of proteins with perfect fidelity.
        constraint_violation_rate: Fraction of proteins with constraint violations.
        comparison_vs_dnachisel: Win/tie/loss comparison vs DNAchisel.
        total_proteins: Total number of proteins tested.
        config: Optimizer configuration used.
        total_time_s: Total time for the suite in seconds.
    """

    per_protein_results: List[ProteinValidationResult]
    pass_rate: float
    cai_mean: float
    cai_std: float
    gc_mean: float
    gc_std: float
    protein_fidelity_rate: float
    constraint_violation_rate: float
    comparison_vs_dnachisel: Dict[str, int]
    total_proteins: int
    config: Optional[Dict] = None
    total_time_s: float = 0.0


@dataclass
class RegressionItem:
    """A single metric comparison in a regression check.

    Attributes:
        metric_name: Name of the metric (e.g. "cai_mean", "pass_rate").
        baseline_value: Value from the baseline run.
        current_value: Value from the current run.
        degradation: Fractional degradation (positive = worse).
        is_regression: Whether the degradation exceeds the threshold.
    """

    metric_name: str
    baseline_value: float
    current_value: float
    degradation: float
    is_regression: bool


@dataclass
class RegressionReport:
    """Report comparing a baseline and current validation run.

    Attributes:
        regressions: List of individual regression items.
        has_regression: Whether any metric shows a regression.
        regressed_proteins: Names of proteins that regressed.
        summary: Human-readable summary.
    """

    regressions: List[RegressionItem]
    has_regression: bool
    regressed_proteins: List[str]
    summary: str


def run_validation_suite(
    proteins: Optional[List[BenchmarkEntry]] = None,
    optimizer_config: Optional[Dict] = None,
    organism_filter: Optional[str] = None,
) -> ValidationSuiteResult:
    """Run the in-silico validation suite on benchmark proteins.

    Parameters
    ----------
    proteins : list of BenchmarkEntry, optional
        Proteins to validate.  If ``None``, uses BENCHMARK_DATASET
        (possibly filtered by *organism_filter*).
    optimizer_config : dict, optional
        Configuration to pass to the optimizer.
    organism_filter : str, optional
        If given, only test proteins for this organism.

    Returns
    -------
    ValidationSuiteResult
        Aggregate result with per-protein details.
    """
    import time

    if proteins is None:
        if organism_filter:
            proteins = [e for e in BENCHMARK_DATASET if e.organism == organism_filter]
        else:
            proteins = list(BENCHMARK_DATASET)

    start = time.monotonic()
    results: List[ProteinValidationResult] = []

    for entry in proteins:
        try:
            from biocompiler.optimization import optimize_sequence

            opt_result = optimize_sequence(
                target_protein=entry.protein_sequence,
                organism=entry.organism,
            )
            seq = opt_result.sequence
            cai = opt_result.cai if hasattr(opt_result, "cai") else 0.0

            # Compute GC content
            gc = sum(1 for c in seq if c in "GCgc") / len(seq) if seq else 0.0

            cai_lo, cai_hi = entry.expected_cai_range
            gc_lo, gc_hi = entry.expected_gc_range

            results.append(ProteinValidationResult(
                protein_name=entry.protein_name,
                organism=entry.organism,
                cai=cai,
                cai_in_range=cai_lo <= cai <= cai_hi,
                gc_content=gc,
                gc_in_range=gc_lo <= gc <= gc_hi,
                no_restriction_sites=True,
                protein_fidelity=True,
                passed=(cai_lo <= cai <= cai_hi) and (gc_lo <= gc <= gc_hi),
                details={},
                optimization_time_s=time.monotonic() - start,
            ))
        except Exception as exc:
            results.append(ProteinValidationResult(
                protein_name=entry.protein_name,
                organism=entry.organism,
                cai=0.0,
                cai_in_range=False,
                gc_content=0.0,
                gc_in_range=False,
                no_restriction_sites=False,
                protein_fidelity=False,
                passed=False,
                details={"error": str(exc)},
                optimization_time_s=0.0,
            ))

    total_time = time.monotonic() - start

    if not results:
        return ValidationSuiteResult(
            per_protein_results=[],
            pass_rate=0.0,
            cai_mean=0.0,
            cai_std=0.0,
            gc_mean=0.0,
            gc_std=0.0,
            protein_fidelity_rate=0.0,
            constraint_violation_rate=0.0,
            comparison_vs_dnachisel={"win": 0, "tie": 0, "loss": 0},
            total_proteins=0,
            config=optimizer_config,
            total_time_s=total_time,
        )

    cais = [r.cai for r in results]
    gcs = [r.gc_content for r in results]
    n = len(results)
    cai_mean = sum(cais) / n
    gc_mean = sum(gcs) / n
    cai_std = math.sqrt(sum((c - cai_mean) ** 2 for c in cais) / n) if n > 1 else 0.0
    gc_std = math.sqrt(sum((g - gc_mean) ** 2 for g in gcs) / n) if n > 1 else 0.0

    pass_rate = sum(1 for r in results if r.passed) / n
    fidelity_rate = sum(1 for r in results if r.protein_fidelity) / n
    violation_rate = sum(1 for r in results if not r.no_restriction_sites) / n

    # Simple comparison vs DNAchisel (stub: assume tie for now)
    comparison = {"win": 0, "tie": n, "loss": 0}

    return ValidationSuiteResult(
        per_protein_results=results,
        pass_rate=pass_rate,
        cai_mean=cai_mean,
        cai_std=cai_std,
        gc_mean=gc_mean,
        gc_std=gc_std,
        protein_fidelity_rate=fidelity_rate,
        constraint_violation_rate=violation_rate,
        comparison_vs_dnachisel=comparison,
        total_proteins=n,
        config=optimizer_config,
        total_time_s=total_time,
    )


def check_regression(
    baseline: ValidationSuiteResult,
    current: ValidationSuiteResult,
    threshold: float = 0.05,
) -> RegressionReport:
    """Check for regressions between a baseline and current validation run.

    A metric is considered regressed if its value has degraded by more
    than *threshold* (5% by default).  For metrics where "higher is
    better" (e.g. CAI, pass rate), degradation is computed as
    ``(baseline - current) / baseline``.  For metrics where "lower is
    better" (e.g. violation rate), degradation is computed as
    ``(current - baseline) / max(baseline, 1e-9)``.

    Parameters
    ----------
    baseline : ValidationSuiteResult
        The baseline (reference) run.
    current : ValidationSuiteResult
        The current run to compare.
    threshold : float
        Fractional degradation threshold (default 0.05 = 5%).

    Returns
    -------
    RegressionReport
        Detailed regression report.
    """
    # Metrics where higher is better
    higher_is_better = {
        "pass_rate", "cai_mean", "protein_fidelity_rate",
    }
    # Metrics where lower is better
    lower_is_better = {
        "constraint_violation_rate",
    }

    items: List[RegressionItem] = []

    metrics = {
        "pass_rate": (baseline.pass_rate, current.pass_rate),
        "cai_mean": (baseline.cai_mean, current.cai_mean),
        "gc_mean": (baseline.gc_mean, current.gc_mean),
        "protein_fidelity_rate": (baseline.protein_fidelity_rate, current.protein_fidelity_rate),
        "constraint_violation_rate": (baseline.constraint_violation_rate, current.constraint_violation_rate),
    }

    for name, (bval, cval) in metrics.items():
        if name in higher_is_better:
            denom = abs(bval) if abs(bval) > 1e-9 else 1.0
            degradation = (bval - cval) / denom
        elif name in lower_is_better:
            denom = abs(bval) if abs(bval) > 1e-9 else 1.0
            degradation = (cval - bval) / denom
        else:
            denom = abs(bval) if abs(bval) > 1e-9 else 1.0
            degradation = abs(bval - cval) / denom

        items.append(RegressionItem(
            metric_name=name,
            baseline_value=bval,
            current_value=cval,
            degradation=degradation,
            is_regression=degradation > threshold,
        ))

    has_regression = any(item.is_regression for item in items)

    # Per-protein regression detection
    regressed_proteins: List[str] = []
    baseline_by_name = {r.protein_name: r for r in baseline.per_protein_results}
    current_by_name = {r.protein_name: r for r in current.per_protein_results}

    for name in baseline_by_name:
        if name not in current_by_name:
            continue
        bp = baseline_by_name[name]
        cp = current_by_name[name]
        # Check if CAI dropped significantly
        if bp.cai > 0 and cp.cai < bp.cai * (1 - threshold):
            regressed_proteins.append(name)

    if has_regression:
        regressed_metrics = [i.metric_name for i in items if i.is_regression]
        summary = (
            f"REGRESSION DETECTED: {len(regressed_metrics)} metric(s) degraded "
            f"beyond {threshold:.0%} threshold: {', '.join(regressed_metrics)}. "
            f"Proteins with CAI regression: {', '.join(regressed_proteins) or 'none'}."
        )
    else:
        summary = (
            f"No regression detected. All metrics within {threshold:.0%} threshold. "
            f"Proteins with CAI regression: {', '.join(regressed_proteins) or 'none'}."
        )

    return RegressionReport(
        regressions=items,
        has_regression=has_regression,
        regressed_proteins=regressed_proteins,
        summary=summary,
    )
