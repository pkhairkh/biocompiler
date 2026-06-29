"""
<<<<<<<< HEAD:papers/demonstrations/shared/diagnostics.py
Shared utilities for BioCompiler paper demonstrations.
========
Shared utilities for BioCompiler demonstrations.
>>>>>>>> 9c54ec0 (fix: comprehensive repo-wide audit — fix 200+ issues across all files):demonstrations/shared/diagnostics.py
Uses the ACTUAL BioCompiler API — 12 canonical predicates from
evaluate_all_predicates plus biosecurity, immunogenicity, and
primer compatibility checks.
"""

import json
import time
from pathlib import Path
from dataclasses import dataclass, asdict

# BioCompiler imports
from biocompiler.optimizer import optimize_sequence
from biocompiler.type_system.predicates import evaluate_all_predicates
from biocompiler.provenance.certificate import (
    generate_certificate,
    compute_certificate,
    format_certificate,
    CertLevel,
)
from biocompiler.biosecurity import screen_hazardous_sequence, check_biosecurity_before_optimize
from biocompiler.immunogenicity import compute_immunogenicity, predict_t_cell_epitopes
from biocompiler.type_system import SLOTMode

DATA_DIR = Path(__file__).parent.parent / "data"

# The 12 canonical predicates from evaluate_all_predicates
CANONICAL_PREDICATES = [
    "NoCrypticSplice",
    "SpliceCorrect",
    "GCInRange",
    "CodonAdapted",
    "NoRestrictionSite",
    "InFrame",
    "NoInstabilityMotif",
    "NoCpGIsland",
    "NoCrypticPromoter",
    "NoUnexpectedTMDomain",
    "mRNASecondaryStructure",
    "CoTranslationalFolding",
]

# Extended diagnostic layers for the paper (12 canonical + 3 extended = 15)
# The paper narrative presents these as "diagnostic layers"
DIAGNOSTIC_LAYERS = CANONICAL_PREDICATES + [
    "BiosecurityScreening",
    "Immunogenicity",
    "PrimerCompatibility",
]


@dataclass
class DiagnosticResult:
    """Full diagnostic of a protein sequence through BioCompiler."""
    protein_name: str
    protein: str
    organism: str

    # Optimization result
    optimized_sequence: str = ""
    gc_content: float = 0.0
    cai: float = 0.0
    optimization_time_s: float = 0.0
    certificate_level: str = ""
    fallback_used: bool = False

    # 12 canonical predicate results
    canonical_predicates: list[dict] = None  # [{predicate, verdict, violation}]

    # Extended diagnostics
    biosecurity_passed: bool = False
    biosecurity_risk: str = ""
    immunogenicity_score: float = 0.0
    immunogenicity_risk: str = ""
    t_cell_epitope_count: int = 0
    primer_compatibility: bool = False

    # Summary
    all_predicates_pass: bool = False
    green_count: int = 0
    red_count: int = 0
    uncertain_count: int = 0

    def __post_init__(self):
        if self.canonical_predicates is None:
            self.canonical_predicates = []

    def to_dict(self) -> dict:
        return asdict(self)

    def save(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, path: Path) -> "DiagnosticResult":
        with open(path) as f:
            return cls(**json.load(f))


def verdict_to_color(verdict_str: str) -> str:
    """Map verdict to traffic-light color."""
    v = verdict_str.upper()
    if v in ("PASS", "LIKELY_PASS"):
        return "green"
    elif v in ("FAIL", "LIKELY_FAIL"):
        return "red"
    else:
        return "amber"


def run_full_diagnostic(
    protein: str,
    protein_name: str = "unknown",
    organism: str = "Homo_sapiens",
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
    cai_threshold: float = 0.5,
    enzymes: list = None,
    slot_mode: SLOTMode = SLOTMode.VERIFIED,
    biosecurity_mode: str = "warn",
    verbose: bool = False,
) -> DiagnosticResult:
    """Run a complete BioCompiler diagnostic on a protein sequence.

    This evaluates all 12 canonical predicates plus biosecurity,
    immunogenicity, and primer compatibility.
    """
    if enzymes is None:
        enzymes = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]

    result = DiagnosticResult(
        protein_name=protein_name,
        protein=protein,
        organism=organism,
    )

    # Step 1: Optimize
    if verbose:
        print(f"  Optimizing {protein_name} ({len(protein)} aa) for {organism}...")

    t0 = time.time()
    try:
        opt = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            cai_threshold=cai_threshold,
            enzymes=enzymes,
            strategy="hybrid",
            strict_mode=False,
            biosecurity_mode=biosecurity_mode,
            track_provenance=True,
        )
        result.optimized_sequence = opt.sequence
        result.gc_content = opt.gc_content
        result.cai = opt.cai
        result.fallback_used = opt.fallback_used
    except Exception as e:
        if verbose:
            print(f"  Optimization failed: {e}")
        result.optimization_time_s = time.time() - t0
        return result

    result.optimization_time_s = time.time() - t0

    # Step 2: Evaluate 12 canonical predicates
    if verbose:
        print(f"  Evaluating 12 canonical predicates...")

    try:
        preds = evaluate_all_predicates(
            seq=result.optimized_sequence,
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
            cai_threshold=cai_threshold,
            enzymes=enzymes,
            slot_mode=slot_mode,
        )

        for p in preds:
            entry = {
                "predicate": p.predicate,
                "verdict": p.verdict.value if hasattr(p.verdict, 'value') else str(p.verdict),
                "violation": p.violation,
            }
            result.canonical_predicates.append(entry)

            color = verdict_to_color(entry["verdict"])
            if color == "green":
                result.green_count += 1
            elif color == "red":
                result.red_count += 1
            else:
                result.uncertain_count += 1
    except Exception as e:
        if verbose:
            print(f"  Predicate evaluation failed: {e}")

    # Step 3: Certificate level
    try:
        cert_level = compute_certificate(
            results=opt.predicate_results if hasattr(opt, 'predicate_results') else [],
            slot_mode=slot_mode,
        )
        result.certificate_level = cert_level.value if hasattr(cert_level, 'value') else str(cert_level)
    except Exception:
        result.certificate_level = "UNKNOWN"

    # Step 4: Biosecurity screening
    if verbose:
        print(f"  Running biosecurity screening...")
    try:
        biosec = screen_hazardous_sequence(
            protein=protein,
            dna=result.optimized_sequence,
        )
        result.biosecurity_passed = not biosec.is_hazardous
        result.biosecurity_risk = str(biosec.risk_level) if hasattr(biosec, 'risk_level') else "unknown"
    except Exception as e:
        if verbose:
            print(f"  Biosecurity check failed: {e}")
        result.biosecurity_passed = False  # Fail-closed: cannot verify safety without screening
        result.biosecurity_risk = f"screening_unavailable: {e}"

    # Step 5: Immunogenicity
    if verbose:
        print(f"  Computing immunogenicity...")
    try:
        immuno = compute_immunogenicity(
            protein=protein,
            organism=organism,
        )
        result.immunogenicity_score = immuno.immunogenicity_score
        result.immunogenicity_risk = immuno.risk_class

        epitopes = predict_t_cell_epitopes(
            protein=protein,
        )
        result.t_cell_epitope_count = len(epitopes) if epitopes else 0
    except Exception as e:
        if verbose:
            print(f"  Immunogenicity check failed: {e}")

    # Step 6: Primer compatibility (basic GC check for PCR)
    gc = result.gc_content
    result.primer_compatibility = (0.35 <= gc <= 0.65)  # Standard PCR GC range

    # Summary
    # Count biosecurity and immunogenicity as extended predicates
    if result.biosecurity_passed:
        result.green_count += 1
    else:
        result.red_count += 1

    if result.immunogenicity_risk in ("low", "moderate"):
        result.green_count += 1
    elif result.immunogenicity_risk == "high":
        result.red_count += 1
    else:
        result.uncertain_count += 1

    if result.primer_compatibility:
        result.green_count += 1
    else:
        result.red_count += 1

    result.all_predicates_pass = (result.red_count == 0)

    if verbose:
        total = result.green_count + result.red_count + result.uncertain_count
        print(f"  Result: {result.green_count}/{total} green, {result.red_count} red, {result.uncertain_count} amber")
        print(f"  Certificate: {result.certificate_level}")

    return result


def print_diagnostic_summary(result: DiagnosticResult):
    """Print a formatted diagnostic summary."""
    print(f"\n{'='*70}")
    print(f"  DIAGNOSTIC: {result.protein_name}")
    print(f"  Protein: {len(result.protein)} aa")
    print(f"  Organism: {result.organism}")
    print(f"{'='*70}")

    print(f"\n  Optimization: {result.optimization_time_s:.2f}s")
    print(f"  GC Content: {result.gc_content:.1%}")
    print(f"  CAI: {result.cai:.3f}")
    print(f"  Certificate: {result.certificate_level}")
    print(f"  Fallback: {'Yes' if result.fallback_used else 'No'}")

    print(f"\n  --- 12 Canonical Predicates ---")
    for p in result.canonical_predicates:
        color = verdict_to_color(p["verdict"])
        icon = {"green": "[OK]", "red": "[FAIL]", "amber": "[?]"}[color]
        violation = f" ({p['violation']})" if p.get("violation") else ""
        print(f"  {icon} {p['predicate']:30s} {p['verdict']:15s}{violation}")

    print(f"\n  --- Extended Diagnostics ---")
    icon = "[OK]" if result.biosecurity_passed else "[FAIL]"
    print(f"  {icon} {'Biosecurity Screening':30s} {'PASS' if result.biosecurity_passed else 'FAIL':15s} risk={result.biosecurity_risk}")

    icon = "[OK]" if result.immunogenicity_risk in ("low", "moderate") else "[FAIL]"
    print(f"  {icon} {'Immunogenicity':30s} {result.immunogenicity_risk:15s} score={result.immunogenicity_score:.3f}")

    icon = "[OK]" if result.primer_compatibility else "[FAIL]"
    print(f"  {icon} {'Primer Compatibility':30s} {'PASS' if result.primer_compatibility else 'FAIL':15s}")

    total = result.green_count + result.red_count + result.uncertain_count
    print(f"\n  Total: {result.green_count}/{total} green  |  {result.red_count} red  |  {result.uncertain_count} amber")
    print(f"  All-pass: {'YES' if result.all_predicates_pass else 'NO'}")
