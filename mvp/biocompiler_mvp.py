"""
BioCompiler MVP — Integrated Pipeline

Combines MaxEntScan splice site scoring and CSP sequence optimization into
a complete end-to-end pipeline with a three-valued (PASS/FAIL/UNCERTAIN)
type system and machine-verifiable certificates.

Pipeline:
  1. CSP optimization: generate a DNA sequence for target protein
  2. MaxEntScan: score all splice sites in the sequence
  3. Type checking: evaluate all 7 type predicates
  4. Certificate generation: if all PASS, generate certificate
  5. Certificate verification: independently re-check
"""

import json
import hashlib
import time
import math
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass, field, asdict
from enum import Enum

from maxentscan import (
    score_donor,
    score_acceptor,
    scan_splice_sites,
    max_donor_score,
    max_acceptor_score,
)
from csp_optimizer import (
    optimize_sequence,
    compute_cai,
    compute_gc_content,
    reverse_complement,
    COMMON_RESTRICTION_SITES,
    CODON_USAGE_TABLES,
)


# ---------------------------------------------------------------------------
# Three-valued type system
# ---------------------------------------------------------------------------
class Verdict(Enum):
    """Three-valued logic for type checking results."""
    PASS = "PASS"
    FAIL = "FAIL"
    UNCERTAIN = "UNCERTAIN"


@dataclass
class TypeCheckResult:
    """Result of a single type predicate check."""
    predicate: str
    verdict: Verdict
    message: str
    details: Dict = field(default_factory=dict)


@dataclass
class Certificate:
    """Machine-verifiable certificate that a sequence passes all type predicates."""
    sequence_hash: str
    protein: str
    organism: str
    sequence_length: int
    gc_content: float
    cai: float
    max_donor_score: float
    max_acceptor_score: float
    predicate_results: List[Dict]
    overall_verdict: str
    timestamp: float
    pipeline_version: str = "biocompiler-mvp-0.1.0"

    def to_json(self) -> str:
        """Serialize certificate to JSON."""
        return json.dumps(asdict(self), indent=2)

    @classmethod
    def from_json(cls, json_str: str) -> "Certificate":
        """Deserialize certificate from JSON."""
        data = json.loads(json_str)
        return cls(**data)


# ---------------------------------------------------------------------------
# Type predicates (the 7 predicates from the BioCompiler specification)
# ---------------------------------------------------------------------------

# Configuration constants
DONOR_THRESHOLD = 3.0
ACCEPTOR_THRESHOLD = 3.0
GC_LO = 0.30
GC_HI = 0.70
CAI_THRESHOLD = 0.2
MAX_CONSECUTIVE_T = 5
INSTABILITY_MOTIF = "ATTTA"


def check_no_cryptic_splice(sequence: str) -> TypeCheckResult:
    """
    NoCrypticSplice: no cryptic donor or acceptor sites above threshold.

    Uses MaxEntScan to score all potential splice sites. A site is 'cryptic'
    if its score exceeds the threshold, suggesting it could be incorrectly
    recognized by the splicing machinery.
    """
    seq = sequence.upper()
    max_d = max_donor_score(seq)
    max_a = max_acceptor_score(seq)

    sites = scan_splice_sites(seq, DONOR_THRESHOLD, ACCEPTOR_THRESHOLD)

    if not sites:
        return TypeCheckResult(
            predicate="NoCrypticSplice",
            verdict=Verdict.PASS,
            message="No cryptic splice sites detected",
            details={"max_donor_score": max_d, "max_acceptor_score": max_a},
        )
    else:
        donor_sites = [s for s in sites if s[1] == "donor"]
        acceptor_sites = [s for s in sites if s[1] == "acceptor"]
        msg_parts = []
        if donor_sites:
            msg_parts.append(f"{len(donor_sites)} donor site(s)")
        if acceptor_sites:
            msg_parts.append(f"{len(acceptor_sites)} acceptor site(s)")

        # If only weak sites (score < 6), mark UNCERTAIN; strong sites = FAIL
        max_score = max(s[2] for s in sites)
        if max_score < 6.0:
            verdict = Verdict.UNCERTAIN
        else:
            verdict = Verdict.FAIL

        return TypeCheckResult(
            predicate="NoCrypticSplice",
            verdict=verdict,
            message=f"Found {', '.join(msg_parts)} above threshold",
            details={
                "max_donor_score": max_d,
                "max_acceptor_score": max_a,
                "donor_sites": [(s[0], round(s[2], 2)) for s in donor_sites],
                "acceptor_sites": [(s[0], round(s[2], 2)) for s in acceptor_sites],
            },
        )


def check_gc_in_range(sequence: str) -> TypeCheckResult:
    """GCInRange: GC content within acceptable range [0.30, 0.70]."""
    gc = compute_gc_content(sequence)
    if GC_LO <= gc <= GC_HI:
        return TypeCheckResult(
            predicate="GCInRange",
            verdict=Verdict.PASS,
            message=f"GC content {gc:.1%} within [{GC_LO:.0%}, {GC_HI:.0%}]",
            details={"gc_content": gc, "lo": GC_LO, "hi": GC_HI},
        )
    else:
        return TypeCheckResult(
            predicate="GCInRange",
            verdict=Verdict.FAIL,
            message=f"GC content {gc:.1%} outside [{GC_LO:.0%}, {GC_HI:.0%}]",
            details={"gc_content": gc, "lo": GC_LO, "hi": GC_HI},
        )


def check_codon_adapted(sequence: str, organism: str = "human") -> TypeCheckResult:
    """CodonAdapted: CAI >= threshold for the given organism."""
    cai = compute_cai(sequence, organism)
    if cai >= CAI_THRESHOLD:
        return TypeCheckResult(
            predicate="CodonAdapted",
            verdict=Verdict.PASS,
            message=f"CAI={cai:.4f} >= {CAI_THRESHOLD}",
            details={"cai": cai, "threshold": CAI_THRESHOLD, "organism": organism},
        )
    else:
        return TypeCheckResult(
            predicate="CodonAdapted",
            verdict=Verdict.FAIL,
            message=f"CAI={cai:.4f} < {CAI_THRESHOLD}",
            details={"cai": cai, "threshold": CAI_THRESHOLD, "organism": organism},
        )


def check_no_restriction_site(sequence: str) -> TypeCheckResult:
    """NoRestrictionSite: no common restriction enzyme recognition sites."""
    seq = sequence.upper()
    found_sites: List[Dict] = []

    for name, site in COMMON_RESTRICTION_SITES.items():
        site_upper = site.upper()
        site_rc = reverse_complement(site_upper)
        if site_upper in seq:
            pos = seq.find(site_upper)
            found_sites.append({"enzyme": name, "site": site_upper, "position": pos})
        if site_rc in seq:
            pos = seq.find(site_rc)
            found_sites.append({"enzyme": name, "site": f"{site_rc}(RC)", "position": pos})

    if not found_sites:
        return TypeCheckResult(
            predicate="NoRestrictionSite",
            verdict=Verdict.PASS,
            message="No restriction enzyme sites found",
            details={"checked_enzymes": list(COMMON_RESTRICTION_SITES.keys())},
        )
    else:
        return TypeCheckResult(
            predicate="NoRestrictionSite",
            verdict=Verdict.FAIL,
            message=f"Found {len(found_sites)} restriction site(s)",
            details={"found": found_sites},
        )


def check_in_frame(sequence: str, target_protein: str) -> TypeCheckResult:
    """InFrame: no stop codons in reading frame, correct translation."""
    seq = sequence.upper()
    target = target_protein.upper()

    # Check length
    if len(seq) % 3 != 0:
        return TypeCheckResult(
            predicate="InFrame",
            verdict=Verdict.FAIL,
            message=f"Sequence length {len(seq)} not a multiple of 3",
            details={},
        )

    # Check for stop codons
    from csp_optimizer import CODON_TABLE
    stop_positions = []
    translated = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, "?")
        if aa == "*":
            stop_positions.append(i)
        translated.append(aa)

    # Check translation matches target
    translated_protein = "".join(translated)
    matches = translated_protein == target

    if stop_positions:
        return TypeCheckResult(
            predicate="InFrame",
            verdict=Verdict.FAIL,
            message=f"Found {len(stop_positions)} stop codon(s) in frame",
            details={"stop_positions": stop_positions},
        )
    elif not matches:
        return TypeCheckResult(
            predicate="InFrame",
            verdict=Verdict.FAIL,
            message="Translation does not match target protein",
            details={"expected": target[:20] + "...", "got": translated_protein[:20] + "..."},
        )
    else:
        return TypeCheckResult(
            predicate="InFrame",
            verdict=Verdict.PASS,
            message="No stop codons in frame, correct translation",
            details={"length_aa": len(target)},
        )


def check_no_instability_motif(sequence: str) -> TypeCheckResult:
    """NoInstabilityMotif: no ATTTA and no 6+ consecutive T."""
    seq = sequence.upper()
    issues = []

    # Check for ATTTA
    atta_positions = []
    pos = seq.find(INSTABILITY_MOTIF)
    while pos != -1:
        atta_positions.append(pos)
        pos = seq.find(INSTABILITY_MOTIF, pos + 1)

    if atta_positions:
        issues.append(f"ATTTA at position(s) {atta_positions}")

    # Check for 6+ consecutive T
    long_t_positions = []
    i = 0
    while i < len(seq):
        if seq[i] == "T":
            j = i
            while j < len(seq) and seq[j] == "T":
                j += 1
            if j - i >= 6:
                long_t_positions.append((i, j - i))
            i = j
        else:
            i += 1

    if long_t_positions:
        issues.append(f"6+ consecutive T at: {long_t_positions}")

    if not issues:
        return TypeCheckResult(
            predicate="NoInstabilityMotif",
            verdict=Verdict.PASS,
            message="No instability motifs found",
            details={},
        )
    else:
        return TypeCheckResult(
            predicate="NoInstabilityMotif",
            verdict=Verdict.FAIL,
            message="; ".join(issues),
            details={"atta_positions": atta_positions, "long_t_runs": long_t_positions},
        )


def check_sequence_integrity(sequence: str, target_protein: str) -> TypeCheckResult:
    """SequenceIntegrity: sequence contains only valid DNA bases and encodes the protein."""
    seq = sequence.upper()

    # Check valid bases
    invalid_chars = set(seq) - set("ACGT")
    if invalid_chars:
        return TypeCheckResult(
            predicate="SequenceIntegrity",
            verdict=Verdict.FAIL,
            message=f"Invalid characters: {invalid_chars}",
            details={"invalid_chars": list(invalid_chars)},
        )

    # Check length matches protein
    expected_len = len(target_protein) * 3
    if len(seq) != expected_len:
        return TypeCheckResult(
            predicate="SequenceIntegrity",
            verdict=Verdict.FAIL,
            message=f"Length mismatch: expected {expected_len}, got {len(seq)}",
            details={},
        )

    return TypeCheckResult(
        predicate="SequenceIntegrity",
        verdict=Verdict.PASS,
        message="Sequence integrity verified",
        details={"length": len(seq), "valid_bases": True},
    )


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def run_type_checks(
    sequence: str,
    target_protein: str,
    organism: str = "human",
) -> List[TypeCheckResult]:
    """Run all 7 type predicates against the sequence."""
    results = [
        check_sequence_integrity(sequence, target_protein),
        check_in_frame(sequence, target_protein),
        check_no_cryptic_splice(sequence),
        check_gc_in_range(sequence),
        check_codon_adapted(sequence, organism),
        check_no_restriction_site(sequence),
        check_no_instability_motif(sequence),
    ]
    return results


def generate_certificate(
    sequence: str,
    target_protein: str,
    organism: str,
    results: List[TypeCheckResult],
) -> Certificate:
    """Generate a certificate from the type check results."""
    seq_hash = hashlib.sha256(sequence.encode()).hexdigest()
    gc = compute_gc_content(sequence)
    cai = compute_cai(sequence, organism)
    max_d = max_donor_score(sequence)
    max_a = max_acceptor_score(sequence)

    # Overall verdict: PASS only if all predicates PASS
    # UNCERTAIN if any is UNCERTAIN (and none FAIL)
    # FAIL if any FAIL
    verdicts = [r.verdict for r in results]
    if Verdict.FAIL in verdicts:
        overall = "FAIL"
    elif Verdict.UNCERTAIN in verdicts:
        overall = "UNCERTAIN"
    else:
        overall = "PASS"

    return Certificate(
        sequence_hash=seq_hash,
        protein=target_protein,
        organism=organism,
        sequence_length=len(sequence),
        gc_content=gc,
        cai=cai,
        max_donor_score=max_d,
        max_acceptor_score=max_a,
        predicate_results=[
            {
                "predicate": r.predicate,
                "verdict": r.verdict.value,
                "message": r.message,
            }
            for r in results
        ],
        overall_verdict=overall,
        timestamp=time.time(),
    )


def verify_certificate(
    certificate: Certificate,
    sequence: str,
    target_protein: str,
    organism: str,
) -> bool:
    """
    Independently verify a certificate by re-running all type checks
    and confirming the results match.
    """
    # Verify hash
    expected_hash = hashlib.sha256(sequence.encode()).hexdigest()
    if certificate.sequence_hash != expected_hash:
        return False

    # Re-run type checks
    results = run_type_checks(sequence, target_protein, organism)

    # Verify each predicate result matches
    if len(results) != len(certificate.predicate_results):
        return False

    for result, cert_entry in zip(results, certificate.predicate_results):
        if result.predicate != cert_entry["predicate"]:
            return False
        if result.verdict.value != cert_entry["verdict"]:
            return False

    # Verify overall verdict
    verdicts = [r.verdict for r in results]
    if Verdict.FAIL in verdicts:
        expected_overall = "FAIL"
    elif Verdict.UNCERTAIN in verdicts:
        expected_overall = "UNCERTAIN"
    else:
        expected_overall = "PASS"

    if certificate.overall_verdict != expected_overall:
        return False

    return True


def run_pipeline(
    target_protein: str,
    organism: str = "human",
    protein_name: str = "unnamed",
    save_certificate: bool = True,
) -> Tuple[str, Certificate, bool]:
    """
    Run the complete BioCompiler pipeline.

    Steps:
      1. CSP optimization: generate a DNA sequence
      2. MaxEntScan: score splice sites
      3. Type checking: evaluate all 7 predicates
      4. Certificate generation
      5. Certificate verification

    Returns:
        (sequence, certificate, verification_passed)
    """
    print("=" * 72)
    print(f"BioCompiler MVP Pipeline — {protein_name}")
    print(f"Target protein: {len(target_protein)} amino acids")
    print(f"Organism: {organism}")
    print("=" * 72)

    # Step 1: CSP Optimization
    print("\n[Step 1] CSP Sequence Optimization (z3)")
    print("-" * 40)
    t0 = time.time()
    result = optimize_sequence(
        target_protein=target_protein,
        organism=organism,
        gc_lo=GC_LO,
        gc_hi=GC_HI,
        restriction_sites=list(COMMON_RESTRICTION_SITES.values()),
        cai_threshold=CAI_THRESHOLD,
    )
    t1 = time.time()
    sequence = result.sequence
    print(f"  Sequence length: {len(sequence)} bp")
    print(f"  CAI: {result.cai:.4f}")
    print(f"  GC content: {result.gc_content:.1%}")
    print(f"  Optimization time: {t1-t0:.2f}s")
    print(f"  Fallback solver used: {result.fallback_used}")

    if result.unsat_core:
        print(f"  UNSAT core: {result.unsat_core}")

    # Step 2: MaxEntScan Scoring
    print("\n[Step 2] MaxEntScan Splice Site Scoring")
    print("-" * 40)
    sites = scan_splice_sites(sequence, DONOR_THRESHOLD, ACCEPTOR_THRESHOLD)
    max_d = max_donor_score(sequence)
    max_a = max_acceptor_score(sequence)
    print(f"  Max donor score: {max_d:.4f} (threshold: {DONOR_THRESHOLD})")
    print(f"  Max acceptor score: {max_a:.4f} (threshold: {ACCEPTOR_THRESHOLD})")
    print(f"  Sites above threshold: {len(sites)}")
    if sites:
        for pos, stype, score in sites[:10]:  # Show first 10
            print(f"    pos {pos}: {stype} score={score:.2f}")
        if len(sites) > 10:
            print(f"    ... and {len(sites)-10} more")

    # Step 3: Type Checking
    print("\n[Step 3] Type Checking (7 Predicates)")
    print("-" * 40)
    type_results = run_type_checks(sequence, target_protein, organism)
    for r in type_results:
        symbol = {"PASS": "✓", "FAIL": "✗", "UNCERTAIN": "?"}[r.verdict.value]
        print(f"  [{symbol}] {r.predicate}: {r.message}")

    # Step 4: Certificate Generation
    print("\n[Step 4] Certificate Generation")
    print("-" * 40)
    certificate = generate_certificate(sequence, target_protein, organism, type_results)
    print(f"  Overall verdict: {certificate.overall_verdict}")
    print(f"  Sequence hash: {certificate.sequence_hash[:16]}...")
    print(f"  Pipeline version: {certificate.pipeline_version}")

    if save_certificate:
        import os
        cert_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificates")
        os.makedirs(cert_dir, exist_ok=True)
        cert_filename = os.path.join(cert_dir, f"certificate_{protein_name.lower().replace(' ', '_')}.json")
        with open(cert_filename, "w") as f:
            f.write(certificate.to_json())
        print(f"  Certificate saved: {cert_filename}")

    # Step 5: Certificate Verification
    print("\n[Step 5] Certificate Verification")
    print("-" * 40)
    verified = verify_certificate(certificate, sequence, target_protein, organism)
    print(f"  Verification: {'PASSED' if verified else 'FAILED'}")

    # Summary
    print("\n" + "=" * 72)
    print("PIPELINE SUMMARY")
    print("=" * 72)
    print(f"  Protein: {protein_name} ({len(target_protein)} aa)")
    print(f"  DNA sequence: {len(sequence)} bp")
    print(f"  CAI: {result.cai:.4f}")
    print(f"  GC: {result.gc_content:.1%}")
    print(f"  Max donor score: {max_d:.2f}")
    print(f"  Max acceptor score: {max_a:.2f}")
    print(f"  Overall verdict: {certificate.overall_verdict}")
    print(f"  Certificate verified: {verified}")
    print(f"  Total time: {time.time()-t0:.2f}s")

    # Print first/last 60bp of sequence
    if len(sequence) > 120:
        print(f"\n  Sequence (first 60bp): {sequence[:60]}")
        print(f"  Sequence (last 60bp):  {sequence[-60:]}")
    else:
        print(f"\n  Sequence: {sequence}")

    print()

    return sequence, certificate, verified


# ---------------------------------------------------------------------------
# Target proteins
# ---------------------------------------------------------------------------
HUMAN_INSULIN = "FVNQHLCGSHLVEALYLVCGERGFFYTPKTGIVEQCCTSICSLYQLENYCN"
# Mature insulin: B-chain (30 aa) + A-chain (21 aa) = 51 aa

EGFP = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTY"
    "GVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKED"
    "GNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHY"
    "LSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)


if __name__ == "__main__":
    # Run both demos
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║                    BioCompiler MVP Demo                             ║")
    print("║  MaxEntScan + CSP Optimization with Three-Valued Type System       ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    print()

    # Demo 1: Human Insulin
    run_pipeline(
        target_protein=HUMAN_INSULIN,
        organism="human",
        protein_name="Human_Insulin",
    )

    print("\n\n")

    # Demo 2: eGFP
    run_pipeline(
        target_protein=EGFP,
        organism="human",
        protein_name="eGFP",
    )

    print("\nDemo complete. Certificates saved as JSON files.")
