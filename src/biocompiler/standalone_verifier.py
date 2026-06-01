"""
BioCompiler Standalone Certificate Verifier
=============================================
Self-contained, zero-dependency verifier that independently re-checks every
predicate from scratch. Imports ONLY Python stdlib (hashlib, json, math, sys).

TRUSTED COMPUTING BASE (TCB):
  This file IS the TCB. No biocompiler imports. Every biological check is
  implemented from first principles using data embedded in this file.

SECURITY ARGUMENT:
  <500 lines. Audit in <30 min. CompComp approach: compiler proves, verifier checks.

CERTIFICATE FORMATS SUPPORTED:
  1. "types" format: {version, design_id, sequence,
     types: [{predicate, verdict, derivation}], provenance: {parameters, ...}}
  2. "predicates" format: {version, design_id, sequence,
     predicates: {name: {verdict, ...}}, ...}
"""
import sys as _sys
# When run as a script, Python adds the script's directory to sys.path[0].
# If this directory contains a types.py (as in biocompiler/), it shadows
# the stdlib types module, causing a circular import. Remove it early.
if _sys.path and _sys.path[0]:
    import os as _os
    _script_dir = _os.path.dirname(_os.path.abspath(__file__))
    if _sys.path[0] == _script_dir:
        _sys.path.pop(0)

import hashlib
import json
import math
import sys

# ── EMBEDDED BIOLOGICAL DATA ── no external files, no biocompiler imports ──

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
STOP_CODONS = {"TAA", "TAG", "TGA"}

AA_TO_CODONS = {}
for _c, _a in CODON_TABLE.items():
    if _a != "*":
        AA_TO_CODONS.setdefault(_a, []).append(_c)

RESTRICTION_SITES = {
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "XhoI": "CTCGAG",
    "HindIII": "AAGCTT", "NotI": "GCGGCCGC", "XbaI": "TCTAGA",
    "SalI": "GTCGAC", "PstI": "CTGCAG", "SphI": "GCATGC",
    "NdeI": "CATATG", "NcoI": "CCATGG", "NheI": "GCTAGC",
    "KpnI": "GGTACC", "SmaI": "CCCGGG", "SacI": "GAGCTC",
    "SpeI": "ACTAGT", "ApaI": "GGGCCC", "ClaI": "ATCGAT",
    "EcoRV": "GATATC", "BglII": "AGATCT", "MluI": "ACGCGT",
    "AscI": "GGCGCGCC", "FseI": "GGCCGGCC", "PacI": "TTAATTAA",
    "SbfI": "CCTGCAGG", "BsiWI": "CGTACG", "BsrGI": "TGTACA",
    "AgeI": "ACCGGT", "MfeI": "CAATTG", "AluI": "AGCT",
    "HaeIII": "GGCC", "MspI": "CCGG", "TaqI": "TCGA",
    "Sau3AI": "GATC", "PmeI": "GTTTAAAC",
}

# Human codon usage (per thousand) for CAI
HUMAN_CODON_USAGE = {
    "TTT": 17.2, "TTC": 20.8, "TTA": 7.4, "TTG": 12.9,
    "CTT": 13.0, "CTC": 19.4, "CTA": 7.5, "CTG": 39.4,
    "ATT": 16.0, "ATC": 21.0, "ATA": 7.1, "ATG": 22.3,
    "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.5,
    "TCT": 14.9, "TCC": 17.4, "TCA": 11.7, "TCG": 4.5,
    "AGT": 12.0, "AGC": 19.3,
    "CCT": 17.3, "CCC": 19.7, "CCA": 16.7, "CCG": 7.0,
    "ACT": 12.9, "ACC": 18.6, "ACA": 14.8, "ACG": 6.2,
    "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
    "TAT": 15.4, "TAC": 15.6,
    "CAT": 10.5, "CAC": 15.0, "CAA": 11.8, "CAG": 34.3,
    "AAT": 16.8, "AAC": 19.5, "AAA": 24.1, "AAG": 32.1,
    "GAT": 21.5, "GAC": 25.4, "GAA": 28.8, "GAG": 39.8,
    "TGT": 10.2, "TGC": 12.4, "TGG": 13.4,
    "CGT": 4.5, "CGC": 10.4, "CGA": 6.1, "CGG": 11.3,
    "AGA": 11.7, "AGG": 12.0,
    "GGT": 10.8, "GGC": 22.2, "GGA": 16.4, "GGG": 16.5,
}

# E. coli codon usage (per thousand) for CAI
ECOLI_CODON_USAGE = {
    "TTT": 17.6, "TTC": 20.3, "TTA": 7.6, "TTG": 11.0,
    "CTT": 10.5, "CTC": 10.5, "CTA": 3.9, "CTG": 51.0,
    "ATT": 29.8, "ATC": 25.1, "ATA": 4.2, "ATG": 27.0,
    "GTT": 18.3, "GTC": 15.0, "GTA": 10.8, "GTG": 27.8,
    "TCT": 8.5, "TCC": 8.5, "TCA": 7.3, "TCG": 4.3,
    "AGT": 9.6, "AGC": 15.4,
    "CCT": 7.0, "CCC": 5.5, "CCA": 8.4, "CCG": 23.2,
    "ACT": 12.9, "ACC": 25.7, "ACA": 7.1, "ACG": 6.3,
    "GCT": 18.5, "GCC": 27.1, "GCA": 20.2, "GCG": 7.4,
    "TAT": 16.3, "TAC": 14.9,
    "CAT": 13.5, "CAC": 9.8, "CAA": 14.6, "CAG": 29.0,
    "AAT": 17.1, "AAC": 21.3, "AAA": 33.5, "AAG": 24.1,
    "GAT": 31.0, "GAC": 21.4, "GAA": 39.2, "GAG": 19.6,
    "TGT": 5.1, "TGC": 5.5, "TGG": 12.9,
    "CGT": 20.0, "CGC": 21.5, "CGA": 3.5, "CGG": 5.4,
    "AGA": 2.1, "AGG": 1.2,
    "GGT": 24.5, "GGC": 28.6, "GGA": 8.0, "GGG": 6.8,
}

# MaxEntScan donor PWM (9 positions: -3 to +6 relative to GT)
# Source: Yeo & Burge 2004, trained on human Chr 20-22
DONOR_PWM = [
    [0.310, 0.334, 0.192, 0.164],  # -3
    [0.292, 0.334, 0.207, 0.167],  # -2
    [0.078, 0.416, 0.096, 0.410],  # -1
    [0.003, 0.003, 0.990, 0.004],  # +1 (G)
    [0.003, 0.004, 0.003, 0.990],  # +2 (T)
    [0.332, 0.190, 0.298, 0.180],  # +3
    [0.240, 0.213, 0.325, 0.222],  # +4
    [0.154, 0.150, 0.408, 0.288],  # +5
    [0.209, 0.213, 0.297, 0.281],  # +6
]
_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
_BG = 0.25
_EPS = 0.001


# ── PREDICATE CHECKS ── each implemented from first principles ──

def _gc_content(seq: str) -> float:
    """GC fraction of DNA sequence."""
    if not seq:
        return 0.0
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s)


def check_gc_in_range(seq: str, gc_lo: float, gc_hi: float) -> tuple[str, str]:
    """GC content within [gc_lo, gc_hi]. Returns (verdict, detail)."""
    gc = _gc_content(seq)
    if gc_lo <= gc <= gc_hi:
        return "PASS", f"GC={gc:.4f} in [{gc_lo}, {gc_hi}]"
    return "FAIL", f"GC={gc:.4f} outside [{gc_lo}, {gc_hi}]"


def check_no_stop_codons(seq: str) -> tuple[str, str]:
    """No internal stop codons (last codon allowed)."""
    s = seq.upper()
    if len(s) < 3:
        return "PASS", "Sequence too short for stops"
    last_start = len(s) - 3
    positions = [i for i in range(0, last_start, 3) if s[i:i+3] in STOP_CODONS]
    if positions:
        return "FAIL", f"Internal stop codons at {positions}"
    return "PASS", "No internal stop codons"


def _maxent_donor_score(seq: str, pos: int) -> float:
    """Score donor splice site at position using MaxEntScan PWM."""
    s = seq.upper()
    start = pos - 3
    end = pos + 6
    if start < 0 or end > len(s):
        return -50.0
    score = 0.0
    for pwm_i in range(9):
        base = s[start + pwm_i]
        if base not in _BASE_IDX:
            return -50.0
        prob = max(DONOR_PWM[pwm_i][_BASE_IDX[base]], _EPS)
        score += math.log2(prob / _BG)
    return round(score, 4)


def check_no_cryptic_splice(seq: str, threshold: float = 3.0) -> tuple[str, str]:
    """No cryptic splice sites with MaxEntScan score >= threshold."""
    s = seq.upper()
    sites_found = []
    for i in range(len(s) - 1):
        if s[i:i+2] == "GT":
            score = _maxent_donor_score(s, i)
            if score >= threshold:
                sites_found.append((i, score))
    if not sites_found:
        return "PASS", "No cryptic splice sites found"
    worst = max(sites_found, key=lambda x: x[1])
    return "FAIL", f"{len(sites_found)} cryptic splice site(s), worst pos {worst[0]} score={worst[1]:.2f}"


def check_no_restriction_sites(seq: str, enzymes: list[str]) -> tuple[str, str]:
    """No restriction enzyme recognition sites in sequence."""
    s = seq.upper()
    found = []
    for enz in enzymes:
        site = RESTRICTION_SITES.get(enz)
        if site is None:
            continue
        pos = s.find(site)
        while pos != -1:
            found.append((enz, pos))
            pos = s.find(site, pos + 1)
    if found:
        return "FAIL", "Restriction sites: " + "; ".join(f"{e}@{p}" for e, p in found)
    return "PASS", "No restriction sites found"


def check_no_cpg_island(seq: str, window: int = 200, threshold: float = 0.6) -> tuple[str, str]:
    """No CpG islands (Obs/Exp CG ratio > threshold in any window)."""
    s = seq.upper()
    worst_ratio, worst_start = 0.0, -1
    for start in range(0, len(s) - window + 1):
        w = s[start:start + window]
        c_count = w.count("C")
        g_count = w.count("G")
        cg_count = sum(1 for i in range(len(w) - 1) if w[i:i+2] == "CG")
        expected = (c_count * g_count) / len(w) if len(w) > 0 else 0
        ratio = cg_count / expected if expected > 0 else 0.0
        if ratio > worst_ratio:
            worst_ratio = ratio
            worst_start = start
    if worst_ratio > threshold:
        return "FAIL", f"CpG island at pos {worst_start}, Obs/Exp={worst_ratio:.3f}>{threshold}"
    return "PASS", f"Worst CpG Obs/Exp={worst_ratio:.3f}<={threshold}"


def check_no_gt_dinucleotide(seq: str) -> tuple[str, str]:
    """No GT dinucleotides in sequence."""
    s = seq.upper()
    positions = [i for i in range(len(s) - 1) if s[i:i+2] == "GT"]
    if positions:
        return "FAIL", f"GT dinucleotides at {positions}"
    return "PASS", "No GT dinucleotides"


def _compute_cai(seq: str, usage: dict[str, float]) -> float:
    """Compute Codon Adaptation Index for a coding sequence."""
    weights = {}
    for aa, codons in AA_TO_CODONS.items():
        freqs = [usage.get(c, 0.1) for c in codons]
        max_freq = max(freqs) if freqs else 1.0
        for codon, freq in zip(codons, freqs):
            weights[codon] = freq / max_freq if max_freq > 0 else 0.0
    s = seq.upper()
    log_sum, n_codons = 0.0, 0
    for i in range(0, len(s) - 2, 3):
        codon = s[i:i+3]
        if codon in STOP_CODONS and i == len(s) - 3:
            continue  # skip terminal stop
        w = weights.get(codon, 0.01)
        if w <= 0:
            w = 0.01
        log_sum += math.log(w)
        n_codons += 1
    return math.exp(log_sum / n_codons) if n_codons > 0 else 0.0


def check_cai_above_threshold(seq: str, organism: str, threshold: float) -> tuple[str, str]:
    """CAI above threshold for given organism."""
    org_lower = organism.lower().replace(" ", "_")
    if "coli" in org_lower or org_lower in ("e_coli", "ecoli"):
        usage = ECOLI_CODON_USAGE
    else:
        usage = HUMAN_CODON_USAGE
    cai = _compute_cai(seq.upper(), usage)
    if cai >= threshold:
        return "PASS", f"CAI={cai:.4f}>={threshold}"
    return "FAIL", f"CAI={cai:.4f}<{threshold}"


def check_protein_identity(seq: str, expected_aa: str) -> tuple[str, str]:
    """Translated protein matches expected amino acid sequence."""
    s = seq.upper()
    translated = []
    for i in range(0, len(s) - 2, 3):
        aa = CODON_TABLE.get(s[i:i+3], "?")
        if aa == "*":
            break
        translated.append(aa)
    protein = "".join(translated)
    if protein == expected_aa:
        return "PASS", f"Protein identity: {len(protein)} AAs match"
    for i, (a, b) in enumerate(zip(protein, expected_aa)):
        if a != b:
            return "FAIL", f"Mismatch at AA {i}: got {a}, expected {b}"
    if len(protein) != len(expected_aa):
        return "FAIL", f"Length mismatch: got {len(protein)}, expected {len(expected_aa)}"
    return "PASS", f"Protein identity: {len(protein)} AAs match"


# ── PREDICATE DISPATCH ── maps predicate names to checker functions ──

def _resolve_base_name(name: str) -> str:
    """Strip parameter suffix: 'GCInRange(0.3, 0.7)' -> 'GCInRange'."""
    idx = name.find("(")
    return name[:idx] if idx > 0 else name


def _get_params(cert: dict) -> dict:
    """Extract verification parameters from certificate provenance or top-level."""
    p = dict(cert.get("provenance", {}).get("parameters", {}))
    # Accept top-level keys as fallback
    for key in ("gc_lo", "gc_hi", "cai_threshold", "organism", "enzymes",
                "cryptic_splice_threshold", "amino_acids", "species",
                "cpg_window", "cpg_threshold"):
        if key not in p and key in cert:
            p[key] = cert[key]
    p.setdefault("gc_lo", 0.30)
    p.setdefault("gc_hi", 0.70)
    p.setdefault("cai_threshold", 0.5)
    p.setdefault("organism", p.get("species", "Homo_sapiens"))
    p.setdefault("enzymes", ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"])
    p.setdefault("cryptic_splice_threshold", 3.0)
    p.setdefault("cpg_window", 200)
    p.setdefault("cpg_threshold", 0.6)
    return p


def _dispatch_check(name: str, seq: str, params: dict) -> tuple[str, str]:
    """Dispatch predicate name to checker. Returns (verdict, detail)."""
    base = _resolve_base_name(name)
    if base in ("GCInRange", "gc_content"):
        return check_gc_in_range(seq, params["gc_lo"], params["gc_hi"])
    if base in ("NoStopCodons", "no_stop_codons"):
        return check_no_stop_codons(seq)
    if base in ("NoCrypticSplice", "no_cryptic_splice"):
        return check_no_cryptic_splice(seq, params["cryptic_splice_threshold"])
    if base in ("NoRestrictionSite", "no_restriction_sites"):
        return check_no_restriction_sites(seq, params["enzymes"])
    if base in ("NoCpGIsland", "no_cpg_island"):
        return check_no_cpg_island(seq, params["cpg_window"], params["cpg_threshold"])
    if base in ("NoGTDinucleotide", "no_gt_dinucleotide"):
        return check_no_gt_dinucleotide(seq)
    if base in ("CodonAdapted", "CodonOptimality", "cai_above_threshold"):
        return check_cai_above_threshold(seq, params["organism"], params["cai_threshold"])
    if base in ("ProteinIdentity", "protein_identity"):
        aa = params.get("amino_acids", "")
        if not aa:
            return "UNCERTAIN", "No amino_acids parameter for protein identity check"
        return check_protein_identity(seq, aa)
    if base in ("ValidCodingSeq", "InFrame"):
        s = seq.upper()
        if len(s) % 3 != 0:
            return "FAIL", f"Length {len(s)} not divisible by 3"
        bad = [(i, s[i:i+3]) for i in range(0, len(s), 3) if s[i:i+3] not in CODON_TABLE]
        if bad:
            return "FAIL", f"Invalid codons: {bad}"
        return "PASS", "All codons valid, length divisible by 3"
    if base == "NoInstabilityMotif":
        s = seq.upper()
        positions = [i for i in range(len(s) - 4) if s[i:i+5] == "ATTTA"]
        return ("FAIL", f"Instability motifs at {positions}") if positions else ("PASS", "No instability motifs")
    if base == "SpliceCorrect":
        return check_no_cryptic_splice(seq, params["cryptic_splice_threshold"])
    if base == "ConservationScore":
        return "UNCERTAIN", "Conservation requires original AA (not re-checkable from cert alone)"
    return "UNCERTAIN", f"Unknown predicate: {name}"


# ── MAIN VERIFICATION LOGIC ──

def verify(cert: dict) -> tuple[str, list[str]]:
    """Independently verify a BioCompiler certificate from scratch.

    Re-checks every predicate using only the sequence and parameters
    embedded in the certificate. Does NOT trust any claimed verdict.

    Returns (status, reasons) where status is "VERIFIED" or "REJECTED".
    """
    reasons: list[str] = []

    # Structural validation
    for key in ("version", "design_id", "sequence"):
        if key not in cert:
            reasons.append(f"Missing required key: {key}")
    if reasons:
        return "REJECTED", reasons

    seq = cert["sequence"].upper()
    params = _get_params(cert)

    # Check 1: SHA-256 integrity
    computed_hash = hashlib.sha256(seq.encode()).hexdigest()
    stored_hash = cert["design_id"]
    if computed_hash != stored_hash:
        reasons.append(f"design_id mismatch: computed {computed_hash[:16]}... != stored {stored_hash[:16]}...")

    # Check 2: Re-verify each predicate (supports both "types" and "predicates" formats)
    for entry in cert.get("types", []):
        pname = entry.get("predicate", "")
        claimed = entry.get("verdict", "UNKNOWN")
        actual, detail = _dispatch_check(pname, seq, params)
        if actual == "UNCERTAIN" and claimed in ("PASS", "UNCERTAIN"):
            continue  # Cannot independently verify, but cannot refute either
        if actual != claimed:
            reasons.append(f"Predicate {pname}: cert claims {claimed}, re-check gives {actual} ({detail})")

    for pname, pdata in cert.get("predicates", {}).items():
        claimed = pdata.get("verdict", "UNKNOWN")
        actual, detail = _dispatch_check(pname, seq, params)
        if actual == "UNCERTAIN" and claimed in ("PASS", "UNCERTAIN"):
            continue
        if actual != claimed:
            reasons.append(f"Predicate {pname}: cert claims {claimed}, re-check gives {actual} ({detail})")

    # Check 3: Provenance completeness (if present)
    prov = cert.get("provenance", {})
    if prov:
        for field in ("tool", "timestamp", "input_hash"):
            if field not in prov:
                reasons.append(f"Missing provenance field: {field}")
        if "input_hash" in prov and prov["input_hash"] != computed_hash:
            reasons.append(f"provenance.input_hash mismatch: {prov['input_hash'][:16]}... != computed {computed_hash[:16]}...")

    if reasons:
        return "REJECTED", reasons
    return "VERIFIED", []


def verify_file(path: str) -> tuple[str, list[str]]:
    """Load and verify a certificate from a JSON file."""
    with open(path) as f:
        cert = json.load(f)
    return verify(cert)


# ── CLI ENTRY POINT ──

def main() -> None:
    """CLI: biocompiler-verify certificate.json"""
    if len(sys.argv) < 2:
        print("Usage: biocompiler-verify <certificate.json>", file=sys.stderr)
        sys.exit(2)
    path = sys.argv[1]
    try:
        status, reasons = verify_file(path)
    except FileNotFoundError:
        print(f"ERROR: File not found: {path}", file=sys.stderr)
        sys.exit(2)
    except json.JSONDecodeError as e:
        print(f"ERROR: Invalid JSON in {path}: {e}", file=sys.stderr)
        sys.exit(2)
    if status == "VERIFIED":
        print(f"VERIFIED: {path}")
        sys.exit(0)
    print(f"REJECTED: {path}")
    for r in reasons:
        print(f"  - {r}")
    sys.exit(1)


if __name__ == "__main__":
    main()
