#!/usr/bin/env python3
"""Validate organism configuration data in the biocompiler repository.

Checks:
  a. For each organism, verifies:
     - All 20 amino acids have codon entries
     - Codon frequencies sum to reasonable values per amino acid
     - Relative adaptiveness values are in [0, 1]
     - The most frequent codon per AA has w=1.0
  b. Checks codon pair bias data against published Buchan et al. (2006) values
  c. Verifies UTR model data has reasonable stability scores
  d. Checks that the Kazusa-derived frequencies match known E. coli codon usage:
     - Most common Leu codon should be CTG (~51 per thousand)
     - Most common Gly codon should be GGC (~28 per thousand)
     - Most common Ala codon should be GCC (~27 per thousand)
"""

from __future__ import annotations

import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from biocompiler.organisms import (
    CODON_USAGE_TABLES,
    CODON_ADAPTIVENESS_TABLES,
    PREFERRED_CODON_TABLES,
    SUPPORTED_ORGANISMS,
    E_COLI_CODON_USAGE,
    E_COLI_CODON_ADAPTIVENESS,
    E_COLI_PREFERRED_CODONS,
    ECOLI_CODON_USAGE,
    HUMAN_CODON_USAGE,
    HUMAN_CODON_ADAPTIVENESS,
    HUMAN_PREFERRED_CODONS,
    MOUSE_CODON_USAGE,
    MOUSE_CODON_ADAPTIVENESS,
    MOUSE_PREFERRED_CODONS,
    CHO_CODON_USAGE,
    CHO_CODON_ADAPTIVENESS,
    CHO_PREFERRED_CODONS,
    YEAST_CODON_USAGE,
    YEAST_CODON_ADAPTIVENESS,
    YEAST_PREFERRED_CODONS,
)
from biocompiler.organisms.e_coli import (
    E_COLI_CODON_PAIR_BIAS,
    E_COLI_EXPRESSION_OPTIMIZATION_PARAMS,
    compute_codon_pair_bias,
)
from biocompiler.organisms.human import (
    HUMAN_CODON_PAIR_BIAS,
    HUMAN_EXPRESSION_OPTIMIZATION_PARAMS,
    HUMAN_UTR_MODELS,
)
from biocompiler.organisms.mouse import (
    MOUSE_CODON_PAIR_BIAS,
    MOUSE_EXPRESSION_OPTIMIZATION_PARAMS,
    MOUSE_UTR_MODELS,
)
from biocompiler.codon_pair_scoring import (
    get_codon_pair_data,
    compute_cpb,
    score_codon_pair,
)
from biocompiler.utr_models import (
    ORGANISM_UTR_CONFIGS,
    AVAILABLE_ORGANISMS as UTR_ORGANISMS,
    score_5utr,
    score_3utr,
    suggest_5utr,
    suggest_3utr,
    UTRConfig,
)
from biocompiler.type_system import AA_TO_CODONS

# Standard 20 amino acids (excluding stop codon)
STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")


def validate_organism_configs() -> dict:
    """Run all organism config validation checks."""
    report: dict = {}
    issues: list[str] = []
    warnings: list[str] = []
    organism_reports: dict = {}

    print("=" * 70)
    print("PART B: Organism Config Validation")
    print("=" * 70)
    print(f"\nSupported organisms: {SUPPORTED_ORGANISMS}")
    print(f"UTR organisms: {UTR_ORGANISMS}")

    # (a) For each organism, verify codon usage tables
    print("\n--- Codon Usage Table Validation ---\n")

    for org_name in SUPPORTED_ORGANISMS:
        usage = CODON_USAGE_TABLES.get(org_name)
        adaptiveness = CODON_ADAPTIVENESS_TABLES.get(org_name)
        preferred = PREFERRED_CODON_TABLES.get(org_name)

        if usage is None:
            issues.append(f"{org_name}: No codon usage table found")
            continue

        org_info: dict = {"organism": org_name}

        # Check all 20 amino acids have codon entries
        aas_with_codons = set()
        for codon, (aa, frac, freq, count) in usage.items():
            if aa != "*":
                aas_with_codons.add(aa)

        missing_aas = STANDARD_AAS - aas_with_codons
        if missing_aas:
            issues.append(f"{org_name}: Missing codon entries for amino acids: {sorted(missing_aas)}")
        org_info["missing_amino_acids"] = sorted(missing_aas) if missing_aas else []
        org_info["amino_acids_with_codons"] = len(aas_with_codons)

        # Check codon frequencies sum to reasonable values per amino acid
        aa_freqs: dict[str, list[tuple[str, float]]] = {}
        for codon, (aa, frac, freq, count) in usage.items():
            if aa == "*":
                continue
            aa_freqs.setdefault(aa, []).append((codon, frac))

        freq_sum_issues = []
        for aa, codon_fracs in sorted(aa_freqs.items()):
            total_frac = sum(f for _, f in codon_fracs)
            if not (0.95 <= total_frac <= 1.05):
                freq_sum_issues.append((aa, total_frac, len(codon_fracs)))

        if freq_sum_issues:
            for aa, total, n in freq_sum_issues:
                warnings.append(
                    f"{org_name}: AA {aa} frequency sum = {total:.4f} (expected ~1.0, {n} codons)"
                )
        org_info["freq_sum_issues"] = len(freq_sum_issues)

        # Check relative adaptiveness values are in [0, 1]
        if adaptiveness:
            invalid_w = [(codon, w) for codon, w in adaptiveness.items() if not (0.0 <= w <= 1.0)]
            if invalid_w:
                issues.append(
                    f"{org_name}: {len(invalid_w)} adaptiveness values outside [0,1]: "
                    f"{invalid_w[:5]}"
                )
        org_info["adaptiveness_range_ok"] = not invalid_w if adaptiveness else "N/A"

        # Check that the most frequent codon per AA has w=1.0
        if adaptiveness and preferred:
            max_w_per_aa = {}
            for aa, codon_list in aa_freqs.items():
                # Find the codon with highest frequency for this AA
                best_codon = max(codon_list, key=lambda x: x[1])[0]
                w_val = adaptiveness.get(best_codon, -1)
                max_w_per_aa[aa] = (best_codon, w_val)

            not_one = [(aa, codon, w) for aa, (codon, w) in max_w_per_aa.items() if abs(w - 1.0) > 0.01]
            if not_one:
                issues.append(
                    f"{org_name}: {len(not_one)} AAs where most frequent codon doesn't have w=1.0: "
                    f"{not_one[:5]}"
                )
            org_info["preferred_codon_w1_ok"] = len(not_one) == 0
        else:
            org_info["preferred_codon_w1_ok"] = "N/A"

        # Print summary
        print(f"  {org_name}:")
        print(f"    Amino acids with codons: {len(aas_with_codons)}/20")
        print(f"    Missing AAs: {sorted(missing_aas) if missing_aas else 'None'}")
        print(f"    Frequency sum issues: {len(freq_sum_issues)}")
        print(f"    Adaptiveness range OK: {not invalid_w if adaptiveness else 'N/A'}")
        if adaptiveness:
            print(f"    Adaptiveness entries: {len(adaptiveness)}")
        print()

        organism_reports[org_name] = org_info

    # (d) Kazusa-derived E. coli frequency checks
    print("--- E. coli Kazusa Frequency Verification ---\n")

    # Check against legacy ECOLI_CODON_USAGE (per-thousand)
    ecoli_legacy = ECOLI_CODON_USAGE
    leu_codons_legacy = {c: ecoli_legacy[c] for c in AA_TO_CODONS.get("L", [])}
    gly_codons_legacy = {c: ecoli_legacy[c] for c in AA_TO_CODONS.get("G", [])}
    ala_codons_legacy = {c: ecoli_legacy[c] for c in AA_TO_CODONS.get("A", [])}

    most_common_leu = max(leu_codons_legacy, key=leu_codons_legacy.get)
    most_common_gly = max(gly_codons_legacy, key=gly_codons_legacy.get)
    most_common_ala = max(ala_codons_legacy, key=ala_codons_legacy.get)

    print(f"  Legacy ECOLI_CODON_USAGE (per-thousand):")
    print(f"    Leu codons: {dict(sorted(leu_codons_legacy.items(), key=lambda x: -x[1]))}")
    print(f"    Most common Leu: {most_common_leu} = {leu_codons_legacy[most_common_leu]}/k (expected CTG ~51/k)")

    leu_ok = most_common_leu == "CTG" and abs(leu_codons_legacy["CTG"] - 51.0) < 5.0
    if not leu_ok:
        issues.append(f"E. coli Leu: Most common codon is {most_common_leu}={leu_codons_legacy[most_common_leu]}/k, expected CTG ~51/k")
    else:
        print(f"    ✓ Leu CTG frequency check passed ({leu_codons_legacy['CTG']}/k ≈ 51/k)")

    print(f"    Gly codons: {dict(sorted(gly_codons_legacy.items(), key=lambda x: -x[1]))}")
    print(f"    Most common Gly: {most_common_gly} = {gly_codons_legacy[most_common_gly]}/k (expected GGC ~28/k)")

    gly_ok = most_common_gly == "GGC" and abs(gly_codons_legacy["GGC"] - 28.0) < 5.0
    if not gly_ok:
        issues.append(f"E. coli Gly: Most common codon is {most_common_gly}={gly_codons_legacy[most_common_gly]}/k, expected GGC ~28/k")
    else:
        print(f"    ✓ Gly GGC frequency check passed ({gly_codons_legacy['GGC']}/k ≈ 28/k)")

    print(f"    Ala codons: {dict(sorted(ala_codons_legacy.items(), key=lambda x: -x[1]))}")
    print(f"    Most common Ala: {most_common_ala} = {ala_codons_legacy[most_common_ala]}/k (expected GCC ~27/k)")

    ala_ok = most_common_ala == "GCC" and abs(ala_codons_legacy["GCC"] - 27.0) < 5.0
    if not ala_ok:
        issues.append(f"E. coli Ala: Most common codon is {most_common_ala}={ala_codons_legacy[most_common_ala]}/k, expected GCC ~27/k")
    else:
        print(f"    ✓ Ala GCC frequency check passed ({ala_codons_legacy['GCC']}/k ≈ 27/k)")

    # Also check the richer tuple format
    ecoli_rich = E_COLI_CODON_USAGE
    leu_codons_rich = {c: ecoli_rich[c] for c in AA_TO_CODONS.get("L", []) if c in ecoli_rich}
    gly_codons_rich = {c: ecoli_rich[c] for c in AA_TO_CODONS.get("G", []) if c in ecoli_rich}
    ala_codons_rich = {c: ecoli_rich[c] for c in AA_TO_CODONS.get("A", []) if c in ecoli_rich}

    most_common_leu_rich = max(leu_codons_rich, key=lambda c: leu_codons_rich[c][2])
    most_common_gly_rich = max(gly_codons_rich, key=lambda c: gly_codons_rich[c][2])
    most_common_ala_rich = max(ala_codons_rich, key=lambda c: ala_codons_rich[c][2])

    print(f"\n  E_COLI_CODON_USAGE (tuple format, per-thousand from Kazusa):")
    print(f"    Most common Leu: {most_common_leu_rich} = {leu_codons_rich[most_common_leu_rich][2]}/k")
    print(f"    Most common Gly: {most_common_gly_rich} = {gly_codons_rich[most_common_gly_rich][2]}/k")
    print(f"    Most common Ala: {most_common_ala_rich} = {ala_codons_rich[most_common_ala_rich][2]}/k")

    # (b) Codon pair bias validation
    print("\n--- Codon Pair Bias Validation ---\n")

    # Check E. coli CPB data
    ecoli_cpb = E_COLI_CODON_PAIR_BIAS
    print(f"  E. coli codon pair bias entries: {len(ecoli_cpb)}")
    positive_cpbs = [(k, v) for k, v in ecoli_cpb.items() if v > 0]
    negative_cpbs = [(k, v) for k, v in ecoli_cpb.items() if v < 0]
    print(f"    Over-represented pairs (positive): {len(positive_cpbs)}")
    print(f"    Under-represented pairs (negative): {len(negative_cpbs)}")

    # Validate against known Buchan et al. patterns
    # CTG-CTG should be the most over-represented pair in E. coli
    if "CTG-CTG" in ecoli_cpb:
        ctg_ctg = ecoli_cpb["CTG-CTG"]
        print(f"    CTG-CTG CPB = {ctg_ctg} (expected positive, most over-represented)")
        if ctg_ctg <= 0:
            issues.append(f"E. coli CTG-CTG CPB = {ctg_ctg}, expected positive (over-represented)")
    else:
        issues.append("E. coli: CTG-CTG pair not found in CPB data")

    # AGG-AGA and AGA-AGG should be under-represented (rare Arg pairs)
    for pair in ["AGG-AGA", "AGA-AGG", "CUA-ATA"]:
        if pair in ecoli_cpb:
            val = ecoli_cpb[pair]
            print(f"    {pair} CPB = {val} (expected negative)")
            if val >= 0:
                warnings.append(f"E. coli {pair} CPB = {val}, expected negative (under-represented)")

    # Test compute_codon_pair_bias function
    test_dna = "ATGCTGGAAGAACTGATG"  # ATG-CTG-GAA-GAA-CTG-ATG
    cpb_score = compute_codon_pair_bias(test_dna, "e_coli")
    print(f"    Test CPB computation for '{test_dna}': {cpb_score:.4f}")

    # Check Human CPB data
    human_cpb = HUMAN_CODON_PAIR_BIAS
    print(f"\n  Human codon pair bias entries: {len(human_cpb)}")
    positive_human = sum(1 for v in human_cpb.values() if v > 0)
    negative_human = sum(1 for v in human_cpb.values() if v < 0)
    print(f"    Over-represented pairs: {positive_human}")
    print(f"    Under-represented pairs: {negative_human}")

    # Human: CTG_CTC should be most over-represented (from the data)
    if "CTG_CTC" in human_cpb:
        print(f"    CTG_CTC CPB = {human_cpb['CTG_CTC']} (expected most over-represented)")

    # Test codon_pair_scoring module
    print("\n  Testing codon_pair_scoring module:")
    for org in ["Escherichia_coli", "Homo_sapiens", "Mus_musculus", "CHO_K1", "Saccharomyces_cerevisiae"]:
        data = get_codon_pair_data(org)
        print(f"    {org}: {len(data)} CPB entries")

    score = score_codon_pair("CTG", "CTG", "Escherichia_coli")
    print(f"    score_codon_pair('CTG','CTG','E_coli'): {score}")

    cpb_e = compute_cpb("ATGCTGCTGGAAGAA", "Escherichia_coli")
    print(f"    compute_cpb test (E. coli): {cpb_e:.4f}")

    # (c) UTR model validation
    print("\n--- UTR Model Validation ---\n")

    for org_name, config in ORGANISM_UTR_CONFIGS.items():
        print(f"  {org_name}:")
        print(f"    5' UTR consensus: {config.utr5_consensus}")
        print(f"    3' UTR consensus: {config.utr3_consensus}")
        print(f"    Kozak: {config.kozak_sequence}")
        print(f"    Shine-Dalgarno: {config.shine_dalgarno}")
        print(f"    PolyA signal: {config.polya_signal}")
        print(f"    Stability motifs: {config.stability_motifs}")
        print(f"    Instability motifs: {config.instability_motifs}")
        print(f"    Splicing signals: {config.splicing_signals}")

        # Verify organism type constraints
        if config.shine_dalgarno is not None:
            # Prokaryotic — should NOT have Kozak or PolyA
            if config.kozak_sequence is not None:
                issues.append(f"{org_name}: Has both Shine-Dalgarno and Kozak — should be prokaryotic only")
            if config.polya_signal is not None:
                issues.append(f"{org_name}: Has both Shine-Dalgarno and PolyA — unexpected for prokaryote")
        else:
            # Eukaryotic — should have Kozak or PolyA
            if config.kozak_sequence is None and config.polya_signal is None:
                if org_name != "Saccharomyces_cerevisiae":
                    issues.append(f"{org_name}: Eukaryote with no Kozak and no PolyA signal")

        # Test scoring functions
        s5 = score_5utr(config.utr5_consensus + "ATG", org_name)
        s3 = score_3utr(config.utr3_consensus, org_name)
        print(f"    5' UTR score (consensus): {s5:.3f}")
        print(f"    3' UTR score (consensus): {s3:.3f}")

        # Verify consensus UTRs score reasonably well
        if s5 < 0.3:
            warnings.append(f"{org_name}: 5' UTR consensus scores low ({s5:.3f})")
        if s3 < 0.3:
            warnings.append(f"{org_name}: 3' UTR consensus scores low ({s3:.3f})")

        # Test suggested UTRs
        suggested_5 = suggest_5utr(org_name)
        suggested_3 = suggest_3utr(org_name)
        s5_suggested = score_5utr(suggested_5, org_name)
        s3_suggested = score_3utr(suggested_3, org_name)
        print(f"    Suggested 5' UTR: {suggested_5} (score: {s5_suggested:.3f})")
        print(f"    Suggested 3' UTR: {suggested_3} (score: {s3_suggested:.3f})")

        if s5_suggested < 0.5:
            warnings.append(f"{org_name}: Suggested 5' UTR scores low ({s5_suggested:.3f})")
        if s3_suggested < 0.5:
            warnings.append(f"{org_name}: Suggested 3' UTR scores low ({s3_suggested:.3f})")

        print()

    # Check Human UTR models detailed structure
    print("--- Human UTR Models Detail ---\n")
    human_utr = HUMAN_UTR_MODELS
    if "5utr" in human_utr:
        h5 = human_utr["5utr"]
        print(f"  Human 5' UTR model:")
        if "kozak_preferences" in h5:
            for gene_type, prefs in h5["kozak_preferences"].items():
                tei = prefs.get("typical_tei_range", "N/A")
                print(f"    {gene_type}: strength={prefs.get('strength')}, TEI={tei}")
        print(f"    Max 5' UTR length: {h5.get('max_5utr_length')}")
        print(f"    Min 5' UTR length: {h5.get('min_5utr_length')}")
        print(f"    Avoid upstream AUG: {h5.get('avoid_upstream_aug')}")
    if "3utr" in human_utr:
        h3 = human_utr["3utr"]
        print(f"  Human 3' UTR model:")
        if "polya_signal" in h3:
            print(f"    PolyA signal: {h3['polya_signal'].get('consensus')}")
            print(f"    Variant signals: {h3['polya_signal'].get('variant_signals')}")
        if "au_rich_elements" in h3:
            are = h3["au_rich_elements"]
            print(f"    ARE consensus: {are.get('consensus')}")
            if "classes" in are:
                for cls_name, cls_info in are["classes"].items():
                    print(f"    ARE {cls_name}: effect={cls_info.get('effect')}, half_life={cls_info.get('half_life_range_min')}")
        print(f"    Max 3' UTR length: {h3.get('max_3utr_length')}")
        print(f"    Min 3' UTR length: {h3.get('min_3utr_length')}")

    # Check expression optimization params
    print("\n--- Expression Optimization Parameters ---\n")
    ecoli_params = E_COLI_EXPRESSION_OPTIMIZATION_PARAMS
    print(f"  E. coli params: {ecoli_params}")
    human_params = HUMAN_EXPRESSION_OPTIMIZATION_PARAMS
    print(f"  Human params: {human_params}")

    # Verify GC content targets are reasonable
    for name, params in [("E. coli", ecoli_params), ("Human", human_params)]:
        gc_target = params.get("gc_content_target")
        gc_min = params.get("gc_content_min")
        gc_max = params.get("gc_content_max")
        if gc_target is not None:
            if not (0.2 <= gc_target <= 0.8):
                issues.append(f"{name}: GC target {gc_target} outside reasonable range [0.2, 0.8]")
            if gc_min is not None and gc_max is not None:
                if not (gc_min < gc_target < gc_max):
                    issues.append(f"{name}: GC range [{gc_min}, {gc_max}] doesn't bracket target {gc_target}")

    # Final summary
    print("\n" + "=" * 70)
    print("ORGANISM CONFIG VALIDATION SUMMARY")
    print("=" * 70)
    print(f"Total organisms checked: {len(SUPPORTED_ORGANISMS)}")
    print(f"Total issues found: {len(issues)}")
    print(f"Total warnings: {len(warnings)}")

    print(f"\nE. coli Kazusa verification:")
    print(f"  Leu: CTG {'PASS' if leu_ok else 'FAIL'} ({leu_codons_legacy.get('CTG', 'N/A')}/k)")
    print(f"  Gly: GGC {'PASS' if gly_ok else 'FAIL'} ({gly_codons_legacy.get('GGC', 'N/A')}/k)")
    print(f"  Ala: GCC {'PASS' if ala_ok else 'FAIL'} ({ala_codons_legacy.get('GCC', 'N/A')}/k)")

    if issues:
        print("\nISSUES:")
        for i, issue in enumerate(issues, 1):
            print(f"  {i}. {issue}")

    if warnings:
        print("\nWARNINGS:")
        for i, w in enumerate(warnings, 1):
            print(f"  {i}. {w}")

    report["organism_reports"] = organism_reports
    report["ecoli_kazusa"] = {"leu_ok": leu_ok, "gly_ok": gly_ok, "ala_ok": ala_ok}
    report["issues"] = issues
    report["warnings"] = warnings
    report["passed"] = len(issues) == 0

    return report


if __name__ == "__main__":
    result = validate_organism_configs()
    sys.exit(0 if result["passed"] else 1)
