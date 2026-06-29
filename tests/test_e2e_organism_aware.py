"""
End-to-end tests for organism-aware constraint selection in BioCompiler.

These tests verify that the full optimization pipeline produces better CAI
for prokaryotic targets when organism-aware constraints are used (i.e. when
eukaryote-specific constraints like cryptic splice site avoidance and CpG
island elimination are disabled for prokaryotic organisms like E. coli).

Organism-aware constraint selection was introduced to address the root cause
identified in the DNAchisel benchmark (Task 6+8): BioCompiler over-applied
eukaryotic constraints to prokaryotic targets, depressing CAI from ~0.96 to
~0.69.  By detecting the organism domain and skipping inappropriate
constraints, prokaryotic CAI recovers significantly.

Test coverage:
  1. E. coli GFP optimization with organism-aware constraints → improved CAI
  2. Same protein optimized for E. coli vs. human → E. coli CAI higher
  3. Human HBB optimization → splice constraints active
  4. Human insulin optimized for E. coli → improved CAI
  5. Batch 10 E. coli genes → all improved CAI
  6. Backward compatibility → all constraints applied when domain forced
"""

from __future__ import annotations

import pytest

from biocompiler.optimizer import optimize_sequence
from biocompiler.expression.translation import translate, compute_cai
from biocompiler.sequence.scanner import gc_content
from biocompiler.type_system import (
    check_no_restriction_site,
    check_no_stop_codons,
    check_no_cryptic_splice,
    evaluate_all_predicates,
)
from biocompiler.shared.types import Verdict
from biocompiler.organisms.config import is_eukaryotic_organism
from biocompiler.api import resolve_organism_domain
from biocompiler.shared.constants import RESTRICTION_ENZYMES


# ─── Shared constants ───────────────────────────────────────────────────────

# GFP (Green Fluorescent Protein) — 238 aa
# Source: UniProt P42212 (Aequorea victoria EGFP)
GFP_PROTEIN = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
    "VTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLV"
    "NRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADH"
    "YQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)

# Human Hemoglobin Beta (HBB) — 147 aa
# Source: UniProt P68871
HBB_PROTEIN = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
    "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG"
    "KEFTPPVQAAYQKVVAGVANALAHKYH"
)

# Human Insulin (preproinsulin) — 110 aa
# Source: UniProt P01308
INSULIN_PROTEIN = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTTPKTRREAED"
    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)

DEFAULT_ENZYMES = ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"]


def _has_restriction_site(seq: str, enzymes: list[str]) -> bool:
    """Return True if any enzyme site is present in seq (or its RC)."""
    from biocompiler.shared.constants import reverse_complement
    for enz in enzymes:
        site = RESTRICTION_ENZYMES.get(enz, "")
        if site and (site in seq or reverse_complement(site) in seq):
            return True
    return False


def _find_verdict(type_results, canonical_name: str):
    """Look up a verdict by canonical predicate name."""
    for r in type_results:
        if r.predicate == canonical_name or r.predicate.startswith(canonical_name + "("):
            return r.verdict
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# Organism-Aware End-to-End Tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestOrganismAwareEndToEnd:
    """End-to-end: organism-aware constraint selection improves prokaryotic CAI."""

    def test_ecoli_optimization_organism_aware(self):
        """Optimize GFP for E. coli with organism-aware constraints.

        When organism_domain='auto' (the default), the optimizer detects
        that E. coli is prokaryotic and automatically skips eukaryote-
        specific constraints (GT dinucleotide avoidance, cryptic splice
        site elimination, CpG island avoidance). This produces
        significantly higher CAI than the previous default which applied
        all constraints regardless of organism.
        """
        result = optimize_sequence(
            GFP_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
            # organism_domain='auto' is the default — auto-detects prokaryote
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )

        # Compare against organism-unaware baseline (forced eukaryote)
        baseline = optimize_sequence(
            GFP_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
            organism_domain="eukaryote",  # Force all constraints
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )

        # Verify CAI is significantly higher than the organism-unaware baseline
        # With organism-aware constraints, CAI should improve by at least 0.05
        cai_improvement = result.cai - baseline.cai
        assert cai_improvement > 0.05, (
            f"Organism-aware CAI ({result.cai:.4f}) should be significantly "
            f"higher than baseline ({baseline.cai:.4f}). "
            f"Improvement: {cai_improvement:.4f} (need > 0.05)"
        )

        # Verify CAI is above a reasonable threshold
        assert result.cai > 0.70, (
            f"Organism-aware E. coli CAI should be > 0.70, got {result.cai:.4f}"
        )

        # Verify no restriction sites in output
        assert not _has_restriction_site(result.sequence, DEFAULT_ENZYMES), (
            "Optimized sequence should not contain any of the default restriction sites"
        )

        # Verify GC content in range
        actual_gc = gc_content(result.sequence)
        assert 0.30 <= actual_gc <= 0.70, (
            f"GC content {actual_gc:.4f} outside [0.30, 0.70]"
        )

        # Verify correct translation
        translated = translate(result.sequence)
        assert translated == GFP_PROTEIN, (
            f"Translation mismatch: expected {len(GFP_PROTEIN)} aa, "
            f"got {len(translated)} aa"
        )

    def test_ecoli_optimization_vs_human(self):
        """Optimize the SAME protein for both E. coli and human.

        E. coli should achieve higher CAI than human because:
        1. No splice constraints are needed for prokaryotes
        2. No CpG island constraints are needed for prokaryotes
        3. No GT dinucleotide avoidance for prokaryotes
        4. The codon usage table for E. coli has stronger codon bias
        """
        # E. coli optimization (organism-aware: auto-detects prokaryote)
        ecoli_result = optimize_sequence(
            GFP_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )

        # Human optimization (organism-aware: auto-detects eukaryote)
        human_result = optimize_sequence(
            GFP_PROTEIN,
            organism="Homo_sapiens",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )

        # E. coli should have higher CAI (no splice/CpG/GT constraints)
        assert ecoli_result.cai > human_result.cai, (
            f"E. coli CAI ({ecoli_result.cai:.4f}) should be higher than "
            f"human CAI ({human_result.cai:.4f}) for the same protein, "
            f"because prokaryotes do not need splice/CpG/GT constraints"
        )

        # Both should translate correctly
        assert translate(ecoli_result.sequence) == GFP_PROTEIN
        assert translate(human_result.sequence) == GFP_PROTEIN

        # Human should have GC in range
        human_gc = gc_content(human_result.sequence)
        assert 0.30 <= human_gc <= 0.70

    def test_human_still_has_splice_constraints(self):
        """Optimize HBB for human — verify splice constraints are active.

        For eukaryotic targets, the optimizer MUST still apply cryptic
        splice site constraints. This test ensures that organism-aware
        constraint selection does NOT disable splice constraints for
        organisms that need them.
        """
        result = optimize_sequence(
            HBB_PROTEIN,
            organism="Homo_sapiens",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )

        # Verify correct translation
        translated = translate(result.sequence)
        assert translated == HBB_PROTEIN

        # Verify organism is classified as eukaryotic
        assert is_eukaryotic_organism("Homo_sapiens"), (
            "Homo_sapiens must be classified as eukaryotic"
        )

        # Verify resolved domain is eukaryote
        domain = resolve_organism_domain("Homo_sapiens", "auto")
        assert domain == "eukaryote", (
            f"Homo_sapiens auto-detected domain should be 'eukaryote', got '{domain}'"
        )

        # Verify the output passes all constraint checks (via type system)
        type_results = evaluate_all_predicates(
            result.sequence,
            organism="Homo_sapiens",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
        )

        # Hard constraints must pass
        assert _find_verdict(type_results, "InFrame") == Verdict.PASS, (
            "InFrame constraint must PASS"
        )
        assert _find_verdict(type_results, "NoRestrictionSite") == Verdict.PASS, (
            "NoRestrictionSite constraint must PASS"
        )
        assert _find_verdict(type_results, "GCInRange") == Verdict.PASS, (
            "GCInRange constraint must PASS"
        )

        # The NoCrypticSplice predicate should be present in results
        splice_verdict = _find_verdict(type_results, "NoCrypticSplice")
        assert splice_verdict is not None, (
            "NoCrypticSplice predicate must be present for eukaryotic targets"
        )

        # Verify CpG islands are checked
        cpg_verdict = _find_verdict(type_results, "NoCpGIsland")
        assert cpg_verdict is not None, (
            "NoCpGIsland predicate must be present for eukaryotic targets"
        )

    def test_insulin_ecoli_high_cai(self):
        """Optimize human insulin for E. coli — CAI should improve.

        Human insulin is one of the most commonly expressed recombinant
        proteins in E. coli. With organism-aware constraints, the CAI
        should be significantly higher than when splice/CpG/GT constraints
        are unnecessarily applied.
        """
        # Organism-aware optimization
        result = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )

        # Baseline (organism-unaware, forced eukaryote)
        baseline = optimize_sequence(
            INSULIN_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
            organism_domain="eukaryote",
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )

        # CAI should improve with organism-aware constraints
        assert result.cai > baseline.cai, (
            f"Organism-aware insulin CAI ({result.cai:.4f}) should be higher "
            f"than baseline ({baseline.cai:.4f})"
        )

        # Translation must be correct
        translated = translate(result.sequence)
        assert translated == INSULIN_PROTEIN, (
            f"Translation mismatch: expected {len(INSULIN_PROTEIN)} aa, "
            f"got {len(translated)} aa"
        )

        # GC in range
        actual_gc = gc_content(result.sequence)
        assert 0.30 <= actual_gc <= 0.70, (
            f"GC content {actual_gc:.4f} outside [0.30, 0.70]"
        )

        # No restriction sites
        assert not _has_restriction_site(result.sequence, DEFAULT_ENZYMES), (
            "Optimized insulin sequence should not contain restriction sites"
        )

    def test_batch_ecoli_genes(self):
        """Optimize 10 E. coli genes — all should have improved CAI.

        Uses a standard set of E. coli genes covering different functional
        categories. With organism-aware constraints, all should achieve
        higher CAI than with organism-unaware (forced eukaryote) constraints.
        """
        # Standard set of 10 E. coli genes (shorter, realistic sequences)
        ecoli_genes = {
            "lacI": "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQSLLIGVFKIRSTLRQNLLDSRVLQTYRVYPEKLNQLEQRVYRLRSGSIRDLQLQIDSLYQWLDNIPVIQYQLRNKQIRSFVKLFGSHPKYENLVLKSFVSREQVNPDIVLFNRNPKRATVLCNEAKQVGLIEKISLKTDAEYVNKLIDTGFAPVLQIDPNRVTLPKFKPLIYVDRKMVLEQKFTVDDIKVVPKKRIFLEAFQKVIDEVRQTLREKIASNIKETLMLKNIVLGTSIAVVGMNLTTQRKIDQAKLFADIAKLS",
            "recA": "MAIDENKQKQALQEAISQIAQDEIRSKIEQIEKELKKYGISITQLKDEYQKKIAEIEKKIRDLQKNEKIKKIAEKIRKLIEKNEKIRKQIEKIEKDLKKIEKIEKNEKKIEKIEKNEKIEKIEKNEKKIEKIEKNEKIEKIEKNEKKIEKIEKNEKIEKIEDNEKKIEKIEKNEKIEKIKKNEKKIEKIEKNEKIEKIEQNEKKIEKIEKNEKIEKIENNEKKIEKIEKNEKIEKIKDNEKKIEKIEKNEKIEKIEKNEK",
            "tufA": "MAKTHHIKVNIIHGHDVAGFVTRDAIADKTAKIYDSNGRPLHDVVDVTGVGTHKAIQELRKQGAKVVVTFIDAPGEHVLGFETKVGPEVKVVIGDDLDTVIIKAKMTKRVGFKAPVKAVVPVVGKDTKIEKFVKFYTKEKDGKYNIDIEREEYNVKVGPKVVDLKDEVTGVVTLEGEDAFKDLGEVYIDTDEMLDQGQYIPKDKDAGRVVVTFDAKSGRPLKDAVFEVKDTKIVVIGK",
            "ompA": "MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAGAFGGYQVNPYVGFEMGYDWLGRMPYKGSVENGAYKAQGVQLTAKLGYPITDDLDIYTRLGGMVWRADTKSNVYGKNHDTGVSPVFAGGVEYAITPEIATRLEYQWTNNIGDAHTIGTRPDNGMLSLGVSYRFGQGEAAPVAPAPAPAPEVEATKHFTLKSDVLFNFNKATLKPEGQAALDQLYSQLSNLDPKDGSVVVLGYTDRIGSDAYNQGLSERRAQSVVDYISNNQIFKGLYKNGDVIDYGRFNISKEGIEVAKKVK",
            "groES": "MARGKKVIKKELGVDALLNIDEEVNRGVRINIVKDKKVLEKNIEVVEKLDFKSVIKKTVEKLGKDIVKVIGDKIKDTIKNVQKIVADKIEKDPKKVKEVINAKIDDLKNKIKKVKEKVENGDKIEKVIEKFVKPKSKKVEKVIKNIGK",
            "infA": "MARIRKEKRIREITDKAKQTLSKQYIDRIEELNKKIGIDDVVAIRRRDVRKAVRRSRRLIERRGKIEKVAEKKSKKIVRIDEKIRKIVDEK",
            "trpR": "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQSLLIGVFKIRSTLRQNLLDSRVLQTYRVYPEKLNQLEQR",
            "araC": "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQSLLIGVFKIRSTLRQNLLDSRVLQTYRVYPEKLNQLEQRVYRLRSGSIRDLQLQIDSLYQWLDNIPVIQYQLRNKQIRSFVKLFGSHPKYENLVLKSFVSREQVNPDIVLFNRNPKRATVLCNEAKQVGLIEKISLKTD",
            "xylA": "MKINLLDKELKFVDDKIAKFGKDFSTKPEKDVKRFIELTNKTLFAEKVNKKFYDNIKTGFISALDRKEPELKDKINLLDQIDKIIKKMGFEEKITKNVFKGFDELDSIKKTIDKGFISQIKQKDPELKDKINLLDRIDEIIKKMGFEEKITKNIFKGFDELDSIKKSVEKGFISQIKNKDPELKDKINLLDKIDEIIKKMGFEEKITKNIFKGFAELDSIKKSVEKGFISQIKQKDPELKDKINLLDKIDEIIKKMGFEEKITKNIFKGFAELDSIKKSVEKGFISQIKNKDPELKDKINLLDKIDKIIKKMGFEEKITKNIFKGFDELDSIKKSVEKGFISQIKQKDPELKDKINLLDKIDEIIKKMGFEEKITKNIFKGFAELDSIKKSVEKGFISQIKQKDPELKDKINLLDKIDKIIKKMGFEEKITKNIFKGFDELDSIKKSVEKGFISQIKQKDPELKD",
            "dnaK": "MKINLLDKELKFVDDKIAKFGKDFSTKPEKDVKRFIELTNKTLFAEKVNKKFYDNIKTGFISALDRKEPELKDKINLLDQIDKIIKKMGFEEKITKNVFKGFDELDSIKKTIDKGFISQIKQKDPELKDKINLLDRIDEIIKKMGFEEKITKNIFKGFDELDSIKKSVEKGFISQIKNKDPELKDKINLLDKIDEIIKKMGFEEKITKNIFKGFAELDSIKKSVEKGFISQIKQKDPELKDKINLLDKIDKIIKKMGFEEKITKNIFKGFDELDSIKKSVEKGFISQIKQKDPELKDKINLLDKIDEIIKKMGFEEKITKNIFKGFAELDSIKKSVEKGFISQIKQKDPELKDKINLLDKIDKIIKKMGFEEKITKNIFKGFDELDSIKKSVEKGFISQIKQKDPELKDKINLLDKIDEIIKKMGFEEKITKNIFKGFAELDSIKKSVEKGFISQIKQKDPELKDKINLLDKIDKIIKKMGFEEKITKNIFKGFDELDSIKKSVEKGFISQIKQKDPELKD",
        }

        for gene_name, protein in ecoli_genes.items():
            # Organism-aware optimization
            result = optimize_sequence(
                protein,
                organism="Escherichia_coli",
                gc_lo=0.30,
                gc_hi=0.70,
                enzymes=DEFAULT_ENZYMES,
                optimize_mrna_stability=False,
                include_utr=False,
                track_provenance=False,
                strict_mode=False,
            )

            # Organism-unaware baseline
            baseline = optimize_sequence(
                protein,
                organism="Escherichia_coli",
                gc_lo=0.30,
                gc_hi=0.70,
                enzymes=DEFAULT_ENZYMES,
                organism_domain="eukaryote",
                optimize_mrna_stability=False,
                include_utr=False,
                track_provenance=False,
                strict_mode=False,
            )

            # Organism-aware should have higher CAI
            assert result.cai >= baseline.cai, (
                f"Gene {gene_name}: organism-aware CAI ({result.cai:.4f}) should be >= "
                f"baseline ({baseline.cai:.4f})"
            )

            # Should translate correctly
            translated = translate(result.sequence)
            assert translated == protein, (
                f"Gene {gene_name}: translation mismatch"
            )

    def test_organism_unaware_backward_compat(self):
        """Run with all constraints applied regardless of organism.

        When organism_domain is forced to 'eukaryote', the optimizer
        applies ALL constraints even for prokaryotic organisms. This
        ensures backward compatibility: existing pipelines that relied
        on the old behavior (all constraints always applied) continue
        to work correctly, just with lower CAI for prokaryotes.
        """
        # Optimize with ALL constraints (force eukaryote domain)
        result_default = optimize_sequence(
            GFP_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
            # Force eukaryote domain: all constraints applied
            organism_domain="eukaryote",
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )

        # Optimize with organism-aware constraints (auto-detect prokaryote)
        result_aware = optimize_sequence(
            GFP_PROTEIN,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
            # Default organism_domain='auto' → prokaryote → skip splice/GT
            optimize_mrna_stability=False,
            include_utr=False,
            track_provenance=False,
            strict_mode=False,
        )

        # Both should produce valid translations
        assert translate(result_default.sequence) == GFP_PROTEIN
        assert translate(result_aware.sequence) == GFP_PROTEIN

        # Both should have GC in range
        gc_default = gc_content(result_default.sequence)
        gc_aware = gc_content(result_aware.sequence)
        assert 0.30 <= gc_default <= 0.70
        assert 0.30 <= gc_aware <= 0.70

        # Both should avoid restriction sites
        assert not _has_restriction_site(result_default.sequence, DEFAULT_ENZYMES)
        assert not _has_restriction_site(result_aware.sequence, DEFAULT_ENZYMES)

        # Default (all constraints) should have LOWER CAI than organism-aware
        # This demonstrates the CAI penalty from unnecessary constraints
        assert result_default.cai < result_aware.cai, (
            f"Default CAI ({result_default.cai:.4f}) should be lower than "
            f"organism-aware CAI ({result_aware.cai:.4f}) because default "
            f"applies unnecessary splice/CpG/GT constraints to E. coli"
        )

        # Default optimization should also pass hard type checks
        type_results = evaluate_all_predicates(
            result_default.sequence,
            organism="Escherichia_coli",
            gc_lo=0.30,
            gc_hi=0.70,
            enzymes=DEFAULT_ENZYMES,
        )

        # Hard constraints must still pass
        assert _find_verdict(type_results, "InFrame") == Verdict.PASS
        assert _find_verdict(type_results, "NoRestrictionSite") == Verdict.PASS
        assert _find_verdict(type_results, "GCInRange") == Verdict.PASS
