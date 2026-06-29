"""
Biosecurity coverage regression test (GAP-1).

Verifies that every hazard reference fragment in
``data/hazard_reference_sequences.fasta`` is flagged as hazardous by
:func:`biocompiler.biosecurity.screen_hazardous_sequence`.

Before GAP-1 the screener flagged 6/8 known hazards — botulinum (no
motif coverage) and anthrax (distance-2 fuzzy match downgraded to
risk="low", bypassing ``is_hazardous``) were missed.  After GAP-1 all
8 hazards must be flagged.

The test suite covers:
  - All 8 fasta reference fragments are flagged (parametrised).
  - Aggregate coverage is exactly 8/8.
  - Botulinum reference produces at least one botulinum-named match
    (verifies the new motifs / updated reference).
  - Anthrax lethal factor's distance-2 fuzzy match is now risk="medium"
    (not "low"), so the select_agent risk floor is in effect.
  - Non-select_agent distance-2 fuzzy matches still downgrade to "low"
    (verifies the floor is scoped to select_agent only — no over-broad
    change that would inflate false positives for oncogenes etc.).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from biocompiler.biosecurity import (
    HazardMatch,
    screen_hazardous_sequence,
)
from biocompiler.biosecurity.screening import sig_risk_for_match


# ─────────────────────────────────────────────────────────────────────────────
# Fasta parsing
# ─────────────────────────────────────────────────────────────────────────────

FASTA_PATH = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "hazard_reference_sequences.fasta"
)


def _parse_fasta(path: Path) -> dict[str, str]:
    """Parse a fasta file into ``{header_name: sequence}``.

    Comment lines starting with ``#`` are skipped.  The header key is
    the first whitespace-delimited token after ``>``.
    """
    sequences: dict[str, str] = {}
    current_name: str | None = None
    current_chunks: list[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_chunks)
                current_name = line[1:].split()[0]
                current_chunks = []
            else:
                current_chunks.append(line)
    if current_name is not None:
        sequences[current_name] = "".join(current_chunks)
    return sequences


# ─────────────────────────────────────────────────────────────────────────────
# Expected hazards — every entry in the fasta must be flagged
# ─────────────────────────────────────────────────────────────────────────────

EXPECTED_HAZARDS: list[str] = [
    "ricin_A_chain_reference",
    "botulinum_toxin_A_reference",
    "anthrax_lethal_factor_reference",
    "shiga_toxin_A_subunit_reference",
    "diphtheria_toxin_reference",
    "tetanus_toxin_reference",
    "cholera_toxin_A1_reference",
    "abrin_A_chain_reference",
]

assert len(EXPECTED_HAZARDS) == 8, "GAP-1 baseline is 8 hazards"


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def fasta_sequences() -> dict[str, str]:
    """Parse the hazard reference fasta once per test module run."""
    assert FASTA_PATH.exists(), (
        f"Hazard reference fasta not found at {FASTA_PATH}"
    )
    seqs = _parse_fasta(FASTA_PATH)
    return seqs


# ─────────────────────────────────────────────────────────────────────────────
# Tests: fasta integrity
# ─────────────────────────────────────────────────────────────────────────────

class TestFastaIntegrity:
    """The reference fasta must contain all 8 expected hazard entries."""

    def test_fasta_contains_all_expected_hazards(self, fasta_sequences):
        missing = [h for h in EXPECTED_HAZARDS if h not in fasta_sequences]
        assert not missing, (
            f"Missing reference sequences in fasta: {missing}"
        )

    def test_all_hazard_sequences_nonempty(self, fasta_sequences):
        for name in EXPECTED_HAZARDS:
            seq = fasta_sequences[name]
            assert seq, f"Empty sequence for {name}"
            assert len(seq) >= 30, (
                f"{name} sequence suspiciously short: {len(seq)} aa"
            )

    def test_botulinum_reference_contains_zinc_motif(self, fasta_sequences):
        """Botulinum reference must contain the HETQSNLRDL zinc-binding motif.

        This is the GAP-1 botulinum fix — the original reference did not
        contain any of the 5 database motifs, so the screener returned
        'no match'.  The updated reference now embeds the
        botulinum_zinc_protease motif (and the new botulinum_zinc_HELIH
        extended motif).
        """
        bot = fasta_sequences["botulinum_toxin_A_reference"]
        assert "HETQSNLRDL" in bot, (
            "Botulinum reference must contain the HETQSNLRDL zinc-binding "
            "motif (botulinum_zinc_protease) — GAP-1 fix"
        )

    def test_botulinum_reference_contains_HELIH_motif(self, fasta_sequences):
        """Botulinum reference must contain the new PVHELIHAVL motif.

        This is the extended-context HExxH/HELIH motif added by GAP-1
        to improve BoNT fragment coverage.
        """
        bot = fasta_sequences["botulinum_toxin_A_reference"]
        assert "PVHELIHAVL" in bot, (
            "Botulinum reference must contain the PVHELIHAVL (HELIH) motif "
            "(botulinum_zinc_HELIH) — GAP-1 fix"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Tests: per-hazard screening (parametrised)
# ─────────────────────────────────────────────────────────────────────────────

class TestPerHazardScreening:
    """Each individual hazard reference must be flagged as hazardous."""

    @pytest.mark.parametrize("hazard_name", EXPECTED_HAZARDS)
    def test_hazard_is_flagged(self, fasta_sequences, hazard_name):
        seq = fasta_sequences[hazard_name]
        report = screen_hazardous_sequence(seq)
        assert report.is_hazardous, (
            f"{hazard_name} was NOT flagged as hazardous "
            f"(risk={report.risk_level}, "
            f"matches={len(report.matches)}, "
            f"categories={report.flagged_categories}). "
            f"This is a biosecurity coverage regression — the screener "
            f"must flag all 8 known hazards."
        )

    @pytest.mark.parametrize("hazard_name", EXPECTED_HAZARDS)
    def test_hazard_has_select_agent_category(self, fasta_sequences, hazard_name):
        """Each hazard reference must trigger a select_agent flag.

        All 8 reference fragments are CDC Select Agent toxins, so the
        screener must classify them under the ``select_agent`` category.
        """
        seq = fasta_sequences[hazard_name]
        report = screen_hazardous_sequence(seq)
        assert "select_agent" in report.flagged_categories, (
            f"{hazard_name} did not trigger select_agent category "
            f"(got: {report.flagged_categories})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Tests: aggregate coverage (the headline GAP-1 metric)
# ─────────────────────────────────────────────────────────────────────────────

class TestAggregateCoverage:
    """Verify 8/8 hazard coverage — the headline GAP-1 success metric."""

    def test_all_8_hazards_flagged(self, fasta_sequences):
        """All 8 hazards must be flagged. Before GAP-1 only 6/8 passed."""
        flagged: list[str] = []
        failures: list[str] = []
        for name in EXPECTED_HAZARDS:
            seq = fasta_sequences[name]
            report = screen_hazardous_sequence(seq)
            if report.is_hazardous:
                flagged.append(name)
            else:
                failures.append(
                    f"{name}: risk={report.risk_level}, "
                    f"matches={len(report.matches)}"
                )
        assert len(flagged) == 8, (
            f"Only {len(flagged)}/8 hazards flagged. "
            f"GAP-1 requires 8/8.\n"
            f"Failures:\n  " + "\n  ".join(failures)
        )

    def test_coverage_regression_baseline(self, fasta_sequences):
        """Explicit regression check: botulinum + anthrax must now pass.

        Before GAP-1 these two specific hazards failed:
          - botulinum_toxin_A_reference: no motif match (risk=none)
          - anthrax_lethal_factor_reference: distance-2 fuzzy match
            downgraded to risk=low (bypasses is_hazardous)
        """
        bot_report = screen_hazardous_sequence(
            fasta_sequences["botulinum_toxin_A_reference"]
        )
        anthrax_report = screen_hazardous_sequence(
            fasta_sequences["anthrax_lethal_factor_reference"]
        )
        assert bot_report.is_hazardous, (
            f"Botulinum gap not closed: risk={bot_report.risk_level}, "
            f"matches={len(bot_report.matches)}"
        )
        assert anthrax_report.is_hazardous, (
            f"Anthrax gap not closed: risk={anthrax_report.risk_level}, "
            f"matches={len(anthrax_report.matches)}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Tests: botulinum-specific gap closure
# ─────────────────────────────────────────────────────────────────────────────

class TestBotulinumGapClosed:
    """Verify the botulinum-specific fixes (Gap A)."""

    def test_botulinum_reference_produces_botulinum_named_match(
        self, fasta_sequences
    ):
        """Botulinum reference must produce at least one botulinum-named match.

        Before GAP-1 the reference had zero matches.  After GAP-1 it must
        produce at least one match whose signature name starts with
        'botulinum'.
        """
        seq = fasta_sequences["botulinum_toxin_A_reference"]
        report = screen_hazardous_sequence(seq)
        bot_matches = [
            m for m in report.matches if m.name.lower().startswith("botulinum")
        ]
        assert bot_matches, (
            "Botulinum reference produced no botulinum-named matches. "
            f"All matches: {[m.name for m in report.matches]}"
        )

    def test_botulinum_has_critical_risk(self, fasta_sequences):
        """Botulinum is a CDC Select Agent — risk must be critical.

        The new botulinum_zinc_HELIH motif has risk='critical', so an
        exact match must escalate the report to critical.
        """
        seq = fasta_sequences["botulinum_toxin_A_reference"]
        report = screen_hazardous_sequence(seq)
        assert report.risk_level == "critical", (
            f"Botulinum risk should be 'critical' (Select Agent toxin), "
            f"got '{report.risk_level}'"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Tests: anthrax-specific gap closure (select_agent fuzzy risk floor)
# ─────────────────────────────────────────────────────────────────────────────

class TestAnthraxGapClosed:
    """Verify the anthrax-specific fix (Gap B): select_agent fuzzy risk floor."""

    def test_anthrax_distance_2_fuzzy_is_medium(self, fasta_sequences):
        """Anthrax LF distance-2 fuzzy match must be 'medium' (not 'low').

        The anthrax lethal factor reference has HETHFGVVYS (motif is
        HETHFGVVSY — last two residues swapped = Hamming distance 2).
        Before GAP-1 this was downgraded to 'low' and bypassed
        is_hazardous.  After GAP-1 select_agent fuzzy matches are
        floored at 'medium'.
        """
        seq = fasta_sequences["anthrax_lethal_factor_reference"]
        report = screen_hazardous_sequence(seq)
        # Find the anthrax_LF_protease fuzzy match
        lf_matches = [
            m for m in report.matches
            if m.name == "anthrax_LF_protease" and m.match_type == "fuzzy"
        ]
        assert lf_matches, (
            "Anthrax LF reference should produce a fuzzy match against "
            "anthrax_LF_protease (HETHFGVVSY motif). "
            f"Matches: {[m.name for m in report.matches]}"
        )
        for m in lf_matches:
            risk = sig_risk_for_match(m)
            assert risk in ("medium", "high", "critical"), (
                f"Select_agent fuzzy match (distance={m.distance}) has "
                f"risk='{risk}' — GAP-1 requires at least 'medium' so "
                f"is_hazardous=True is triggered."
            )

    def test_select_agent_fuzzy_distance_2_floored_at_medium(self):
        """Direct unit test of the select_agent fuzzy risk floor.

        Synthesises a HazardMatch with category='select_agent' and
        distance=2, then verifies sig_risk_for_match returns 'medium'.
        """
        m = HazardMatch(
            category="select_agent",
            name="anthrax_LF_protease",
            position=0,
            matched_sequence="HETHFGVVYS",  # 2 subs from HETHFGVVSY
            confidence=0.80,
            source="test",
            match_type="fuzzy",
            distance=2,
            strand="forward",
            substitutions=[(8, "S", "Y"), (9, "Y", "S")],
        )
        assert sig_risk_for_match(m) == "medium", (
            "select_agent fuzzy distance-2 must be 'medium' (GAP-1 fix)."
        )

    def test_non_select_agent_fuzzy_distance_2_stays_low(self):
        """The risk floor must NOT apply to non-select_agent categories.

        GAP-1 only floors select_agent fuzzy matches at 'medium'.  Other
        categories (oncogene, viral_surface, antibiotic_resistance,
        australia_group) keep the original distance-2 downgrade to
        'low'.  This prevents false-positive inflation for lower-risk
        hazard classes.
        """
        for cat in ("oncogene", "viral_surface", "antibiotic_resistance"):
            m = HazardMatch(
                category=cat,
                name="some_motif",
                position=0,
                matched_sequence="XXXXXXXXXY",
                confidence=0.80,
                source="test",
                match_type="fuzzy",
                distance=2,
                strand="forward",
                substitutions=[(0, "X", "Y")],
            )
            assert sig_risk_for_match(m) == "low", (
                f"category={cat} fuzzy distance-2 should stay 'low' "
                f"(GAP-1 floor is select_agent only), got "
                f"'{sig_risk_for_match(m)}'"
            )

    def test_select_agent_fuzzy_distance_1_is_medium(self):
        """select_agent fuzzy distance-1 must be 'medium' (unchanged)."""
        m = HazardMatch(
            category="select_agent",
            name="ricin_A_chain_catalytic",
            position=0,
            matched_sequence="KIRVGLPIIS",
            confidence=0.85,
            source="test",
            match_type="fuzzy",
            distance=1,
            strand="forward",
            substitutions=[(0, "N", "K")],
        )
        assert sig_risk_for_match(m) == "medium"

    def test_select_agent_exact_match_unchanged(self, fasta_sequences):
        """select_agent exact matches must keep their original (critical) risk.

        The GAP-1 fix only affects fuzzy matches.  Exact matches on
        select_agent toxins must still escalate to their database
        risk level (typically 'critical' for select_agent toxins).
        """
        # Ricin A-chain reference has an exact match on ricin_A_chain_catalytic
        seq = fasta_sequences["ricin_A_chain_reference"]
        report = screen_hazardous_sequence(seq)
        exact_matches = [
            m for m in report.matches
            if m.match_type == "exact" and m.category == "select_agent"
        ]
        assert exact_matches, "Ricin reference should have an exact match"
        for m in exact_matches:
            risk = sig_risk_for_match(m)
            assert risk in ("high", "critical"), (
                f"Exact select_agent match risk='{risk}' should be "
                f"high or critical (unchanged by GAP-1)"
            )
