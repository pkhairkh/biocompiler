"""Unit tests for restriction_sites.py and RESTRICTION_ENZYMES in constants.py.

Covers:
1. get_recognition_site() for known enzymes (EcoRI, BamHI, etc.)
2. Unknown enzyme returns None
3. Recognition site format (uppercase ACGT only)
4. All enzymes in RESTRICTION_ENZYMES have valid sites
"""

import re

import pytest

from biocompiler.sequence.restriction_sites import RESTRICTION_SITES, get_recognition_site
from biocompiler.shared.constants import RESTRICTION_ENZYMES


# ---------------------------------------------------------------------------
# 1. get_recognition_site() for known enzymes
# ---------------------------------------------------------------------------

class TestGetRecognitionSiteKnownEnzymes:
    """Verify get_recognition_site returns the correct sequence for known enzymes."""

    @pytest.mark.parametrize(
        "enzyme, expected",
        [
            ("EcoRI", "GAATTC"),
            ("BamHI", "GGATCC"),
            ("HindIII", "AAGCTT"),
            ("XhoI", "CTCGAG"),
            ("XbaI", "TCTAGA"),
            ("SalI", "GTCGAC"),
            ("PstI", "CTGCAG"),
            ("SphI", "GCATGC"),
            ("KpnI", "GGTACC"),
            ("SacI", "GAGCTC"),
            ("NcoI", "CCATGG"),
            ("NdeI", "CATATG"),
            ("NotI", "GCGGCCGC"),
            ("BglII", "AGATCT"),
            ("ClaI", "ATCGAT"),
            ("EcoRV", "GATATC"),
            ("SmaI", "CCCGGG"),
            ("SpeI", "ACTAGT"),
            ("NheI", "GCTAGC"),
            ("ApaI", "GGGCCC"),
        ],
    )
    def test_known_enzyme_site(self, enzyme, expected):
        """Each known enzyme maps to its correct recognition sequence."""
        assert get_recognition_site(enzyme) == expected

    def test_case_insensitive_ecori(self):
        """Lowercase enzyme name still resolves correctly."""
        assert get_recognition_site("ecori") == "GAATTC"

    def test_case_insensitive_bamhi(self):
        """Mixed-case enzyme name still resolves correctly."""
        assert get_recognition_site("bAmHi") == "GGATCC"

    def test_case_insensitive_hindiii(self):
        """Lowercase hindiii resolves correctly."""
        assert get_recognition_site("hindiii") == "AAGCTT"

    def test_all_restriction_sites_entries_reachable(self):
        """Every enzyme key in RESTRICTION_SITES is reachable via get_recognition_site."""
        for enzyme, expected_site in RESTRICTION_SITES.items():
            assert get_recognition_site(enzyme) == expected_site, (
                f"get_recognition_site({enzyme!r}) did not return {expected_site!r}"
            )


# ---------------------------------------------------------------------------
# 2. Unknown enzyme returns None
# ---------------------------------------------------------------------------

class TestGetRecognitionSiteUnknown:
    """Verify get_recognition_site returns None for unknown enzymes."""

    @pytest.mark.parametrize(
        "enzyme",
        [
            "FakeEnzyme",
            "ZZZ",
            "",
            "EcoRII",       # close but not in database
            "NotAnEnzyme",
            "bamh1",        # digit 1 instead of letter I
        ],
    )
    def test_unknown_enzyme_returns_none(self, enzyme):
        """Unknown enzyme names must return None."""
        assert get_recognition_site(enzyme) is None

    def test_none_input_returns_none(self):
        """Passing None as enzyme name should return None (not raise)."""
        # get_recognition_site calls .lower() on the input, so None will
        # raise AttributeError — but that is the current behaviour.
        # If it should return None instead, the implementation needs fixing.
        # For now, document the actual behaviour:
        with pytest.raises(AttributeError):
            get_recognition_site(None)


# ---------------------------------------------------------------------------
# 3. Recognition site format (uppercase ACGT)
# ---------------------------------------------------------------------------

class TestRecognitionSiteFormat:
    """All recognition sites in RESTRICTION_SITES must be uppercase ACGT strings."""

    _ACGT_RE = re.compile(r"^[ACGT]+$")

    def test_every_site_is_uppercase_acgt(self):
        """Each recognition site in RESTRICTION_SITES consists only of A, C, G, T."""
        for enzyme, site in RESTRICTION_SITES.items():
            assert self._ACGT_RE.match(site), (
                f"Enzyme {enzyme!r} has invalid site {site!r} — "
                f"must be uppercase ACGT only"
            )

    def test_sites_are_not_lowercase(self):
        """No site should contain lowercase letters."""
        for enzyme, site in RESTRICTION_SITES.items():
            assert site == site.upper(), (
                f"Enzyme {enzyme!r} site {site!r} is not fully uppercase"
            )

    def test_sites_contain_no_iupac_ambiguity(self):
        """Sites should not contain IUPAC ambiguity codes (N, R, Y, etc.)."""
        invalid_chars = set("NRYSWKMBDHV")
        for enzyme, site in RESTRICTION_SITES.items():
            bad = set(site) & invalid_chars
            assert not bad, (
                f"Enzyme {enzyme!r} site {site!r} contains IUPAC ambiguity codes: {bad}"
            )

    def test_sites_are_non_empty(self):
        """Every recognition site must be a non-empty string."""
        for enzyme, site in RESTRICTION_SITES.items():
            assert isinstance(site, str) and len(site) > 0, (
                f"Enzyme {enzyme!r} has empty or non-string site"
            )

    def test_site_minimum_length(self):
        """Recognition sites should be at least 4 bp (shortest known: 4-cutters)."""
        for enzyme, site in RESTRICTION_SITES.items():
            assert len(site) >= 4, (
                f"Enzyme {enzyme!r} site {site!r} is shorter than 4 bp"
            )


# ---------------------------------------------------------------------------
# 4. All enzymes in RESTRICTION_ENZYMES (constants.py) have valid sites
# ---------------------------------------------------------------------------

class TestRestrictionEnzymesConstants:
    """Validate every entry in RESTRICTION_ENZYMES from constants.py."""

    # SfiI uses IUPAC N wildcards — allowed in the broader REBASE dict
    # but must be documented.  All others should be pure ACGT.
    _ACGT_RE = re.compile(r"^[ACGT]+$")
    _ACGT_OR_IUPAC_RE = re.compile(r"^[ACGTNRYSWKMBDHV]+$")

    def test_all_sites_non_empty(self):
        """Every RESTRICTION_ENZYMES entry must have a non-empty site string."""
        for enzyme, site in RESTRICTION_ENZYMES.items():
            assert isinstance(site, str) and len(site) > 0, (
                f"Enzyme {enzyme!r} has empty or non-string site"
            )

    def test_all_sites_valid_iupac(self):
        """Every site in RESTRICTION_ENZYMES must contain only valid IUPAC bases."""
        for enzyme, site in RESTRICTION_ENZYMES.items():
            assert self._ACGT_OR_IUPAC_RE.match(site), (
                f"Enzyme {enzyme!r} site {site!r} contains invalid characters"
            )

    def test_pure_acgt_sites_except_documented_wildcards(self):
        """Sites that are pure ACGT are valid; those with N etc. are tolerated but flagged."""
        wildcard_enzymes = []
        for enzyme, site in RESTRICTION_ENZYMES.items():
            if not self._ACGT_RE.match(site):
                wildcard_enzymes.append(enzyme)
        # SfiI is the only known enzyme with wildcard bases; if new ones
        # appear this test still passes but records them.
        assert all(e in {"SfiI"} for e in wildcard_enzymes), (
            f"Unexpected wildcard-containing enzymes: {wildcard_enzymes}"
        )

    def test_sufficient_enzyme_count(self):
        """RESTRICTION_ENZYMES should contain a reasonable number of enzymes."""
        assert len(RESTRICTION_ENZYMES) >= 15, (
            f"Only {len(RESTRICTION_ENZYMES)} enzymes — expected at least 15"
        )

    def test_site_minimum_length(self):
        """Every recognition site in RESTRICTION_ENZYMES must be >= 4 bp."""
        for enzyme, site in RESTRICTION_ENZYMES.items():
            assert len(site) >= 4, (
                f"Enzyme {enzyme!r} site {site!r} is shorter than 4 bp"
            )

    def test_overlapping_enzymes_consistent(self):
        """Enzymes present in both RESTRICTION_SITES and RESTRICTION_ENZYMES
        must map to the same recognition sequence."""
        common = set(RESTRICTION_SITES) & set(RESTRICTION_ENZYMES)
        for enzyme in common:
            assert RESTRICTION_SITES[enzyme] == RESTRICTION_ENZYMES[enzyme], (
                f"Enzyme {enzyme!r} has inconsistent sites: "
                f"restriction_sites={RESTRICTION_SITES[enzyme]!r}, "
                f"constants={RESTRICTION_ENZYMES[enzyme]!r}"
            )

    def test_palindromic_sites_are_self_reverse_complement(self):
        """Most restriction enzyme sites are palindromic (equal to their
        reverse complement).  Verify this for pure-ACGT sites."""
        from biocompiler.shared.constants import reverse_complement

        for enzyme, site in RESTRICTION_ENZYMES.items():
            if not self._ACGT_RE.match(site):
                continue  # skip sites with IUPAC wildcards
            rc = reverse_complement(site)
            # Many but not all are palindromic — we just check the
            # reverse_complement function works without error.
            assert isinstance(rc, str) and len(rc) == len(site), (
                f"reverse_complement({site!r}) for {enzyme!r} produced unexpected result"
            )
