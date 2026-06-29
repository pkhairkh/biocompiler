"""Test BioCompiler Restriction Site Database — recognition site lookup."""

import pytest
from biocompiler.sequence.restriction_sites import get_recognition_site, RESTRICTION_SITES


class TestRecognitionSiteLookup:
    """Tests for restriction enzyme recognition site lookup."""

    def test_ecori_lookup(self):
        """EcoRI recognition site is GAATTC."""
        assert get_recognition_site("EcoRI") == "GAATTC"

    def test_bamhi_lookup(self):
        """BamHI recognition site is GGATCC."""
        assert get_recognition_site("BamHI") == "GGATCC"

    def test_case_insensitive_lookup(self):
        """Lookup is case-insensitive — lowercase and PascalCase both work."""
        # Exact PascalCase key works
        assert get_recognition_site("EcoRI") == "GAATTC"
        # Lowercase also works (case-insensitive lookup)
        assert get_recognition_site("ecori") == "GAATTC"

    def test_unknown_enzyme_returns_none(self):
        """Unknown enzyme should return None."""
        assert get_recognition_site("FakeEnzyme") is None

    def test_all_sites_non_empty(self):
        """All entries in RESTRICTION_SITES should have non-empty recognition sequences."""
        for enzyme, site in RESTRICTION_SITES.items():
            assert len(site) > 0, f"Enzyme {enzyme} has empty recognition site"
            assert all(b in "ATGCU" for b in site.upper()), (
                f"Enzyme {enzyme} site {site!r} contains invalid bases"
            )

    def test_sufficient_enzymes_registered(self):
        """There should be a reasonable number of enzymes in the database."""
        assert len(RESTRICTION_SITES) >= 15
