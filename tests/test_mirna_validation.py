"""Tests for miRNA validation against published literature."""
from biocompiler.validation.mirna_validation import validate_mirna_predictions


def test_mirna_predictions_against_literature():
    """Verify miRNA predictions match published experimental data."""
    errors = validate_mirna_predictions()
    assert errors == [], f"miRNA validation failures: {'; '.join(errors)}"
