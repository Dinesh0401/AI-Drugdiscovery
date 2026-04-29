"""Tests for chemscreen.reference_data formula -> SMILES lookup."""

from __future__ import annotations

from chemscreen.descriptors import is_valid_smiles
from chemscreen.reference_data import FORMULA_TO_SMILES, lookup_formula


def test_lookup_known_formulas() -> None:
    assert lookup_formula("H2O") == "O"
    assert lookup_formula("CO2") == "O=C=O"


def test_lookup_unknown_returns_none() -> None:
    assert lookup_formula("UNOBTAINIUM") is None


def test_lookup_strips_whitespace() -> None:
    assert lookup_formula("  H2O  ") == "O"


def test_every_reference_smiles_is_valid() -> None:
    """Every entry in the table must round-trip through RDKit."""
    for formula, smiles in FORMULA_TO_SMILES.items():
        assert is_valid_smiles(smiles), f"Invalid SMILES for {formula}: {smiles!r}"
