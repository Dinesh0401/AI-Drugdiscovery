"""Tests for chemscreen.druglikeness rules and QED."""

from __future__ import annotations

import pytest

from chemscreen.descriptors import compute_descriptors
from chemscreen.druglikeness import (
    compute_qed,
    evaluate_druglikeness,
    evaluate_ghose,
    evaluate_lipinski,
    evaluate_veber,
)


def test_aspirin_passes_lipinski(aspirin: str) -> None:
    result = evaluate_druglikeness(aspirin)
    assert result is not None
    assert result.lipinski.passed is True
    assert result.lipinski.violations == []


def test_caffeine_passes_lipinski(caffeine: str) -> None:
    result = evaluate_druglikeness(caffeine)
    assert result is not None
    assert result.lipinski.passed is True


def test_ibuprofen_passes_lipinski(ibuprofen: str) -> None:
    result = evaluate_druglikeness(ibuprofen)
    assert result is not None
    assert result.lipinski.passed is True


def test_atorvastatin_fails_lipinski(atorvastatin: str) -> None:
    """Lipinski must fire on atorvastatin and report the violating MW number."""
    result = evaluate_druglikeness(atorvastatin)
    assert result is not None
    assert result.lipinski.passed is False
    # Phase 5 explainability hook: violation strings carry the actual numbers.
    assert any("MW=" in v and "> 500" in v for v in result.lipinski.violations)


def test_qed_caffeine_in_range(caffeine: str) -> None:
    qed = compute_qed(caffeine)
    assert qed is not None
    assert 0.0 < qed < 1.0
    # Published QED for caffeine sits around 0.54; allow generous tolerance
    # for RDKit-version drift in QED's underlying ADS coefficients.
    assert qed == pytest.approx(0.54, abs=0.1)


def test_qed_invalid_returns_none() -> None:
    assert compute_qed("") is None
    assert compute_qed("not_a_smiles") is None


def test_lipinski_violation_messages_have_numbers(atorvastatin: str) -> None:
    d = compute_descriptors(atorvastatin)
    assert d is not None
    result = evaluate_lipinski(d)
    assert not result.passed
    # Every violation string must be feature-grounded (contain "=" and ">")
    for violation in result.violations:
        assert "=" in violation
        assert ">" in violation


def test_veber_passes_aspirin(aspirin: str) -> None:
    d = compute_descriptors(aspirin)
    assert d is not None
    assert evaluate_veber(d).passed is True


def test_ghose_filter_runs_without_error(caffeine: str, atorvastatin: str) -> None:
    for smiles in (caffeine, atorvastatin):
        d = compute_descriptors(smiles)
        assert d is not None
        result = evaluate_ghose(d)
        assert isinstance(result.passed, bool)
        assert isinstance(result.violations, list)


def test_druglikeness_bundles_all_rules(aspirin: str) -> None:
    result = evaluate_druglikeness(aspirin)
    assert result is not None
    assert result.lipinski.name == "Lipinski Rule of Five"
    assert result.veber.name == "Veber Rule"
    assert result.ghose.name == "Ghose Filter"


def test_druglikeness_invalid_returns_none() -> None:
    assert evaluate_druglikeness("") is None
    assert evaluate_druglikeness("not_a_smiles") is None
