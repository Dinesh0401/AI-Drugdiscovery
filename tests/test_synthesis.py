"""Tests for chemscreen.synthesis (SA score)."""

from __future__ import annotations

from chemscreen.synthesis import is_synthetically_feasible, sa_score


def test_aspirin_easy_to_synthesize(aspirin: str) -> None:
    """Small simple drugs should land in the easy half of the SA scale."""
    score = sa_score(aspirin)
    assert score is not None
    assert 1.0 <= score <= 4.0


def test_ibuprofen_in_drug_range(ibuprofen: str) -> None:
    score = sa_score(ibuprofen)
    assert score is not None
    assert 1.0 <= score <= 5.0


def test_caffeine_in_drug_range(caffeine: str) -> None:
    score = sa_score(caffeine)
    assert score is not None
    assert 1.0 <= score <= 5.0


def test_score_is_finite_for_atorvastatin(atorvastatin: str) -> None:
    """Even larger molecules should produce a finite score in [1,10]."""
    score = sa_score(atorvastatin)
    assert score is not None
    assert 1.0 <= score <= 10.0


def test_feasibility_threshold() -> None:
    assert is_synthetically_feasible(3.0) is True
    assert is_synthetically_feasible(6.0) is True       # boundary
    assert is_synthetically_feasible(7.5) is False
    assert is_synthetically_feasible(5.0, threshold=4.0) is False


def test_sa_score_invalid_returns_none() -> None:
    assert sa_score("") is None
    assert sa_score("not_a_smiles") is None
    assert sa_score(None) is None
