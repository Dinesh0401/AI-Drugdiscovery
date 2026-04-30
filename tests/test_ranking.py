"""Tests for chemscreen.ranking (Phase 3)."""

from __future__ import annotations

import pytest

from chemscreen.ranking import filter_top_k, rank_candidates, screen_smiles
from chemscreen.scoring import (
    CandidateScore,
    ScoreComponents,
    ScoreWeights,
    ScoringMethod,
)


def _make(smiles: str, score: float) -> CandidateScore:
    return CandidateScore(
        smiles=smiles,
        components=ScoreComponents(score, score, score, score),
        final_score=score,
        method="linear",
        weights=ScoreWeights(),
    )


def test_rank_candidates_sorts_descending() -> None:
    a = _make("a", 0.5)
    b = _make("b", 0.9)
    c = _make("c", 0.1)

    ranked = rank_candidates([a, b, c])

    assert [r.smiles for r in ranked] == ["b", "a", "c"]
    assert [r.rank for r in ranked] == [1, 2, 3]


def test_rank_candidates_preserves_original_objects() -> None:
    """Ranking should not mutate the input list."""
    a = _make("a", 0.5)
    b = _make("b", 0.9)
    inputs = [a, b]

    ranked = rank_candidates(inputs)

    assert inputs == [a, b]  # original list unchanged
    assert ranked[0].rank == 1
    assert a.rank is None  # original object unchanged


def test_rank_candidates_handles_ties_stably() -> None:
    a = _make("a", 0.5)
    b = _make("b", 0.5)
    ranked = rank_candidates([a, b])
    assert [r.smiles for r in ranked] == ["a", "b"]  # original order kept on tie


def test_rank_candidates_empty() -> None:
    assert rank_candidates([]) == []


def test_screen_smiles_drops_invalid(aspirin: str, caffeine: str) -> None:
    ranked = screen_smiles(
        [aspirin, "not_a_smiles", caffeine, ""],
        weights=ScoreWeights(),
    )
    assert len(ranked) == 2
    assert {r.smiles for r in ranked} == {aspirin, caffeine}
    assert {r.rank for r in ranked} == {1, 2}


def test_screen_smiles_uses_method(aspirin: str, caffeine: str) -> None:
    linear = screen_smiles([aspirin, caffeine], method=ScoringMethod.LINEAR)
    desirability = screen_smiles([aspirin, caffeine], method=ScoringMethod.DESIRABILITY)
    assert all(r.method == "linear" for r in linear)
    assert all(r.method == "desirability" for r in desirability)


def test_filter_top_k_returns_k() -> None:
    candidates = [_make(f"m{i}", i / 10) for i in range(10)]
    top3 = filter_top_k(candidates, k=3)
    assert len(top3) == 3
    assert top3[0].rank == 1
    assert top3[0].final_score == pytest.approx(0.9)


def test_filter_top_k_zero_returns_empty() -> None:
    candidates = [_make("a", 0.5)]
    assert filter_top_k(candidates, k=0) == []


def test_filter_top_k_more_than_available() -> None:
    """Asking for more than exists returns everything ranked."""
    candidates = [_make("a", 0.5), _make("b", 0.7)]
    result = filter_top_k(candidates, k=10)
    assert len(result) == 2
    assert result[0].smiles == "b"


def test_filter_top_k_negative_raises() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        filter_top_k([], k=-1)
