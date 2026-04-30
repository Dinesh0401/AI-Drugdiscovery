"""Candidate ranking and end-to-end batch screening.

Ranking is intentionally separate from scoring: ``score_candidate`` produces
one ``CandidateScore`` (no rank field set); ``rank_candidates`` takes a list
of those and returns a sorted copy with ranks assigned. ``screen_smiles``
is the user-facing convenience that wires scoring + ranking together for a
list of SMILES, dropping invalid inputs along the way.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Sequence

from chemscreen.scoring import (
    CandidateScore,
    ScoreWeights,
    ScoringMethod,
    score_candidate,
)


def rank_candidates(candidates: Sequence[CandidateScore]) -> list[CandidateScore]:
    """Sort by ``final_score`` descending and assign ``rank`` 1..N.

    The input order is preserved for ties (Python's sort is stable).
    """
    ordered = sorted(candidates, key=lambda c: c.final_score, reverse=True)
    return [replace(c, rank=i + 1) for i, c in enumerate(ordered)]


def screen_smiles(
    smiles_list: Sequence[str | None],
    weights: ScoreWeights | None = None,
    method: ScoringMethod = ScoringMethod.LINEAR,
) -> list[CandidateScore]:
    """End-to-end: score every SMILES, drop invalid ones, return ranked results."""
    scored: list[CandidateScore] = []
    for smiles in smiles_list:
        result = score_candidate(smiles, weights=weights, method=method)
        if result is not None:
            scored.append(result)
    return rank_candidates(scored)


def filter_top_k(candidates: Sequence[CandidateScore], k: int) -> list[CandidateScore]:
    """Return the top-k candidates by score (ranked)."""
    if k < 0:
        raise ValueError("k must be non-negative")
    return rank_candidates(candidates)[:k]
