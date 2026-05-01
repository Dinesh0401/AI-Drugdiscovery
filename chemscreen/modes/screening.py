"""Screening mode: batch score + rank a list of SMILES with rejection accounting.

Wraps Phase 3's ``screen_smiles`` to produce a structured ``ScreeningReport``
that also tracks how many inputs were dropped as invalid SMILES — useful for
real-world usage where the input may contain typos or rendering errors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from chemscreen.descriptors import is_valid_smiles
from chemscreen.ranking import rank_candidates
from chemscreen.scoring import (
    CandidateScore,
    ScoreWeights,
    ScoringMethod,
    score_candidate,
)


@dataclass(frozen=True)
class ScreeningReport:
    """Structured outcome of a batch screening run."""

    candidates: list[CandidateScore]
    n_input: int
    n_scored: int
    n_rejected: int
    method: str
    weights: ScoreWeights
    require_lipinski: bool


def screen_batch(
    smiles_list: Sequence[str | None],
    weights: ScoreWeights | None = None,
    method: ScoringMethod = ScoringMethod.LINEAR,
    top_k: int | None = None,
    require_lipinski: bool = False,
) -> ScreeningReport:
    """Score and rank a batch of SMILES, returning a ``ScreeningReport``.

    Args:
        smiles_list: candidate SMILES strings.
        weights: scoring weights; defaults to ``ScoreWeights()``.
        method: linear (default) or desirability.
        top_k: keep only the top ``k`` after ranking.
        require_lipinski: when True, drop any candidate that fails the
            Lipinski Rule of Five before ranking. Useful for early-stage
            screening where Lipinski-violators aren't worth carrying.
    """
    weights = weights or ScoreWeights()

    n_input = len(smiles_list)
    scored: list[CandidateScore] = []
    n_rejected = 0

    for smiles in smiles_list:
        if not is_valid_smiles(smiles):
            n_rejected += 1
            continue
        result = score_candidate(smiles, weights=weights, method=method)
        if result is None:
            n_rejected += 1
            continue
        if require_lipinski and result.components.lipinski < 1.0:
            # any Lipinski violation drops lipinski_component below 1.0
            continue
        scored.append(result)

    ranked = rank_candidates(scored)
    if top_k is not None:
        ranked = ranked[: max(0, top_k)]

    return ScreeningReport(
        candidates=ranked,
        n_input=n_input,
        n_scored=len(scored),
        n_rejected=n_rejected,
        method=method.value,
        weights=weights,
        require_lipinski=require_lipinski,
    )
