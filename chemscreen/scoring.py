"""Multi-objective scoring engine.

Two aggregation methods, both ship in the same module:

  - LINEAR (default): the standard weighted sum
        final = w_qed * QED + w_lip * LipinskiCompliance
              + w_sa * SyntheticEase + w_tox * ToxicityPenalty

  - DESIRABILITY (Derringer & Suich, 1980): weighted geometric mean of the
    same components. Penalizes molecules that fail any single component
    instead of letting strong components mask weak ones — the standard
    multi-objective method in process chemistry and computational drug
    design literature.

Each component is a function that maps a raw measurement to ``[0, 1]``
where 1 is best. A ``CandidateScore`` carries the components, the
aggregated score, the method used, and (after ranking) the rank.

All public functions return ``None`` on invalid input rather than raising,
matching the contract used everywhere else in chemscreen.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Sequence

from chemscreen.descriptors import compute_descriptors
from chemscreen.druglikeness import evaluate_druglikeness
from chemscreen.synthesis import sa_score


class ScoringMethod(str, Enum):
    """How to combine the four components into a single score."""

    LINEAR = "linear"
    DESIRABILITY = "desirability"


@dataclass(frozen=True)
class ScoreWeights:
    """Weights for the four scoring components. Must sum to 1.0."""

    qed: float = 0.3
    lipinski: float = 0.2
    synthesis: float = 0.2
    toxicity: float = 0.3

    def __post_init__(self) -> None:
        if min(self.qed, self.lipinski, self.synthesis, self.toxicity) < 0:
            raise ValueError("ScoreWeights cannot be negative")
        total = self.qed + self.lipinski + self.synthesis + self.toxicity
        if abs(total - 1.0) > 1e-6:
            raise ValueError(
                f"ScoreWeights must sum to 1.0, got {total:.6f}"
            )


@dataclass(frozen=True)
class ScoreComponents:
    """Per-component scores, each in [0, 1] where 1 is best."""

    qed: float
    lipinski: float
    synthesis: float
    toxicity: float


@dataclass(frozen=True)
class CandidateScore:
    """Outcome of scoring one molecule. ``rank`` is filled by ``rank_candidates``."""

    smiles: str
    components: ScoreComponents
    final_score: float
    method: str
    weights: ScoreWeights
    rank: int | None = None


# --- Component normalizers --------------------------------------------------


def lipinski_component(num_violations: int) -> float:
    """Map number of Lipinski violations (0..4) to ``[0, 1]`` (1 = no violations)."""
    return max(0.0, 1.0 - num_violations / 4.0)


def synthesis_component(sa: float) -> float:
    """Map SA score (1=easy..10=hard) to ``[0, 1]`` (1 = easy to make)."""
    return max(0.0, min(1.0, (10.0 - sa) / 9.0))


def toxicity_component(risk_scores: Sequence[float] | None) -> float:
    """Convert mean toxicity probability to a desirability ``[0, 1]``.

    Empty list (no toxicity data) returns 0.5 — neither rewarded nor
    punished, since we have no evidence either way.
    """
    if not risk_scores:
        return 0.5
    mean = sum(risk_scores) / len(risk_scores)
    return max(0.0, min(1.0, 1.0 - mean))


# --- Aggregators ------------------------------------------------------------


def _linear(c: ScoreComponents, w: ScoreWeights) -> float:
    return (
        w.qed * c.qed
        + w.lipinski * c.lipinski
        + w.synthesis * c.synthesis
        + w.toxicity * c.toxicity
    )


def _desirability(c: ScoreComponents, w: ScoreWeights) -> float:
    """Weighted geometric mean (Derringer–Suich). Any zero component drives total to ~0."""
    eps = 1e-9
    log_d = (
        w.qed * math.log(max(c.qed, eps))
        + w.lipinski * math.log(max(c.lipinski, eps))
        + w.synthesis * math.log(max(c.synthesis, eps))
        + w.toxicity * math.log(max(c.toxicity, eps))
    )
    return math.exp(log_d)


def aggregate(
    components: ScoreComponents,
    weights: ScoreWeights,
    method: ScoringMethod = ScoringMethod.LINEAR,
) -> float:
    """Combine per-component scores into a single ``[0, 1]`` final score."""
    if method == ScoringMethod.LINEAR:
        return _linear(components, weights)
    if method == ScoringMethod.DESIRABILITY:
        return _desirability(components, weights)
    raise ValueError(f"Unknown ScoringMethod: {method}")


# --- End-to-end -------------------------------------------------------------


def _resolve_tox_risks(smiles: str) -> list[float]:
    """Best-effort toxicity lookup; returns empty list if models aren't available."""
    try:
        from chemscreen.toxicity import predict_all_toxicity

        result = predict_all_toxicity(smiles)
    except (ImportError, FileNotFoundError):
        return []
    if not result:
        return []
    return [prediction.risk_score for prediction in result.values()]


def score_candidate(
    smiles: str | None,
    weights: ScoreWeights | None = None,
    method: ScoringMethod = ScoringMethod.LINEAR,
    tox_risk_scores: Sequence[float] | None = None,
) -> CandidateScore | None:
    """Compute the multi-objective score for one SMILES.

    If ``tox_risk_scores`` is None, the function tries to compute toxicity
    via the trained endpoint models. If models aren't present, the
    toxicity component falls back to neutral (0.5) so scoring still
    produces a usable number for fresh checkouts.

    Returns None if the SMILES is invalid.
    """
    weights = weights or ScoreWeights()

    descriptors = compute_descriptors(smiles)
    if descriptors is None:
        return None

    druglikeness = evaluate_druglikeness(smiles)
    sa = sa_score(smiles)
    if druglikeness is None or sa is None:
        return None

    if tox_risk_scores is None:
        assert smiles is not None  # validated above
        tox_risk_scores = _resolve_tox_risks(smiles)

    components = ScoreComponents(
        qed=druglikeness.qed,
        lipinski=lipinski_component(len(druglikeness.lipinski.violations)),
        synthesis=synthesis_component(sa),
        toxicity=toxicity_component(tox_risk_scores),
    )
    final = aggregate(components, weights, method)

    return CandidateScore(
        smiles=smiles,  # type: ignore[arg-type]
        components=components,
        final_score=final,
        method=method.value,
        weights=weights,
    )
