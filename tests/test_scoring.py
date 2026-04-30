"""Tests for chemscreen.scoring (Phase 3)."""

from __future__ import annotations

import math

import pytest

from chemscreen.scoring import (
    CandidateScore,
    ScoreComponents,
    ScoreWeights,
    ScoringMethod,
    aggregate,
    lipinski_component,
    score_candidate,
    synthesis_component,
    toxicity_component,
)


# --- ScoreWeights ---------------------------------------------------------


def test_default_weights_sum_to_one() -> None:
    w = ScoreWeights()
    assert w.qed + w.lipinski + w.synthesis + w.toxicity == pytest.approx(1.0)


def test_weights_must_sum_to_one() -> None:
    with pytest.raises(ValueError, match="must sum to 1.0"):
        ScoreWeights(qed=0.5, lipinski=0.5, synthesis=0.5, toxicity=0.5)


def test_negative_weights_rejected() -> None:
    with pytest.raises(ValueError, match="cannot be negative"):
        ScoreWeights(qed=-0.1, lipinski=0.4, synthesis=0.3, toxicity=0.4)


def test_custom_weights_accepted() -> None:
    w = ScoreWeights(qed=0.4, lipinski=0.1, synthesis=0.1, toxicity=0.4)
    assert w.qed == 0.4


# --- Component normalizers ------------------------------------------------


def test_lipinski_component_zero_violations() -> None:
    assert lipinski_component(0) == 1.0


def test_lipinski_component_max_violations() -> None:
    assert lipinski_component(4) == 0.0


def test_lipinski_component_partial() -> None:
    assert lipinski_component(2) == 0.5


def test_lipinski_component_clamps_above_4() -> None:
    assert lipinski_component(5) == 0.0


def test_synthesis_component_easy() -> None:
    assert synthesis_component(1.0) == pytest.approx(1.0)


def test_synthesis_component_hard() -> None:
    assert synthesis_component(10.0) == pytest.approx(0.0)


def test_synthesis_component_mid() -> None:
    # SA=5.5 -> (10-5.5)/9 = 0.5
    assert synthesis_component(5.5) == pytest.approx(0.5)


def test_synthesis_component_clamps() -> None:
    assert synthesis_component(11.0) == 0.0
    assert synthesis_component(0.0) == 1.0


def test_toxicity_component_no_data_returns_neutral() -> None:
    assert toxicity_component(None) == 0.5
    assert toxicity_component([]) == 0.5


def test_toxicity_component_zero_risk() -> None:
    assert toxicity_component([0.0, 0.0, 0.0]) == 1.0


def test_toxicity_component_max_risk() -> None:
    assert toxicity_component([1.0, 1.0, 1.0]) == 0.0


def test_toxicity_component_mean_50_percent() -> None:
    assert toxicity_component([0.5]) == pytest.approx(0.5)


# --- Aggregators ----------------------------------------------------------


def test_linear_aggregate_matches_formula() -> None:
    components = ScoreComponents(qed=0.8, lipinski=1.0, synthesis=0.7, toxicity=0.6)
    weights = ScoreWeights(qed=0.3, lipinski=0.2, synthesis=0.2, toxicity=0.3)
    expected = 0.3 * 0.8 + 0.2 * 1.0 + 0.2 * 0.7 + 0.3 * 0.6
    assert aggregate(components, weights, ScoringMethod.LINEAR) == pytest.approx(expected)


def test_desirability_geometric_mean() -> None:
    """Equal weights, equal components -> desirability equals the component."""
    components = ScoreComponents(qed=0.5, lipinski=0.5, synthesis=0.5, toxicity=0.5)
    weights = ScoreWeights(qed=0.25, lipinski=0.25, synthesis=0.25, toxicity=0.25)
    assert aggregate(components, weights, ScoringMethod.DESIRABILITY) == pytest.approx(0.5)


def test_desirability_penalizes_zero_component() -> None:
    """The defining property: any single zero pulls geometric mean to ~0."""
    components = ScoreComponents(qed=0.0, lipinski=1.0, synthesis=1.0, toxicity=1.0)
    weights = ScoreWeights(qed=0.25, lipinski=0.25, synthesis=0.25, toxicity=0.25)

    linear = aggregate(components, weights, ScoringMethod.LINEAR)
    desirability = aggregate(components, weights, ScoringMethod.DESIRABILITY)

    assert linear == pytest.approx(0.75)
    # geometric mean with one zero (eps=1e-9, weight=0.25) ~= (1e-9)^0.25 ~= 0.006
    assert desirability < 0.01
    assert desirability < linear / 50  # qualitative: dramatically penalized


def test_aggregate_unknown_method_raises() -> None:
    components = ScoreComponents(0.5, 0.5, 0.5, 0.5)
    weights = ScoreWeights()
    with pytest.raises(ValueError, match="Unknown ScoringMethod"):
        aggregate(components, weights, "not_a_method")  # type: ignore[arg-type]


# --- score_candidate end-to-end ------------------------------------------


def test_score_aspirin_is_drug_like(aspirin: str) -> None:
    """Aspirin should score reasonably well — passes Lipinski, simple, low tox."""
    result = score_candidate(aspirin, tox_risk_scores=[0.1, 0.1, 0.1])
    assert result is not None
    assert isinstance(result, CandidateScore)
    assert 0.0 < result.final_score <= 1.0
    assert result.method == "linear"
    assert result.smiles == aspirin


def test_score_invalid_returns_none() -> None:
    assert score_candidate("") is None
    assert score_candidate("not_a_smiles") is None
    assert score_candidate(None) is None


def test_score_atorvastatin_lower_than_aspirin(aspirin: str, atorvastatin: str) -> None:
    """Atorvastatin fails Lipinski; aspirin doesn't. Aspirin should score higher."""
    aspirin_score = score_candidate(aspirin, tox_risk_scores=[0.1, 0.1, 0.1])
    ator_score = score_candidate(atorvastatin, tox_risk_scores=[0.1, 0.1, 0.1])
    assert aspirin_score is not None and ator_score is not None
    assert aspirin_score.final_score > ator_score.final_score


def test_score_components_recorded(aspirin: str) -> None:
    result = score_candidate(aspirin, tox_risk_scores=[0.0, 0.0, 0.0])
    assert result is not None
    assert 0.0 <= result.components.qed <= 1.0
    assert 0.0 <= result.components.lipinski <= 1.0
    assert 0.0 <= result.components.synthesis <= 1.0
    assert result.components.toxicity == 1.0  # zero risks -> max tox component


def test_score_with_custom_weights(aspirin: str) -> None:
    weights = ScoreWeights(qed=0.7, lipinski=0.1, synthesis=0.1, toxicity=0.1)
    result = score_candidate(aspirin, weights=weights, tox_risk_scores=[0.1])
    assert result is not None
    assert result.weights == weights


def test_score_desirability_method(aspirin: str) -> None:
    result = score_candidate(
        aspirin, method=ScoringMethod.DESIRABILITY, tox_risk_scores=[0.1, 0.1, 0.1]
    )
    assert result is not None
    assert result.method == "desirability"
    assert 0.0 < result.final_score <= 1.0


def test_score_falls_back_when_tox_unavailable(aspirin: str, monkeypatch) -> None:
    """If predict_all_toxicity is unavailable (e.g. models missing), tox falls to neutral."""

    def _raise_filenotfound(_smiles: str) -> None:
        raise FileNotFoundError("Test: simulated missing model")

    monkeypatch.setattr(
        "chemscreen.toxicity.predict_all_toxicity", _raise_filenotfound
    )
    result = score_candidate(aspirin)  # tox_risk_scores=None -> tries lookup
    assert result is not None
    assert result.components.toxicity == 0.5
    assert math.isfinite(result.final_score)
