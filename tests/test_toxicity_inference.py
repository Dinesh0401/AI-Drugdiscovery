"""Inference tests for the toxicity scaffold using a tiny synthetic RF.

These tests deliberately avoid depending on real trained models — they verify
the SMILES->fingerprint->predict_proba->classification pipeline in isolation.
The training script ``scripts/train_toxicity_models.py`` produces the actual
endpoint models against published datasets; ``test_toxicity_models.py``
covers those when the joblib files are present.
"""

from __future__ import annotations

import numpy as np
import pytest
from sklearn.ensemble import RandomForestClassifier

from chemscreen.toxicity.base import RFToxicityModel, ToxicityPrediction


@pytest.fixture
def tiny_rf_model() -> RFToxicityModel:
    """Train a small RF on synthetic Morgan-shaped data so we can call predict."""
    rng = np.random.default_rng(seed=7)
    X = rng.integers(0, 2, size=(50, 2048), dtype=np.uint8)
    y = rng.integers(0, 2, size=50)
    clf = RandomForestClassifier(n_estimators=20, random_state=7)
    clf.fit(X, y)
    return RFToxicityModel(endpoint="synthetic", estimator=clf)


def test_predict_returns_prediction(tiny_rf_model: RFToxicityModel, aspirin: str) -> None:
    result = tiny_rf_model.predict(aspirin)
    assert result is not None
    assert isinstance(result, ToxicityPrediction)
    assert result.endpoint == "synthetic"
    assert 0.0 <= result.risk_score <= 1.0
    assert result.classification in {"low", "moderate", "high"}


def test_predict_invalid_returns_none(
    tiny_rf_model: RFToxicityModel, bad_input: str | None
) -> None:
    assert tiny_rf_model.predict(bad_input) is None


def test_classification_thresholds() -> None:
    """Internal classification mapping is monotone and respects bounds."""
    model = RFToxicityModel(
        endpoint="t",
        estimator=_DummyEstimator(0.5),
        threshold_high=0.7,
        threshold_moderate=0.4,
    )
    assert model._classify(0.39) == "low"
    assert model._classify(0.40) == "moderate"
    assert model._classify(0.69) == "moderate"
    assert model._classify(0.70) == "high"
    assert model._classify(0.99) == "high"


def test_thresholds_propagate_to_prediction(aspirin: str) -> None:
    model = RFToxicityModel(
        endpoint="t",
        estimator=_DummyEstimator(0.95),
        threshold_high=0.5,
        threshold_moderate=0.2,
    )
    result = model.predict(aspirin)
    assert result is not None
    assert result.threshold_high == 0.5
    assert result.threshold_moderate == 0.2
    assert result.classification == "high"
    assert result.risk_score == pytest.approx(0.95)


class _DummyEstimator:
    """Sklearn-compatible stub that always returns a fixed positive probability."""

    def __init__(self, positive_proba: float) -> None:
        self._p = positive_proba

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        n = X.shape[0]
        return np.tile([1.0 - self._p, self._p], (n, 1))
