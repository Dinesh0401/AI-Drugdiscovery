"""Integration tests against the actual trained endpoint models.

Skipped automatically when ``models/`` is empty so unit-test runs don't
require running the training script first. Once
``scripts/train_toxicity_models.py`` has been run, these tests exercise the
end-to-end pipeline (lazy load -> featurize -> predict -> classify -> alerts).
"""

from __future__ import annotations

from pathlib import Path

import pytest

MODELS_DIR = Path(__file__).resolve().parent.parent / "models"

_REQUIRED = ["ames_rf.joblib", "herg_rf.joblib", "dili_rf.joblib"]
_MISSING = [name for name in _REQUIRED if not (MODELS_DIR / name).exists()]

pytestmark = pytest.mark.skipif(
    bool(_MISSING),
    reason=f"Trained models not present (missing: {_MISSING}); "
    "run scripts/train_toxicity_models.py to produce them.",
)


def test_predict_ames_returns_prediction(aspirin: str) -> None:
    from chemscreen.toxicity import predict_ames

    result = predict_ames(aspirin)
    assert result is not None
    assert result.endpoint == "ames"
    assert 0.0 <= result.risk_score <= 1.0


def test_predict_herg_returns_prediction(aspirin: str) -> None:
    from chemscreen.toxicity import predict_herg

    result = predict_herg(aspirin)
    assert result is not None
    assert result.endpoint == "herg"


def test_predict_dili_returns_prediction(aspirin: str) -> None:
    from chemscreen.toxicity import predict_dili

    result = predict_dili(aspirin)
    assert result is not None
    assert result.endpoint == "dili"


def test_predict_all_runs_and_attaches_alerts(aspirin: str) -> None:
    from chemscreen.toxicity import predict_all_toxicity

    result = predict_all_toxicity(aspirin)
    assert result is not None
    assert set(result.keys()) == {"ames", "herg", "dili"}
    for prediction in result.values():
        assert 0.0 <= prediction.risk_score <= 1.0
        assert prediction.classification in {"low", "moderate", "high"}
        # contributing_alerts must be a list (possibly empty for clean molecules)
        assert isinstance(prediction.contributing_alerts, list)


def test_predict_all_invalid_returns_none() -> None:
    from chemscreen.toxicity import predict_all_toxicity

    assert predict_all_toxicity("not_a_smiles") is None
