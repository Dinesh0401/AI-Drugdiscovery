"""Ames mutagenicity prediction.

The Ames test detects compounds that cause mutations in bacterial DNA.
Mutagenic hits are early-stage red flags for genotoxicity in humans.

Trained on the Hansen 2009 Ames dataset (~6,500 compounds) via TDC.
"""

from __future__ import annotations

from pathlib import Path

from chemscreen.toxicity.base import RFToxicityModel, ToxicityPrediction

_MODEL_PATH = (
    Path(__file__).resolve().parent.parent.parent / "models" / "ames_rf.joblib"
)
_THRESHOLD_HIGH = 0.7
_THRESHOLD_MODERATE = 0.4

_MODEL: RFToxicityModel | None = None


def _get_model() -> RFToxicityModel:
    global _MODEL
    if _MODEL is None:
        if not _MODEL_PATH.exists():
            raise FileNotFoundError(
                f"Ames model not found at {_MODEL_PATH}. "
                "Run scripts/train_toxicity_models.py to produce it."
            )
        _MODEL = RFToxicityModel.load(
            _MODEL_PATH,
            endpoint="ames",
            threshold_high=_THRESHOLD_HIGH,
            threshold_moderate=_THRESHOLD_MODERATE,
        )
    return _MODEL


def predict_ames(smiles: str | None) -> ToxicityPrediction | None:
    """Predict Ames mutagenicity risk for a SMILES. None on invalid input."""
    return _get_model().predict(smiles)
