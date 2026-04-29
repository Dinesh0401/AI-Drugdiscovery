"""DILI (Drug-Induced Liver Injury) hepatotoxicity prediction.

Hepatotoxicity is the single most common cause of post-market drug withdrawal.
The DILI dataset (Xu et al.) is small (~475 compounds) and class-imbalanced,
so models on it are noisier than Ames/hERG — we surface the prediction with
a wider "moderate" band to reflect that uncertainty.

Trained on the TDC DILI benchmark.
"""

from __future__ import annotations

from pathlib import Path

from chemscreen.toxicity.base import RFToxicityModel, ToxicityPrediction

_MODEL_PATH = (
    Path(__file__).resolve().parent.parent.parent / "models" / "dili_rf.joblib"
)
_THRESHOLD_HIGH = 0.7
_THRESHOLD_MODERATE = 0.35

_MODEL: RFToxicityModel | None = None


def _get_model() -> RFToxicityModel:
    global _MODEL
    if _MODEL is None:
        if not _MODEL_PATH.exists():
            raise FileNotFoundError(
                f"DILI model not found at {_MODEL_PATH}. "
                "Run scripts/train_toxicity_models.py to produce it."
            )
        _MODEL = RFToxicityModel.load(
            _MODEL_PATH,
            endpoint="dili",
            threshold_high=_THRESHOLD_HIGH,
            threshold_moderate=_THRESHOLD_MODERATE,
        )
    return _MODEL


def predict_dili(smiles: str | None) -> ToxicityPrediction | None:
    """Predict DILI (hepatotoxicity) risk for a SMILES. None on invalid input."""
    return _get_model().predict(smiles)
