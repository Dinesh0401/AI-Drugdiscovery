"""hERG cardiac risk prediction.

The hERG potassium channel (encoded by KCNH2) controls cardiac repolarization.
Inhibitors prolong the QT interval and can trigger torsades de pointes — a
fatal arrhythmia. hERG screening is mandatory in modern drug development.

Trained on the Karim et al. hERG inhibition dataset via TDC. We use a
slightly tighter "high" threshold than Ames since hERG hits are more
attrition-determining in late-stage development.
"""

from __future__ import annotations

from pathlib import Path

from chemscreen.toxicity.base import RFToxicityModel, ToxicityPrediction

_MODEL_PATH = (
    Path(__file__).resolve().parent.parent.parent / "models" / "herg_rf.joblib"
)
_THRESHOLD_HIGH = 0.6
_THRESHOLD_MODERATE = 0.3

_MODEL: RFToxicityModel | None = None


def _get_model() -> RFToxicityModel:
    global _MODEL
    if _MODEL is None:
        if not _MODEL_PATH.exists():
            raise FileNotFoundError(
                f"hERG model not found at {_MODEL_PATH}. "
                "Run scripts/train_toxicity_models.py to produce it."
            )
        _MODEL = RFToxicityModel.load(
            _MODEL_PATH,
            endpoint="herg",
            threshold_high=_THRESHOLD_HIGH,
            threshold_moderate=_THRESHOLD_MODERATE,
        )
    return _MODEL


def predict_herg(smiles: str | None) -> ToxicityPrediction | None:
    """Predict hERG inhibition risk for a SMILES. None on invalid input."""
    return _get_model().predict(smiles)
