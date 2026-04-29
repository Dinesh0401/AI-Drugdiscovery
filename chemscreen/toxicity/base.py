"""Toxicity prediction interfaces and shared inference scaffold.

The ``RFToxicityModel`` wraps any scikit-learn classifier with ``predict_proba``
and exposes a SMILES-in / ``ToxicityPrediction``-out API. Endpoint-specific
modules (``ames``, ``herg``, ``dili``) instantiate this with their own
serialized estimator and risk-classification thresholds.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol

import joblib
from sklearn.base import BaseEstimator

from chemscreen.descriptors import morgan_fingerprint


@dataclass(frozen=True)
class ToxicityPrediction:
    """Outcome of evaluating one toxicity endpoint for one molecule.

    ``risk_score`` is the model's predicted probability of the toxic class
    (0..1). ``classification`` bins it into low / moderate / high using the
    endpoint-specific thresholds the model was constructed with.
    """

    endpoint: str
    risk_score: float
    classification: str  # "low" | "moderate" | "high"
    threshold_high: float
    threshold_moderate: float
    contributing_alerts: list[str] = field(default_factory=list)


class ToxicityModel(Protocol):
    """Structural type for any object that can predict one toxicity endpoint."""

    endpoint: str

    def predict(self, smiles: str | None) -> ToxicityPrediction | None: ...


class RFToxicityModel:
    """Random Forest binary classifier on 2048-bit Morgan fingerprints."""

    def __init__(
        self,
        endpoint: str,
        estimator: BaseEstimator,
        threshold_high: float = 0.7,
        threshold_moderate: float = 0.4,
        fp_radius: int = 2,
        fp_bits: int = 2048,
    ) -> None:
        self.endpoint = endpoint
        self._estimator = estimator
        self.threshold_high = threshold_high
        self.threshold_moderate = threshold_moderate
        self._fp_radius = fp_radius
        self._fp_bits = fp_bits

    @classmethod
    def load(
        cls,
        path: Path | str,
        endpoint: str,
        threshold_high: float = 0.7,
        threshold_moderate: float = 0.4,
    ) -> RFToxicityModel:
        """Load a serialized estimator from disk."""
        estimator = joblib.load(path)
        return cls(
            endpoint=endpoint,
            estimator=estimator,
            threshold_high=threshold_high,
            threshold_moderate=threshold_moderate,
        )

    def predict(self, smiles: str | None) -> ToxicityPrediction | None:
        """Predict toxicity for one SMILES. Returns None on invalid input."""
        fp = morgan_fingerprint(smiles, radius=self._fp_radius, n_bits=self._fp_bits)
        if fp is None:
            return None

        proba = float(
            self._estimator.predict_proba(fp.reshape(1, -1))[0, 1]
        )
        return ToxicityPrediction(
            endpoint=self.endpoint,
            risk_score=proba,
            classification=self._classify(proba),
            threshold_high=self.threshold_high,
            threshold_moderate=self.threshold_moderate,
        )

    def _classify(self, proba: float) -> str:
        if proba >= self.threshold_high:
            return "high"
        if proba >= self.threshold_moderate:
            return "moderate"
        return "low"
