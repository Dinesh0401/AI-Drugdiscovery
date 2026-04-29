"""Run all three toxicity endpoints + structural alerts in one call."""

from __future__ import annotations

from dataclasses import replace

from chemscreen.descriptors import is_valid_smiles
from chemscreen.toxicity.alerts import get_structural_alerts
from chemscreen.toxicity.ames import predict_ames
from chemscreen.toxicity.base import ToxicityPrediction
from chemscreen.toxicity.dili import predict_dili
from chemscreen.toxicity.herg import predict_herg

_ENDPOINTS = (
    ("ames", predict_ames),
    ("herg", predict_herg),
    ("dili", predict_dili),
)


def predict_all_toxicity(
    smiles: str | None,
) -> dict[str, ToxicityPrediction] | None:
    """Run Ames, hERG, and DILI predictions plus structural-alert annotation.

    Each prediction's ``contributing_alerts`` field carries the union of
    Brenk + PAINS hits, so a single risk-analysis call gives both the
    probabilistic and the rule-based view of what's wrong with the molecule.

    Returns None if the SMILES is invalid.
    """
    if not is_valid_smiles(smiles):
        return None

    alerts = get_structural_alerts(smiles) or []
    out: dict[str, ToxicityPrediction] = {}
    for endpoint, predict_fn in _ENDPOINTS:
        prediction = predict_fn(smiles)
        if prediction is None:
            continue
        out[endpoint] = replace(prediction, contributing_alerts=alerts)
    return out
