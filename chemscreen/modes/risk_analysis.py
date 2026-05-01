"""Risk analysis mode: focused toxicity + structural-alerts report.

Same data Phase 5's explainer pulls, packaged for the "is this molecule
safe enough to pursue?" workflow. Returns a ``RiskReport`` with:

  - per-endpoint toxicity classifications and probabilities
  - structural alerts (Brenk + PAINS) by name
  - the most-similar approved drugs (helpful context — if the query
    resembles a known toxic drug, that's a red flag itself)
  - a single risk_level verdict and a chemist-actionable recommendation
"""

from __future__ import annotations

from dataclasses import dataclass, field

from chemscreen.descriptors import is_valid_smiles
from chemscreen.similarity import SimilarMolecule, find_similar
from chemscreen.toxicity.alerts import get_structural_alerts
from chemscreen.toxicity.base import ToxicityPrediction


@dataclass(frozen=True)
class RiskReport:
    """Structured risk assessment for one molecule."""

    smiles: str
    risk_level: str  # "low" | "moderate" | "high"
    tox_predictions: dict[str, ToxicityPrediction]
    structural_alerts: list[str]
    high_risk_endpoints: list[str]
    moderate_risk_endpoints: list[str]
    similar_drugs: list[SimilarMolecule] = field(default_factory=list)
    summary: str = ""
    recommendation: str = ""


_ENDPOINT_FULL = {
    "ames": "Ames mutagenicity",
    "herg": "hERG cardiac inhibition",
    "dili": "DILI hepatotoxicity",
}


def _classify_overall_risk(
    tox: dict[str, ToxicityPrediction], num_alerts: int
) -> str:
    """Roll up per-endpoint classifications + alerts into one risk_level."""
    high = [k for k, p in tox.items() if p.classification == "high"]
    if len(high) >= 2:
        return "high"
    if len(high) == 1:
        return "high"  # single HIGH still counts as high overall
    moderate = [k for k, p in tox.items() if p.classification == "moderate"]
    if num_alerts >= 2:
        return "moderate"
    if len(moderate) >= 2:
        return "moderate"
    return "low"


def _summarize(
    risk_level: str,
    high: list[str],
    moderate: list[str],
    num_alerts: int,
) -> str:
    if risk_level == "high":
        endpoints = ", ".join(_ENDPOINT_FULL.get(e, e) for e in high)
        return f"HIGH risk: predicted to be toxic for {endpoints}."
    if risk_level == "moderate":
        notes = []
        if moderate:
            endpoints = ", ".join(_ENDPOINT_FULL.get(e, e) for e in moderate)
            notes.append(f"moderate risk on {endpoints}")
        if num_alerts >= 2:
            notes.append(f"{num_alerts} structural alerts")
        return "MODERATE risk: " + "; ".join(notes) + "."
    return "LOW risk: predictions and structural filters are clean."


def _recommend(
    risk_level: str, high: list[str], num_alerts: int
) -> str:
    if "herg" in high:
        return (
            "Cardiac safety dominates. Reduce hERG affinity (commonly: "
            "remove basic amines, lower lipophilicity) before progressing."
        )
    if "ames" in high:
        return (
            "Mutagenicity flag — investigate the structural feature driving "
            "the prediction (aromatic amines, nitro groups, Michael acceptors)."
        )
    if "dili" in high:
        return (
            "Hepatotoxicity flag — schedule in vitro hepatocyte testing "
            "before further optimization."
        )
    if num_alerts >= 2:
        return (
            "Multiple Brenk/PAINS alerts. Replace the flagged substructures "
            "(rhodanines, quinones, catechols, reactive Michael acceptors) "
            "with cleaner bioisosteres."
        )
    if risk_level == "moderate":
        return (
            "Borderline risk. Acceptable for early-stage exploration; revisit "
            "after lead optimization."
        )
    return (
        "No tox flags. Standard in vitro safety panel still recommended "
        "before in vivo work."
    )


def analyze_risk(smiles: str | None) -> RiskReport | None:
    """Produce a focused risk report for one SMILES.

    Returns None if the SMILES is invalid. If toxicity models aren't
    trained, returns a report with an empty ``tox_predictions`` dict and
    a clear note in ``summary``.
    """
    if not is_valid_smiles(smiles):
        return None
    assert smiles is not None  # is_valid_smiles guards None

    try:
        from chemscreen.toxicity import predict_all_toxicity

        tox = predict_all_toxicity(smiles) or {}
    except (ImportError, FileNotFoundError):
        tox = {}

    structural_alerts = get_structural_alerts(smiles) or []
    high = [k for k, p in tox.items() if p.classification == "high"]
    moderate = [k for k, p in tox.items() if p.classification == "moderate"]

    if not tox:
        return RiskReport(
            smiles=smiles,
            risk_level="unknown",
            tox_predictions={},
            structural_alerts=structural_alerts,
            high_risk_endpoints=[],
            moderate_risk_endpoints=[],
            similar_drugs=find_similar(smiles, top_k=3) or [],
            summary=(
                "Toxicity models not trained — run "
                "scripts/train_toxicity_models.py to enable tox predictions."
            ),
            recommendation=(
                "Train toxicity models before relying on this report; "
                "structural alerts alone are not a complete risk picture."
            ),
        )

    risk_level = _classify_overall_risk(tox, len(structural_alerts))

    return RiskReport(
        smiles=smiles,
        risk_level=risk_level,
        tox_predictions=tox,
        structural_alerts=structural_alerts,
        high_risk_endpoints=high,
        moderate_risk_endpoints=moderate,
        similar_drugs=find_similar(smiles, top_k=3) or [],
        summary=_summarize(risk_level, high, moderate, len(structural_alerts)),
        recommendation=_recommend(risk_level, high, len(structural_alerts)),
    )
