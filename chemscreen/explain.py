"""Phase 5 — explainability layer.

Turns the raw outputs of Phases 1–4 (descriptors, drug-likeness rule
results, SA score, toxicity predictions, structural alerts, similarity
hits, multi-objective score) into a structured, human-readable
``Explanation`` object.

Every bullet in the explanation is feature-grounded — it references the
specific number that produced the verdict. There is no LLM in this
module: the templates are deterministic, reproducible, and cannot
hallucinate. A future polish layer can paraphrase these sentences with
an LLM if desired, but the *content* is computed here.

Output structure:
  - severity         "favorable" | "marginal" | "unfavorable"
  - headline         single-line verdict with score
  - bullets          detailed findings, one per finding
  - recommendation   what a chemist should do next

Methods to_text() and to_dict() format the same content for CLI and
API/JSON respectively.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from chemscreen.descriptors import compute_descriptors
from chemscreen.druglikeness import DrugLikeness, evaluate_druglikeness
from chemscreen.scoring import CandidateScore, ScoreWeights, score_candidate
from chemscreen.similarity import SimilarMolecule, find_similar
from chemscreen.synthesis import sa_score
from chemscreen.toxicity.alerts import get_structural_alerts
from chemscreen.toxicity.base import ToxicityPrediction


_ENDPOINT_FULL = {
    "ames": "Ames mutagenicity",
    "herg": "hERG cardiac inhibition",
    "dili": "DILI hepatotoxicity",
}

_SEVERITY_WORD = {
    "favorable": "Favorable",
    "marginal": "Marginal",
    "unfavorable": "Unfavorable",
}


@dataclass(frozen=True)
class Explanation:
    """Structured explanation of a screening verdict for one molecule."""

    smiles: str
    severity: str
    headline: str
    bullets: list[str]
    recommendation: str
    score: float | None = None

    def to_text(self) -> str:
        """Pretty multi-line CLI rendering."""
        lines: list[str] = [self.headline, ""]
        for bullet in self.bullets:
            lines.append(f"  - {bullet}")
        lines.append("")
        lines.append(f"Recommendation: {self.recommendation}")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """JSON-serializable view (suitable for API responses)."""
        return asdict(self)


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------


def _explain_druglikeness(descriptors, dl: DrugLikeness) -> list[str]:
    bullets: list[str] = []

    if dl.lipinski.passed:
        bullets.append("Lipinski Rule of Five: PASS (no violations)")
    else:
        violations = "; ".join(dl.lipinski.violations)
        bullets.append(f"Lipinski Rule of Five: FAIL ({violations})")

    qed_label = (
        "drug-like" if dl.qed >= 0.5
        else "marginal QED" if dl.qed >= 0.3
        else "low QED"
    )
    bullets.append(f"QED = {dl.qed:.2f} ({qed_label}; range 0..1)")

    bullets.append(
        f"MW = {descriptors.molecular_weight:.1f} Da, "
        f"LogP = {descriptors.logp:.2f}, "
        f"TPSA = {descriptors.tpsa:.1f} A^2"
    )

    if not dl.veber.passed:
        bullets.append(
            f"Veber rule: FAIL ({'; '.join(dl.veber.violations)})"
        )

    return bullets


def _explain_synthesis(sa: float) -> list[str]:
    if sa <= 3.0:
        label = "easy to synthesize"
    elif sa <= 6.0:
        label = "moderate synthetic complexity"
    else:
        label = "synthetically challenging"
    return [f"SA score = {sa:.2f} ({label}; 1=trivial, 10=very hard)"]


def _explain_toxicity(
    tox: dict[str, ToxicityPrediction],
) -> list[str]:
    if not tox:
        return [
            "Toxicity: not computed (run scripts/train_toxicity_models.py "
            "to enable Ames/hERG/DILI predictions)"
        ]
    bullets: list[str] = []
    for endpoint, prediction in tox.items():
        full = _ENDPOINT_FULL.get(endpoint, endpoint)
        bullets.append(
            f"{full}: {prediction.classification.upper()} "
            f"risk (probability {prediction.risk_score:.2f}; "
            f"RF on Morgan FP, TDC AUC 0.85+)"
        )
    return bullets


def _explain_alerts(smiles: str) -> list[str]:
    alerts = get_structural_alerts(smiles)
    if alerts is None or not alerts:
        return ["Structural alerts (Brenk + PAINS): none"]
    return [
        f"Structural alerts ({len(alerts)}): {', '.join(alerts)}"
    ]


def _explain_similarity(similar: list[SimilarMolecule]) -> list[str]:
    if not similar:
        return []
    if similar[0].similarity >= 0.999:
        # Exact reference match — show next 2 non-self hits
        non_self = similar[1:3]
        prefix = "Identical to a reference drug; closest neighbors"
    else:
        non_self = similar[:3]
        prefix = "Most similar approved drugs"
    if not non_self:
        return []
    formatted = "; ".join(
        f"{hit.name} ({hit.similarity:.2f}, {hit.drug_class or '-'})"
        for hit in non_self
    )
    return [f"{prefix}: {formatted}"]


# ---------------------------------------------------------------------------
# Verdict synthesis
# ---------------------------------------------------------------------------


def _high_tox_endpoints(tox: dict[str, ToxicityPrediction]) -> list[str]:
    return [k for k, p in tox.items() if p.classification == "high"]


def _classify_severity(
    score: CandidateScore | None,
    dl: DrugLikeness,
    tox: dict[str, ToxicityPrediction],
    num_alerts: int,
) -> str:
    if not dl.lipinski.passed:
        return "unfavorable"
    high = _high_tox_endpoints(tox)
    if len(high) >= 2:
        return "unfavorable"
    if len(high) == 1:
        return "marginal"
    # Multiple structural alerts (PAINS+Brenk) demote favorable to marginal
    if num_alerts >= 2:
        return "marginal"
    if score is not None and score.final_score >= 0.7:
        return "favorable"
    if score is not None and score.final_score >= 0.5:
        return "marginal"
    return "unfavorable"


def _make_headline(
    severity: str,
    score: CandidateScore | None,
    dl: DrugLikeness,
    tox: dict[str, ToxicityPrediction],
    num_alerts: int,
) -> str:
    score_str = f"{score.final_score:.2f}" if score else "n/a"
    word = _SEVERITY_WORD[severity]

    flags: list[str] = []
    if not dl.lipinski.passed:
        flags.append("Lipinski violation")
    high = _high_tox_endpoints(tox)
    if high:
        flags.append(
            ", ".join(_ENDPOINT_FULL.get(e, e) for e in high) + " risk"
        )
    if num_alerts >= 2 and not flags:
        flags.append(f"{num_alerts} structural alerts")

    if flags:
        return f"{word} candidate (score {score_str}) - flagged for {' + '.join(flags)}"
    if severity == "favorable":
        return f"{word} candidate (score {score_str}) - passes all primary filters"
    return f"{word} candidate (score {score_str})"


def _make_recommendation(
    severity: str,
    dl: DrugLikeness,
    tox: dict[str, ToxicityPrediction],
    num_alerts: int,
) -> str:
    if not dl.lipinski.passed:
        return (
            "Reject for oral-drug development. Lipinski rule failure "
            "(typically MW or LogP) predicts poor oral absorption. "
            "Consider scaffold simplification or alternative delivery routes."
        )
    high = _high_tox_endpoints(tox)
    if "herg" in high:
        return (
            "Cardiac safety concern dominates. Explore SAR variants that "
            "reduce hERG affinity (typically removing basic amines or "
            "lipophilic aromatic rings) before progressing."
        )
    if "ames" in high:
        return (
            "Mutagenicity flag - investigate the structural feature driving "
            "the prediction; aromatic amines, nitro groups, and Michael "
            "acceptors are common culprits."
        )
    if "dili" in high:
        return (
            "Hepatotoxicity flag. DILI is the leading cause of post-market "
            "drug withdrawal; warrants careful in vitro hepatocyte review "
            "before progressing."
        )
    if num_alerts >= 2:
        return (
            "Multiple structural alerts (Brenk/PAINS) detected. Investigate "
            "which substructures are flagged before progressing - common "
            "culprits are reactive Michael acceptors, quinones, rhodanines, "
            "and catechols."
        )
    if severity == "favorable":
        return (
            "Strong candidate. Worth advancing to in vitro validation "
            "(metabolic stability, permeability, target binding)."
        )
    if severity == "marginal":
        return (
            "Borderline candidate. Prioritize improving the weakest "
            "component (QED, SA, or a moderate-risk toxicity endpoint) "
            "before progressing."
        )
    return (
        "Insufficient signal to advance. Rebuild scaffold or pick a "
        "different starting point."
    )


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


def explain_molecule(
    smiles: str | None,
    weights: ScoreWeights | None = None,
) -> Explanation | None:
    """Run all chemscreen analyses and produce a structured explanation.

    Returns None if the SMILES is invalid. If the toxicity models aren't
    trained yet, the explanation still produces (with a note that
    toxicity wasn't computed) so the explainer is useful from a fresh
    checkout.
    """
    descriptors = compute_descriptors(smiles)
    if descriptors is None:
        return None

    druglikeness = evaluate_druglikeness(smiles)
    sa = sa_score(smiles)
    if druglikeness is None or sa is None:
        return None

    # Toxicity is best-effort: missing models shouldn't break explanation.
    try:
        from chemscreen.toxicity import predict_all_toxicity

        tox = predict_all_toxicity(smiles) or {}
    except (ImportError, FileNotFoundError):
        tox = {}

    similar = find_similar(smiles, top_k=5) or []

    tox_risks = [p.risk_score for p in tox.values()] if tox else None
    score = score_candidate(
        smiles, weights=weights, tox_risk_scores=tox_risks
    )

    assert smiles is not None  # validated by compute_descriptors
    structural_alerts = get_structural_alerts(smiles) or []

    bullets: list[str] = []
    bullets.extend(_explain_druglikeness(descriptors, druglikeness))
    bullets.extend(_explain_synthesis(sa))
    bullets.extend(_explain_toxicity(tox))
    bullets.extend(_explain_alerts(smiles))
    bullets.extend(_explain_similarity(similar))

    severity = _classify_severity(score, druglikeness, tox, len(structural_alerts))
    headline = _make_headline(
        severity, score, druglikeness, tox, len(structural_alerts)
    )
    recommendation = _make_recommendation(
        severity, druglikeness, tox, len(structural_alerts)
    )

    return Explanation(
        smiles=smiles,
        severity=severity,
        headline=headline,
        bullets=bullets,
        recommendation=recommendation,
        score=score.final_score if score else None,
    )
