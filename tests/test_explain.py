"""Tests for chemscreen.explain (Phase 5).

Verifies the explanation:
  - assigns the right severity given the inputs (Lipinski violator, hERG
    high-risk, multiple structural alerts);
  - carries the actual computed numbers in its bullets (every claim
    must be feature-grounded — that's the whole point of the
    deterministic-template approach);
  - serializes correctly via to_text() and to_dict();
  - returns None on invalid input.

Test molecules picked so each severity branch is exercised.
"""

from __future__ import annotations

import re

import pytest

from chemscreen.explain import Explanation, explain_molecule

RHODANINE = "O=C1CSC(=S)N1"  # 2 structural alerts (PAINS + Brenk)


def _all_text(exp: Explanation) -> str:
    return exp.headline + "\n" + "\n".join(exp.bullets) + "\n" + exp.recommendation


# --- Type / interface -----------------------------------------------------


def test_returns_explanation_dataclass(aspirin: str) -> None:
    exp = explain_molecule(aspirin)
    assert exp is not None
    assert isinstance(exp, Explanation)
    assert exp.smiles == aspirin
    assert isinstance(exp.bullets, list)
    assert all(isinstance(b, str) for b in exp.bullets)


def test_invalid_returns_none() -> None:
    assert explain_molecule("") is None
    assert explain_molecule("not_a_smiles") is None
    assert explain_molecule(None) is None


# --- Severity classification ---------------------------------------------


def test_aspirin_is_favorable(aspirin: str) -> None:
    """Aspirin: passes Lipinski, low tox, ~1 alert — should be favorable."""
    exp = explain_molecule(aspirin)
    assert exp is not None
    assert exp.severity == "favorable"


def test_ibuprofen_is_favorable(ibuprofen: str) -> None:
    exp = explain_molecule(ibuprofen)
    assert exp is not None
    assert exp.severity == "favorable"


def test_atorvastatin_is_unfavorable(atorvastatin: str) -> None:
    """Lipinski violator (MW > 500) — must be unfavorable."""
    exp = explain_molecule(atorvastatin)
    assert exp is not None
    assert exp.severity == "unfavorable"


def test_rhodanine_is_at_least_marginal() -> None:
    """Rhodanine has 2 structural alerts (PAINS + Brenk) — must NOT be favorable."""
    exp = explain_molecule(RHODANINE)
    assert exp is not None
    assert exp.severity in {"marginal", "unfavorable"}


# --- Headline content -----------------------------------------------------


def test_headline_contains_score(aspirin: str) -> None:
    exp = explain_molecule(aspirin)
    assert exp is not None
    assert exp.score is not None
    assert f"{exp.score:.2f}" in exp.headline


def test_atorvastatin_headline_flags_lipinski(atorvastatin: str) -> None:
    """The headline must spell out the failure reason (Lipinski violation)."""
    exp = explain_molecule(atorvastatin)
    assert exp is not None
    assert "Lipinski" in exp.headline


def test_rhodanine_headline_mentions_alerts() -> None:
    exp = explain_molecule(RHODANINE)
    assert exp is not None
    assert "alert" in exp.headline.lower()


# --- Bullets are feature-grounded ----------------------------------------


def test_bullets_carry_qed_number(aspirin: str) -> None:
    """Every QED claim must reference the actual numeric value."""
    exp = explain_molecule(aspirin)
    assert exp is not None
    qed_bullets = [b for b in exp.bullets if "QED" in b]
    assert qed_bullets, "QED must be reported"
    # Match a number like 0.55, 0.5, 1.0 etc.
    assert re.search(r"\d+\.\d+", qed_bullets[0]), \
        f"QED bullet must contain a number: {qed_bullets[0]!r}"


def test_bullets_carry_lipinski_number_on_violation(atorvastatin: str) -> None:
    """A Lipinski failure bullet must include the violating numeric value."""
    exp = explain_molecule(atorvastatin)
    assert exp is not None
    lipinski_bullets = [b for b in exp.bullets if "Lipinski" in b and "FAIL" in b]
    assert lipinski_bullets, "Lipinski FAIL must appear for atorvastatin"
    # Must contain "MW=" with a number > 500
    assert re.search(r"MW=\d+\.\d+\s*>\s*500", lipinski_bullets[0]), \
        f"Lipinski violation bullet must show MW number: {lipinski_bullets[0]!r}"


def test_bullets_carry_sa_number(aspirin: str) -> None:
    exp = explain_molecule(aspirin)
    assert exp is not None
    sa_bullets = [b for b in exp.bullets if "SA score" in b]
    assert sa_bullets
    assert re.search(r"\d+\.\d+", sa_bullets[0])


def test_atorvastatin_bullets_show_dili_high_risk(atorvastatin: str) -> None:
    """When DILI is high, a bullet must say so with the probability."""
    exp = explain_molecule(atorvastatin)
    assert exp is not None
    dili_bullets = [b for b in exp.bullets if "DILI" in b]
    assert dili_bullets
    assert "HIGH" in dili_bullets[0]
    assert re.search(r"probability\s+\d+\.\d+", dili_bullets[0])


def test_rhodanine_alerts_listed_by_name() -> None:
    exp = explain_molecule(RHODANINE)
    assert exp is not None
    alerts_bullets = [b for b in exp.bullets if "alerts" in b.lower()]
    assert alerts_bullets
    text = " ".join(alerts_bullets).lower()
    # Expect at least one of the rhodanine-related catalog names
    assert "rhod" in text or "thiocarbonyl" in text


# --- Similarity bullets --------------------------------------------------


def test_aspirin_similarity_bullet_lists_neighbors(aspirin: str) -> None:
    """An exact reference match should produce a 'closest neighbors' bullet."""
    exp = explain_molecule(aspirin)
    assert exp is not None
    sim_bullets = [b for b in exp.bullets if "neighbors" in b.lower() or "similar" in b.lower()]
    assert sim_bullets


def test_off_reference_similarity_bullet_uses_most_similar() -> None:
    """A non-reference query gets 'Most similar' phrasing, not 'Identical'."""
    exp = explain_molecule(RHODANINE)
    assert exp is not None
    sim_bullets = [b for b in exp.bullets if "similar" in b.lower()]
    assert sim_bullets
    # "Most similar" — not "Identical"
    assert any("Most similar" in b for b in sim_bullets)


# --- Recommendation -------------------------------------------------------


def test_lipinski_violator_recommendation_says_reject(atorvastatin: str) -> None:
    exp = explain_molecule(atorvastatin)
    assert exp is not None
    assert "Reject" in exp.recommendation


def test_favorable_recommendation_advances(aspirin: str) -> None:
    exp = explain_molecule(aspirin)
    assert exp is not None
    assert exp.severity == "favorable"
    assert "validation" in exp.recommendation.lower() or "advance" in exp.recommendation.lower()


def test_alert_heavy_recommendation_calls_out_alerts() -> None:
    exp = explain_molecule(RHODANINE)
    assert exp is not None
    assert "alert" in exp.recommendation.lower()


# --- Serialization --------------------------------------------------------


def test_to_text_is_multiline(aspirin: str) -> None:
    exp = explain_molecule(aspirin)
    assert exp is not None
    text = exp.to_text()
    assert "\n" in text
    assert exp.headline in text
    assert exp.recommendation in text


def test_to_dict_has_required_keys(aspirin: str) -> None:
    exp = explain_molecule(aspirin)
    assert exp is not None
    d = exp.to_dict()
    required = {"smiles", "severity", "headline", "bullets", "recommendation", "score"}
    assert required <= set(d.keys())
    assert isinstance(d["bullets"], list)


def test_to_dict_round_trip_text_present(aspirin: str) -> None:
    exp = explain_molecule(aspirin)
    assert exp is not None
    d = exp.to_dict()
    assert d["smiles"] == aspirin
    assert d["severity"] in {"favorable", "marginal", "unfavorable"}


# --- Robustness against missing toxicity models --------------------------


def test_runs_when_tox_models_unavailable(aspirin: str, monkeypatch) -> None:
    """If predict_all_toxicity raises FileNotFoundError, explanation still produces."""

    def _raise(_smiles: str) -> None:
        raise FileNotFoundError("simulated missing model")

    monkeypatch.setattr("chemscreen.toxicity.predict_all_toxicity", _raise)

    exp = explain_molecule(aspirin)
    assert exp is not None
    # The toxicity-not-computed bullet must appear
    assert any("not computed" in b for b in exp.bullets)


# --- Targeted severity / recommendation branches (monkeypatched tox) -----


def _fake_tox(classifications: dict[str, str]):
    """Build a predict_all_toxicity stub returning HIGH/MOD/LOW per endpoint."""
    from chemscreen.toxicity.base import ToxicityPrediction

    def _predict(_smiles: str):
        risk = {"low": 0.1, "moderate": 0.5, "high": 0.85}
        return {
            ep: ToxicityPrediction(
                endpoint=ep,
                risk_score=risk[cls],
                classification=cls,
                threshold_high=0.7,
                threshold_moderate=0.4,
                contributing_alerts=[],
            )
            for ep, cls in classifications.items()
        }

    return _predict


def test_two_high_tox_is_unfavorable(aspirin: str, monkeypatch) -> None:
    """Two HIGH tox endpoints must drive severity to unfavorable."""
    monkeypatch.setattr(
        "chemscreen.toxicity.predict_all_toxicity",
        _fake_tox({"ames": "high", "herg": "high", "dili": "low"}),
    )
    exp = explain_molecule(aspirin)
    assert exp is not None
    assert exp.severity == "unfavorable"


def test_herg_high_recommendation_calls_out_cardiac(aspirin: str, monkeypatch) -> None:
    monkeypatch.setattr(
        "chemscreen.toxicity.predict_all_toxicity",
        _fake_tox({"ames": "low", "herg": "high", "dili": "low"}),
    )
    exp = explain_molecule(aspirin)
    assert exp is not None
    assert "Cardiac" in exp.recommendation or "hERG" in exp.recommendation


def test_ames_high_recommendation_calls_out_mutagenicity(aspirin: str, monkeypatch) -> None:
    monkeypatch.setattr(
        "chemscreen.toxicity.predict_all_toxicity",
        _fake_tox({"ames": "high", "herg": "low", "dili": "low"}),
    )
    exp = explain_molecule(aspirin)
    assert exp is not None
    assert "Mutagenicity" in exp.recommendation
