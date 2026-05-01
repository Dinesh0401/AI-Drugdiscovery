"""Tests for chemscreen.modes (Phase 6 — operational workflows)."""

from __future__ import annotations

import pytest

from chemscreen.modes import (
    LeadOptimizationReport,
    RiskReport,
    ScreeningReport,
    StructuralVariant,
    analyze_risk,
    screen_batch,
    suggest_variants,
)
from chemscreen.scoring import ScoreWeights, ScoringMethod

CURCUMIN = "OC1=CC=C(C=C1)/C=C/C(=O)CC(=O)/C=C/C1=CC=C(O)C(=C1)OC"


# --- screening mode ------------------------------------------------------


def test_screen_batch_drops_invalid_smiles(aspirin: str, ibuprofen: str) -> None:
    report = screen_batch([aspirin, "not_a_smiles", ibuprofen, ""])
    assert isinstance(report, ScreeningReport)
    assert report.n_input == 4
    assert report.n_scored == 2
    assert report.n_rejected == 2


def test_screen_batch_returns_ranked_candidates(aspirin: str, ibuprofen: str) -> None:
    report = screen_batch([aspirin, ibuprofen])
    ranks = [c.rank for c in report.candidates]
    assert ranks == sorted(ranks)
    scores = [c.final_score for c in report.candidates]
    assert scores == sorted(scores, reverse=True)


def test_screen_batch_top_k_truncates(aspirin: str, ibuprofen: str, caffeine: str) -> None:
    report = screen_batch([aspirin, ibuprofen, caffeine], top_k=2)
    assert len(report.candidates) == 2


def test_screen_batch_require_lipinski_filters_violators(
    aspirin: str, atorvastatin: str
) -> None:
    """Atorvastatin (Lipinski violator) must be dropped when require_lipinski=True."""
    report = screen_batch([aspirin, atorvastatin], require_lipinski=True)
    smiles_kept = {c.smiles for c in report.candidates}
    assert aspirin in smiles_kept
    assert atorvastatin not in smiles_kept


def test_screen_batch_records_method_and_weights(aspirin: str) -> None:
    weights = ScoreWeights(qed=0.4, lipinski=0.3, synthesis=0.2, toxicity=0.1)
    report = screen_batch(
        [aspirin], weights=weights, method=ScoringMethod.DESIRABILITY
    )
    assert report.method == "desirability"
    assert report.weights == weights


def test_screen_batch_empty_input() -> None:
    report = screen_batch([])
    assert report.n_input == 0
    assert report.n_scored == 0
    assert report.candidates == []


# --- risk analysis mode --------------------------------------------------


def test_analyze_risk_returns_report(aspirin: str) -> None:
    report = analyze_risk(aspirin)
    assert isinstance(report, RiskReport)
    assert report.smiles == aspirin
    assert report.risk_level in {"low", "moderate", "high", "unknown"}


def test_analyze_risk_invalid_returns_none() -> None:
    assert analyze_risk("not_a_smiles") is None
    assert analyze_risk("") is None
    assert analyze_risk(None) is None


def test_analyze_risk_atorvastatin_high(atorvastatin: str) -> None:
    """Atorvastatin trains as DILI-positive in our model — risk_level should be high."""
    report = analyze_risk(atorvastatin)
    assert report is not None
    assert report.risk_level == "high"
    assert "dili" in report.high_risk_endpoints
    assert "DILI" in report.summary or "hepatotoxicity" in report.summary.lower()


def test_analyze_risk_recommendation_matches_high_endpoint(atorvastatin: str) -> None:
    report = analyze_risk(atorvastatin)
    assert report is not None
    # DILI flag -> recommendation should reference hepatotoxicity
    assert "hepatocyte" in report.recommendation.lower() \
           or "hepatotoxicity" in report.recommendation.lower()


def test_analyze_risk_includes_alerts_for_curcumin() -> None:
    """Curcumin has Brenk alerts (Michael acceptor)."""
    report = analyze_risk(CURCUMIN)
    assert report is not None
    assert len(report.structural_alerts) >= 1


def test_analyze_risk_lists_similar_drugs(aspirin: str) -> None:
    report = analyze_risk(aspirin)
    assert report is not None
    assert len(report.similar_drugs) > 0


def test_analyze_risk_when_models_missing(aspirin: str, monkeypatch) -> None:
    """If toxicity models aren't trained, returns 'unknown' risk_level with note."""

    def _raise(_smi: str) -> None:
        raise FileNotFoundError("simulated missing model")

    monkeypatch.setattr("chemscreen.toxicity.predict_all_toxicity", _raise)
    report = analyze_risk(aspirin)
    assert report is not None
    assert report.risk_level == "unknown"
    assert report.tox_predictions == {}
    assert "not trained" in report.summary or "not computed" in report.summary


# --- lead optimization mode ----------------------------------------------


def test_suggest_variants_returns_report(aspirin: str) -> None:
    report = suggest_variants(aspirin)
    assert isinstance(report, LeadOptimizationReport)
    assert report.parent_score is not None
    assert isinstance(report.variants, list)


def test_suggest_variants_invalid_returns_none() -> None:
    assert suggest_variants("not_a_smiles") is None
    assert suggest_variants("") is None
    assert suggest_variants(None) is None


def test_variants_are_valid_smiles(aspirin: str) -> None:
    """Every reported variant must be parseable by RDKit."""
    from chemscreen.descriptors import is_valid_smiles

    report = suggest_variants(aspirin)
    assert report is not None
    for variant in report.variants:
        assert isinstance(variant, StructuralVariant)
        assert is_valid_smiles(variant.smiles), \
            f"Invalid variant SMILES: {variant.smiles}"


def test_variants_are_distinct_from_parent(aspirin: str) -> None:
    """Lead-opt should never report the parent itself as a 'variant'."""
    report = suggest_variants(aspirin)
    assert report is not None
    parent_canonical = report.parent_smiles
    variant_smiles = {v.smiles for v in report.variants}
    assert parent_canonical not in variant_smiles


def test_variants_have_no_atom_maps(aspirin: str) -> None:
    """Atom map artifacts ([X:1]) must not leak into variant SMILES."""
    report = suggest_variants(aspirin)
    assert report is not None
    for variant in report.variants:
        assert ":" not in variant.smiles or variant.smiles.count(":") < 1


def test_variants_are_not_fragmented(aspirin: str) -> None:
    """Bond-breaking transformations must be filtered (no '.' in SMILES)."""
    report = suggest_variants(aspirin)
    assert report is not None
    for variant in report.variants:
        assert "." not in variant.smiles


def test_variants_include_transformation_description(aspirin: str) -> None:
    report = suggest_variants(aspirin)
    assert report is not None
    assert len(report.variants) > 0
    for variant in report.variants:
        assert len(variant.transformation) > 0
        assert "->" in variant.transformation


def test_max_variants_limits_output(aspirin: str) -> None:
    report = suggest_variants(aspirin, max_variants=1)
    assert report is not None
    assert len(report.variants) <= 1


def test_improved_variants_scored_higher_than_parent(atorvastatin: str) -> None:
    """If improved_variants is non-empty, every entry must beat the parent score."""
    report = suggest_variants(atorvastatin)
    assert report is not None
    if report.parent_score is None:
        pytest.skip("parent could not be scored")
    for variant in report.improved_variants:
        assert variant.score is not None
        assert variant.score.final_score > report.parent_score.final_score


def test_improved_variants_sorted_descending(atorvastatin: str) -> None:
    report = suggest_variants(atorvastatin)
    assert report is not None
    scores = [
        v.score.final_score
        for v in report.improved_variants
        if v.score is not None
    ]
    assert scores == sorted(scores, reverse=True)
