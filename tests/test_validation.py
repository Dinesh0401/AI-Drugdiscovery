"""Tests for chemscreen.validation (Phase 7 — drugs-vs-decoys benchmark)."""

from __future__ import annotations

import pytest

from chemscreen.data.approved_drugs import APPROVED_DRUGS
from chemscreen.data.decoys import DECOYS
from chemscreen.descriptors import canonical_smiles, is_valid_smiles
from chemscreen.scoring import ScoreWeights, ScoringMethod
from chemscreen.validation import (
    BENCHMARK_WEIGHTS,
    BenchmarkResult,
    run_benchmark,
)


# --- Decoy set integrity --------------------------------------------------


def test_every_decoy_smiles_is_valid() -> None:
    """Regression: a typo would silently shrink the benchmark set."""
    for name, smiles in DECOYS:
        assert is_valid_smiles(smiles), f"Invalid SMILES for {name}: {smiles!r}"


def test_decoys_do_not_overlap_approved_drugs() -> None:
    """If a decoy collides with the drug set, the benchmark is bogus."""
    drug_canon = {canonical_smiles(s) for _, s, _ in APPROVED_DRUGS}
    decoy_canon = {canonical_smiles(s) for _, s in DECOYS}
    overlap = decoy_canon & drug_canon
    assert overlap == set(), f"Decoy/drug overlap: {overlap}"


def test_decoy_set_minimum_size() -> None:
    """Need a reasonable sample for AUC to be meaningful."""
    assert len(DECOYS) >= 50


# --- BENCHMARK_WEIGHTS contract ------------------------------------------


def test_benchmark_weights_exclude_toxicity() -> None:
    """Toxicity weight must be 0 for the default benchmark to avoid leakage."""
    assert BENCHMARK_WEIGHTS.toxicity == 0.0


def test_benchmark_weights_sum_to_one() -> None:
    w = BENCHMARK_WEIGHTS
    assert w.qed + w.lipinski + w.synthesis + w.toxicity == pytest.approx(1.0)


# --- run_benchmark behavior ----------------------------------------------


def test_benchmark_returns_result_dataclass() -> None:
    drugs = [smi for _, smi, _ in APPROVED_DRUGS[:20]]
    decoys = [smi for _, smi in DECOYS[:20]]
    result = run_benchmark(drugs, decoys)
    assert isinstance(result, BenchmarkResult)
    assert result.n_positive > 0
    assert result.n_negative > 0


def test_benchmark_full_set_auc_above_random() -> None:
    """Full drugs-vs-decoys benchmark must beat random chance (AUC > 0.5)."""
    drugs = [smi for _, smi, _ in APPROVED_DRUGS]
    decoys = [smi for _, smi in DECOYS]
    result = run_benchmark(drugs, decoys)
    assert result.auc > 0.5, \
        f"Benchmark AUC {result.auc:.3f} is no better than random"


def test_benchmark_metrics_in_valid_ranges() -> None:
    drugs = [smi for _, smi, _ in APPROVED_DRUGS]
    decoys = [smi for _, smi in DECOYS]
    result = run_benchmark(drugs, decoys)
    assert 0.0 <= result.auc <= 1.0
    assert 0.0 <= result.accuracy <= 1.0
    assert 0.0 <= result.sensitivity <= 1.0
    assert 0.0 <= result.specificity <= 1.0
    assert 0.0 <= result.f1 <= 1.0


def test_benchmark_confusion_matrix_sums_to_total() -> None:
    drugs = [smi for _, smi, _ in APPROVED_DRUGS]
    decoys = [smi for _, smi in DECOYS]
    result = run_benchmark(drugs, decoys)
    cm = result.confusion_matrix
    assert cm["tp"] + cm["fn"] == result.n_positive
    assert cm["tn"] + cm["fp"] == result.n_negative


def test_benchmark_records_method_and_weights() -> None:
    drugs = [smi for _, smi, _ in APPROVED_DRUGS[:10]]
    decoys = [smi for _, smi in DECOYS[:10]]
    weights = ScoreWeights(qed=0.5, lipinski=0.5, synthesis=0.0, toxicity=0.0)
    result = run_benchmark(drugs, decoys, weights=weights, method=ScoringMethod.DESIRABILITY)
    assert result.method == "desirability"
    assert result.weights == weights


def test_benchmark_skips_invalid_smiles() -> None:
    drugs = [smi for _, smi, _ in APPROVED_DRUGS[:5]]
    decoys = [smi for _, smi in DECOYS[:5]]
    drugs_with_junk = drugs + ["not_a_smiles", ""]
    result = run_benchmark(drugs_with_junk, decoys)
    assert result.n_skipped == 2
    assert result.n_positive == 5


def test_benchmark_summary_is_multiline() -> None:
    drugs = [smi for _, smi, _ in APPROVED_DRUGS[:10]]
    decoys = [smi for _, smi in DECOYS[:10]]
    summary = run_benchmark(drugs, decoys).summary()
    assert "ROC-AUC" in summary
    assert "Confusion matrix" in summary
    assert "\n" in summary


def test_benchmark_empty_class_raises() -> None:
    """Both classes must be non-empty for AUC to be defined."""
    drugs = [smi for _, smi, _ in APPROVED_DRUGS[:3]]
    with pytest.raises(ValueError, match="at least one valid"):
        run_benchmark(drugs, [])
    with pytest.raises(ValueError, match="at least one valid"):
        run_benchmark([], [smi for _, smi in DECOYS[:3]])


def test_benchmark_drug_score_mean_higher_than_decoy() -> None:
    """Sanity: on average drugs should score higher than non-drugs."""
    drugs = [smi for _, smi, _ in APPROVED_DRUGS]
    decoys = [smi for _, smi in DECOYS]
    result = run_benchmark(drugs, decoys)
    assert result.positive_score_mean > result.negative_score_mean
