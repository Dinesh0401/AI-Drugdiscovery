"""Phase 7 — validation harness.

Treats the chemscreen scoring pipeline as a binary classifier — "is this
molecule drug-like enough to advance?" — and benchmarks it on:

  - positive class: known FDA-approved drugs
  - negative class: random non-drug commercial compounds (decoys)

Reports ROC-AUC, accuracy, sensitivity, specificity, F1, and the optimal
classification threshold (Youden's J statistic), plus the score
distributions for each class.

By default, toxicity is **excluded** from the benchmark scoring (weights
0.4/0.3/0.3/0.0 for QED/Lipinski/SA/Tox). Reason: the toxicity models
were trained on TDC datasets that may share molecules with random
non-drug compounds, so including toxicity would risk indirect data
leakage into this benchmark. The toxicity models' own held-out AUCs are
already reported separately in ``models/metrics.json`` (0.85 / 0.86 /
0.85 on AMES / hERG / DILI).

This benchmark answers a different question:

    "Does the deterministic scoring (QED + Lipinski + SA) distinguish
     approved drugs from random commercial compounds?"

Re-runs are deterministic given a fixed input set, so the AUC is a
reproducible claim — not a one-off estimate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    roc_auc_score,
    roc_curve,
)

from chemscreen.descriptors import is_valid_smiles
from chemscreen.scoring import (
    ScoreWeights,
    ScoringMethod,
    score_candidate,
)

# Default benchmark weights: toxicity disabled to prevent indirect data leakage.
BENCHMARK_WEIGHTS = ScoreWeights(qed=0.4, lipinski=0.3, synthesis=0.3, toxicity=0.0)


@dataclass(frozen=True)
class BenchmarkResult:
    """Outcome of running the chemscreen scoring on a labeled dataset."""

    n_positive: int
    n_negative: int
    n_skipped: int
    auc: float
    accuracy: float
    sensitivity: float
    specificity: float
    f1: float
    optimal_threshold: float
    confusion_matrix: dict[str, int]  # tp, fp, tn, fn
    positive_score_mean: float
    positive_score_std: float
    negative_score_mean: float
    negative_score_std: float
    method: str
    weights: ScoreWeights

    def summary(self) -> str:
        """Human-readable one-shot summary."""
        cm = self.confusion_matrix
        return (
            f"Drugs vs decoys benchmark "
            f"({self.n_positive} drugs, {self.n_negative} decoys, "
            f"{self.n_skipped} skipped)\n"
            f"  ROC-AUC          {self.auc:.3f}\n"
            f"  Accuracy         {self.accuracy:.3f}\n"
            f"  Sensitivity      {self.sensitivity:.3f}\n"
            f"  Specificity      {self.specificity:.3f}\n"
            f"  F1               {self.f1:.3f}\n"
            f"  Optimal threshold {self.optimal_threshold:.3f}\n"
            f"  Confusion matrix  TP={cm['tp']} FN={cm['fn']} "
            f"FP={cm['fp']} TN={cm['tn']}\n"
            f"  Score (drug)     mean={self.positive_score_mean:.3f} "
            f"std={self.positive_score_std:.3f}\n"
            f"  Score (decoy)    mean={self.negative_score_mean:.3f} "
            f"std={self.negative_score_std:.3f}"
        )


def _score_set(
    smiles_list: Sequence[str],
    weights: ScoreWeights,
    method: ScoringMethod,
) -> tuple[list[float], int]:
    """Score every SMILES; return (scores, n_skipped) for invalid inputs."""
    scores: list[float] = []
    skipped = 0
    for smiles in smiles_list:
        if not is_valid_smiles(smiles):
            skipped += 1
            continue
        # tox_risk_scores=[] -> toxicity component contributes 0.5 uniformly,
        # but with weight=0 in BENCHMARK_WEIGHTS it doesn't affect score anyway.
        result = score_candidate(
            smiles, weights=weights, method=method, tox_risk_scores=[]
        )
        if result is None:
            skipped += 1
            continue
        scores.append(result.final_score)
    return scores, skipped


def run_benchmark(
    positive_smiles: Sequence[str],
    negative_smiles: Sequence[str],
    weights: ScoreWeights | None = None,
    method: ScoringMethod = ScoringMethod.LINEAR,
) -> BenchmarkResult:
    """Run the chemscreen scoring as a binary classifier on labeled inputs.

    Positives (label=1) score high if the pipeline thinks they are
    drug-like; negatives (label=0) should score low. ROC-AUC measures how
    well the score separates the two classes.

    The optimal threshold is chosen by maximizing Youden's J statistic
    (sensitivity + specificity - 1), which is the standard choice when
    sensitivity and specificity matter equally.
    """
    weights = weights or BENCHMARK_WEIGHTS

    pos_scores, pos_skipped = _score_set(positive_smiles, weights, method)
    neg_scores, neg_skipped = _score_set(negative_smiles, weights, method)
    if not pos_scores or not neg_scores:
        raise ValueError(
            "Both positive and negative classes must have at least one "
            "valid SMILES after filtering."
        )

    y_true = np.array([1] * len(pos_scores) + [0] * len(neg_scores))
    y_score = np.array(pos_scores + neg_scores, dtype=float)

    auc = float(roc_auc_score(y_true, y_score))

    # Youden's J -> optimal threshold from the ROC curve
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    youden = tpr - fpr
    optimal_idx = int(np.argmax(youden))
    optimal_threshold = float(thresholds[optimal_idx])
    # roc_curve can pad threshold[0] to inf; clamp to score range
    if not np.isfinite(optimal_threshold):
        optimal_threshold = float(np.median(y_score))

    y_pred = (y_score >= optimal_threshold).astype(int)
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) else 0.0
    specificity = tn / (tn + fp) if (tn + fp) else 0.0

    return BenchmarkResult(
        n_positive=len(pos_scores),
        n_negative=len(neg_scores),
        n_skipped=pos_skipped + neg_skipped,
        auc=auc,
        accuracy=float(accuracy_score(y_true, y_pred)),
        sensitivity=float(sensitivity),
        specificity=float(specificity),
        f1=float(f1_score(y_true, y_pred, zero_division=0)),
        optimal_threshold=optimal_threshold,
        confusion_matrix={"tp": int(tp), "fp": int(fp), "tn": int(tn), "fn": int(fn)},
        positive_score_mean=float(np.mean(pos_scores)),
        positive_score_std=float(np.std(pos_scores)),
        negative_score_mean=float(np.mean(neg_scores)),
        negative_score_std=float(np.std(neg_scores)),
        method=method.value,
        weights=weights,
    )
