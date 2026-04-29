"""Train Random Forest toxicity classifiers from TDC datasets.

Run once after installing the training extras:

    uv pip install -e .[training]
    uv run python scripts/train_toxicity_models.py

Pipeline per endpoint:
  1. Pull dataset from Therapeutics Data Commons (TDC).
  2. Scaffold-split 70/10/20 (val ignored — held-out test only).
  3. Featurize SMILES as 2048-bit Morgan FP (radius=2).
  4. Train RandomForest (500 trees, balanced class weights).
  5. Report AUC / accuracy / sensitivity / specificity / F1 on the test set.
  6. Persist the estimator to ``models/<endpoint>_rf.joblib``.

Outputs are deterministic given a fixed random_state. Re-running reproduces
the same numbers, which is what makes "AUC = 0.83" a defensible claim.
"""

from __future__ import annotations

import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    roc_auc_score,
)

try:
    from tdc.single_pred import Tox
except ImportError as exc:
    raise SystemExit(
        "PyTDC is required to run this script. "
        "Install with: uv pip install -e .[training]"
    ) from exc

from chemscreen.descriptors import morgan_fingerprint

ROOT = Path(__file__).resolve().parent.parent
MODELS_DIR = ROOT / "models"
DATA_CACHE = ROOT / "data" / "tdc_cache"

ENDPOINTS: dict[str, str] = {
    "ames": "AMES",
    "herg": "hERG_Karim",
    "dili": "DILI",
}

RANDOM_STATE = 42
N_ESTIMATORS = 500


def featurize(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Convert (Drug, Y) dataframe into Morgan FP matrix + label vector."""
    rows: list[np.ndarray] = []
    labels: list[int] = []
    skipped = 0
    for smiles, label in zip(df["Drug"], df["Y"], strict=False):
        fp = morgan_fingerprint(smiles)
        if fp is None:
            skipped += 1
            continue
        rows.append(fp)
        labels.append(int(label))
    if skipped:
        print(f"  (skipped {skipped} unparseable SMILES)")
    return np.vstack(rows), np.asarray(labels, dtype=int)


def evaluate(y_true: np.ndarray, y_pred: np.ndarray, y_proba: np.ndarray) -> dict:
    """Compute the standard binary-classification metric set."""
    auc = roc_auc_score(y_true, y_proba)
    acc = accuracy_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred, zero_division=0)

    cm = confusion_matrix(y_true, y_pred)
    if cm.shape == (2, 2):
        tn, fp, fn, tp = cm.ravel()
    else:
        tn = fp = fn = tp = 0
    sensitivity = tp / (tp + fn) if (tp + fn) else 0.0
    specificity = tn / (tn + fp) if (tn + fp) else 0.0

    return {
        "auc": auc,
        "accuracy": acc,
        "f1": f1,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    }


def train_endpoint(endpoint: str, tdc_name: str) -> dict:
    print(f"\n=== {endpoint.upper()} ({tdc_name}) ===")

    DATA_CACHE.mkdir(parents=True, exist_ok=True)
    data = Tox(name=tdc_name, path=str(DATA_CACHE))
    split = data.get_split(method="scaffold", seed=RANDOM_STATE, frac=[0.7, 0.1, 0.2])
    train_df = split["train"]
    test_df = split["test"]
    print(f"  train: {len(train_df):>5}   test: {len(test_df):>5}")

    X_train, y_train = featurize(train_df)
    X_test, y_test = featurize(test_df)
    pos_rate = y_train.mean()
    print(f"  positive class rate (train): {pos_rate:.3f}")

    clf = RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        n_jobs=-1,
        class_weight="balanced",
        random_state=RANDOM_STATE,
    )
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    y_proba = clf.predict_proba(X_test)[:, 1]
    metrics = evaluate(y_test, y_pred, y_proba)

    print(f"  AUC          {metrics['auc']:.3f}")
    print(f"  accuracy     {metrics['accuracy']:.3f}")
    print(f"  sensitivity  {metrics['sensitivity']:.3f}")
    print(f"  specificity  {metrics['specificity']:.3f}")
    print(f"  F1           {metrics['f1']:.3f}")
    print(
        f"  confusion    TN={metrics['tn']} FP={metrics['fp']} "
        f"FN={metrics['fn']} TP={metrics['tp']}"
    )

    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    out = MODELS_DIR / f"{endpoint}_rf.joblib"
    joblib.dump(clf, out)
    print(f"  saved -> {out.relative_to(ROOT)}")

    return {"endpoint": endpoint, "tdc_name": tdc_name, **metrics}


def main() -> int:
    summary: list[dict] = []
    for endpoint, tdc_name in ENDPOINTS.items():
        try:
            summary.append(train_endpoint(endpoint, tdc_name))
        except Exception as exc:  # noqa: BLE001 — surface any TDC/sklearn error
            print(f"FAILED {endpoint}: {exc}", file=sys.stderr)
            return 1

    print("\n=== Summary ===")
    df = pd.DataFrame(summary)[
        ["endpoint", "auc", "accuracy", "sensitivity", "specificity", "f1"]
    ]
    print(df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    metrics_path = MODELS_DIR / "metrics.json"
    df.to_json(metrics_path, orient="records", indent=2)
    print(f"\nMetrics written to {metrics_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
