"""Run the drugs-vs-decoys validation benchmark.

Loads the curated 77 approved drugs (positives) and 88 decoys (negatives)
that ship with chemscreen, scores each through the deterministic pipeline
(QED + Lipinski + SA, toxicity weight=0 to avoid indirect data leakage),
and reports ROC-AUC + classifier metrics.

Persists results to ``data/benchmark_results.json`` for citation.

Run from the repo root:

    uv run python scripts/run_benchmark.py
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from chemscreen.data.approved_drugs import APPROVED_DRUGS
from chemscreen.data.decoys import DECOYS
from chemscreen.validation import run_benchmark

ROOT = Path(__file__).resolve().parent.parent
RESULTS_PATH = ROOT / "data" / "benchmark_results.json"


def main() -> int:
    drugs = [smi for _, smi, _ in APPROVED_DRUGS]
    decoys = [smi for _, smi in DECOYS]

    print(f"Running drugs-vs-decoys benchmark:")
    print(f"  positives (drugs):   {len(drugs)}")
    print(f"  negatives (decoys):  {len(decoys)}")

    result = run_benchmark(drugs, decoys)

    print()
    print(result.summary())

    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(result)
    payload["weights"] = asdict(result.weights)
    RESULTS_PATH.write_text(json.dumps(payload, indent=2))
    print(f"\nResults written to {RESULTS_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
