"""Score canonical examples to demonstrate the Phase 3 scoring engine.

Run from the repo root:

    uv run python scripts/score_examples.py

Prints two ranked tables — one per scoring method — over a hand-picked set
of molecules: well-known drugs, a Lipinski violator, and a couple of
PAINS-positive structures. Useful as a sanity check and as a reference
output for slides / reports.
"""

from __future__ import annotations

from chemscreen.ranking import screen_smiles
from chemscreen.scoring import ScoreWeights, ScoringMethod

EXAMPLES: dict[str, str] = {
    "caffeine":       "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "aspirin":        "CC(=O)OC1=CC=CC=C1C(=O)O",
    "ibuprofen":      "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "atorvastatin":   "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)"
                      "c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
    "rhodanine":      "O=C1CSC(=S)N1",
    "benzoquinone":   "O=C1C=CC(=O)C=C1",
    "curcumin":       "OC1=CC=C(C=C1)/C=C/C(=O)CC(=O)/C=C/C1=CC=C(O)C(=C1)OC",
}


def _print_table(method: ScoringMethod, weights: ScoreWeights) -> None:
    smiles_list = list(EXAMPLES.values())
    name_by_smiles = {smi: name for name, smi in EXAMPLES.items()}

    ranked = screen_smiles(smiles_list, weights=weights, method=method)

    print(f"\n=== {method.value.upper()} "
          f"(weights: qed={weights.qed} lip={weights.lipinski} "
          f"sa={weights.synthesis} tox={weights.toxicity}) ===")
    header = f"{'rank':<5}{'name':<14}{'qed':>7}{'lip':>7}{'sa':>7}{'tox':>7}{'final':>9}"
    print(header)
    print("-" * len(header))
    for r in ranked:
        c = r.components
        name = name_by_smiles.get(r.smiles, r.smiles[:14])
        print(
            f"{r.rank:<5}{name:<14}"
            f"{c.qed:>7.3f}{c.lipinski:>7.3f}"
            f"{c.synthesis:>7.3f}{c.toxicity:>7.3f}"
            f"{r.final_score:>9.3f}"
        )


def main() -> int:
    weights = ScoreWeights()  # 0.3 / 0.2 / 0.2 / 0.3 default
    _print_table(ScoringMethod.LINEAR, weights)
    _print_table(ScoringMethod.DESIRABILITY, weights)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
