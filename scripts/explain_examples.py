"""Demo: explain_molecule on canonical molecules across severity classes.

Run from the repo root:

    uv run python scripts/explain_examples.py

Picks one molecule per severity class (favorable, marginal, unfavorable)
plus a Lipinski violator and a structural-alert hit, and prints the full
``Explanation.to_text()`` rendering. This is the single best
demonstration of what chemscreen produces — every other phase feeds into
this output.
"""

from __future__ import annotations

from chemscreen.explain import explain_molecule

EXAMPLES: list[tuple[str, str]] = [
    ("Ibuprofen (clean drug-like)",
     "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ("Aspirin (drug-like, 1 alert)",
     "CC(=O)Oc1ccccc1C(=O)O"),
    ("Atorvastatin (Lipinski violator + DILI risk)",
     "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)"
     "c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O"),
    ("Rhodanine (PAINS-positive)",
     "O=C1CSC(=S)N1"),
    ("Curcumin (drug-like, 2 alerts)",
     "OC1=CC=C(C=C1)/C=C/C(=O)CC(=O)/C=C/C1=CC=C(O)C(=C1)OC"),
]


def main() -> int:
    for label, smiles in EXAMPLES:
        print("=" * 78)
        print(f"QUERY: {label}")
        print("SMILES: " + smiles)
        print("=" * 78)
        explanation = explain_molecule(smiles)
        if explanation is None:
            print("  (invalid SMILES)\n")
            continue
        print(explanation.to_text())
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
