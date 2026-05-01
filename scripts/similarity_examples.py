"""Demo: Tanimoto similarity search vs the approved-drugs reference set.

Run from the repo root:

    uv run python scripts/similarity_examples.py

Searches the bundled 77-drug reference index for known canonical molecules
plus a few non-drug structures, prints top-5 hits per query. Useful as a
sanity check that the index returns chemically sensible neighbors --
statins retrieve other statins, NSAIDs cluster, etc.
"""

from __future__ import annotations

from chemscreen.similarity import find_similar, reference_size

QUERIES: dict[str, str] = {
    "aspirin (NSAID)":           "CC(=O)Oc1ccccc1C(=O)O",
    "atorvastatin (statin)":     "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)"
                                 "c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
    "metoprolol (beta blocker)": "COCCc1ccc(OCC(O)CNC(C)C)cc1",
    "fluoxetine (SSRI)":         "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",
    "nicotine (non-drug)":       "CN1CCCC1c1cccnc1",
    "rhodanine (PAINS hit)":     "O=C1CSC(=S)N1",
}


def _print_query(label: str, query_smiles: str) -> None:
    results = find_similar(query_smiles, top_k=5)
    print(f"\n--- {label} ---")
    if results is None:
        print("  (invalid SMILES)")
        return
    if not results:
        print("  (no hits)")
        return
    print(f"  {'sim':>5}  {'name':<18} {'class'}")
    for hit in results:
        print(
            f"  {hit.similarity:>5.3f}  {hit.name:<18} "
            f"({hit.drug_class or '-'})"
        )


def main() -> int:
    print(f"Reference set: {reference_size()} approved drugs")
    for label, smiles in QUERIES.items():
        _print_query(label, smiles)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
