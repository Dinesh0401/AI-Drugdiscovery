"""Demo: Phase 6 operational modes (screening / risk / lead optimization).

Runs all three modes on canonical molecules and prints the structured
output. This is what the Phase 8 Flask UI will surface as separate pages.

Run from the repo root:

    uv run python scripts/modes_examples.py
"""

from __future__ import annotations

from chemscreen.modes import analyze_risk, screen_batch, suggest_variants

ATORVASTATIN = (
    "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)"
    "c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O"
)
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"
CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
CURCUMIN = "OC1=CC=C(C=C1)/C=C/C(=O)CC(=O)/C=C/C1=CC=C(O)C(=C1)OC"


def _section(title: str) -> None:
    bar = "=" * 78
    print(f"\n{bar}\n{title}\n{bar}")


def _demo_screening() -> None:
    _section("MODE 1: SCREENING — batch score + rank with rejection accounting")
    report = screen_batch(
        [ASPIRIN, IBUPROFEN, CAFFEINE, ATORVASTATIN, CURCUMIN, "not_a_smiles"],
        top_k=10,
    )
    print(f"\n  inputs: {report.n_input}, scored: {report.n_scored}, "
          f"rejected: {report.n_rejected}")
    print(f"  method: {report.method}\n")
    print(f"  {'rank':<5}{'score':<8}{'qed':<7}{'lip':<7}{'sa':<7}{'tox':<7}smiles")
    for c in report.candidates:
        co = c.components
        print(
            f"  {c.rank:<5}"
            f"{c.final_score:<8.3f}"
            f"{co.qed:<7.3f}{co.lipinski:<7.3f}"
            f"{co.synthesis:<7.3f}{co.toxicity:<7.3f}"
            f"{c.smiles[:40]}"
        )


def _demo_risk(label: str, smiles: str) -> None:
    _section(f"MODE 2: RISK ANALYSIS — {label}")
    report = analyze_risk(smiles)
    if report is None:
        print("  (invalid SMILES)")
        return
    print(f"\n  RISK LEVEL: {report.risk_level.upper()}")
    print(f"  Summary:    {report.summary}")
    if report.tox_predictions:
        print("  Toxicity predictions:")
        for endpoint, prediction in report.tox_predictions.items():
            print(
                f"    - {endpoint:<5}: {prediction.classification.upper():<8} "
                f"(p={prediction.risk_score:.2f})"
            )
    if report.structural_alerts:
        print(f"  Structural alerts: {', '.join(report.structural_alerts)}")
    if report.similar_drugs:
        print("  Closest approved drugs:")
        for hit in report.similar_drugs:
            print(
                f"    - {hit.similarity:.2f}  {hit.name:<18} "
                f"({hit.drug_class or '-'})"
            )
    print(f"  Recommendation: {report.recommendation}")


def _demo_lead_opt(label: str, smiles: str) -> None:
    _section(f"MODE 3: LEAD OPTIMIZATION — {label}")
    report = suggest_variants(smiles, max_variants=20)
    if report is None:
        print("  (invalid SMILES)")
        return
    print(f"\n  Parent: {report.parent_smiles}")
    if report.parent_score is not None:
        print(f"  Parent score: {report.parent_score.final_score:.3f}")
    print(f"  Variants generated: {len(report.variants)}")
    print(f"  Improved variants: {len(report.improved_variants)}")
    if report.improved_variants:
        print("\n  Top improvements:")
        for variant in report.improved_variants[:3]:
            assert variant.score is not None and report.parent_score is not None
            delta = variant.score.final_score - report.parent_score.final_score
            print(
                f"    +{delta:.3f}  {variant.transformation}\n"
                f"           -> {variant.smiles}"
            )
    elif report.variants:
        print("\n  All variants (none improved):")
        for variant in report.variants[:5]:
            if variant.score is None or report.parent_score is None:
                continue
            delta = variant.score.final_score - report.parent_score.final_score
            print(
                f"    {delta:+.3f}  {variant.transformation}\n"
                f"           -> {variant.smiles}"
            )


def main() -> int:
    _demo_screening()
    _demo_risk("aspirin", ASPIRIN)
    _demo_risk("atorvastatin", ATORVASTATIN)
    _demo_lead_opt("aspirin", ASPIRIN)
    _demo_lead_opt("atorvastatin", ATORVASTATIN)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
