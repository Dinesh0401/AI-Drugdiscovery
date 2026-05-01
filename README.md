# chemscreen

**AI-Assisted Drug Candidate Screening and Ranking System.**

A computational drug-candidate screening engine that evaluates, ranks, and
explains molecular viability using cheminformatics rules and machine-learning
toxicity predictors validated on public benchmark datasets.

This is the upgrade path from the legacy `ChemDesigner/` Flask demo into a
defensible scientific tool. Every claim is backed by code, computation, and
held-out test-set numbers.

---

## Status

| Phase | Scope | Status |
|---|---|---|
| **1** | Cheminformatics core (descriptors, drug-likeness, SA score) + tests | ✅ done |
| **2** | Toxicity ML (Ames, hERG, DILI) + structural alerts | ✅ done |
| **3** | Multi-objective scoring engine + ranking | ✅ done |
| **4** | Similarity search vs approved-drug reference set | ✅ done |
| **5** | Explainability layer (deterministic, feature-grounded) | ✅ done |
| **6** | Operational modes (screening / risk / lead-opt) | ✅ done |
| 7 | Validation harness | planned |
| 8 | Flask UI integration | planned |

---

## Toxicity model performance

Trained on Therapeutics Data Commons (Huang et al., NeurIPS 2021) with
scaffold-split 70/10/20, RandomForest (500 trees, balanced class weights) on
2048-bit Morgan fingerprints (radius=2). Metrics on held-out test set:

| Endpoint | Dataset | n_train | n_test | AUC | Accuracy | Sensitivity | Specificity | F1 |
|---|---|---|---|---|---|---|---|---|
| **AMES** mutagenicity | Hansen et al. | 5,091 | 1,456 | **0.850** | 0.784 | 0.832 | 0.716 | 0.817 |
| **hERG** cardiac risk | Karim et al. | 9,411 | 2,690 | **0.864** | 0.793 | 0.764 | 0.822 | 0.785 |
| **DILI** hepatotoxicity | Xu et al. | 332 | 96 | **0.851** | 0.740 | 0.820 | 0.652 | 0.766 |

These numbers are reproducible — re-running `scripts/train_toxicity_models.py`
produces the same results given a fixed random seed (42).

---

## Multi-objective scoring (Phase 3)

Combine QED, Lipinski compliance, synthetic accessibility, and toxicity into
a single ranked score. Two methods ship in the same module:

| Method | Formula | When to use |
|---|---|---|
| **Linear** *(default)* | `final = w_qed·QED + w_lip·LipFrac + w_sa·SAease + w_tox·ToxPenalty` | Quick screening; what most reviewers expect |
| **Desirability** | weighted geometric mean of the same components (Derringer & Suich, 1980) | Multi-objective optimization where any single weak component should disqualify a candidate |

Default weights `0.3 / 0.2 / 0.2 / 0.3` are configurable. Sample output from
`uv run python scripts/score_examples.py`:

```
=== LINEAR (weights: qed=0.3 lip=0.2 sa=0.2 tox=0.3) ===
rank name              qed    lip     sa    tox    final
1    ibuprofen       0.822  1.000  0.868  0.818    0.865
2    aspirin         0.550  1.000  0.936  0.902    0.823
3    caffeine        0.538  1.000  0.856  0.722    0.749
...
7    atorvastatin    0.163  0.500  0.744  0.504    0.449   <- Lipinski violator
```

---

## Similarity search (Phase 4)

Tanimoto fingerprint search against a curated 77-drug reference set
covering 25+ therapeutic classes (NSAIDs, statins, SSRIs, beta blockers,
antibiotics, etc.). Vectorized over the whole reference matrix in numpy
for ~50x speedup vs per-row Python loops.

```python
from chemscreen import find_similar

results = find_similar("CC(=O)Oc1ccccc1C(=O)O", top_k=5)  # aspirin
for hit in results:
    print(f"{hit.similarity:.3f}  {hit.name}  ({hit.drug_class})")
```

Sample output validates same-class clustering:

```
--- metoprolol query ---
1.000  metoprolol      (Beta blocker)
0.622  atenolol        (Beta blocker)
0.617  bisoprolol      (Beta blocker)
0.449  propranolol     (Beta blocker)
0.273  carvedilol      (Beta blocker)
```

Swap the bundled set for a full DrugBank or ChEMBL approved-drug export
via `set_reference_drugs(...)`.

---

## Explainability (Phase 5)

`explain_molecule(smiles)` runs every analysis (descriptors, drug-likeness
rules, SA score, toxicity ML, structural alerts, similarity search,
multi-objective score) and produces a structured `Explanation` with:

  - **severity** — `favorable` / `marginal` / `unfavorable`
  - **headline** — single-line verdict including the score
  - **bullets** — feature-grounded findings (every claim references the
    actual computed number)
  - **recommendation** — concrete next step for a chemist

There is no LLM in this layer. The templates are deterministic, so
explanations are reproducible and cannot hallucinate. A future polish
step can paraphrase these sentences with a Hugging Face model if
desired, but the *content* is computed.

Example output (`scripts/explain_examples.py` for atorvastatin):

```
Unfavorable candidate (score 0.45) - flagged for Lipinski violation + DILI hepatotoxicity risk

  - Lipinski Rule of Five: FAIL (MW=558.65 > 500; LogP=6.31 > 5)
  - QED = 0.16 (low QED; range 0..1)
  - MW = 558.7 Da, LogP = 6.31, TPSA = 111.8 A^2
  - Veber rule: FAIL (RotatableBonds=12 > 10)
  - SA score = 3.31 (moderate synthetic complexity; 1=trivial, 10=very hard)
  - Ames mutagenicity: LOW risk (probability 0.39; RF on Morgan FP, TDC AUC 0.85+)
  - hERG cardiac inhibition: MODERATE risk (probability 0.34; ...)
  - DILI hepatotoxicity: HIGH risk (probability 0.75; ...)
  - Structural alerts (Brenk + PAINS): none
  - Identical to a reference drug; closest neighbors: rosuvastatin (0.31, Statin); pravastatin (0.24, Statin)

Recommendation: Reject for oral-drug development. Lipinski rule failure
(typically MW or LogP) predicts poor oral absorption. Consider scaffold
simplification or alternative delivery routes.
```

---

## Operational modes (Phase 6)

Three named workflows that map to real medicinal-chemistry tasks. Each
returns a structured report (dataclass) — easy to consume from CLI,
Python, or a future Flask UI.

```python
from chemscreen.modes import screen_batch, analyze_risk, suggest_variants

# 1. Screening — score and rank a list, with rejection accounting
report = screen_batch(["CC(=O)Oc1ccccc1C(=O)O", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"],
                      top_k=10, require_lipinski=True)

# 2. Risk analysis — focused tox + alerts + similar-drug context
risk = analyze_risk(atorvastatin)
# -> RiskReport(risk_level="high", high_risk_endpoints=["dili"], ...)

# 3. Lead optimization — bioisostere swaps + re-scored variants
opt = suggest_variants(atorvastatin)
# -> LeadOptimizationReport(parent_score=0.45, variants=[...], improved_variants=[...])
```

The lead-optimization rules are real medchem bioisostere transformations
(Meanwell 2011): methyl→trifluoromethyl (metabolic stability),
hydroxyl→fluoro (H-bond donor swap), chloro↔fluoro (halogen size),
primary amine→hydroxyl (basicity reduction), methoxy→hydroxyl
(demethylation), and ethyl→methyl (chain trimming). Not generative —
deterministic SMARTS replacements, atom maps stripped, fragmented
products filtered.

---

## Module layout

```
chemscreen/
├── descriptors.py       Morgan fingerprints + RDKit descriptors
├── druglikeness.py      Lipinski, Veber, Ghose, QED
├── synthesis.py         Synthetic Accessibility (Ertl & Schuffenhauer)
├── reference_data.py    Common formula -> SMILES lookup
├── scoring.py           Multi-objective scoring (linear + desirability)
├── ranking.py           rank_candidates / screen_smiles / filter_top_k
├── similarity.py        Tanimoto search vs approved-drug reference
├── explain.py           Deterministic feature-grounded explanations
├── data/
│   └── approved_drugs.py   Curated 77-drug reference (25+ classes)
├── toxicity/
│   ├── base.py          ToxicityPrediction + RFToxicityModel
│   ├── ames.py          Ames mutagenicity endpoint
│   ├── herg.py          hERG cardiac-risk endpoint
│   ├── dili.py          DILI hepatotoxicity endpoint
│   ├── alerts.py        Brenk + PAINS structural alerts
│   └── ensemble.py      predict_all_toxicity convenience
└── modes/
    ├── screening.py        Batch screen + rejection accounting
    ├── risk_analysis.py    Focused toxicity report
    └── lead_optimization.py  Bioisostere-swap variant generation

scripts/
├── train_toxicity_models.py    Train and persist all 3 RFs from TDC
├── score_examples.py            Demo scoring on canonical molecules
├── similarity_examples.py       Demo similarity search
├── explain_examples.py          Demo full explanations across severity classes
└── modes_examples.py            Demo all three operational modes

models/                  Persisted joblib estimators (gitignored)
tests/                   205 tests, 94% coverage
ChemDesigner/            Legacy Flask demo (kept for reference)
```

Public functions return `None` on invalid input rather than raising.
Drug-likeness rule failures and toxicity alerts carry feature-grounded
strings (`"MW=558.65 > 500"`, `"phenol_ester"`) so explainability in
later phases reads off the same data structure.

---

## Setup and verification

```bash
# Install runtime + dev deps
uv sync --extra dev

# Run the test suite
uv run pytest --cov=chemscreen --cov-report=term-missing
# 97 passed, coverage 97%

# (One-time) train the toxicity models
uv sync --extra dev --extra training
uv run python scripts/train_toxicity_models.py
# produces models/{ames,herg,dili}_rf.joblib + models/metrics.json
```

## Smoke check

```python
from chemscreen import compute_descriptors, evaluate_druglikeness, sa_score
from chemscreen.toxicity import predict_all_toxicity

aspirin = "CC(=O)OC1=CC=CC=C1C(=O)O"

print(compute_descriptors(aspirin).molecular_weight)              # ~180.16
print(evaluate_druglikeness(aspirin).lipinski.passed)             # True
print(sa_score(aspirin))                                          # ~1.5-2.5

risk = predict_all_toxicity(aspirin)
for endpoint, prediction in risk.items():
    print(f"{endpoint}: risk={prediction.risk_score:.2f}, "
          f"class={prediction.classification}, "
          f"alerts={prediction.contributing_alerts}")
```

---

## Citations

- Lipinski et al. *Adv. Drug Deliv. Rev.* (1997) — Rule of Five
- Bickerton et al. *Nat. Chem.* (2012) — QED
- Ertl & Schuffenhauer *J. Cheminform.* (2009) — SA score
- Hansen et al. *J. Chem. Inf. Model.* (2009) — Ames benchmark
- Baell & Holloway *J. Med. Chem.* (2010) — PAINS
- Brenk et al. *ChemMedChem* (2008) — Brenk filter
- Huang et al. *NeurIPS* (2021) — Therapeutics Data Commons
