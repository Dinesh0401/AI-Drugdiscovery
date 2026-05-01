"""Lead optimization mode: rule-based bioisostere swaps + re-scoring.

Bioisosteric replacement is a foundational medicinal-chemistry strategy
(Meanwell, *J. Med. Chem.* 2011) — replacing a substructure with a
chemically distinct but biologically similar fragment to tune properties
like metabolic stability, lipophilicity, or hERG affinity, without
disrupting the pharmacophore.

This module exposes a curated set of well-known bioisostere transformations
as RDKit SMARTS-replacement rules. Given a parent molecule, it applies each
rule, deduplicates the resulting variants, scores them via Phase 3, and
returns the parent + variants + the subset that improved over the parent.

Not a generative model. Not a latent walk. Each variant is a deterministic
single-atom or single-fragment swap with a citable rationale.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import AllChem

from chemscreen.descriptors import canonical_smiles, parse_smiles
from chemscreen.scoring import (
    CandidateScore,
    ScoreWeights,
    ScoringMethod,
    score_candidate,
)


@dataclass(frozen=True)
class StructuralVariant:
    """One transformation product."""

    smiles: str
    transformation: str  # human-readable description
    score: CandidateScore | None


@dataclass(frozen=True)
class LeadOptimizationReport:
    """Outcome of running bioisostere swaps on a parent molecule."""

    parent_smiles: str
    parent_score: CandidateScore | None
    variants: list[StructuralVariant]
    improved_variants: list[StructuralVariant] = field(default_factory=list)


# --- Bioisostere catalog --------------------------------------------------
# Each entry: (description, source_smarts, replacement_smiles, rationale).
# Sources: Meanwell 2011 (bioisosterism review), Ritchie et al. 2011
# (aromatic-ring drug-likeness), and standard medchem texts.

_TRANSFORMATIONS: list[tuple[str, str, str, str]] = [
    (
        "methyl -> trifluoromethyl",
        "[CH3]",
        "C(F)(F)F",
        "improves metabolic stability and lipophilicity",
    ),
    (
        "hydroxyl -> fluoro",
        "[OX2H1]",
        "F",
        "classic H-bond donor -> halogen swap; reduces metabolism at OH",
    ),
    (
        "chloro -> fluoro",
        "[Cl]",
        "F",
        "smaller halogen reduces metabolic liability",
    ),
    (
        "primary amine -> hydroxyl",
        "[NX3H2]",
        "O",
        "reduces basicity; commonly lowers hERG and CNS exposure",
    ),
    (
        "methoxy -> hydroxyl",
        "[OX2][CH3]",
        "O",
        "demethylation; changes H-bonding pattern",
    ),
    (
        "fluoro -> chloro",
        "[F]",
        "Cl",
        "increases lipophilicity (reverse halogen swap)",
    ),
    (
        "ethyl -> methyl",
        "[CH2][CH3]",
        "C",
        "shrinks aliphatic chain (reduces logP)",
    ),
]


def _strip_atom_maps(mol: Chem.Mol) -> None:
    """RDKit's ReplaceSubstructs can carry atom maps into products; strip them."""
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)


def _apply_transformation(
    parent_mol: Chem.Mol,
    source_smarts: str,
    replacement_smiles: str,
) -> list[str]:
    """Apply one SMARTS replacement rule, return list of unique product SMILES.

    Filters out:
      - products that fail to canonicalize
      - fragmented products (containing '.', meaning the swap broke a bond)
      - products identical to the parent (the rule didn't really change anything)
    """
    source = Chem.MolFromSmarts(source_smarts)
    replacement = Chem.MolFromSmiles(replacement_smiles)
    if source is None or replacement is None:
        return []
    if not parent_mol.HasSubstructMatch(source):
        return []

    parent_smi = Chem.MolToSmiles(parent_mol)

    products = AllChem.ReplaceSubstructs(parent_mol, source, replacement)
    unique: set[str] = set()
    for product in products:
        try:
            _strip_atom_maps(product)
            smi = Chem.MolToSmiles(product)
        except Exception:  # noqa: BLE001 — RDKit can throw for weird products
            continue
        if not smi or "." in smi:
            continue
        if smi == parent_smi:
            continue
        if Chem.MolFromSmiles(smi) is None:
            continue
        unique.add(smi)
    return sorted(unique)


def suggest_variants(
    smiles: str | None,
    weights: ScoreWeights | None = None,
    method: ScoringMethod = ScoringMethod.LINEAR,
    max_variants: int = 20,
) -> LeadOptimizationReport | None:
    """Apply bioisostere swaps and return scored variants.

    Returns None if the parent SMILES is invalid. Otherwise returns a
    ``LeadOptimizationReport`` with the parent score (if computable),
    every successful variant scored, and a separate list of variants
    whose score improved over the parent.
    """
    parent_mol = parse_smiles(smiles)
    if parent_mol is None:
        return None
    assert smiles is not None
    parent_canonical = canonical_smiles(smiles) or smiles

    parent_score = score_candidate(smiles, weights=weights, method=method)

    seen: set[str] = {parent_canonical}
    variants: list[StructuralVariant] = []

    for description, src_smarts, repl_smiles, rationale in _TRANSFORMATIONS:
        product_smiles = _apply_transformation(
            parent_mol, src_smarts, repl_smiles
        )
        for product_smi in product_smiles:
            canonical = canonical_smiles(product_smi)
            if canonical is None or canonical in seen:
                continue
            seen.add(canonical)

            variant_score = score_candidate(
                canonical, weights=weights, method=method
            )
            variants.append(
                StructuralVariant(
                    smiles=canonical,
                    transformation=f"{description} ({rationale})",
                    score=variant_score,
                )
            )
            if len(variants) >= max_variants:
                break
        if len(variants) >= max_variants:
            break

    improved: list[StructuralVariant] = []
    if parent_score is not None:
        for variant in variants:
            if variant.score is None:
                continue
            if variant.score.final_score > parent_score.final_score:
                improved.append(variant)
        improved.sort(
            key=lambda v: (v.score.final_score if v.score else 0),  # type: ignore[union-attr]
            reverse=True,
        )

    return LeadOptimizationReport(
        parent_smiles=parent_canonical,
        parent_score=parent_score,
        variants=variants,
        improved_variants=improved,
    )
