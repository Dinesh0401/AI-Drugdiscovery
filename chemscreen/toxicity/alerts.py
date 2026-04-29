"""Structural alert filters: Brenk + PAINS via RDKit FilterCatalog.

These are published filter sets used in real medicinal-chemistry triage:

  - PAINS (Baell & Holloway 2010): pan-assay interference compounds — react
    promiscuously and produce false positives in screening assays.
  - Brenk (Brenk et al. 2008): 105 structural alerts for unwanted/reactive
    fragments (Michael acceptors, reactive halogens, etc.).

Hits do NOT mean "this molecule is toxic"; they flag liabilities a chemist
should review. Phase 2 surfaces them as ``contributing_alerts`` on each
toxicity prediction so reviewers can see *why* a molecule looks risky beyond
the model's probability.
"""

from __future__ import annotations

from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from chemscreen.descriptors import parse_smiles


def _build_catalog() -> FilterCatalog.FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)


def _build_pains_catalog() -> FilterCatalog.FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog.FilterCatalog(params)


_FULL_CATALOG = _build_catalog()
_PAINS_CATALOG = _build_pains_catalog()


def get_structural_alerts(smiles: str | None) -> list[str] | None:
    """Return a list of triggered Brenk + PAINS alert names.

    Empty list means the molecule is alert-free. None means the SMILES is
    invalid.
    """
    mol = parse_smiles(smiles)
    if mol is None:
        return None
    return [match.GetDescription() for match in _FULL_CATALOG.GetMatches(mol)]


def has_pains(smiles: str | None) -> bool | None:
    """True iff the molecule triggers any PAINS alert. None if invalid SMILES."""
    mol = parse_smiles(smiles)
    if mol is None:
        return None
    return _PAINS_CATALOG.HasMatch(mol)
