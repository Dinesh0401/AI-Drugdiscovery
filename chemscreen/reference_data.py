"""Reference lookups (chemical formula -> canonical SMILES).

Lifted from ChemDesigner/services/molecular.py so the same human-friendly
formulas (H2O, CO2, caffeine, ...) work uniformly across the new package.
"""

from __future__ import annotations

FORMULA_TO_SMILES: dict[str, str] = {
    "H2O": "O",
    "CO2": "O=C=O",
    "CH4": "C",
    "NH3": "N",
    "HCl": "Cl",
    "H2SO4": "OS(=O)(=O)O",
    "HNO3": "O[N+](=O)[O-]",
    "C2H4": "C=C",
    "C2H6": "CC",
    "C2H5OH": "CCO",
    "CH3OH": "CO",
    "C6H6": "c1ccccc1",
    "C6H12O6": "OCC(O)C(O)C(O)C(O)C=O",
    "CH3COOH": "CC(=O)O",
    "NaCl": "[Na+].[Cl-]",
    "CaCO3": "[Ca+2].[O-]C([O-])=O",
    "H2O2": "OO",
    "NH4OH": "[NH4+].[OH-]",
    "C8H10N4O2": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # caffeine
    "C9H8O4": "CC(=O)OC1=CC=CC=C1C(=O)O",          # aspirin
    "C13H18O2": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",   # ibuprofen
}


def lookup_formula(formula: str) -> str | None:
    """Resolve a common chemical formula to its canonical SMILES, or None."""
    return FORMULA_TO_SMILES.get(formula.strip())
