"""Molecular descriptors and Morgan fingerprints (RDKit).

Public functions all accept SMILES strings and return None on invalid input
rather than raising — validation happens once at the boundary.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

try:
    from rdkit.Chem import rdFingerprintGenerator

    _MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    _HAS_FP_GENERATOR = True
except (ImportError, AttributeError):
    from rdkit.Chem import AllChem  # noqa: F401  (used in fallback)

    _HAS_FP_GENERATOR = False


@dataclass(frozen=True)
class MolecularDescriptors:
    """Numerical descriptors for one molecule."""

    smiles: str
    canonical_smiles: str
    molecular_weight: float
    logp: float
    hbd: int
    hba: int
    tpsa: float
    rotatable_bonds: int
    num_heavy_atoms: int
    num_rings: int
    num_aromatic_rings: int
    formula: str


def parse_smiles(smiles: str | None) -> Chem.Mol | None:
    """Parse SMILES into an RDKit Mol, or return None if invalid."""
    if not smiles or not isinstance(smiles, str):
        return None
    return Chem.MolFromSmiles(smiles)


def is_valid_smiles(smiles: str | None) -> bool:
    """True iff the input parses as a valid molecule."""
    return parse_smiles(smiles) is not None


def canonical_smiles(smiles: str | None) -> str | None:
    """Return RDKit canonical SMILES for the input, or None if invalid."""
    mol = parse_smiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def compute_descriptors(smiles: str | None) -> MolecularDescriptors | None:
    """Compute the full descriptor set for one molecule, or None if invalid."""
    mol = parse_smiles(smiles)
    if mol is None:
        return None

    return MolecularDescriptors(
        smiles=smiles,  # type: ignore[arg-type]  # parse_smiles guarantees non-None
        canonical_smiles=Chem.MolToSmiles(mol),
        molecular_weight=float(Descriptors.MolWt(mol)),
        logp=float(Crippen.MolLogP(mol)),
        hbd=int(Lipinski.NumHDonors(mol)),
        hba=int(Lipinski.NumHAcceptors(mol)),
        tpsa=float(Descriptors.TPSA(mol)),
        rotatable_bonds=int(Descriptors.NumRotatableBonds(mol)),
        num_heavy_atoms=int(mol.GetNumHeavyAtoms()),
        num_rings=int(Descriptors.RingCount(mol)),
        num_aromatic_rings=int(Descriptors.NumAromaticRings(mol)),
        formula=rdMolDescriptors.CalcMolFormula(mol),
    )


def morgan_fingerprint(
    smiles: str | None,
    radius: int = 2,
    n_bits: int = 2048,
) -> np.ndarray | None:
    """Morgan circular fingerprint as a uint8 numpy array of length ``n_bits``.

    Returns None on invalid input. Uses the modern ``rdFingerprintGenerator``
    when available and falls back to ``AllChem.GetMorganFingerprintAsBitVect``
    on older RDKit versions.
    """
    mol = parse_smiles(smiles)
    if mol is None:
        return None

    if _HAS_FP_GENERATOR and radius == 2 and n_bits == 2048:
        fp = _MORGAN_GEN.GetFingerprint(mol)
    elif _HAS_FP_GENERATOR:
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
        fp = gen.GetFingerprint(mol)
    else:
        from rdkit.Chem import AllChem

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)

    arr = np.zeros((n_bits,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr
