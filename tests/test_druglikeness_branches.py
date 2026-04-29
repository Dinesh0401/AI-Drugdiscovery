"""Targeted tests that drive each rule's individual violation branches.

The reference drugs (caffeine, aspirin, ibuprofen, atorvastatin) cover the
common cases but not every threshold, so here we construct synthetic
``MolecularDescriptors`` to exercise the remaining branches.
"""

from __future__ import annotations

from typing import Any

from chemscreen.descriptors import MolecularDescriptors
from chemscreen.druglikeness import evaluate_ghose, evaluate_lipinski, evaluate_veber


def _descriptors(**overrides: Any) -> MolecularDescriptors:
    """Build a baseline drug-like descriptor set, with field overrides."""
    base = dict(
        smiles="x",
        canonical_smiles="x",
        molecular_weight=300.0,
        logp=2.0,
        hbd=2,
        hba=4,
        tpsa=60.0,
        rotatable_bonds=4,
        num_heavy_atoms=25,
        num_rings=2,
        num_aromatic_rings=1,
        formula="X",
    )
    base.update(overrides)
    return MolecularDescriptors(**base)  # type: ignore[arg-type]


def test_lipinski_high_logp_violates() -> None:
    result = evaluate_lipinski(_descriptors(logp=6.0))
    assert not result.passed
    assert any("LogP=" in v for v in result.violations)


def test_lipinski_high_hbd_violates() -> None:
    result = evaluate_lipinski(_descriptors(hbd=6))
    assert not result.passed
    assert any("HBD=" in v for v in result.violations)


def test_lipinski_high_hba_violates() -> None:
    result = evaluate_lipinski(_descriptors(hba=11))
    assert not result.passed
    assert any("HBA=" in v for v in result.violations)


def test_veber_high_rotatable_bonds_violates() -> None:
    result = evaluate_veber(_descriptors(rotatable_bonds=15))
    assert not result.passed
    assert any("RotatableBonds=" in v for v in result.violations)


def test_veber_high_tpsa_violates() -> None:
    result = evaluate_veber(_descriptors(tpsa=200.0))
    assert not result.passed
    assert any("TPSA=" in v for v in result.violations)


def test_ghose_low_mw_violates() -> None:
    result = evaluate_ghose(_descriptors(molecular_weight=100.0))
    assert not result.passed
    assert any("< 160" in v for v in result.violations)


def test_ghose_low_logp_violates() -> None:
    result = evaluate_ghose(_descriptors(logp=-1.0))
    assert not result.passed
    assert any("< -0.4" in v for v in result.violations)


def test_ghose_high_logp_violates() -> None:
    result = evaluate_ghose(_descriptors(logp=7.0))
    assert not result.passed
    assert any("> 5.6" in v for v in result.violations)


def test_ghose_low_heavy_atoms_violates() -> None:
    result = evaluate_ghose(_descriptors(num_heavy_atoms=10))
    assert not result.passed
    assert any("HeavyAtoms=" in v and "< 20" in v for v in result.violations)


def test_ghose_high_heavy_atoms_violates() -> None:
    result = evaluate_ghose(_descriptors(num_heavy_atoms=80))
    assert not result.passed
    assert any("HeavyAtoms=" in v and "> 70" in v for v in result.violations)
