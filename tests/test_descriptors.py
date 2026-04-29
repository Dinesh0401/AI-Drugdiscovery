"""Tests for chemscreen.descriptors against known reference values."""

from __future__ import annotations

import numpy as np
import pytest

from chemscreen.descriptors import (
    canonical_smiles,
    compute_descriptors,
    is_valid_smiles,
    morgan_fingerprint,
    parse_smiles,
)


def test_parse_valid_returns_mol(caffeine: str) -> None:
    assert parse_smiles(caffeine) is not None


def test_parse_invalid_returns_none() -> None:
    assert parse_smiles("not_a_smiles") is None
    assert parse_smiles("") is None
    assert parse_smiles(None) is None


def test_is_valid_smiles(caffeine: str) -> None:
    assert is_valid_smiles(caffeine) is True
    assert is_valid_smiles("not_a_smiles") is False
    assert is_valid_smiles("") is False


def test_canonical_smiles_idempotent(aspirin: str) -> None:
    once = canonical_smiles(aspirin)
    twice = canonical_smiles(once)
    assert once is not None
    assert once == twice


def test_caffeine_descriptors(caffeine: str) -> None:
    d = compute_descriptors(caffeine)
    assert d is not None
    assert d.molecular_weight == pytest.approx(194.19, abs=0.1)
    assert d.hbd == 0
    assert d.hba >= 3
    assert d.formula == "C8H10N4O2"
    assert d.num_aromatic_rings == 2
    assert d.num_heavy_atoms == 14


def test_aspirin_descriptors(aspirin: str) -> None:
    d = compute_descriptors(aspirin)
    assert d is not None
    assert d.molecular_weight == pytest.approx(180.16, abs=0.1)
    assert d.formula == "C9H8O4"
    assert d.hbd == 1
    assert d.num_aromatic_rings == 1


def test_ibuprofen_descriptors(ibuprofen: str) -> None:
    d = compute_descriptors(ibuprofen)
    assert d is not None
    assert d.molecular_weight == pytest.approx(206.28, abs=0.1)
    assert d.logp == pytest.approx(3.5, abs=0.5)
    assert d.hbd == 1


def test_atorvastatin_is_heavy(atorvastatin: str) -> None:
    """Atorvastatin must report MW > 500 (this is what triggers Lipinski)."""
    d = compute_descriptors(atorvastatin)
    assert d is not None
    assert d.molecular_weight > 500


def test_morgan_fingerprint_shape(caffeine: str) -> None:
    fp = morgan_fingerprint(caffeine)
    assert fp is not None
    assert fp.shape == (2048,)
    assert fp.dtype == np.uint8
    assert fp.sum() > 0  # non-trivial molecule has bits set
    assert set(np.unique(fp).tolist()).issubset({0, 1})


def test_morgan_fingerprint_custom_size(aspirin: str) -> None:
    fp = morgan_fingerprint(aspirin, radius=3, n_bits=1024)
    assert fp is not None
    assert fp.shape == (1024,)


def test_morgan_fingerprint_distinguishes_molecules(caffeine: str, aspirin: str) -> None:
    """Different molecules should produce different fingerprints."""
    fp_a = morgan_fingerprint(caffeine)
    fp_b = morgan_fingerprint(aspirin)
    assert fp_a is not None and fp_b is not None
    assert not np.array_equal(fp_a, fp_b)


def test_compute_descriptors_invalid_returns_none() -> None:
    assert compute_descriptors("") is None
    assert compute_descriptors("not_a_smiles") is None


def test_morgan_fingerprint_invalid_returns_none() -> None:
    assert morgan_fingerprint("") is None
    assert morgan_fingerprint("not_a_smiles") is None
