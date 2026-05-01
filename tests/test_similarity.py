"""Tests for chemscreen.similarity (Phase 4)."""

from __future__ import annotations

import numpy as np
import pytest

from chemscreen.data.approved_drugs import APPROVED_DRUGS
from chemscreen.descriptors import is_valid_smiles, morgan_fingerprint
from chemscreen.similarity import (
    SimilarMolecule,
    filter_by_class,
    find_similar,
    reference_size,
    set_reference_drugs,
    tanimoto,
)

_BUILTIN_REFERENCE = list(APPROVED_DRUGS)


@pytest.fixture(autouse=True)
def restore_reference_after_test():
    """Tests that override the reference set must not leak into other tests."""
    yield
    set_reference_drugs(_BUILTIN_REFERENCE)


# --- Reference data integrity --------------------------------------------


def test_every_approved_drug_smiles_is_valid() -> None:
    """Regression: a typo in approved_drugs.py would silently shrink the index."""
    for name, smiles, _drug_class in APPROVED_DRUGS:
        assert is_valid_smiles(smiles), f"Invalid SMILES for {name}: {smiles!r}"


def test_reference_size_matches_data() -> None:
    set_reference_drugs(_BUILTIN_REFERENCE)
    assert reference_size() == len(APPROVED_DRUGS)


# --- Tanimoto math --------------------------------------------------------


def test_tanimoto_self_similarity_is_one(aspirin: str) -> None:
    fp = morgan_fingerprint(aspirin)
    assert fp is not None
    assert tanimoto(fp, fp) == pytest.approx(1.0)


def test_tanimoto_disjoint_is_zero() -> None:
    a = np.array([1, 1, 1, 0, 0, 0], dtype=np.uint8)
    b = np.array([0, 0, 0, 1, 1, 1], dtype=np.uint8)
    assert tanimoto(a, b) == 0.0


def test_tanimoto_partial_overlap() -> None:
    a = np.array([1, 1, 1, 0], dtype=np.uint8)
    b = np.array([1, 1, 0, 1], dtype=np.uint8)
    # intersect=2, union=4 -> 0.5
    assert tanimoto(a, b) == pytest.approx(0.5)


def test_tanimoto_both_empty_returns_zero() -> None:
    a = np.zeros(8, dtype=np.uint8)
    b = np.zeros(8, dtype=np.uint8)
    assert tanimoto(a, b) == 0.0


def test_tanimoto_in_range_for_real_molecules(aspirin: str, caffeine: str) -> None:
    fp_a = morgan_fingerprint(aspirin)
    fp_c = morgan_fingerprint(caffeine)
    assert fp_a is not None and fp_c is not None
    sim = tanimoto(fp_a, fp_c)
    assert 0.0 <= sim <= 1.0


# --- find_similar end-to-end ---------------------------------------------


def test_find_similar_self_match_is_top(aspirin: str) -> None:
    """An exact match in the reference set should score similarity ~1.0."""
    results = find_similar(aspirin, top_k=1)
    assert results is not None
    assert len(results) == 1
    assert results[0].name == "aspirin"
    assert results[0].similarity == pytest.approx(1.0)


def test_find_similar_returns_top_k(aspirin: str) -> None:
    results = find_similar(aspirin, top_k=3)
    assert results is not None
    assert len(results) == 3
    # Sorted descending
    sims = [r.similarity for r in results]
    assert sims == sorted(sims, reverse=True)


def test_find_similar_all_return_smilarmolecule(aspirin: str) -> None:
    results = find_similar(aspirin, top_k=5)
    assert results is not None
    for hit in results:
        assert isinstance(hit, SimilarMolecule)
        assert isinstance(hit.name, str)
        assert isinstance(hit.smiles, str)
        assert 0.0 <= hit.similarity <= 1.0


def test_find_similar_threshold_filters(aspirin: str) -> None:
    """High threshold should filter out weak matches."""
    high_thresh = find_similar(aspirin, top_k=100, threshold=0.99)
    assert high_thresh is not None
    # Only aspirin itself (sim=1.0) and any ~identical analogs survive
    for hit in high_thresh:
        assert hit.similarity >= 0.99


def test_find_similar_atorvastatin_retrieves_other_statins() -> None:
    """A statin should retrieve other statins as its top non-self hits."""
    atorvastatin = (
        "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)"
        "c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O"
    )
    results = find_similar(atorvastatin, top_k=3)
    assert results is not None
    # At least one of the other top-3 hits is also a statin
    non_self = [r for r in results if r.name != "atorvastatin"]
    assert any(r.drug_class == "Statin" for r in non_self)


def test_find_similar_invalid_returns_none() -> None:
    assert find_similar("") is None
    assert find_similar("not_a_smiles") is None
    assert find_similar(None) is None


def test_find_similar_negative_topk_raises(aspirin: str) -> None:
    with pytest.raises(ValueError, match="top_k"):
        find_similar(aspirin, top_k=-1)


def test_find_similar_threshold_out_of_range_raises(aspirin: str) -> None:
    with pytest.raises(ValueError, match="threshold"):
        find_similar(aspirin, threshold=1.5)
    with pytest.raises(ValueError, match="threshold"):
        find_similar(aspirin, threshold=-0.1)


def test_find_similar_topk_zero_returns_empty(aspirin: str) -> None:
    assert find_similar(aspirin, top_k=0) == []


# --- filter_by_class -----------------------------------------------------


def test_filter_by_class_returns_only_that_class(aspirin: str) -> None:
    results = filter_by_class(aspirin, drug_class="NSAID", top_k=10)
    assert results is not None
    assert len(results) > 0
    assert all(r.drug_class == "NSAID" for r in results)


def test_filter_by_class_unknown_returns_empty(aspirin: str) -> None:
    assert filter_by_class(aspirin, drug_class="UNKNOWN_CLASS") == []


def test_filter_by_class_invalid_returns_none() -> None:
    assert filter_by_class("not_a_smiles", drug_class="NSAID") is None


# --- set_reference_drugs --------------------------------------------------


def test_set_reference_drugs_replaces_index(aspirin: str) -> None:
    """Overriding the reference set must take effect immediately."""
    custom: list[tuple[str, str, str | None]] = [
        ("test_aspirin", aspirin, "TestClass"),
    ]
    set_reference_drugs(custom)
    assert reference_size() == 1

    results = find_similar(aspirin, top_k=5)
    assert results is not None
    assert len(results) == 1
    assert results[0].name == "test_aspirin"
    assert results[0].drug_class == "TestClass"


def test_set_reference_drugs_drops_invalid_entries(aspirin: str) -> None:
    custom: list[tuple[str, str, str | None]] = [
        ("good", aspirin, None),
        ("bad", "not_a_smiles", None),
    ]
    set_reference_drugs(custom)
    assert reference_size() == 1  # invalid entry dropped


def test_empty_reference_returns_empty_list(aspirin: str) -> None:
    set_reference_drugs([])
    assert reference_size() == 0
    assert find_similar(aspirin) == []
