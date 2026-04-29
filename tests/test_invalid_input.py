"""Cross-cutting: every public function must handle invalid input gracefully.

Phase 1 contract: invalid SMILES never raise; functions return None or False.
This is what lets downstream pipeline code (Phase 2+) treat invalid inputs as
filterable rather than fatal.
"""

from __future__ import annotations

from chemscreen.descriptors import (
    canonical_smiles,
    compute_descriptors,
    is_valid_smiles,
    morgan_fingerprint,
    parse_smiles,
)
from chemscreen.druglikeness import compute_qed, evaluate_druglikeness
from chemscreen.synthesis import sa_score


def test_parse_smiles_handles_bad_input(bad_input: str | None) -> None:
    assert parse_smiles(bad_input) is None


def test_is_valid_smiles_handles_bad_input(bad_input: str | None) -> None:
    assert is_valid_smiles(bad_input) is False


def test_canonical_smiles_handles_bad_input(bad_input: str | None) -> None:
    assert canonical_smiles(bad_input) is None


def test_compute_descriptors_handles_bad_input(bad_input: str | None) -> None:
    assert compute_descriptors(bad_input) is None


def test_morgan_fingerprint_handles_bad_input(bad_input: str | None) -> None:
    assert morgan_fingerprint(bad_input) is None


def test_compute_qed_handles_bad_input(bad_input: str | None) -> None:
    assert compute_qed(bad_input) is None


def test_evaluate_druglikeness_handles_bad_input(bad_input: str | None) -> None:
    assert evaluate_druglikeness(bad_input) is None


def test_sa_score_handles_bad_input(bad_input: str | None) -> None:
    assert sa_score(bad_input) is None
