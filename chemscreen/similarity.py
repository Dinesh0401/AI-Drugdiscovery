"""Tanimoto similarity search vs a reference set of approved drugs.

Phase 4 of chemscreen. The reference set ships in
``chemscreen.data.approved_drugs`` and can be replaced at runtime via
``set_reference_drugs`` (e.g. with a full DrugBank or ChEMBL export).

Index structure:
  - parallel arrays of names, SMILES, drug classes, and a stacked uint8
    fingerprint matrix (N x 2048). Built lazily on first call so importing
    chemscreen stays fast.

Search uses vectorized Tanimoto over the matrix:
    sim_i = popcount(query AND ref_i) / popcount(query OR ref_i)
which is ~50x faster than per-row Python loops.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from chemscreen.data.approved_drugs import APPROVED_DRUGS
from chemscreen.descriptors import morgan_fingerprint


@dataclass(frozen=True)
class SimilarMolecule:
    """One hit from the similarity search."""

    name: str
    smiles: str
    similarity: float  # Tanimoto, 0..1
    drug_class: str | None


# Mutable module-level reference index. Initialized lazily.
_REF_NAMES: list[str] | None = None
_REF_SMILES: list[str] | None = None
_REF_CLASSES: list[str | None] | None = None
_REF_MATRIX: np.ndarray | None = None
_REF_COUNTS: np.ndarray | None = None
_REF_SOURCE: list[tuple[str, str, str | None]] = list(APPROVED_DRUGS)


def _build_index() -> None:
    """Compute fingerprints for every reference drug and stack them.

    Drops entries whose SMILES fail to parse so the index stays clean.
    """
    global _REF_NAMES, _REF_SMILES, _REF_CLASSES, _REF_MATRIX, _REF_COUNTS

    names: list[str] = []
    smiles_list: list[str] = []
    classes: list[str | None] = []
    fps: list[np.ndarray] = []

    for entry in _REF_SOURCE:
        name = entry[0]
        smiles = entry[1]
        drug_class = entry[2] if len(entry) > 2 else None
        fp = morgan_fingerprint(smiles)
        if fp is None:
            continue
        names.append(name)
        smiles_list.append(smiles)
        classes.append(drug_class)
        fps.append(fp)

    if fps:
        matrix = np.vstack(fps)
        counts = matrix.sum(axis=1)
    else:
        matrix = np.zeros((0, 2048), dtype=np.uint8)
        counts = np.zeros((0,), dtype=np.int64)

    _REF_NAMES = names
    _REF_SMILES = smiles_list
    _REF_CLASSES = classes
    _REF_MATRIX = matrix
    _REF_COUNTS = counts


def _ensure_index() -> None:
    if _REF_MATRIX is None:
        _build_index()


def set_reference_drugs(drugs: Sequence[tuple[str, str, str | None]]) -> None:
    """Replace the reference drug set and rebuild the fingerprint index.

    Each entry is ``(name, smiles, drug_class | None)``. Useful for users
    who want to point similarity search at a full DrugBank or ChEMBL
    approved-drugs export instead of the small bundled set.
    """
    global _REF_SOURCE, _REF_MATRIX
    _REF_SOURCE = list(drugs)
    _REF_MATRIX = None  # force rebuild on next call
    _ensure_index()


def reference_size() -> int:
    """Number of valid molecules in the current reference index."""
    _ensure_index()
    assert _REF_NAMES is not None
    return len(_REF_NAMES)


def tanimoto(fp_a: np.ndarray, fp_b: np.ndarray) -> float:
    """Tanimoto similarity between two binary fingerprint vectors."""
    a = fp_a.astype(bool)
    b = fp_b.astype(bool)
    intersect = int(np.logical_and(a, b).sum())
    union = int(np.logical_or(a, b).sum())
    return intersect / union if union else 0.0


def find_similar(
    query_smiles: str | None,
    top_k: int = 5,
    threshold: float = 0.0,
) -> list[SimilarMolecule] | None:
    """Find Tanimoto-similar approved drugs for a query SMILES.

    Returns up to ``top_k`` drugs with similarity >= ``threshold``, sorted
    descending. Returns None if the query SMILES is invalid.
    """
    if top_k < 0:
        raise ValueError("top_k must be non-negative")
    if not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must be in [0, 1]")

    query_fp = morgan_fingerprint(query_smiles)
    if query_fp is None:
        return None

    _ensure_index()
    assert _REF_MATRIX is not None
    assert _REF_COUNTS is not None
    assert _REF_NAMES is not None
    assert _REF_SMILES is not None
    assert _REF_CLASSES is not None

    if _REF_MATRIX.shape[0] == 0:
        return []

    # Vectorized Tanimoto: query-vs-all in one shot
    query_count = int(query_fp.sum())
    intersect = (_REF_MATRIX & query_fp).sum(axis=1)
    union = _REF_COUNTS + query_count - intersect
    sims = np.where(union > 0, intersect / np.maximum(union, 1), 0.0)

    # Filter by threshold, sort descending, take top_k
    keep = np.where(sims >= threshold)[0]
    ordered = keep[np.argsort(-sims[keep])][:top_k]

    return [
        SimilarMolecule(
            name=_REF_NAMES[i],
            smiles=_REF_SMILES[i],
            similarity=float(sims[i]),
            drug_class=_REF_CLASSES[i],
        )
        for i in ordered
    ]


def filter_by_class(
    query_smiles: str | None,
    drug_class: str,
    top_k: int = 5,
) -> list[SimilarMolecule] | None:
    """Same as ``find_similar`` but restricted to one drug class."""
    results = find_similar(query_smiles, top_k=10_000, threshold=0.0)
    if results is None:
        return None
    return [r for r in results if r.drug_class == drug_class][:top_k]
