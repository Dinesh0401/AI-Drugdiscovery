"""Shared pytest fixtures: canonical reference molecules used across tests.

These SMILES are the same strings RDKit canonicalizes to, so canonical-form
assertions in the test files will be stable across RDKit minor versions.
"""

from __future__ import annotations

import pytest

CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"
IBUPROFEN = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
ATORVASTATIN = (
    "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)"
    "c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O"
)


@pytest.fixture
def caffeine() -> str:
    return CAFFEINE


@pytest.fixture
def aspirin() -> str:
    return ASPIRIN


@pytest.fixture
def ibuprofen() -> str:
    return IBUPROFEN


@pytest.fixture
def atorvastatin() -> str:
    """Known Lipinski violator: MW > 500."""
    return ATORVASTATIN


@pytest.fixture(params=["", "not_a_smiles", "C(C(C", None])
def bad_input(request: pytest.FixtureRequest) -> str | None:
    return request.param
