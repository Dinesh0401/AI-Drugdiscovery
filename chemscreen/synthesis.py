"""Synthetic Accessibility scoring (Ertl & Schuffenhauer, 2009).

Wraps RDKit's contributed ``sascorer`` so it integrates with the rest of the
package. Score is on a 1 (easy to make) to 10 (very hard) scale.
"""

from __future__ import annotations

import os
import sys

from rdkit.Chem import RDConfig

from chemscreen.descriptors import parse_smiles

_SA_PATH = os.path.join(RDConfig.RDContribDir, "SA_Score")
if _SA_PATH not in sys.path:
    sys.path.insert(0, _SA_PATH)

import sascorer  # noqa: E402  (path append must happen before import)


def sa_score(smiles: str | None) -> float | None:
    """Synthetic accessibility score, 1 (easy) -> 10 (hard). None on invalid input."""
    mol = parse_smiles(smiles)
    if mol is None:
        return None
    return float(sascorer.calculateScore(mol))


def is_synthetically_feasible(score: float, threshold: float = 6.0) -> bool:
    """True iff the SA score is at or below the feasibility threshold."""
    return score <= threshold
