"""Drug-likeness rule evaluation: Lipinski, Veber, Ghose, plus QED.

Each rule produces a ``RuleResult`` with feature-grounded violation strings
(e.g. ``"MW=558.65 > 500"``). Phase 5 explainability consumes these directly:
no separate explanation computation needed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit.Chem import QED

from chemscreen.descriptors import MolecularDescriptors, compute_descriptors, parse_smiles


@dataclass(frozen=True)
class RuleResult:
    """Outcome of evaluating one drug-likeness rule against a molecule."""

    name: str
    passed: bool
    violations: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class DrugLikeness:
    """Bundled drug-likeness verdict for one molecule."""

    qed: float
    lipinski: RuleResult
    veber: RuleResult
    ghose: RuleResult


def evaluate_lipinski(d: MolecularDescriptors) -> RuleResult:
    """Lipinski's Rule of Five (Lipinski et al., 1997).

    Thresholds: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10.
    """
    violations: list[str] = []
    if d.molecular_weight > 500:
        violations.append(f"MW={d.molecular_weight:.2f} > 500")
    if d.logp > 5:
        violations.append(f"LogP={d.logp:.2f} > 5")
    if d.hbd > 5:
        violations.append(f"HBD={d.hbd} > 5")
    if d.hba > 10:
        violations.append(f"HBA={d.hba} > 10")

    return RuleResult(name="Lipinski Rule of Five", passed=not violations, violations=violations)


def evaluate_veber(d: MolecularDescriptors) -> RuleResult:
    """Veber's rule for oral bioavailability (Veber et al., 2002).

    Thresholds: rotatable bonds <= 10, TPSA <= 140 A^2.
    """
    violations: list[str] = []
    if d.rotatable_bonds > 10:
        violations.append(f"RotatableBonds={d.rotatable_bonds} > 10")
    if d.tpsa > 140:
        violations.append(f"TPSA={d.tpsa:.2f} > 140")

    return RuleResult(name="Veber Rule", passed=not violations, violations=violations)


def evaluate_ghose(d: MolecularDescriptors) -> RuleResult:
    """Ghose filter (Ghose et al., 1999).

    Thresholds: 160 <= MW <= 480, -0.4 <= LogP <= 5.6, 20 <= heavy atoms <= 70.
    """
    violations: list[str] = []
    if d.molecular_weight < 160:
        violations.append(f"MW={d.molecular_weight:.2f} < 160")
    if d.molecular_weight > 480:
        violations.append(f"MW={d.molecular_weight:.2f} > 480")
    if d.logp < -0.4:
        violations.append(f"LogP={d.logp:.2f} < -0.4")
    if d.logp > 5.6:
        violations.append(f"LogP={d.logp:.2f} > 5.6")
    if d.num_heavy_atoms < 20:
        violations.append(f"HeavyAtoms={d.num_heavy_atoms} < 20")
    if d.num_heavy_atoms > 70:
        violations.append(f"HeavyAtoms={d.num_heavy_atoms} > 70")

    return RuleResult(name="Ghose Filter", passed=not violations, violations=violations)


def compute_qed(smiles: str | None) -> float | None:
    """Quantitative Estimate of Drug-likeness (Bickerton et al., 2012). Range 0..1."""
    mol = parse_smiles(smiles)
    if mol is None:
        return None
    return float(QED.qed(mol))


def evaluate_druglikeness(smiles: str | None) -> DrugLikeness | None:
    """Evaluate QED + all three rules for one molecule, or None if invalid."""
    descriptors = compute_descriptors(smiles)
    if descriptors is None:
        return None

    qed = compute_qed(smiles)
    assert qed is not None  # compute_descriptors already validated the SMILES

    return DrugLikeness(
        qed=qed,
        lipinski=evaluate_lipinski(descriptors),
        veber=evaluate_veber(descriptors),
        ghose=evaluate_ghose(descriptors),
    )
