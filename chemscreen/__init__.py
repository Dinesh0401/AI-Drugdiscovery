"""chemscreen: drug-candidate screening and ranking core."""

from chemscreen.descriptors import (
    MolecularDescriptors,
    canonical_smiles,
    compute_descriptors,
    is_valid_smiles,
    morgan_fingerprint,
    parse_smiles,
)
from chemscreen.druglikeness import (
    DrugLikeness,
    RuleResult,
    compute_qed,
    evaluate_druglikeness,
    evaluate_ghose,
    evaluate_lipinski,
    evaluate_veber,
)
from chemscreen.explain import Explanation, explain_molecule
from chemscreen.modes import (
    LeadOptimizationReport,
    RiskReport,
    ScreeningReport,
    StructuralVariant,
    analyze_risk,
    screen_batch,
    suggest_variants,
)
from chemscreen.ranking import filter_top_k, rank_candidates, screen_smiles
from chemscreen.similarity import (
    SimilarMolecule,
    filter_by_class,
    find_similar,
    reference_size,
    set_reference_drugs,
    tanimoto,
)
from chemscreen.scoring import (
    CandidateScore,
    ScoreComponents,
    ScoreWeights,
    ScoringMethod,
    aggregate,
    lipinski_component,
    score_candidate,
    synthesis_component,
    toxicity_component,
)
from chemscreen.synthesis import is_synthetically_feasible, sa_score

__all__ = [
    # descriptors
    "MolecularDescriptors",
    "parse_smiles",
    "is_valid_smiles",
    "canonical_smiles",
    "compute_descriptors",
    "morgan_fingerprint",
    # drug-likeness
    "DrugLikeness",
    "RuleResult",
    "evaluate_lipinski",
    "evaluate_veber",
    "evaluate_ghose",
    "compute_qed",
    "evaluate_druglikeness",
    # synthesis
    "sa_score",
    "is_synthetically_feasible",
    # scoring
    "ScoreWeights",
    "ScoreComponents",
    "ScoringMethod",
    "CandidateScore",
    "lipinski_component",
    "synthesis_component",
    "toxicity_component",
    "aggregate",
    "score_candidate",
    # ranking
    "rank_candidates",
    "screen_smiles",
    "filter_top_k",
    # similarity
    "SimilarMolecule",
    "tanimoto",
    "find_similar",
    "filter_by_class",
    "set_reference_drugs",
    "reference_size",
    # explain
    "Explanation",
    "explain_molecule",
    # operational modes (Phase 6)
    "ScreeningReport",
    "screen_batch",
    "RiskReport",
    "analyze_risk",
    "StructuralVariant",
    "LeadOptimizationReport",
    "suggest_variants",
]
