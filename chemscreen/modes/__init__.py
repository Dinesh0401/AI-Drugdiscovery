"""Phase 6 — operational modes.

Three named workflows that wrap the Phase 1–5 building blocks:

  - ``screening`` — batch score + rank with rejection accounting
  - ``risk_analysis`` — focused toxicity + structural-alerts report
  - ``lead_optimization`` — bioisostere-swap variant generation + scoring
"""

from chemscreen.modes.lead_optimization import (
    LeadOptimizationReport,
    StructuralVariant,
    suggest_variants,
)
from chemscreen.modes.risk_analysis import RiskReport, analyze_risk
from chemscreen.modes.screening import ScreeningReport, screen_batch

__all__ = [
    "ScreeningReport",
    "screen_batch",
    "RiskReport",
    "analyze_risk",
    "StructuralVariant",
    "LeadOptimizationReport",
    "suggest_variants",
]
