"""Tests for chemscreen.toxicity.alerts (Brenk + PAINS structural filters).

Test molecules picked empirically against RDKit's FilterCatalog so the
assertions are stable across minor RDKit versions:

  - Ibuprofen, caffeine: alert-clean (used as negative controls).
  - Aspirin: a phenol ester — Brenk flags it (real medchem behavior).
  - Rhodanine: textbook PAINS positive (rhod_sat_A).
  - para-Benzoquinone: textbook PAINS positive (quinone_A).
"""

from __future__ import annotations

from chemscreen.toxicity.alerts import get_structural_alerts, has_pains

RHODANINE = "O=C1CSC(=S)N1"
BENZOQUINONE = "O=C1C=CC(=O)C=C1"


def test_caffeine_is_alert_clean(caffeine: str) -> None:
    """Caffeine has no PAINS or Brenk hits — true negative for the catalog."""
    assert get_structural_alerts(caffeine) == []
    assert has_pains(caffeine) is False


def test_ibuprofen_is_alert_clean(ibuprofen: str) -> None:
    assert get_structural_alerts(ibuprofen) == []
    assert has_pains(ibuprofen) is False


def test_aspirin_triggers_brenk_phenol_ester(aspirin: str) -> None:
    """Aspirin really is a phenol ester — Brenk catalog correctly flags it.

    This is one of the cases where a marketed drug carries a structural
    alert. The alert system flags chemistry, not approvability.
    """
    alerts = get_structural_alerts(aspirin)
    assert alerts is not None
    assert "phenol_ester" in alerts


def test_rhodanine_is_pains_positive() -> None:
    """Rhodanine cores are a canonical PAINS family (Baell & Holloway 2010)."""
    assert has_pains(RHODANINE) is True
    alerts = get_structural_alerts(RHODANINE)
    assert alerts is not None
    assert any("rhod" in a.lower() for a in alerts)


def test_benzoquinone_is_pains_positive() -> None:
    """Quinones are a flagged PAINS class (redox-active / nonspecific)."""
    assert has_pains(BENZOQUINONE) is True
    alerts = get_structural_alerts(BENZOQUINONE)
    assert alerts is not None
    assert any("quinone" in a.lower() or "chinone" in a.lower() for a in alerts)


def test_alerts_invalid_returns_none(bad_input: str | None) -> None:
    assert get_structural_alerts(bad_input) is None
    assert has_pains(bad_input) is None
