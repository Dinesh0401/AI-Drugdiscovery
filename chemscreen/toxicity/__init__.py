"""Toxicity prediction package: Ames mutagenicity, hERG cardiac risk, DILI hepatotoxicity.

Each endpoint is a Random Forest classifier trained on 2048-bit Morgan
fingerprints. Models are persisted as joblib files under ``models/`` and loaded
lazily on first call.

Datasets come from Therapeutics Data Commons (Huang et al., NeurIPS 2021):
  - AMES: Hansen et al. mutagenicity benchmark (Hansen 2009)
  - hERG: Karim et al. hERG inhibition
  - DILI: Drug-induced liver injury (Xu et al.)

The training script lives at ``scripts/train_toxicity_models.py`` — run once
to populate ``models/`` with reproducible joblib artifacts.
"""

from chemscreen.toxicity.alerts import get_structural_alerts, has_pains
from chemscreen.toxicity.ames import predict_ames
from chemscreen.toxicity.base import (
    RFToxicityModel,
    ToxicityModel,
    ToxicityPrediction,
)
from chemscreen.toxicity.dili import predict_dili
from chemscreen.toxicity.ensemble import predict_all_toxicity
from chemscreen.toxicity.herg import predict_herg

__all__ = [
    "RFToxicityModel",
    "ToxicityModel",
    "ToxicityPrediction",
    "get_structural_alerts",
    "has_pains",
    "predict_ames",
    "predict_herg",
    "predict_dili",
    "predict_all_toxicity",
]
