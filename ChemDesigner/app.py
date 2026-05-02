"""chemscreen Flask UI (Phase 8).

Wires the legacy ChemDesigner shell to the chemscreen package. Replaces
the fake "latent space generator" / "LSTM state" pages with real,
deterministic, scientifically-grounded views:

  - /analyze   single-molecule full explanation
  - /screen    batch screening with rejection accounting
  - /risk      toxicity + structural-alerts focused report
  - /optimize  lead optimization via bioisostere swaps

Run from the repo root:

    uv run --extra web python ChemDesigner/app.py
"""

from __future__ import annotations

import logging
import os
from dataclasses import asdict

from flask import Flask, jsonify, render_template, request, Response
from rdkit.Chem import Draw

from chemscreen import (
    canonical_smiles,
    explain_molecule,
    is_valid_smiles,
    parse_smiles,
)
from chemscreen.modes import analyze_risk, screen_batch, suggest_variants
from chemscreen.reference_data import lookup_formula

logging.basicConfig(level=logging.INFO)

app = Flask(__name__)
app.secret_key = os.environ.get("SESSION_SECRET", "chemscreen-flask-key")


def _normalize_input(text: str | None) -> str | None:
    """Accept either SMILES or a common chemical formula. Return canonical SMILES."""
    if not text:
        return None
    text = text.strip()
    if not text:
        return None
    # Try formula lookup first ("H2O", "CO2", etc.)
    formula_hit = lookup_formula(text)
    if formula_hit:
        return canonical_smiles(formula_hit)
    if is_valid_smiles(text):
        return canonical_smiles(text)
    return None


# --- Pages ----------------------------------------------------------------


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/analyze")
def analyze_page():
    return render_template("analyze.html")


@app.route("/screen")
def screen_page():
    return render_template("screen.html")


@app.route("/risk")
def risk_page():
    return render_template("risk.html")


@app.route("/optimize")
def optimize_page():
    return render_template("optimize.html")


# --- API ------------------------------------------------------------------


@app.post("/api/analyze")
def api_analyze():
    """Full Phase-1-through-5 explanation for one molecule."""
    payload = request.get_json(silent=True) or request.form
    raw = payload.get("input", "")
    smiles = _normalize_input(raw)
    if smiles is None:
        return jsonify({"success": False, "error": "Invalid SMILES or unknown formula"}), 400

    explanation = explain_molecule(smiles)
    if explanation is None:
        return jsonify({"success": False, "error": "Could not analyze"}), 500

    return jsonify({"success": True, "input": raw, **explanation.to_dict()})


@app.post("/api/screen")
def api_screen():
    """Batch screen a list of SMILES (one per line, comma, or whitespace separated)."""
    payload = request.get_json(silent=True) or request.form
    raw = payload.get("smiles_list", "") or payload.get("smiles", "")
    smiles_list = [
        s.strip()
        for s in raw.replace(",", "\n").splitlines()
        if s.strip()
    ]
    if not smiles_list:
        return jsonify({"success": False, "error": "Provide at least one SMILES"}), 400

    require_lipinski = str(payload.get("require_lipinski", "false")).lower() in {"true", "on", "1"}
    try:
        top_k = int(payload.get("top_k", 25))
    except (TypeError, ValueError):
        top_k = 25
    top_k = max(1, min(top_k, 200))

    report = screen_batch(
        smiles_list,
        top_k=top_k,
        require_lipinski=require_lipinski,
    )
    return jsonify({
        "success": True,
        "n_input": report.n_input,
        "n_scored": report.n_scored,
        "n_rejected": report.n_rejected,
        "method": report.method,
        "candidates": [
            {
                "rank": c.rank,
                "smiles": c.smiles,
                "final_score": round(c.final_score, 4),
                "qed": round(c.components.qed, 3),
                "lipinski": round(c.components.lipinski, 3),
                "synthesis": round(c.components.synthesis, 3),
                "toxicity": round(c.components.toxicity, 3),
            }
            for c in report.candidates
        ],
    })


@app.post("/api/risk")
def api_risk():
    """Focused toxicity + structural-alerts report for one molecule."""
    payload = request.get_json(silent=True) or request.form
    smiles = _normalize_input(payload.get("input", ""))
    if smiles is None:
        return jsonify({"success": False, "error": "Invalid SMILES or unknown formula"}), 400

    report = analyze_risk(smiles)
    if report is None:
        return jsonify({"success": False, "error": "Could not analyze"}), 500

    return jsonify({
        "success": True,
        "smiles": report.smiles,
        "risk_level": report.risk_level,
        "summary": report.summary,
        "recommendation": report.recommendation,
        "tox_predictions": {
            ep: {
                "classification": p.classification,
                "risk_score": round(p.risk_score, 3),
            }
            for ep, p in report.tox_predictions.items()
        },
        "structural_alerts": report.structural_alerts,
        "high_risk_endpoints": report.high_risk_endpoints,
        "moderate_risk_endpoints": report.moderate_risk_endpoints,
        "similar_drugs": [
            {
                "name": d.name,
                "similarity": round(d.similarity, 3),
                "drug_class": d.drug_class,
            }
            for d in report.similar_drugs
        ],
    })


@app.post("/api/optimize")
def api_optimize():
    """Bioisostere-swap variant generation + re-scoring."""
    payload = request.get_json(silent=True) or request.form
    smiles = _normalize_input(payload.get("input", ""))
    if smiles is None:
        return jsonify({"success": False, "error": "Invalid SMILES or unknown formula"}), 400

    report = suggest_variants(smiles)
    if report is None:
        return jsonify({"success": False, "error": "Could not optimize"}), 500

    parent_score = report.parent_score.final_score if report.parent_score else None

    def _serialize(variant):
        return {
            "smiles": variant.smiles,
            "transformation": variant.transformation,
            "final_score": (round(variant.score.final_score, 4)
                            if variant.score else None),
            "delta": (round(variant.score.final_score - parent_score, 4)
                      if (variant.score and parent_score is not None) else None),
        }

    return jsonify({
        "success": True,
        "parent_smiles": report.parent_smiles,
        "parent_score": round(parent_score, 4) if parent_score is not None else None,
        "variants": [_serialize(v) for v in report.variants],
        "improved_variants": [_serialize(v) for v in report.improved_variants],
    })


@app.get("/api/svg")
def api_svg():
    """Render a SMILES as an SVG image. ?smiles=<encoded SMILES>&w=300&h=200."""
    smiles = request.args.get("smiles", "")
    try:
        width = max(100, min(int(request.args.get("w", 300)), 800))
        height = max(80, min(int(request.args.get("h", 200)), 600))
    except (TypeError, ValueError):
        width, height = 300, 200

    mol = parse_smiles(smiles)
    if mol is None:
        return Response("invalid SMILES", status=400)

    drawer = Draw.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return Response(drawer.GetDrawingText(), mimetype="image/svg+xml")


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=False)
