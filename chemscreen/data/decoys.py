"""Curated set of ~80 non-drug molecules used as benchmark decoys.

Categories: industrial intermediates, agrochemicals, organic solvents,
simple aromatics, fragrances, and other commercially-available compounds
that are NOT FDA-approved drugs.

Used by ``chemscreen.validation`` to benchmark whether the deterministic
scoring pipeline distinguishes drug-like molecules from random non-drugs.
A test in tests/test_validation.py verifies every entry parses with RDKit
and that no decoy SMILES collides with the approved-drug reference set.

Sourced from PubChem; canonicalized by RDKit on indexing.
"""

from __future__ import annotations

# fmt: off
DECOYS: list[tuple[str, str]] = [
    # (name, SMILES)

    # Solvents / industrial intermediates
    ("methanol",            "CO"),
    ("ethanol",             "CCO"),
    ("isopropanol",         "CC(C)O"),
    ("acetone",             "CC(=O)C"),
    ("methyl-ethyl-ketone", "CCC(=O)C"),
    ("cyclohexanone",       "O=C1CCCCC1"),
    ("ethyl-acetate",       "CCOC(C)=O"),
    ("tetrahydrofuran",     "C1CCOC1"),
    ("dioxane",             "C1COCCO1"),
    ("formamide",           "O=CN"),
    ("dmf",                 "CN(C)C=O"),
    ("dmso",                "CS(C)=O"),
    ("acetic-acid",         "CC(=O)O"),
    ("propionic-acid",      "CCC(=O)O"),
    ("formic-acid",         "OC=O"),
    ("triethylamine",       "CCN(CC)CC"),
    ("pyridine",            "c1ccncc1"),
    ("imidazole",           "c1cnc[nH]1"),
    ("pyrrolidine",         "C1CCNC1"),
    ("piperidine",          "C1CCNCC1"),
    ("morpholine",          "C1COCCN1"),
    ("furan",               "c1ccoc1"),
    ("thiophene",           "c1ccsc1"),

    # Simple aromatics
    ("benzene",             "c1ccccc1"),
    ("toluene",             "Cc1ccccc1"),
    ("o-xylene",            "Cc1ccccc1C"),
    ("m-xylene",            "Cc1cccc(C)c1"),
    ("p-xylene",            "Cc1ccc(C)cc1"),
    ("ethylbenzene",        "CCc1ccccc1"),
    ("styrene",             "C=Cc1ccccc1"),
    ("mesitylene",          "Cc1cc(C)cc(C)c1"),
    ("biphenyl",            "c1ccc(-c2ccccc2)cc1"),
    ("naphthalene",         "c1ccc2ccccc2c1"),
    ("anthracene",          "c1ccc2cc3ccccc3cc2c1"),
    ("phenanthrene",        "c1ccc2ccc3ccccc3c2c1"),
    ("indene",              "C1=CC2=CC=CC=C2C1"),
    ("tetralin",            "C1CCC2=CC=CC=C2C1"),

    # Phenols / catechols / quinones
    ("phenol",              "Oc1ccccc1"),
    ("4-chlorophenol",      "Oc1ccc(Cl)cc1"),
    ("hydroquinone",        "Oc1ccc(O)cc1"),
    ("resorcinol",          "Oc1cccc(O)c1"),
    ("catechol",            "Oc1ccccc1O"),
    ("anisole",             "COc1ccccc1"),
    ("guaiacol",            "COc1ccccc1O"),
    ("eugenol",             "C=CCc1ccc(O)c(OC)c1"),
    ("benzoquinone",        "O=C1C=CC(=O)C=C1"),

    # Aldehydes / ketones / esters
    ("benzaldehyde",        "O=Cc1ccccc1"),
    ("vanillin",            "COc1cc(C=O)ccc1O"),
    ("cinnamaldehyde",      "O=C/C=C/c1ccccc1"),
    ("acetaldehyde",        "CC=O"),
    ("formaldehyde",        "C=O"),
    ("benzoic-acid",        "O=C(O)c1ccccc1"),
    ("salicylic-acid",      "O=C(O)c1ccccc1O"),
    ("methyl-salicylate",   "COC(=O)c1ccccc1O"),
    ("benzyl-acetate",      "CC(=O)OCc1ccccc1"),

    # Halogenated / agrochemicals
    ("chlorobenzene",       "Clc1ccccc1"),
    ("dichlorobenzene",     "Clc1ccc(Cl)cc1"),
    ("trichloroethylene",   "ClC=C(Cl)Cl"),
    ("ddt",                 "ClC(c1ccc(Cl)cc1)(c1ccc(Cl)cc1)C(Cl)(Cl)Cl"),
    ("atrazine",            "CCNc1nc(Cl)nc(NC(C)C)n1"),
    ("2,4-d",               "OC(=O)COc1ccc(Cl)cc1Cl"),
    ("malathion",           "CCOC(=O)CC(SP(=S)(OC)OC)C(=O)OCC"),
    ("alachlor",            "CCc1cccc(CC)c1N(COC)C(=O)CCl"),

    # Amines / basic
    ("aniline",             "Nc1ccccc1"),
    ("4-nitroaniline",      "Nc1ccc([N+](=O)[O-])cc1"),
    ("dimethylamine",       "CNC"),
    ("ethanolamine",        "NCCO"),
    ("benzylamine",         "NCc1ccccc1"),

    # Polymers / monomers
    ("ethylene-glycol",     "OCCO"),
    ("propylene-glycol",    "CC(O)CO"),
    ("acrylic-acid",        "C=CC(=O)O"),
    ("methacrylic-acid",    "C=C(C)C(=O)O"),
    ("vinyl-acetate",       "C=COC(C)=O"),
    ("epichlorohydrin",     "ClCC1CO1"),
    ("caprolactam",         "O=C1CCCCCN1"),
    ("melamine",            "Nc1nc(N)nc(N)n1"),

    # Heterocycles / dyes / misc
    ("indole",              "c1ccc2[nH]ccc2c1"),
    ("benzofuran",          "c1ccc2occc2c1"),
    ("benzothiophene",      "c1ccc2sccc2c1"),
    ("quinoline",           "c1ccc2ncccc2c1"),
    ("isoquinoline",        "c1ccc2cnccc2c1"),
    ("acridine",            "c1ccc2nc3ccccc3cc2c1"),
    ("urea",                "NC(=O)N"),
    ("biuret",              "NC(=O)NC(=O)N"),
    ("nicotine",            "CN1CCCC1c1cccnc1"),  # alkaloid, not an FDA drug
    ("caffeic-acid",        "O=C(O)/C=C/c1ccc(O)c(O)c1"),
    ("ferulic-acid",        "COc1cc(/C=C/C(=O)O)ccc1O"),
    ("vanillyl-alcohol",    "COc1cc(CO)ccc1O"),
]
# fmt: on


def smiles_only() -> list[str]:
    """Just the SMILES strings, in catalog order."""
    return [smi for _, smi in DECOYS]
