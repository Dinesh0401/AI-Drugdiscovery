"""Curated reference set of FDA-approved drugs for similarity search.

Names and SMILES sourced from PubChem and DrugBank Open. This is a small
representative set covering broad therapeutic classes -- not a complete
drug database. For production use, swap with a full DrugBank-Open or
ChEMBL approved-drugs export via ``similarity.set_reference_drugs()``.

Each entry: (name, SMILES, drug_class). RDKit-canonicalizes on indexing
so the SMILES strings here can be in any valid form.

A test in tests/test_similarity.py verifies every entry parses.
"""

from __future__ import annotations

# fmt: off
APPROVED_DRUGS: list[tuple[str, str, str]] = [
    # NSAIDs / analgesics
    ("aspirin",        "CC(=O)Oc1ccccc1C(=O)O",                              "NSAID"),
    ("ibuprofen",      "CC(C)Cc1ccc(C(C)C(=O)O)cc1",                         "NSAID"),
    ("naproxen",       "COc1ccc2cc(C(C)C(=O)O)ccc2c1",                       "NSAID"),
    ("diclofenac",     "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl",                     "NSAID"),
    ("celecoxib",      "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1", "COX-2 inhibitor"),
    ("ketoprofen",     "CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1",                  "NSAID"),
    ("indomethacin",   "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",     "NSAID"),

    # Acetaminophen / opioids
    ("acetaminophen",  "CC(=O)Nc1ccc(O)cc1",                                 "Analgesic"),
    ("tramadol",       "COc1cccc(C2(O)CCCCC2CN(C)C)c1",                      "Opioid"),
    ("morphine",       "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "Opioid"),
    ("codeine",        "CN1CC[C@]23c4c5ccc(OC)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "Opioid"),

    # Antibiotics
    ("amoxicillin",    "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O", "Antibiotic"),
    ("ampicillin",     "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O",   "Antibiotic"),
    ("ciprofloxacin",  "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",         "Antibiotic"),
    ("levofloxacin",   "CC1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23",   "Antibiotic"),
    ("metronidazole",  "Cc1ncc([N+](=O)[O-])n1CCO",                          "Antibiotic"),
    ("doxycycline",    "CC(O)C1=C(O)C2=C(C(=O)c3c(O)cccc3C2(C)O)C(=O)C1C(=O)C(N)=O", "Antibiotic"),
    ("trimethoprim",   "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC",                    "Antibiotic"),
    ("sulfamethoxazole","Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1",                    "Antibiotic"),
    ("azithromycin",   "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@](C)(O)C[C@@H](C)CN(C)[C@H](C)[C@@H]1O", "Antibiotic"),

    # Antifungals
    ("fluconazole",    "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F",                  "Antifungal"),
    ("ketoconazole",   "CC(=O)N1CCN(c2ccc(OCC3COC(Cn4ccnc4)(c4ccc(Cl)cc4Cl)O3)cc2)CC1", "Antifungal"),
    ("itraconazole",   "CCC(C)N1N=CN(c2ccc(N3CCN(c4ccc(OCC5COC(Cn6ccnc6)(c6ccc(Cl)cc6Cl)O5)cc4)CC3)cc2)C1=O", "Antifungal"),

    # Antivirals
    ("acyclovir",      "Nc1nc2c(ncn2COCCO)c(=O)[nH]1",                       "Antiviral"),
    ("oseltamivir",    "CCOC(=O)C1=C[C@H](OC(CC)CC)[C@@H](NC(C)=O)[C@H](N)C1", "Antiviral"),
    ("ribavirin",      "Nc1ncn(C2OC(CO)C(O)C2O)n1",                          "Antiviral"),

    # Statins
    ("atorvastatin",   "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O", "Statin"),
    ("simvastatin",    "CCC(C)(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12", "Statin"),
    ("rosuvastatin",   "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1/C=C/C(O)CC(O)CC(=O)O", "Statin"),
    ("lovastatin",     "CCC(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12", "Statin"),
    ("pravastatin",    "CCC(C)C(=O)OC1CC(O)C=C2C=CC(C)C(CCC(O)CC(O)CC(=O)O)C12", "Statin"),

    # ACE inhibitors / ARBs / antihypertensives
    ("lisinopril",     "NCCCC[C@@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O", "ACE inhibitor"),
    ("captopril",      "C[C@@H](C(=O)N1CCC[C@H]1C(=O)O)CS",                  "ACE inhibitor"),
    ("enalapril",      "CCOC(=O)[C@@H](CCc1ccccc1)N[C@@H](C)C(=O)N1CCC[C@H]1C(=O)O", "ACE inhibitor"),
    ("losartan",       "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1", "ARB"),
    ("valsartan",      "CCCCC(=O)N(Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1)C(C(C)C)C(=O)O", "ARB"),
    ("amlodipine",     "CCOC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccccc1Cl",        "Calcium blocker"),

    # Beta blockers
    ("propranolol",    "CC(C)NCC(O)COc1cccc2ccccc12",                        "Beta blocker"),
    ("atenolol",       "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",                     "Beta blocker"),
    ("metoprolol",     "COCCc1ccc(OCC(O)CNC(C)C)cc1",                        "Beta blocker"),
    ("carvedilol",     "COc1ccccc1OCCNCC(O)COc1cccc2[nH]c3ccccc3c12",        "Beta blocker"),
    ("bisoprolol",     "CC(C)NCC(O)COc1ccc(COCCOC(C)C)cc1",                  "Beta blocker"),

    # Diabetes
    ("metformin",      "CN(C)C(=N)NC(N)=N",                                  "Antidiabetic"),
    ("glipizide",      "Cc1ncc(C(=O)NCCc2ccc(S(=O)(=O)NC(=O)NC3CCCCC3)cc2)cn1", "Antidiabetic"),
    ("sitagliptin",    "Nc1nnc(C(F)(F)F)n1CC(=O)N1CCN2CC(N)c3ncccc3CC12",    "Antidiabetic"),
    ("pioglitazone",   "CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1",            "Antidiabetic"),

    # Antihistamines
    ("diphenhydramine","CN(C)CCOC(c1ccccc1)c1ccccc1",                         "Antihistamine"),
    ("loratadine",     "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc23)CC1",       "Antihistamine"),
    ("cetirizine",     "OC(=O)COCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1",        "Antihistamine"),
    ("fexofenadine",   "CC(C)(C(=O)O)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1", "Antihistamine"),

    # SSRIs / antidepressants
    ("fluoxetine",     "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",                  "SSRI"),
    ("sertraline",     "CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccccc21",        "SSRI"),
    ("paroxetine",     "Fc1ccc(C2CCNCC2COc2ccc3OCOc3c2)cc1",                  "SSRI"),
    ("citalopram",     "N#CC1=CC2=C(C=C1)C(C1=CC=C(F)C=C1)(CCCN(C)C)OC2",     "SSRI"),
    ("venlafaxine",    "COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1",                   "SNRI"),
    ("bupropion",      "CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1",                     "Antidepressant"),

    # Anxiolytics / benzodiazepines
    ("diazepam",       "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",                "Benzodiazepine"),
    ("lorazepam",      "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O",                "Benzodiazepine"),
    ("alprazolam",     "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",              "Benzodiazepine"),

    # Antipsychotics
    ("haloperidol",    "OC1(c2ccc(Cl)cc2)CCN(CCCC(=O)c2ccc(F)cc2)CC1",       "Antipsychotic"),
    ("olanzapine",     "CN1CCN(C2=Nc3cc(C)sc3Nc3ccccc32)CC1",                "Antipsychotic"),
    ("risperidone",    "Cc1noc2cc(F)ccc12.O=C1c2ccccc2N(CCN2CCCN(C)C2=O)CC1", "Antipsychotic"),
    ("quetiapine",     "OCCOCCN1CCN(C2=Nc3ccccc3Sc3ccccc32)CC1",             "Antipsychotic"),

    # Stimulants / CNS
    ("caffeine",       "Cn1cnc2c1c(=O)n(C)c(=O)n2C",                         "Stimulant"),
    ("amphetamine",    "CC(N)Cc1ccccc1",                                     "Stimulant"),
    ("methylphenidate","COC(=O)C(c1ccccc1)C1CCCCN1",                         "Stimulant"),

    # PPIs
    ("omeprazole",     "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",         "PPI"),
    ("pantoprazole",   "COc1ccnc(CS(=O)c2nc3cc(OC(F)F)ccc3[nH]2)c1OC",       "PPI"),

    # Anticoagulants
    ("warfarin",       "CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O",             "Anticoagulant"),
    ("apixaban",       "COc1ccc(-n2nc(C(N)=O)c3CCN(c4ccc(N5CCCCC5=O)cc4)C(=O)c23)cc1", "Anticoagulant"),

    # Bronchodilators
    ("albuterol",      "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",                      "Bronchodilator"),
    ("theophylline",   "Cn1c(=O)c2[nH]cnc2n(C)c1=O",                         "Bronchodilator"),

    # Hormones / steroids
    ("prednisone",     "O=C1C=CC2(C)C(CCC3(C)C2C(=O)CC3(O)C(=O)CO)C1",       "Corticosteroid"),
    ("hydrocortisone", "O=C(CO)C1(O)CCC2C3CCC4=CC(=O)CCC4(C)C3C(O)CC21C",    "Corticosteroid"),

    # Misc
    ("sildenafil",     "CCCc1nn(C)c2c1nc([nH]2)C(=O)c1ccc(OCC)c(S(=O)(=O)N2CCN(C)CC2)c1", "PDE5 inhibitor"),
    ("ranitidine",     "CNC(=C[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1",            "H2 blocker"),
    ("levothyroxine",  "Nc1c(I)cc(Oc2cc(I)c(O)c(I)c2)cc1I.OC(=O)CN",         "Thyroid hormone"),
]
# fmt: on
