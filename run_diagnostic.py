"""Diagnostic : pourquoi les pics des echantillons sont-ils 'douteux' ?

Pour chaque composé, compare le temps de rétention de référence (GAM-6) aux
temps mesurés sur les échantillons BLPC (sans les d20).

Usage (depuis le dossier md4-GCMSproject) :
    python run_diagnostic.py
"""
from pathlib import Path

import pandas as pd

from pipeline import analyser, executer

BATCH = Path("outputs/20251103")
META = "metadata.xlsx"
CONCENTRATIONS = [0.025, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]

ref, detail, points, calib, sf = executer(BATCH, META, CONCENTRATIONS, 1.0)

morceaux = []
for d in sorted(p for p in BATCH.iterdir() if p.is_dir() and "BLPC" in p.name and " d20" not in p.name):
    m = analyser(d, ref)
    m["echantillon"] = d.name
    morceaux.append(m)
x = pd.concat(morceaux, ignore_index=True)
x = x[x["type"].astype(str).str.lower() == "target"]

rt_ref = ref.set_index("compose")["rt_ref"]
g = x.groupby("compose")
res = pd.DataFrame({
    "rt_reference_min": rt_ref,
    "rt_median_echantillons_min": g["rt"].median(),
    "ecart_rt_median_%": g["ecart_rt"].median() * 100,
    "ecart_rt_max_%": g["ecart_rt"].max() * 100,
    "ecart_ratio_median_%": g["ecart_ratio"].median() * 100,
    "nb_confirme": g["statut"].apply(lambda s: int((s == "confirmé").sum())),
    "nb_douteux": g["statut"].apply(lambda s: int((s == "douteux").sum())),
    "nb_non_detecte": g["statut"].apply(lambda s: int((s == "non détecté").sum())),
}).dropna(subset=["rt_median_echantillons_min"])

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)
print("\n=== DIAGNOSTIC PAR COMPOSE (tolerance RT = 0,2 %, ratio = 20 %) ===")
print(res.round(3).to_string())
