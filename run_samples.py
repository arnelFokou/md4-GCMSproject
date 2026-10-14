"""Calcule les concentrations des echantillons BLPC avec la calibration du batch 20251103.

Usage (depuis le dossier md4-GCMSproject) :
    python run_samples.py

Regle pour les echantillons dilues (d20) : on garde la mesure normale si elle est
dans la gamme de calibration (<= 5 ppm). Sinon on prend la mesure d20, dont le
multiplier contient deja le facteur de dilution x20.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline import analyser, executer, inverser, multiplicateur

BATCH = Path("outputs/20251103")
META = "metadata.xlsx"
CONCENTRATIONS = [0.025, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]  # GAM-1 a GAM-8, en ppm
CIBLE_SF = 1.0
LIMITE_HAUTE = 5.0  # plus haut point de la gamme, en ppm

ref, detail, points, calib, sf = executer(BATCH, META, CONCENTRATIONS, CIBLE_SF)
coefs = calib.set_index("compose")[["a", "b", "c"]]
istd = calib.set_index("compose")["istd"].to_dict()

# 1. Mesure brute de chaque injection d'echantillon
lignes = []
for d in sorted(p for p in BATCH.iterdir() if p.is_dir() and "BLPC" in p.name):
    m = analyser(d, ref).set_index("compose")
    mult = multiplicateur(d)
    for c in coefs.index:
        a_i = m["aire_quant"].get(istd.get(c), np.nan)
        a_c = m.loc[c, "aire_quant"]
        ratio = a_c / a_i if a_i > 0 else np.nan
        conc = inverser(coefs.loc[c].to_numpy(), ratio) * mult if ratio == ratio else np.nan
        lignes.append({
            "echantillon": d.name,
            "base": d.name.replace(" d20", ""),
            "dilue_d20": " d20" in d.name,
            "compose": c,
            "conc_ppm": conc,
            "statut": m.loc[c, "statut"],
        })
brut = pd.DataFrame(lignes)

# 2. Choix de la mesure finale pour chaque echantillon et chaque composé
finals = []
for (base, c), g in brut.groupby(["base", "compose"], sort=False):
    normal = g[~g["dilue_d20"]]
    dilue = g[g["dilue_d20"]]
    conc_n = normal["conc_ppm"].iloc[0] if len(normal) else np.nan
    conc_d = dilue["conc_ppm"].iloc[0] if len(dilue) else np.nan
    statut = normal["statut"].iloc[0] if len(normal) else ""
    if conc_n == conc_n and conc_n <= LIMITE_HAUTE:
        valeur, source = conc_n, "normale"
    elif conc_d == conc_d:
        valeur, source = conc_d, "d20"
    elif conc_n == conc_n:
        valeur, source = conc_n, "au-dessus de la gamme"
    else:
        valeur, source = np.nan, "non calculable"
    finals.append({"echantillon": base, "compose": c, "conc_finale_ppm": valeur,
                   "source": source, "statut_pic": statut})
finals = pd.DataFrame(finals)

Path("outputs").mkdir(exist_ok=True)
brut.to_csv("outputs/mesures_brutes_echantillons.csv", index=False)
finals.to_csv("outputs/concentrations_echantillons.csv", index=False)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
tableau = finals.pivot(index="echantillon", columns="compose", values="conc_finale_ppm").round(3)
print("\n=== CONCENTRATIONS (ppm) PAR ECHANTILLON ===")
print(tableau.to_string())

print("\n=== ORIGINE DES VALEURS ===")
print(finals["source"].value_counts().to_string())

print("\n=== PICS NON CONFIRMES (douteux ou non detectes) ===")
print(finals[finals["statut_pic"] != "confirmé"][["echantillon", "compose", "statut_pic"]].to_string(index=False))

print("\nFichiers ecrits : outputs/concentrations_echantillons.csv et outputs/mesures_brutes_echantillons.csv")
