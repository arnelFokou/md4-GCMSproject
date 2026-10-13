"""Resultats finaux du batch 20251103 avec les criteres du sujet.

Criteres (issus du sujet) :
- R2 de calibration minimal : 0,99 (sinon le composé n'est pas valide)
- Contrôle SF : |concentration - 1 ppm| <= 0,23 ppm (tolérance absolue)
- Ratio ionique qualifiant/quantifiant : écart <= 20 %
- Temps de rétention : écart <= 0,2 %
- Gamme de calibration : de 0,025 à 5 ppm

Usage (depuis le dossier md4-GCMSproject) :
    python run_final.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline import analyser, executer, inverser, multiplicateur

BATCH = Path("outputs/20251103")
META = "metadata.xlsx"
CONCENTRATIONS = [0.025, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]  # GAM-1 a GAM-8, en ppm
CIBLE_SF = 1.0

R2_MIN = 0.99
TOL_RATIO = 0.20
TOL_SF_PPM = 0.23
TOL_RT = 0.002
GAMME_BASSE = 0.025
GAMME_HAUTE = 5.0

ref, detail, points, calib, sf = executer(BATCH, META, CONCENTRATIONS, CIBLE_SF)
coefs = calib.set_index("compose")[["a", "b", "c"]]
istd = calib.set_index("compose")["istd"].to_dict()
r2 = calib.set_index("compose")["r2"]
calib_valide = (r2 >= R2_MIN).to_dict()

# 1. Contrôle SF avec la tolerance absolue du sujet
sf = sf.copy()
sf["valide_sujet"] = (sf["conc_ppm"] - CIBLE_SF).abs() <= TOL_SF_PPM

# 2. Mesures brutes des echantillons BLPC
lignes = []
for d in sorted(p for p in BATCH.iterdir() if p.is_dir() and "BLPC" in p.name):
    m = analyser(d, ref).set_index("compose")
    mult = multiplicateur(d)
    for c in coefs.index:
        a_i = m["aire_quant"].get(istd.get(c), np.nan)
        a_c = m.loc[c, "aire_quant"]
        ratio = a_c / a_i if a_i > 0 else np.nan
        conc = inverser(coefs.loc[c].to_numpy(), ratio) * mult if ratio == ratio else np.nan
        e_rt = m.loc[c, "ecart_rt"]
        e_ratio = m.loc[c, "ecart_ratio"]
        rt_ok = bool(e_rt <= TOL_RT)
        ratio_ok = bool(e_ratio <= TOL_RATIO)
        if rt_ok and ratio_ok:
            ident = "identifié"
        elif rt_ok:
            ident = "ratio ionique hors tolérance"
        else:
            ident = "pic non trouvé"
        lignes.append({
            "base": d.name.replace(" d20", ""),
            "dilue": " d20" in d.name,
            "compose": c,
            "conc_ppm": conc,
            "identification": ident,
        })
brut = pd.DataFrame(lignes)

# 3. Choix de la mesure finale pour chaque echantillon et chaque composé
finals = []
for (base, c), g in brut.groupby(["base", "compose"], sort=False):
    choisi = None
    for dilue in [False, True]:
        r = g[g["dilue"] == dilue]
        if len(r) and r.iloc[0]["conc_ppm"] == r.iloc[0]["conc_ppm"]:
            conc = r.iloc[0]["conc_ppm"]
            if GAMME_BASSE <= conc <= GAMME_HAUTE:
                choisi = r.iloc[0]
                break
    normale = g[~g["dilue"]].iloc[0]

    if choisi is not None:
        valeur = round(float(choisi["conc_ppm"]), 4)
        texte = f"{valeur}"
        ident = choisi["identification"]
        source = "d20" if choisi["dilue"] else "normale"
    else:
        valeur = np.nan
        ident = normale["identification"]
        source = "aucune"
        if normale["conc_ppm"] == normale["conc_ppm"] and normale["conc_ppm"] < GAMME_BASSE:
            texte = f"< {GAMME_BASSE} ppm (sous la gamme)"
        else:
            texte = f"> {GAMME_HAUTE} ppm (hors gamme)"

    if not calib_valide.get(c, False):
        decision = "non valide : R² de calibration < 0,99"
    elif choisi is not None and ident == "identifié":
        decision = "retenu"
    elif choisi is not None:
        decision = f"indicatif : {ident}"
    else:
        decision = "non quantifiable"

    finals.append({
        "echantillon": base,
        "compose": c,
        "resultat_ppm": valeur if choisi is not None else texte,
        "source": source,
        "identification": ident,
        "decision": decision,
    })
finals = pd.DataFrame(finals)

Path("outputs").mkdir(exist_ok=True)
brut.to_csv("outputs/mesures_brutes_finales.csv", index=False)
finals.to_csv("outputs/resultats_finaux.csv", index=False)
sf.to_csv("outputs/controle_sf_sujet.csv", index=False)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)

print("\n=== CALIBRATION (R² minimal 0,99) ===")
print(calib.assign(conforme=lambda t: t["r2"] >= R2_MIN)[["compose", "r2", "conforme"]].to_string(index=False))

print("\n=== CONTROLE SF (tolérance ±0,23 ppm) ===")
resume_sf = sf.groupby("compose")["valide_sujet"].agg(["sum", "count"])
resume_sf.columns = ["valides", "total"]
print(resume_sf.to_string())

print("\n=== DECISIONS ===")
print(finals["decision"].value_counts().to_string())

print("\n=== RESULTATS RETENUS (ppm), sans les composés non valides ===")
retenus = finals[finals["decision"] == "retenu"]
tableau = retenus.pivot(index="echantillon", columns="compose", values="resultat_ppm")
print(tableau.to_string())

print("\nFichiers écrits : outputs/resultats_finaux.csv, outputs/controle_sf_sujet.csv, outputs/mesures_brutes_finales.csv")
