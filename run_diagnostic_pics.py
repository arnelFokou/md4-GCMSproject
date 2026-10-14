"""Diagnostic : pour quelques échantillons, compare nos pics (temps de rétention, aire)
avec ceux de MassHunter, composé par composé.

Si les temps de rétention concordent mais pas les aires : problème d'intégration.
Si les temps de rétention ne concordent pas : on ne prend pas le même pic.

Usage (depuis le dossier md4-GCMSproject) :
    python run_diagnostic_pics.py
"""
import csv
from pathlib import Path

import pandas as pd

from pipeline import analyser, executer

BATCH = Path("outputs/20251103")
QUANT = Path("data/20251103/QuantReports/20251103")
META = "metadata.xlsx"
CONCENTRATIONS = [0.025, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]
ECHANTILLONS = ["25BLPC003045-010", "25BLPC003018-001", "25BLPC003045-014"]

ref, detail, points, calib, sf = executer(BATCH, META, CONCENTRATIONS, 1.0)


def lire_masshunter(fichier):
    """Renvoie {composé: (aire, temps de rétention, concentration)} pour les cibles."""
    lignes = Path(fichier).read_text(encoding="utf-8", errors="ignore").splitlines()
    if "[Results Signal 1]" not in lignes:
        return {}
    debut = lignes.index("[Results Signal 1]") + 1
    res = {}
    for ligne in csv.reader(lignes[debut:]):
        if len(ligne) < 7 or not ligne[0].strip() or ligne[6].strip() != "T":
            continue
        try:
            res[ligne[0].strip().lower()] = (float(ligne[1]), float(ligne[2]), float(ligne[3]))
        except ValueError:
            continue
    return res


lignes = []
for s in ECHANTILLONS:
    m = analyser(BATCH / s, ref).set_index("compose")
    mh = lire_masshunter(QUANT / f"Env_Results_{s}.csv")
    for compose in m.index:
        if compose not in [c.lower() for c in mh] and compose not in mh:
            pass
        cle = compose.lower()
        if cle not in mh:
            continue
        aire_mh, rt_mh, conc_mh = mh[cle]
        lignes.append({
            "echantillon": s,
            "compose": compose,
            "rt_nous": round(m.loc[compose, "rt"], 3) if m.loc[compose, "rt"] == m.loc[compose, "rt"] else None,
            "rt_masshunter": round(rt_mh, 3),
            "aire_nous": round(m.loc[compose, "aire_quant"], 1),
            "aire_masshunter": round(aire_mh, 1),
            "conc_masshunter": round(conc_mh, 3),
        })

tab = pd.DataFrame(lignes)
tab["ecart_rt_min"] = (tab["rt_nous"] - tab["rt_masshunter"]).round(3)
tab["aire_nous_sur_mh"] = (tab["aire_nous"] / tab["aire_masshunter"]).round(2)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)
pd.set_option("display.max_rows", 200)
print("\n=== PICS : NOUS vs MASSHUNTER ===")
print(tab.to_string(index=False))
tab.to_csv("outputs/diagnostic_pics.csv", index=False)
print("\nFichier écrit : outputs/diagnostic_pics.csv")
