"""Diagnostic : sur des solutions de concentration connue, qui retrouve la bonne valeur ?

Compare, pour chaque composé :
- la valeur nominale (connue : GAM-2 = 0,05 ppm, GAM-8 = 5 ppm, SF = 1 ppm)
- ce que trouve MassHunter
- ce que trouve notre programme (contrôle SF seulement)

Si MassHunter retrouve la valeur nominale mais pas nous, le problème vient de notre calcul.
Si MassHunter ne la retrouve pas non plus, le problème vient de la méthode ou de la référence.

Usage (depuis le dossier md4-GCMSproject, après run_final.py) :
    python run_diagnostic_mh.py
"""
import csv
from pathlib import Path

import pandas as pd

QUANT = Path("data/20251103/QuantReports")
SF_NOUS = Path("outputs/controle_sf_sujet.csv")

NOMINAL = {"GAM-25-385-2": 0.05, "GAM-25-385-8": 5.0}
for i in range(1, 5):
    NOMINAL[f"SF-25-1038-{i}"] = 1.0


def lire_masshunter(fichier):
    lignes = Path(fichier).read_text(encoding="utf-8", errors="ignore").splitlines()
    if "[Results Signal 1]" not in lignes:
        return {}
    debut = lignes.index("[Results Signal 1]") + 1
    res = {}
    for ligne in csv.reader(lignes[debut:]):
        if len(ligne) < 7 or not ligne[0].strip():
            continue
        try:
            conc = float(ligne[3])
        except ValueError:
            continue
        if ligne[6].strip() == "T":
            res[ligne[0].strip().lower()] = conc
    return res


# 1. MassHunter sur les solutions connues
mh = {}
for f in sorted(QUANT.glob("*/Env_Results_*.csv"), key=lambda p: len(p.parent.name)):
    stem = f.name[len("Env_Results_"):-len(".csv")]
    if stem.replace("20251103-", "") in NOMINAL and stem not in mh:
        mh[stem] = lire_masshunter(f)

lignes = []
for stem, valeurs in mh.items():
    cle = stem.replace("20251103-", "")
    for compose, conc in valeurs.items():
        lignes.append({"solution": cle, "compose": compose, "nominal": NOMINAL[cle],
                       "masshunter": conc})
mh_df = pd.DataFrame(lignes)

# 2. Notre programme sur les SF
nous = pd.read_csv(SF_NOUS)
nous["compose"] = nous["compose"].str.lower()
nous["solution"] = nous["echantillon"].str.replace("20251103-", "", regex=False)
nous = nous[["solution", "compose", "conc_ppm"]].rename(columns={"conc_ppm": "nous"})

tab = mh_df.merge(nous, on=["solution", "compose"], how="left")
tab["ecart_MH_%"] = ((tab["masshunter"] - tab["nominal"]) / tab["nominal"] * 100).round(1)
tab["ecart_nous_%"] = ((tab["nous"] - tab["nominal"]) / tab["nominal"] * 100).round(1)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)
pd.set_option("display.max_rows", 200)

print("\n=== RESULTATS SUR SOLUTIONS CONNUES ===")
print(tab.round(3).to_string(index=False))

resume = tab.groupby("solution").agg(
    nominal=("nominal", "first"),
    mh_median=("masshunter", "median"),
    ecart_MH_median_pct=("ecart_MH_%", "median"),
    ecart_nous_median_pct=("ecart_nous_%", "median"),
)
print("\n=== RESUME PAR SOLUTION (médiane des composés) ===")
print(resume.round(2).to_string())

tab.to_csv("outputs/diagnostic_solutions_connues.csv", index=False)
print("\nFichier écrit : outputs/diagnostic_solutions_connues.csv")
