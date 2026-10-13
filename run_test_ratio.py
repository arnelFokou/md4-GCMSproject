"""Test de stabilite du ratio ionique (qualifiant / quantifiant) sur les gammes GAM-1 a GAM-8.

Question : le ratio change-t-il deja beaucoup sur des solutions propres, sans matrice de sol ?
- Si oui : le probleme vient de la methode de calcul (ou du critere trop strict).
- Si non : le probleme vient de la matrice des echantillons.

Usage (depuis le dossier md4-GCMSproject) :
    python run_test_ratio.py
"""
from pathlib import Path

import pandas as pd

from pipeline import analyser, executer

BATCH = Path("outputs/20251103")
META = "metadata.xlsx"
CONCENTRATIONS = [0.025, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]

ref, detail, points, calib, sf = executer(BATCH, META, CONCENTRATIONS, 1.0)
cibles = calib["compose"].tolist()

gams = sorted(
    (p for p in BATCH.iterdir() if p.is_dir() and "GAM" in p.name),
    key=lambda p: int(p.name.split("-")[-1]),
)
colonnes = [f"GAM{i}" for i in range(1, len(gams) + 1)]

tableau = {}
for nom, d in zip(colonnes, gams):
    m = analyser(d, ref).set_index("compose")
    tableau[nom] = m.loc[cibles, "ratio"]
ratios = pd.DataFrame(tableau)

ref_ratio = ref.set_index("compose").loc[cibles, "ratio_ref"]
ratios.insert(0, "ratio_ref_GAM6", ref_ratio)
ratios["mediane"] = ratios[colonnes].median(axis=1)
ratios["ecart_type_%_sur_moyenne"] = ratios[colonnes].std(axis=1) / ratios[colonnes].mean(axis=1) * 100

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)

print("\n=== RATIO QUALIFIANT / QUANTIFIANT SUR LES GAMMES PROPRES ===")
print(ratios.round(3).to_string())

print("\nLecture : ecart_type_%_sur_moyenne = variabilite du ratio d'une gamme a l'autre.")
print("Plus ce nombre est grand, moins le ratio est stable, meme sur des solutions propres.")

ratios.round(4).to_csv("outputs/test_ratio_gammes.csv")
print("\nFichier ecrit : outputs/test_ratio_gammes.csv")
