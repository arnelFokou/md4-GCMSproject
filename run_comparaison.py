"""Compare nos concentrations avec celles de MassHunter (Agilent) pour le batch 20251103.

Critères du sujet :
- écart relatif maximal : 10 %
- concentration MassHunter minimale pour qu'un échantillon soit positif : 0,5 ppm

Prérequis : avoir lancé run_final.py (qui crée outputs/mesures_brutes_finales.csv).

Usage (depuis le dossier md4-GCMSproject) :
    python run_comparaison.py
"""
import csv
from pathlib import Path

import numpy as np
import pandas as pd

QUANT = Path("data/20251103/QuantReports")
NOUS = Path("outputs/mesures_brutes_finales.csv")
TOL_ECART = 0.10
SEUIL_POSITIF = 0.5


def lire_masshunter(fichier):
    """Renvoie {composé: concentration en ppm} lu dans la section [Results Signal 1]."""
    lignes = Path(fichier).read_text(encoding="utf-8", errors="ignore").splitlines()
    if "[Results Signal 1]" not in lignes:
        return {}
    debut = lignes.index("[Results Signal 1]") + 1
    resultats = {}
    for ligne in csv.reader(lignes[debut:]):
        if len(ligne) < 7 or not ligne[0].strip():
            continue
        try:
            conc = float(ligne[3])
        except ValueError:
            continue
        if ligne[6].strip() != "T":  # on ne garde que les composés cibles
            continue
        resultats[ligne[0].strip().lower()] = conc
    return resultats


# 1. Lire tous les exports MassHunter (un fichier par échantillon)
# Le dossier "20251103" est le rapport complet ; on le privilégie.
mh = {}
for f in sorted(QUANT.glob("*/Env_Results_*.csv"), key=lambda p: len(p.parent.name)):
    stem = f.name[len("Env_Results_"):-len(".csv")]
    if stem in mh:
        continue
    mh[stem] = lire_masshunter(f)

# 2. Nos résultats
nous = pd.read_csv(NOUS)
nous["compose"] = nous["compose"].str.lower()
nous["stem"] = nous.apply(
    lambda r: r["base"] + (" d20" if r["dilue"] else ""), axis=1
)

# 3. Comparaison composé par composé, échantillon par échantillon
lignes = []
for _, r in nous.iterrows():
    if r["stem"] not in mh:
        continue
    cible = mh[r["stem"]].get(r["compose"])
    if cible is None:
        continue
    ours = r["conc_ppm"]
    if ours != ours or cible == 0:
        continue
    ecart = abs(ours - cible) / abs(cible)
    lignes.append({
        "echantillon": r["stem"],
        "compose": r["compose"],
        "nous_ppm": round(ours, 4),
        "masshunter_ppm": round(cible, 4),
        "ecart_relatif_%": round(ecart * 100, 1),
        "dans_10_%": ecart <= TOL_ECART,
        "positif_nous": ours >= SEUIL_POSITIF,
        "positif_masshunter": cible >= SEUIL_POSITIF,
    })

comp = pd.DataFrame(lignes)
Path("outputs").mkdir(exist_ok=True)
comp.to_csv("outputs/comparaison_masshunter.csv", index=False)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)

print(f"\nEchantillons MassHunter lus : {len(mh)}")
print(f"Couples (échantillon, composé) comparés : {len(comp)}")

if comp.empty:
    print("Aucune comparaison possible : vérifie les noms des fichiers et que run_final.py a bien tourné.")
else:
    print("\n=== RESUME PAR COMPOSE ===")
    resume = comp.groupby("compose").agg(
        n=("dans_10_%", "size"),
        dans_10_pct=("dans_10_%", "sum"),
        ecart_median_pct=("ecart_relatif_%", "median"),
    )
    print(resume.to_string())

    print("\n=== ACCORD POSITIF / NEGATIF (seuil 0,5 ppm) ===")
    accord = np.where(comp["positif_nous"] == comp["positif_masshunter"], "accord", "desaccord")
    print(pd.Series(accord).value_counts().to_string())

    print("\n=== TAUX GLOBAL D'ACCORD A 10 % ===")
    print(f"{comp['dans_10_%'].sum()} / {len(comp)} couples ({comp['dans_10_%'].mean() * 100:.0f} %)")

print("\nFichier écrit : outputs/comparaison_masshunter.csv")
