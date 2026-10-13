"""Lance la calibration et le controle SF sur le batch 20251103.

Usage (depuis le dossier md4-GCMSproject) :
    python run_calibration.py
"""
from pathlib import Path

import pandas as pd

from pipeline import executer

BATCH = Path("outputs/20251103")
META = "metadata.xlsx"
CONCENTRATIONS = [0.025, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]  # GAM-1 a GAM-8, en ppm
CIBLE_SF = 1.0  # concentration attendue du controle SF, en ppm

ref, detail, points, calib, sf = executer(BATCH, META, CONCENTRATIONS, CIBLE_SF)

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 20)

print("\n=== CALIBRATION (R2 proche de 1 = bonne courbe) ===")
print(calib[["compose", "r2", "n_points"]].to_string(index=False))

print("\n=== CONTROLE SF (recuperation proche de 100 % = valide) ===")
print(sf[["echantillon", "compose", "conc_ppm", "recuperation_%", "valide"]].to_string(index=False))

Path("outputs").mkdir(exist_ok=True)
calib.to_csv("outputs/calibration.csv", index=False)
sf.to_csv("outputs/controle_sf.csv", index=False)
print("\nFichiers ecrits : outputs/calibration.csv et outputs/controle_sf.csv")
