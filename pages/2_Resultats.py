from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(page_title="Résultats HAP", page_icon="📋", layout="wide")
st.title("Résultats finaux")
st.caption("Concentrations par échantillon et décision, selon les critères du sujet.")

with st.sidebar:
    dossier = Path(st.text_input("Dossier des résultats", "outputs")).expanduser()

fichier = dossier / "resultats_finaux.csv"
try:
    df = pd.read_csv(fichier)
except (OSError, ValueError) as e:
    st.error(f"Impossible de lire {fichier} : {e}")
    st.info("Lance d'abord run_final.py pour créer ce fichier.")
    st.stop()

df["decision"] = df["decision"].fillna("")
decisions_possibles = sorted(df["decision"].unique())
composes = sorted(df["compose"].unique())
echantillons = sorted(df["echantillon"].unique())

with st.sidebar:
    choix_decision = st.multiselect("Décision", decisions_possibles, default=decisions_possibles)
    choix_compose = st.multiselect("Composés", composes, default=composes)
    choix_echantillon = st.multiselect("Échantillons", echantillons, default=echantillons)

filtre = df[
    df["decision"].isin(choix_decision)
    & df["compose"].isin(choix_compose)
    & df["echantillon"].isin(choix_echantillon)
]

c1, c2, c3, c4 = st.columns(4)
c1.metric("Retenus", int((df["decision"] == "retenu").sum()))
c2.metric("Calibration non valide", int(df["decision"].str.startswith("non valide").sum()))
c3.metric("Indicatifs", int(df["decision"].str.startswith("indicatif").sum()))
c4.metric("Non quantifiables", int((df["decision"] == "non quantifiable").sum()))

onglet1, onglet2 = st.tabs(["Tableau détaillé", "Retenus (par échantillon)"])

with onglet1:
    st.write(f"{len(filtre)} ligne(s) affichée(s) sur {len(df)}")
    st.dataframe(filtre, use_container_width=True, hide_index=True)
    st.download_button("Télécharger ce tableau (CSV)", filtre.to_csv(index=False),
                       "resultats_filtres.csv")

with onglet2:
    retenus = df[df["decision"] == "retenu"]
    if retenus.empty:
        st.info("Aucun résultat retenu.")
    else:
        tableau = retenus.pivot_table(index="echantillon", columns="compose",
                                      values="resultat_ppm", aggfunc="first")
        st.write("Concentrations retenues, en ppm. Les cases vides ne sont pas quantifiables.")
        st.dataframe(tableau, use_container_width=True)
