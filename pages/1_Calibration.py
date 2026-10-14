from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import streamlit as st

import pipeline as pl

st.set_page_config(page_title="Calibration GC/MS", page_icon="🧪", layout="wide")
st.title("Calibration et contrôle des SF")
st.caption("16 HAP cibles, calibrés par rapport à leur étalon interne deutéré.")


def lister_batchs(root: Path):
    if any(root.glob("*GAM*")):
        return [root]
    return sorted(p for p in root.iterdir() if p.is_dir() and any(p.glob("*GAM*")))


@st.cache_data(show_spinner="Calcul en cours...")
def calculer(batch, meta, concs, cible, n_ref):
    return pl.executer(batch, meta, list(concs), cible, n_ref)


with st.sidebar:
    root = Path(st.text_input("Dossier des CSV (batch)", "outputs/20251103")).expanduser()
    meta = st.text_input("metadata.xlsx", "metadata.xlsx")
    concs_txt = st.text_input("Concentrations GAM-1 à GAM-8 (ppm)",
                              "0.025, 0.05, 0.1, 0.2, 0.5, 1, 2, 5")
    cible = st.number_input("Concentration cible des SF (ppm)", value=1.0)
    n_ref = st.number_input("GAM de référence (n°)", 1, 8, 6)
    if st.button("Recalculer", use_container_width=True):
        st.cache_data.clear()

try:
    batchs = lister_batchs(root)
    concs = tuple(float(x) for x in concs_txt.split(","))
except (OSError, ValueError) as e:
    st.error(e)
    st.stop()
if not batchs:
    st.warning("Aucun batch (dossiers GAM) trouvé.")
    st.stop()

batch = st.selectbox("Batch", batchs, format_func=lambda p: p.name)
try:
    ref, detail, points, calib, sf = calculer(str(batch), meta, concs, cible, int(n_ref))
except Exception as e:
    st.error(f"{type(e).__name__} : {e}")
    st.stop()

ok = int(sf.groupby("compose")["valide"].all().sum()) if len(sf) else 0
c1, c2, c3 = st.columns(3)
c1.metric("Composés calibrés", f"{len(calib)} / 16")
c2.metric("Composés validés par les SF", f"{ok} / {len(calib)}")
c3.metric("R² minimum", f"{calib.r2.min():.4f}" if len(calib) else "-")

t1, t2, t3, t4 = st.tabs(["Référence GAM", "Calibration", "Contrôle SF", "Voir un pic"])

with t1:
    st.dataframe(ref, use_container_width=True, hide_index=True)

with t2:
    nom = st.selectbox("Composé", calib.compose)
    p = points[points.compose == nom].dropna()
    a, b, c = calib.set_index("compose").loc[nom, ["a", "b", "c"]]
    xs = np.linspace(0, p.conc.max() * 1.05, 200)
    fig = go.Figure()
    fig.add_scatter(x=p.conc, y=p.ratio, mode="markers", name="GAM")
    fig.add_scatter(x=xs, y=a * xs**2 + b * xs + c, mode="lines", name="ratio = aC² + bC + c")
    s = sf[sf.compose == nom]
    fig.add_scatter(x=[cible] * len(s), y=s.ratio, mode="markers", name="SF",
                    marker={"symbol": "x", "size": 11})
    fig.update_layout(xaxis_title="Concentration (ppm)", yaxis_title="Aire composé / aire ISTD")
    st.plotly_chart(fig, use_container_width=True)
    st.dataframe(calib, use_container_width=True, hide_index=True)

with t3:
    st.dataframe(sf, use_container_width=True, hide_index=True)
    st.download_button("Télécharger le contrôle SF (CSV)", sf.to_csv(index=False), "controle_sf.csv")

with t4:
    cc1, cc2 = st.columns(2)
    ech = cc1.selectbox("Échantillon", sorted(detail.echantillon.unique()))
    cmp_ = cc2.selectbox("Composé ", ref.compose, key="pic")
    r = detail[(detail.echantillon == ech) & (detail.compose == cmp_)].iloc[0]
    st.write(f"Statut : **{r.statut}**")
    t, y = pl.lire_signal(batch / ech, r.mz_quant)
    tq, yq = pl.lire_signal(batch / ech, r.mz_qual)
    rt_ref = float(ref.set_index("compose").loc[cmp_, "rt_ref"])
    rt = r.rt if r.rt == r.rt else rt_ref
    st.write(f"Temps de rétention de référence : **{rt_ref:.3f} min**"
             + (f" | temps trouvé dans l'échantillon : **{r.rt:.3f} min**" if r.rt == r.rt else ""))
    zone = np.abs(t - rt) < 0.4
    fig = go.Figure()
    fig.add_scatter(x=t[zone], y=y[zone], name=f"quantifiant {r.mz_quant}")
    fig.add_scatter(x=tq[zone], y=yq[zone], name=f"qualifiant {r.mz_qual}")
    # Zone où le programme cherche le sommet du pic (±0,05 min autour de la référence)
    fig.add_vrect(x0=rt_ref - pl.DEMI_FENETRE, x1=rt_ref + pl.DEMI_FENETRE,
                  fillcolor="lightblue", opacity=0.25, line_width=0,
                  annotation_text="zone de recherche")
    # Tolérance d'acceptation : ±0,2 % du temps de référence
    fig.add_vrect(x0=rt_ref * (1 - pl.TOL_RT), x1=rt_ref * (1 + pl.TOL_RT),
                  fillcolor="orange", opacity=0.5, line_width=0,
                  annotation_text="tolérance ±0,2 %")
    fig.add_vline(x=rt_ref, line_dash="dot", line_color="black")
    if r.debut == r.debut:
        fig.add_vrect(x0=r.debut, x1=r.fin, fillcolor="grey", opacity=0.2, line_width=0,
                      annotation_text="intégration")
    fig.update_layout(xaxis_title="Temps (min)", yaxis_title="Intensité")
    st.plotly_chart(fig, use_container_width=True)