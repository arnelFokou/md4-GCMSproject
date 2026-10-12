"""Interface Streamlit de visualisation des CSV GC/MS."""

from __future__ import annotations

from pathlib import Path

import plotly.graph_objects as go
import streamlit as st

from visualiser_csv import CsvCatalog, Sample, area_under_curve


st.set_page_config(
    page_title="Visualiseur GC/MS",
    page_icon="📈",
    layout="wide",
)


@st.cache_data(show_spinner=False)
def load_channel(root: str, sample_name: str, channel: str):
    """Charge un canal et met en cache le résultat tant que les choix sont identiques."""
    catalog = CsvCatalog(Path(root))
    catalog.scan()
    sample = next(sample for sample in catalog.samples if sample.name == sample_name)
    return catalog.read_channel(sample, channel)


def scan_catalog(root: Path) -> CsvCatalog:
    catalog = CsvCatalog(root)
    catalog.scan()
    return catalog


st.title("Visualiseur des intensités GC/MS")
st.caption("Comparez plusieurs fichiers CSV `mz_<canal>.csv` sur un même graphique.")

with st.sidebar:
    st.header("Filtres")
    root_text = st.text_input("Dossier cache", value="cache")
    root = Path(root_text).expanduser()
    scan_requested = st.button("Actualiser le cache", use_container_width=True)

    if scan_requested:
        st.cache_data.clear()

    try:
        catalog = scan_catalog(root)
    except (FileNotFoundError, OSError) as error:
        st.error(str(error))
        st.stop()

    sample_type = st.selectbox("Type d'échantillon", ("Tous", "GAM", "SF"))
    channels = catalog.channels_for(sample_type)
    if not channels:
        st.warning("Aucun canal trouvé pour ce type d'échantillon.")
        st.stop()
    channel = st.selectbox("Canal (m/z)", channels)

    available_samples = catalog.samples_for(sample_type, channel)
    sample_names = [sample.name for sample in available_samples]
    sample_types = {sample.name: sample.sample_type for sample in available_samples}
    selected_names = st.multiselect(
        "Échantillons à afficher",
        options=sample_names,
        default=sample_names[:1],
        format_func=lambda name: f"{sample_types[name]} | {name}",
        help="Vous pouvez comparer plusieurs échantillons sur le même graphique.",
    )
    logarithmic = st.checkbox("Échelle logarithmique", value=False)
    show_legend = st.checkbox("Afficher la légende", value=True)

st.info(
    f"{len(catalog.samples)} échantillon(s) détecté(s) | "
    f"{len(channels)} canal(aux) disponible(s) pour « {sample_type} »"
)

selected_samples: list[Sample] = [
    sample for sample in available_samples if sample.name in selected_names
]
if not selected_samples:
    st.warning("Sélectionnez au moins un échantillon dans la barre latérale.")
    st.stop()

figure = go.Figure()
failures: list[str] = []
area_rows: list[dict[str, float | int | str]] = []
plotted = 0
for sample in selected_samples:
    try:
        data = load_channel(str(root.resolve()), sample.name, channel)
    except (OSError, ValueError, StopIteration) as error:
        failures.append(f"{sample.name}: {error}")
        continue
    figure.add_trace(
        go.Scatter(
            x=data["time_min"],
            y=data["intensity"],
            mode="lines",
            name=sample.name,
            line={"width": 1},
        )
    )
    area_rows.append(
        {
            "Échantillon": sample.name,
            "Aire (intensité·min)": area_under_curve(data),
            "Points": len(data),
        }
    )
    plotted += 1

figure.update_layout(
    title=f"Canal m/z {channel} — {plotted} courbe(s)",
    xaxis_title="Temps (min)",
    yaxis_title="Intensité",
    yaxis_type="log" if logarithmic else "linear",
    showlegend=show_legend,
    hovermode="x unified",
    height=650,
    margin={"l": 50, "r": 30, "t": 70, "b": 50},
)
figure.update_xaxes(showgrid=True, gridcolor="rgba(128,128,128,0.25)")
figure.update_yaxes(showgrid=True, gridcolor="rgba(128,128,128,0.25)")
visual_col, area_col = st.columns([3, 1], gap="medium")
with visual_col:
    st.plotly_chart(
        figure,
        use_container_width=True,
        config={
            "scrollZoom": True,
            "displaylogo": False,
            "modeBarButtonsToAdd": ["resetScale2d"],
        },
    )

with area_col:
    st.subheader("Aire")
    if area_rows:
        st.dataframe(
            area_rows,
            use_container_width=True,
            hide_index=True,
            column_config={
                "Aire (intensité·min)": st.column_config.NumberColumn(
                    format="%.6g"
                ),
                "Points": st.column_config.NumberColumn(format="%d"),
            },
        )
    else:
        st.info("Aucune courbe disponible.")

if failures:
    st.warning("Certains fichiers n'ont pas pu être lus :\n\n" + "\n".join(failures))
