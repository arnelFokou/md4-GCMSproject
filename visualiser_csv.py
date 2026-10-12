#!/usr/bin/env python3
"""Application graphique pour visualiser les CSV exportés par ``extraire_csv.py``.

L'application parcourt un dossier cache, détecte les échantillons GAM/SF et
leurs canaux ``mz_<canal>.csv``, puis superpose les intensités en fonction du
temps de rétention.

Lancement:
    python visualiser_csv.py
    python visualiser_csv.py --root cache
"""

from __future__ import annotations

import argparse
import re
from scipy.signal import find_peaks
import tkinter as tk
from dataclasses import dataclass
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Final

import matplotlib
import numpy as np

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import pandas as pd


CSV_COLUMNS: Final = ("time_min", "intensity")
SAMPLE_TYPES: Final = ("Tous", "GAM", "SF")
CHANNEL_PATTERN: Final = re.compile(r"^mz_(\d+(?:\.\d+)?)\.csv$", re.IGNORECASE)


@dataclass(frozen=True)
class Sample:
    """Un échantillon détecté dans le cache."""

    name: str
    sample_type: str
    directory: Path


class CsvCatalog:
    """Indexe les échantillons et canaux disponibles dans un cache."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.samples: list[Sample] = []
        self.channels: dict[str, set[str]] = {}

    def scan(self) -> None:
        if not self.root.is_dir():
            raise FileNotFoundError(f"Le dossier n'existe pas : {self.root}")

        samples: list[Sample] = []
        channels: dict[str, set[str]] = {}
        for directory in sorted(path for path in self.root.rglob("*") if path.is_dir()):
            upper_name = directory.name.upper()
            sample_type = next(
                (value for value in ("GAM", "SF") if value in upper_name), None
            )
            if sample_type is None:
                continue

            channel_names = {
                match.group(1)
                for csv_path in directory.glob("*.csv")
                if (match := CHANNEL_PATTERN.match(csv_path.name)) is not None
            }
            if not channel_names:
                continue
            samples.append(Sample(directory.name, sample_type, directory))
            channels[directory.name] = channel_names

        self.samples = samples
        self.channels = channels

    def channels_for(self, sample_type: str) -> list[str]:
        selected = (
            self.samples
            if sample_type == "Tous"
            else [sample for sample in self.samples if sample.sample_type == sample_type]
        )
        return sorted(
            set().union(*(self.channels[sample.name] for sample in selected)),
            key=lambda value: float(value),
        )

    def samples_for(self, sample_type: str, channel: str) -> list[Sample]:
        return [
            sample
            for sample in self.samples
            if (sample_type == "Tous" or sample.sample_type == sample_type)
            and channel in self.channels.get(sample.name, set())
        ]

    def read_channel(self, sample: Sample, channel: str) -> pd.DataFrame:
        path = sample.directory / f"mz_{channel}.csv"
        data = pd.read_csv(path)
        missing = [column for column in CSV_COLUMNS if column not in data.columns]
        if missing:
            raise ValueError(
                f"{path.name} ne contient pas les colonnes attendues : "
                f"{', '.join(missing)}"
            )
        result = data.loc[:, CSV_COLUMNS].copy()
        result["time_min"] = pd.to_numeric(result["time_min"], errors="coerce")
        result["intensity"] = pd.to_numeric(result["intensity"], errors="coerce")
        return result.dropna(subset=list(CSV_COLUMNS)).sort_values("time_min")


def area_under_curve(data: pd.DataFrame) -> float:
    """Calcule l'aire du pic principal d'une série intensité-temps, par trapèzes.
 
    La fonction repère le pic le plus haut, cherche son début et sa fin,
    retire le niveau de base, puis intègre.
    Résultat en minutes x intensité, comme votre fonction d'origine
    (multiplier par 60 pour obtenir des secondes).
    Renvoie 0.0 si aucun pic n'est trouvé.
    """
    if len(data) < 2:
        return 0.0
    t = data["time_min"].to_numpy(dtype=float)
    y = data["intensity"].to_numpy(dtype=float)

    # 1. Repérer le pic le plus haut sur la courbe lissée
    lisse = np.convolve(y, np.ones(3) / 3, mode="same")
    pics, _ = find_peaks(lisse, prominence=0.03 * lisse.max())
    if len(pics) == 0:
        return 0.0
    sommet = pics[np.argmax(lisse[pics])]

    # 2. Prendre le vrai maximum sur la courbe brute, tout près
    proche = np.where(np.abs(t - t[sommet]) <= 0.03)[0]
    sommet = proche[np.argmax(y[proche])]

    # 3. Niveau de base : mesuré juste avant le pic
    avant = (t > t[sommet] - 0.6) & (t < t[sommet] - 0.25)
    base = float(np.median(y[avant])) if avant.any() else float(np.min(y))
    bruit = max(float(np.std(y[avant])) if avant.any() else 1.0, 1.0)
    seuil = base + 0.02 * (y[sommet] - base) + 2 * bruit

    # 4. Début et fin : on descend jusqu'au retour près de la base,
    #    ou jusqu'à ce que le signal remonte (un autre pic commence)
    def descendre(sens):
        i, meilleur = sommet, sommet
        while 0 <= i + sens < len(y):
            i += sens
            if y[i] < y[meilleur]:
                meilleur = i
            if y[i] <= seuil:
                return i
            if y[i] > y[meilleur] + 3 * bruit:
                return meilleur
        return i

    debut, fin = descendre(-1), descendre(+1)

    # 5. Aire par trapèzes, après retrait de la base
    return float(np.trapezoid(y[debut:fin + 1] - base, t[debut:fin + 1]))


class Viewer(tk.Tk):
    """Fenêtre principale de l'application."""

    def __init__(self, initial_root: Path) -> None:
        super().__init__()
        self.title("Visualiseur des intensités GC/MS")
        self.geometry("1280x780")
        self.minsize(950, 600)
        self.catalog = CsvCatalog(initial_root)
        self._build_ui()
        self._scan()

    def _build_ui(self) -> None:
        controls = ttk.Frame(self, padding=10)
        controls.pack(fill=tk.X)
        controls.columnconfigure(1, weight=1)

        ttk.Label(controls, text="Dossier cache :").grid(row=0, column=0, sticky="w")
        self.root_var = tk.StringVar(value=str(self.catalog.root))
        ttk.Entry(controls, textvariable=self.root_var).grid(
            row=0, column=1, sticky="ew", padx=6
        )
        ttk.Button(controls, text="Parcourir...", command=self._choose_root).grid(
            row=0, column=2
        )
        ttk.Button(controls, text="Actualiser", command=self._scan).grid(
            row=0, column=3, padx=(6, 0)
        )

        filters = ttk.Frame(self, padding=(10, 0, 10, 8))
        filters.pack(fill=tk.X)
        ttk.Label(filters, text="Type :").pack(side=tk.LEFT)
        self.type_var = tk.StringVar(value="Tous")
        self.type_combo = ttk.Combobox(
            filters, textvariable=self.type_var, values=SAMPLE_TYPES, state="readonly", width=9
        )
        self.type_combo.pack(side=tk.LEFT, padx=(5, 16))
        self.type_combo.bind("<<ComboboxSelected>>", self._on_type_changed)

        ttk.Label(filters, text="Canal :").pack(side=tk.LEFT)
        self.channel_var = tk.StringVar()
        self.channel_combo = ttk.Combobox(
            filters, textvariable=self.channel_var, state="readonly", width=9
        )
        self.channel_combo.pack(side=tk.LEFT, padx=(5, 16))
        self.channel_combo.bind("<<ComboboxSelected>>", self._on_channel_changed)

        self.log_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(filters, text="Échelle logarithmique", variable=self.log_var).pack(
            side=tk.LEFT
        )
        self.legend_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(filters, text="Légende", variable=self.legend_var).pack(
            side=tk.LEFT, padx=12
        )
        ttk.Button(filters, text="Tracer la sélection", command=self._plot).pack(side=tk.RIGHT)

        body = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        body.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0, 10))
        selection = ttk.Frame(body, padding=6)
        body.add(selection, weight=1)
        chart = ttk.Frame(body, padding=6)
        body.add(chart, weight=4)

        ttk.Label(selection, text="Échantillons (Ctrl/Clic pour plusieurs) :").pack(
            anchor="w"
        )
        list_frame = ttk.Frame(selection)
        list_frame.pack(fill=tk.BOTH, expand=True, pady=6)
        self.sample_list = tk.Listbox(
            list_frame, selectmode=tk.EXTENDED, exportselection=False
        )
        self.sample_list.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar = ttk.Scrollbar(list_frame, orient=tk.VERTICAL, command=self.sample_list.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.sample_list.configure(yscrollcommand=scrollbar.set)
        buttons = ttk.Frame(selection)
        buttons.pack(fill=tk.X)
        ttk.Button(buttons, text="Tout sélectionner", command=self._select_all).pack(
            side=tk.LEFT, expand=True, fill=tk.X
        )
        ttk.Button(buttons, text="Effacer", command=self._clear_selection).pack(
            side=tk.LEFT, expand=True, fill=tk.X, padx=(5, 0)
        )

        self.figure = Figure(figsize=(8, 5), dpi=100)
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=chart)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.status_var = tk.StringVar()
        ttk.Label(self, textvariable=self.status_var, relief=tk.SUNKEN, anchor="w").pack(
            fill=tk.X, side=tk.BOTTOM
        )

    def _choose_root(self) -> None:
        selected = filedialog.askdirectory(initialdir=self.root_var.get())
        if selected:
            self.root_var.set(selected)
            self._scan()

    def _on_type_changed(self, event: tk.Event[tk.Widget]) -> None:
        self.type_var.set(event.widget.get())
        self._update_filters()

    def _on_channel_changed(self, event: tk.Event[tk.Widget]) -> None:
        self.channel_var.set(event.widget.get())
        self._update_samples()

    def _scan(self) -> None:
        self.catalog.root = Path(self.root_var.get()).expanduser()
        try:
            self.catalog.scan()
        except (FileNotFoundError, OSError) as error:
            self.status_var.set(str(error))
            messagebox.showerror("Dossier cache invalide", str(error))
            return
        self._update_filters()
        self.status_var.set(
            f"{len(self.catalog.samples)} échantillon(s) trouvé(s) dans {self.catalog.root}"
        )

    def _update_filters(self) -> None:
        channels = self.catalog.channels_for(self.type_var.get())
        self.channel_combo["values"] = channels
        if channels:
            self.channel_var.set(channels[0])
            self._update_samples()
        else:
            self.channel_var.set("")
            self.sample_list.delete(0, tk.END)

    def _update_samples(self) -> None:
        samples = self.catalog.samples_for(self.type_var.get(), self.channel_var.get())
        self.sample_list.delete(0, tk.END)
        for sample in samples:
            self.sample_list.insert(tk.END, f"{sample.sample_type} | {sample.name}")
        self._select_all()

    def _select_all(self) -> None:
        if self.sample_list.size():
            self.sample_list.selection_set(0, tk.END)

    def _clear_selection(self) -> None:
        self.sample_list.selection_clear(0, tk.END)

    def _plot(self) -> None:
        channel = self.channel_var.get()
        samples = self.catalog.samples_for(self.type_var.get(), channel)
        selected_indices = set(self.sample_list.curselection())
        selected = [sample for index, sample in enumerate(samples) if index in selected_indices]
        if not selected:
            messagebox.showwarning("Aucune sélection", "Sélectionnez au moins un échantillon.")
            return

        self.axes.clear()
        failures: list[str] = []
        areas: list[str] = []
        for sample in selected:
            try:
                data = self.catalog.read_channel(sample, channel)
            except (OSError, ValueError, pd.errors.ParserError) as error:
                failures.append(f"{sample.name}: {error}")
                continue
            area = area_under_curve(data)
            areas.append(f"{sample.name}: {area:.3g}")
            self.axes.plot(
                data["time_min"],
                data["intensity"],
                linewidth=1.0,
                label=f"{sample.name} (aire: {area:.3g})",
            )

        if not self.axes.lines:
            messagebox.showerror("Lecture impossible", "\n".join(failures))
            return
        self.axes.set_title(f"Canal m/z {channel} - {len(self.axes.lines)} courbe(s)")
        self.axes.set_xlabel("Temps (min)")
        self.axes.set_ylabel("Intensité")
        self.axes.set_yscale("log" if self.log_var.get() else "linear")
        self.axes.grid(True, alpha=0.25)
        if self.legend_var.get():
            self.axes.legend(fontsize="small")
        self.figure.tight_layout()
        self.canvas.draw()
        self.status_var.set(
            f"{len(self.axes.lines)} courbe(s) tracée(s) pour le canal {channel} | "
            f"Aire (intensité·min) : {' ; '.join(areas)}"
        )
        if failures:
            messagebox.showwarning("Certains fichiers n'ont pas été lus", "\n".join(failures))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default="cache", help="dossier contenant les caches CSV")
    args = parser.parse_args()
    Viewer(Path(args.root)).mainloop()


if __name__ == "__main__":
    main()
