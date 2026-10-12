"""pipeline.py : référence GAM -> calibration quadratique -> contrôle SF.

Périmètre : les 16 HAP cibles, calibrés par rapport à leur étalon interne deutéré.
Les surrogates sont lus pour séparer correctement les pics co-élués (ex. m/z 240),
puis retirés de la référence finale.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, medfilt

TOL_RT = 0.002        # ±0,2 % du RT (relatif)
TOL_RATIO = 0.20      # ±20 % sur le ratio qualifiant/quantifiant
DEMI_FENETRE = 0.05   # min : zone de recherche autour du RT de référence
FIN_PIC = 0.02        # le pic s'arrête quand il retombe à 2 % de sa hauteur
SEUIL_SNR = 3         # sous 3 x bruit : non détecté


# ---------- Signal ----------
def lire_signal(dossier, mz):
    df = pd.read_csv(Path(dossier) / f"mz_{int(mz)}.csv")
    return df["time_min"].to_numpy(float), df["intensity"].to_numpy(float)


def lisser(y, n=5):
    return np.convolve(y, np.ones(n) / n, mode="same")


def signal_lisse(y, despike=False):
    """Filtre médian (anti-pics parasites) si demandé, puis lissage."""
    if despike:
        y = medfilt(y, 5)
    return lisser(y)


def bruit(y):
    return float(np.std(np.diff(y)) / np.sqrt(2)) or 1.0


def bornes(ys, apex):
    """Descend de chaque côté de l'apex jusqu'au bruit de fond ou jusqu'à une vallée."""
    base = float(np.percentile(ys, 10))
    haut = ys[apex] - base
    seuil = base + FIN_PIC * haut

    def descendre(sens):
        i, mini = apex, apex
        while 0 <= i + sens < len(ys):
            i += sens
            if ys[i] <= seuil:
                return i
            if ys[i] < ys[mini]:
                mini = i
            elif ys[i] > ys[mini] + 0.05 * haut:   # le signal remonte : autre pic
                return mini
        return i

    return descendre(-1), descendre(+1)


def delimiter(t, ys, apex, double=False):
    """Bornes du pic ; si 'double', on réunit les bornes de tous les sommets proches."""
    d, f = bornes(ys, apex)
    if double:
        base = float(np.percentile(ys, 10))
        idx, _ = find_peaks(ys, prominence=0.2 * (ys[apex] - base))
        for i in idx:
            if abs(t[i] - t[apex]) <= 0.1:
                di, fi = bornes(ys, i)
                d, f = min(d, di), max(f, fi)
    return d, f


def aire(t, y, d, f):
    """Somme de Riemann au-dessus de la droite reliant les deux bornes."""
    seg = y[d:f + 1]
    droite = np.linspace(y[d], y[f], len(seg))
    return float(np.clip(seg - droite, 0, None).sum() * np.median(np.diff(t)))


# ---------- Métadonnées ----------
def _propre(df):
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    for c in df.columns:
        df[c] = df[c].map(lambda v: v.strip() if isinstance(v, str) else v)
    return df


def charger_meta(chemin):
    comp = _propre(pd.read_excel(chemin, sheet_name="compound_informations"))
    deut = _propre(pd.read_excel(chemin, sheet_name="compound_deuterated"))
    return comp, deut


# ---------- 1. Référence (GAM-6) ----------
def construire_reference(dossier, comp):
    app = comp["applicable_to"].fillna("").astype(str)
    garde = ((app == "") | app.str.contains("GAM")) & (comp["use"].astype(str) != "Inutilisé")
    comp = comp[garde]
    lignes = []
    # composés partageant le même m/z : le k-ième pic (par ordre de sortie) -> k-ième composé
    for mz_q, g in comp.groupby("quantifier", sort=False):
        g = g.sort_values("coelution_order", na_position="last", kind="stable")
        t, y = lire_signal(dossier, mz_q)
        ys0 = signal_lisse(y, bool((g["median_despike"] == 1).any()))
        idx, p = find_peaks(ys0, prominence=0)
        top = sorted(idx[np.argsort(p["prominences"])[::-1][:len(g)]])
        for (_, c), apex in zip(g.iterrows(), top):
            desp, dbl = c["median_despike"] == 1, c["double_peak"] == 1
            ys = signal_lisse(y, desp)
            d, f = delimiter(t, ys, apex, dbl)
            _, yq = lire_signal(dossier, c["qualifier"])
            aq, aqual = aire(t, y, d, f), aire(t, yq, d, f)
            lignes.append({
                "compose": c["compound_name"], "type": c["type"],
                "mz_quant": int(mz_q), "mz_qual": int(c["qualifier"]),
                "despike": bool(desp), "double": bool(dbl),
                "rt_ref": t[apex], "debut": t[d], "fin": t[f],
                "aire_quant": aq, "ratio_ref": aqual / aq if aq else np.nan,
            })
    ref = pd.DataFrame(lignes)
    # On garde les 16 HAP cibles + les 4 étalons internes (utiles aux ratios).
    # Les surrogates ont servi à séparer les pics co-élués, ils sortent ici.
    return ref[ref["type"].astype(str).str.lower().isin(["target", "istd"])].reset_index(drop=True)


# ---------- 2. Analyse d'un échantillon avec la référence ----------
def analyser(dossier, ref):
    lignes = []
    for _, r in ref.iterrows():
        t, y = lire_signal(dossier, r["mz_quant"])
        ys = signal_lisse(y, r["despike"])
        zone = np.where(np.abs(t - r["rt_ref"]) <= DEMI_FENETRE)[0]
        base = float(np.percentile(ys, 10))
        apex = zone[np.argmax(ys[zone])] if len(zone) else None
        ligne = {"compose": r["compose"], "type": r["type"], "mz_quant": r["mz_quant"],
                 "mz_qual": r["mz_qual"], "rt": np.nan, "debut": np.nan, "fin": np.nan,
                 "aire_quant": np.nan, "ratio": np.nan, "ecart_rt": np.nan,
                 "ecart_ratio": np.nan, "statut": "non détecté"}
        if apex is not None and ys[apex] - base > SEUIL_SNR * bruit(y):
            d, f = delimiter(t, ys, apex, r["double"])
            _, yq = lire_signal(dossier, r["mz_qual"])
            aq, aqual = aire(t, y, d, f), aire(t, yq, d, f)
            ratio = aqual / aq if aq else np.nan
            e_rt = abs(t[apex] - r["rt_ref"]) / r["rt_ref"]
            e_ratio = abs(ratio - r["ratio_ref"]) / r["ratio_ref"]
            ligne.update(rt=t[apex], debut=t[d], fin=t[f], aire_quant=aq, ratio=ratio,
                         ecart_rt=e_rt, ecart_ratio=e_ratio,
                         statut="confirmé" if (e_rt <= TOL_RT and e_ratio <= TOL_RATIO) else "douteux")
        lignes.append(ligne)
    return pd.DataFrame(lignes)


# ---------- 3. Calibration + 4. Contrôle SF ----------
def inverser(coef, ratio):
    """Résout a·C² + b·C + c = ratio ; garde la racine proche de l'estimation linéaire."""
    a, b, c = coef
    lin = (ratio - c) / b if b else np.nan
    if abs(a) < 1e-15:
        return lin
    racines = np.roots([a, b, c - ratio])
    reelles = racines[np.isreal(racines)].real
    if len(reelles) == 0:
        return np.nan
    return float(reelles[np.argmin(np.abs(reelles - lin))])


def multiplicateur(dossier):
    try:
        return float((Path(dossier) / "multiplier.txt").read_text().strip())
    except (OSError, ValueError):
        return 1.0


def numero(d):
    return int(d.name.split("-")[-1])


def executer(batch, meta_path, concs, cible_sf, n_ref=6):
    batch = Path(batch)
    comp, deut = charger_meta(meta_path)

    # cible -> étalon deutéré, insensible à la casse ("Pyrene" == "pyrene")
    noms = {str(n).lower(): n for n in comp["compound_name"]}
    istd = {}
    for cible, d in zip(deut.iloc[:, 0], deut.iloc[:, 1]):
        istd[noms.get(str(cible).lower(), cible)] = noms.get(str(d).lower(), d)

    gams = sorted((d for d in batch.iterdir() if d.is_dir() and "GAM" in d.name), key=numero)
    sfs = sorted(d for d in batch.iterdir() if d.is_dir() and "SF" in d.name)
    if len(gams) != len(concs):
        raise ValueError(f"{len(gams)} GAM trouvés mais {len(concs)} concentrations données.")

    ref = construire_reference(next(d for d in gams if numero(d) == n_ref), comp)
    cibles = ref.loc[ref["type"].astype(str).str.lower() == "target", "compose"].tolist()

    detail, points = [], []
    for d, conc in zip(gams, concs):
        m = analyser(d, ref).assign(echantillon=d.name)
        detail.append(m)
        m = m.set_index("compose")
        for c in cibles:
            a_i = m["aire_quant"].get(istd.get(c), np.nan)
            points.append({"compose": c, "conc": conc,
                           "ratio": m.loc[c, "aire_quant"] / a_i if a_i > 0 else np.nan})
    points = pd.DataFrame(points)

    calib = []
    for c in cibles:
        p = points[points.compose == c].dropna()
        if len(p) < 3:
            continue
        coef = np.polyfit(p.conc, p.ratio, 2)
        res = p.ratio - np.polyval(coef, p.conc)
        r2 = 1 - (res ** 2).sum() / ((p.ratio - p.ratio.mean()) ** 2).sum()
        calib.append({"compose": c, "istd": istd.get(c), "a": coef[0], "b": coef[1],
                      "c": coef[2], "r2": r2, "n_points": len(p)})
    calib = pd.DataFrame(calib)
    if calib.empty:
        raise ValueError("Aucun composé n'a pu être calibré (pics ou ISTD non détectés sur les GAM).")
    coefs = calib.set_index("compose")[["a", "b", "c"]]

    sf = []
    for d in sfs:
        m = analyser(d, ref).assign(echantillon=d.name)
        detail.append(m)
        mult = multiplicateur(d)
        m = m.set_index("compose")
        for c in coefs.index:
            a_i = m["aire_quant"].get(istd.get(c), np.nan)
            ratio = m.loc[c, "aire_quant"] / a_i if a_i > 0 else np.nan
            conc = inverser(coefs.loc[c].to_numpy(), ratio) * mult if ratio == ratio else np.nan
            rec = conc / cible_sf * 100
            sf.append({"echantillon": d.name, "compose": c, "ratio": ratio,
                       "multiplier": mult, "conc_ppm": conc, "recuperation_%": rec,
                       "valide": bool(abs(rec - 100) <= 20) if rec == rec else False})
    return ref, pd.concat(detail, ignore_index=True), points, calib, pd.DataFrame(sf)