#!/usr/bin/env python3
"""
extraire_csv.py : transforme les dossiers bruts de la machine (.D) en fichiers CSV.

Pour chaque échantillon (ex. 20251103-GAM-25-385-6.D) on écrit :
    <sortie>/<batch>/<échantillon>/mz_128.csv   (colonnes : time_min, intensity)
    <sortie>/<batch>/<échantillon>/mz_136.csv   ... un fichier par m/z
    <sortie>/<batch>/<échantillon>/multiplier.txt   (facteur de dilution)

Utilisation :
    python extraire_csv.py --batch 20251103 --sortie cache
    python extraire_csv.py --batch 20251103 --types GAM SF Hexane --zip
    python extraire_csv.py --batch 20251103 --verifier     (compare avec 'rainbow')

Le seul fichier lu dans chaque dossier .D est data.ms. Il contient, pour chaque
instant de mesure (un "scan"), une liste de couples (m/z, intensité).
Ce script le décode lui-même, étape par étape (voir lire_data_ms).
Dépendances : numpy, pandas  (rainbow-api seulement pour --verifier)
"""
import argparse
import shutil
import struct
from pathlib import Path

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# 1. Décoder data.ms
# --------------------------------------------------------------------------
def lire_data_ms(chemin):
    """Lit un fichier data.ms de type 'GC / MS Data File' (Agilent ChemStation).

    Renvoie : (temps_en_minutes, liste_des_m/z, matrice[temps x m/z]).

    Structure du fichier (déduite du code de la bibliothèque 'rainbow') :
      - début : signature 0x0132 puis le texte "GC / MS Data File"
      - à l'octet 0x142 : le nombre de scans (entier 2 octets, little-endian)
      - à l'octet 0x10A : la position du début des mesures
      - puis, pour chaque scan :  [2 octets ignorés] [temps : 4 octets, en ms]
                                  [6 ignorés] [nb de couples : 2 octets] [4 ignorés]
                                  [nb de couples x 4 octets : m/z et intensité]
                                  [10 ignorés]
    """
    d = Path(chemin).read_bytes()
    if struct.unpack_from(">I", d, 0)[0] != 0x01320000:
        raise ValueError("Signature inconnue : ce n'est pas un data.ms standard.")
    longueur = d[4]
    type_fichier = d[5:5 + longueur].decode("ascii", errors="replace")
    if type_fichier != "GC / MS Data File":
        raise ValueError(f"Type de fichier non géré par ce script : {type_fichier!r}")

    nb_scans = struct.unpack_from("<H", d, 0x142)[0]
    pos = struct.unpack_from(">H", d, 0x10A)[0] * 2 - 2   # début des mesures

    temps_ms = np.empty(nb_scans, dtype=np.uint32)
    nb_couples = np.empty(nb_scans, dtype=np.int64)
    blocs = []                                            # octets des couples
    for i in range(nb_scans):
        pos += 2
        temps_ms[i] = struct.unpack_from(">I", d, pos)[0]; pos += 4
        pos += 6
        n = struct.unpack_from(">H", d, pos)[0]; pos += 2
        pos += 4
        blocs.append(d[pos:pos + 4 * n]); pos += 4 * n
        pos += 10
        nb_couples[i] = n

    brut = b"".join(blocs)
    total = int(nb_couples.sum())
    # m/z : entier 2 octets divisé par 20 (la machine enregistre sur une grille de 0,05)
    mz = np.ndarray(total, ">H", brut, 0, 4) / 20
    # intensité : 2 bits d'"exposant" + 14 bits de valeur, intensité = valeur x 8^exposant
    codee = np.ndarray(total, ">H", brut, 2, 4)
    intensite = (codee & 0x3FFF).astype(np.uint32) * (8 ** (codee >> 14).astype(np.uint32))

    # On range les couples dans un tableau [scan x m/z], en arrondissant le m/z à
    # l'entier le plus proche (127,1 -> 127). Les couples qui tombent dans le même
    # entier seraient additionnés ; sur ces fichiers (33 m/z enregistrés), ça n'arrive pas.
    colonne_mz = np.rint(mz).astype(int)
    liste_mz = np.unique(colonne_mz)
    colonne = np.searchsorted(liste_mz, colonne_mz)
    ligne = np.repeat(np.arange(nb_scans), nb_couples)
    matrice = np.zeros((nb_scans, len(liste_mz)), dtype=np.uint32)
    np.add.at(matrice, (ligne, colonne), intensite)

    return temps_ms / 60000, liste_mz, matrice   # ms -> minutes


# --------------------------------------------------------------------------
# 2. Facteur de dilution : lu dans les journaux de séquence (.TSV) du batch
# --------------------------------------------------------------------------
def lire_multiplicateurs(dossier_batch):
    """Renvoie {nom_échantillon: facteur}. Un échantillon peut être dans l'un
    OU l'autre des journaux (ex. un SF passé le lendemain)."""
    res = {}
    for tsv in Path(dossier_batch).glob("*.TSV"):
        lignes = tsv.read_text(encoding="utf-8-sig").splitlines()
        entetes = lignes[1].split("\t")
        i_nom, i_mult = entetes.index("_dataname$"), entetes.index("_multiplr")
        for l in lignes[2:]:
            c = l.split("\t")
            if len(c) > i_mult:
                res[c[i_nom]] = c[i_mult]
    return res


# --------------------------------------------------------------------------
# 3. Export
# --------------------------------------------------------------------------
def exporter(batch, sortie, types, zipper):
    batch = Path(batch)
    cible_batch = Path(sortie) / batch.name
    mult = lire_multiplicateurs(batch)
    dossiers = sorted(d for d in batch.glob("*.D") if any(t in d.name for t in types))
    if not dossiers:
        raise SystemExit(f"Aucun dossier .D contenant {types} dans {batch}")

    for d in dossiers:
        nom = d.name.removesuffix(".D")
        t, liste_mz, matrice = lire_data_ms(d / "data.ms")
        dest = cible_batch / nom
        dest.mkdir(parents=True, exist_ok=True)
        for j, mz in enumerate(liste_mz):
            pd.DataFrame({"time_min": t, "intensity": matrice[:, j].astype(float)}) \
              .to_csv(dest / f"mz_{mz}.csv", index=False)
        (dest / "multiplier.txt").write_text(f"{float(mult.get(nom, '1')):.1f}\n")
        print(f"{nom}: {len(t)} instants, {len(liste_mz)} m/z")

    if zipper:
        archive = shutil.make_archive(str(Path(sortie) / f"cache_{batch.name}"),
                                      "zip", root_dir=sortie, base_dir=batch.name)
        print("Archive :", archive)


# --------------------------------------------------------------------------
# 4. Vérification contre la bibliothèque 'rainbow' (outil de référence)
# --------------------------------------------------------------------------
def verifier(batch, types):
    import rainbow as rb
    dossiers = sorted(d for d in Path(batch).glob("*.D") if any(t in d.name for t in types))
    for d in dossiers:
        t, mz, m = lire_data_ms(d / "data.ms")
        ref = rb.read(str(d)).datafiles[0]
        ok = (np.array_equal(t, ref.xlabels) and np.array_equal(mz, ref.ylabels)
              and np.array_equal(m, ref.data))
        print(f"{d.name}: {'identique à rainbow' if ok else 'DIFFÉRENT'}")


def main():
    p = argparse.ArgumentParser(description="Dossiers .D -> fichiers CSV par m/z")
    p.add_argument("--batch", required=True, help="dossier du batch (ex. 20251103)")
    p.add_argument("--sortie", default="cache", help="dossier de destination")
    p.add_argument("--types", nargs="+", default=["GAM", "SF"],
                   help="morceaux de nom à garder (défaut : GAM SF)")
    p.add_argument("--zip", action="store_true", help="créer aussi une archive .zip")
    p.add_argument("--verifier", action="store_true", help="comparer avec rainbow, sans écrire")
    a = p.parse_args()
    if a.verifier:
        verifier(a.batch, a.types)
    else:
        exporter(a.batch, a.sortie, a.types, a.zip)


if __name__ == "__main__":
    main()
