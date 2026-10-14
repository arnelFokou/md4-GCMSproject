# md4-GCMSproject

Projet GC/MS : calibration et contrôle qualité des HAP

Automatisation de la calibration et du contrôle qualité de 16 HAP (hydrocarbures
aromatiques polycycliques) dans des échantillons de sol agricole, analysés par GC-MS.
Projet réalisé dans le cadre du cours d'analytics engineering (HETIC).


# Données (non versionnées)

Les données brutes ne sont pas sur GitHub (voir `.gitignore`). Il faut placer :

- le dossier du batch dans `data/20251103/` (fichiers `.D`) ;
- le fichier `metadata.xlsx` à la racine du projet, fourni par l'enseignant sur le drive partagé.

# Ordre d'exécution

1. Extraire les données brutes en CSV :
   `python extraire_csv.py --batch data\20251103 --types GAM SF BLPC --sortie outputs`
2. Calibration et contrôle SF : `python run_calibration.py`
3. Résultats finaux selon les critères du sujet : `python run_final.py`
4. Comparaison avec MassHunter : `python run_comparaison.py`

Les autres scripts servent à vérifier et expliquer les résultats :

- `run_test_ratio.py` : stabilité du ratio ionique sur les gammes propres ;
- `run_diagnostic.py` : temps de rétention et ratios sur les échantillons ;
- `run_diagnostic_mh.py` : contrôle sur les solutions de concentration connue, comparé à MassHunter ;
- `run_diagnostic_pics.py` : comparaison pic par pic (temps de rétention et aire) avec MassHunter ;
- `run_samples.py` : première version de la quantification, remplacée par `run_final.py`.

# Visualisation (Streamlit)

Lancer avec `streamlit run streamlit_app.py`. Le chemin par défaut est `outputs/20251103`.

- `streamlit_app.py` : visualisation des courbes `mz_XXX.csv`, comparaison d'échantillons.
- Page Calibration : référence GAM, courbes de calibration, contrôle SF, visualisation d'un pic
  (zone de recherche et tolérance du temps de rétention).
- Page Résultats : décisions par composé et concentrations retenues.

# Limites connues

- Dibenz(a,h)anthracène : calibration non valide (R² 0,986 < 0,99).
- Ratio ionique instable en dessous de 0,5 ppm.
- Phénanthrène/anthracène, benzo(a)anthracène/chrysène et benzo(b)/(k)fluoranthène : co-élués,
  notre méthode leur attribue la même aire, MassHunter les sépare.
- L'intégration du pic de Streamlit (`visualiser_csv.py`) utilise une méthode distincte de
  `pipeline.py` : ses aires ne doivent pas être comparées aux résultats du rapport.

# Fichiers

- `pipeline.py` : fonctions de calcul (pics, aires, calibration, inversion).
- `extraire_csv.py` : extraction des données brutes en CSV.
- `run_*.py` : scripts de chaque étape.
- `streamlit_app.py`, `pages/` : interface Streamlit.
- `visualiser_csv.py` : visionneuse et calcul d'aire utilisés par Streamlit.


# Contributions

- Extraction des données (`extraire_csv.py`), application Streamlit (`streamlit_app.py`, `visualiser_csv.py`, `pages/`) et `pipeline.py` : Arnel Fokou.
- Scripts de comparaison et de diagnostic (`run_comparaison.py`, `run_diagnostic*.py`) : Biolohf Tchatchoua.
- Calibration, résultats finaux, tests et quantification (`run_calibration.py`, `run_final.py`, `run_test_ratio.py`, `run_samples.py`), ajustements de l'application Streamlit (paramètres par défaut, cadre de recherche du pic, page Résultats) et rédaction du README : Franck Kevin Tambou.
- Une partie du code et de la documentation a été rédigée avec l'aide d'un assistant IA. Les scripts ont été lancés et les résultats vérifiés par les auteurs.

