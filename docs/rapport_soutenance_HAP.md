# Soutenance : dosage des HAP dans les enrobés

Ce document prépare les questions que le jury peut poser. Pour chaque question, on donne d'abord la réponse en une phrase, puis l'explication, puis la limite.

Les chiffres ont été relus et vérifiés à partir des fichiers produits par nos scripts. Ce document a été rédigé avec l'aide d'un assistant IA. Les calculs et les résultats ont été contrôlés par le groupe en lançant les programmes.

---

## Thème 1 : le processus métier

### T1-Q1 : Pourquoi mesure-t-on les HAP dans les enrobés ?
**Réponse :** certains HAP sont cancérigènes, donc on vérifie leur teneur avant d'utiliser le matériau.

**Explication :** les HAP se forment lors de la fabrication et de la combustion. Une teneur trop élevée peut rendre le matériau impropre à son usage. Le résultat dit donc si l'enrobé est conforme, mais la décision finale revient au laboratoire.

### T1-Q2 : Comment fonctionne la chromatographie en phase gazeuse ?
**Réponse :** on fait passer le mélange dans une colonne, et chaque composé en sort à un moment différent, son temps de rétention.

**Explication :** les composés les plus légers sortent en premier. Le détecteur de masse enregistre l'intensité de plusieurs ions, et chaque composé se reconnaît à son temps de rétention et à ses ions caractéristiques.

**Limite :** deux composés de masses proches peuvent sortir presque en même temps. C'est le cas de certains HAP, traité au thème 3.

### T1-Q3 : À quoi sert l'étalon interne deutéré ?
**Réponse :** il corrige les variations de l'appareil d'une injection à l'autre.

**Explication :** on ajoute à chaque échantillon une quantité connue d'une molécule presque identique, mais marquée au deutérium. On ne regarde pas l'aire du HAP seul, mais son rapport à l'aire de l'étalon. Si le signal baisse un peu, il baisse pour les deux.

### T1-Q4 : À quoi servent les références et les contrôles ?
**Réponse :** ils permettent de savoir si la mesure est juste avant de faire confiance au résultat.

**Explication :** la solution GAM-6 à 1 ppm sert de référence : on y mesure le temps de rétention et les ions de chaque composé. Les solutions SF, également à 1 ppm, servent de contrôle : on les calcule comme des inconnues et on vérifie qu'on retrouve la bonne valeur.

### T1-Q5 : Quelles sont les limites du processus ?
**Réponse :** il ne donne des résultats fiables qu'entre 0,025 et 5 ppm, et il sépare mal certains composés.

**Explication :** en dehors de la gamme, on ne conclut pas. Certains composés sortent presque ensemble, et notre méthode ne les sépare pas toujours.

---

## Thème 2 : données d'entrée et courbe de calibration

### T2-Q1 : Quelles sont les données d'entrée ?
**Réponse :** pour chaque injection, un fichier par ion, avec le temps et l'intensité.

**Explication :** l'extraction transforme les fichiers de l'appareil en fichiers `mz_XXX.csv`. Chaque injection donne donc 33 fichiers, et on en tire un signal par composé.

**Défauts possibles :** du bruit, une ligne de base qui dérive, des pics parasites. Pour les réduire, on lisse le signal sur 5 points et, si besoin, on applique un filtre médian qui supprime les pics isolés.

### T2-Q2 : Pourquoi une courbe quadratique ?
**Réponse :** parce que le modèle imposé par le cours est une courbe de degré 2.

**Explication :** on ajuste rapport = a × C² + b × C + c, où C est la concentration et le rapport est l'aire du composé divisée par l'aire de son étalon. Le terme en C² permet de suivre une légère courbure.

**Limite :** nous n'avons pas comparé systématiquement avec une droite simple. Si le jury le demande, reconnais-le : nous avons appliqué le modèle demandé.

### T2-Q3 : Quels critères valident la courbe ?
**Réponse :** un coefficient de détermination R² d'au moins 0,99.

**Explication :** le R² mesure à quel point les points sont proches de la courbe. Sur nos données, 15 composés sur 16 dépassent ce seuil, avec des valeurs entre 0,9996 et 0,99999 pour ces composés. Le dibenz(a,h)anthracène a un R² de 0,986 : il est rejeté.

**Ce qu'on fait quand la courbe ne passe pas :** on ne l'utilise pas. Ses résultats sont marqués « non valide », et on ne les présente pas comme des résultats. Le benzo(b)fluoranthène et le benzo(k)fluoranthène ont exactement le même R², ce qui s'explique par le problème de co-élution.

### T2-Q4 : Comment passe-t-on d'une mesure sur un pic à une concentration ?
**Réponse :** on calcule le rapport entre l'aire du composé et l'aire de son étalon, puis on résout la courbe pour trouver la concentration.

**Explication :** résoudre la courbe revient à trouver la concentration C qui donne ce rapport. L'équation peut avoir deux solutions : on garde celle qui est la plus proche d'une estimation par droite simple.

**Si l'échantillon est trop concentré :** au-dessus de 5 ppm, on utilise la solution diluée au 20e, puis on multiplie le résultat par 20. Le facteur de dilution est lu dans le fichier `multiplier.txt`.

### T2-Q5 : Quel est le domaine de validité ?
**Réponse :** de 0,025 à 5 ppm.

**Explication :** c'est la gamme couverte par les huit solutions GAM. En dehors, on n'extrapole pas : le résultat est marqué « sous la gamme » ou « hors gamme ».

**Limite :** sous 0,5 ppm, le rapport entre ion qualifiant et ion quantifiant varie beaucoup d'une solution propre à l'autre. Le test de stabilité du ratio le montre.

---

## Thème 3 : l'algorithme de détermination des temps d'intégration

### T3-Q1 : Comment l'algorithme trouve-t-il un pic ?
**Réponse :** il cherche le point le plus haut dans une petite fenêtre autour du temps de rétention de référence.

**Explication :** la fenêtre fait ±0,05 minute autour du temps de référence. Le point le plus haut du signal lissé dans cette fenêtre est le sommet du pic. Le pic n'est retenu que si ce sommet dépasse trois fois le bruit de fond.

**Pourquoi une fenêtre et pas tout le signal :** sur tout le signal, il y a beaucoup de pics. La fenêtre permet de ne chercher que le composé attendu, au temps où il doit sortir.

### T3-Q2 : Comment trouve-t-il le début et la fin du pic ?
**Réponse :** il descend de chaque côté du sommet jusqu'à ce que le signal retombe à 2 % de la hauteur du pic.

**Explication :** si le signal remonte avant, cela veut dire qu'un autre pic commence, et l'algorithme s'arrête au creux. On obtient ainsi deux bornes : le début et la fin.

### T3-Q3 : Comment calcule-t-il l'aire ?
**Réponse :** il additionne les intensités entre les deux bornes, au-dessus d'une droite qui relie le début et la fin.

**Explication :** cette droite représente la ligne de base. On retire donc le fond, et il ne reste que la surface du pic.

### T3-Q4 : Quels paramètres utilises-tu, et d'où viennent-ils ?
**Réponse :** une fenêtre de ±0,05 minute, une détection à trois fois le bruit, une fermeture du pic à 2 % de sa hauteur, et un lissage sur cinq points.

**Explication :** ces valeurs sont fixées dans le programme `pipeline.py`. Elles n'ont pas été imposées par le cahier des charges. Nous les avons vérifiées sur les solutions connues, et elles donnent des résultats proches de MassHunter pour la plupart des composés.

### T3-Q5 : Que se passe-t-il quand deux pics se chevauchent ?
**Réponse :** notre méthode leur attribue la même aire, donc l'un des deux résultats est faux.

**Explication :** le phénanthrène et l'anthracène ont presque le même temps de rétention et utilisent le même ion. Notre algorithme prend le même sommet pour les deux, donc la même aire. MassHunter, lui, les sépare. Il en va de même pour le benzo(a)anthracène et la chrysène, et pour le benzo(b)fluoranthène et le benzo(k)fluoranthène.

**Ce qu'il faudrait faire :** attribuer à chaque composé le sommet le plus proche de son propre temps de rétention, sans reprendre un sommet déjà pris. Nous ne l'avons pas fait faute de temps, et une telle modification changerait tous les résultats.

### T3-Q6 : Quel est l'effet d'une borne mal placée ?
**Réponse :** une borne trop courte fait perdre de la surface, une borne trop large ajoute du bruit, et les deux modifient la concentration finale.

**Explication :** la concentration vient directement de l'aire. Une aire fausse de 10 % donne une concentration fausse d'environ 10 %, après passage dans la courbe.

### T3-Q7 : Comment as-tu vérifié que l'algorithme fonctionne ?
**Réponse :** sur les solutions de concentration connue, nos résultats sont très proches de ceux de MassHunter.

**Explication :**
- Sur les solutions GAM, MassHunter retrouve la valeur nominale : à 0,05 ppm les écarts vont de −4 % à +2 %, et à 5 ppm de −0,2 % à +0,5 %.
- Sur les solutions SF à 1 ppm, nos valeurs et celles de MassHunter sont proches pour la plupart des composés, à 1 à 3 % près. Les deux sous-estiment la solution d'environ 10 %, probablement à cause de la préparation de la solution. Ce n'est pas une erreur de calcul, mais il faut le signaler.
- Le benzo(b)/(k)fluoranthène s'écarte davantage, d'environ 3 à 4 %. Le dibenz(a,h)anthracène s'écarte beaucoup, de 25 % à 37 %, ce qui correspond à sa calibration non valide.

**Limite :** sur les échantillons réels, les résultats ne concordent pas aussi bien. Sur 405 comparaisons avec MassHunter, 92 sont à moins de 10 % (soit 23 %). Le taux d'accord sur le caractère positif ou négatif, au seuil de 0,5 ppm, est de 85 % (344 couples sur 405). Les écarts sont grands pour le naphtalène, le fluorène et les benzo(b)/(k)fluoranthène. Une partie vient de la sélection des pics dans des échantillons plus complexes, et une autre des très faibles concentrations : un petit écart en valeur absolue donne un grand écart en pourcentage.

---

## Questions pièges

- **Pourquoi n'avez-vous pas corrigé le problème des co-élués ?** Nous l'avons identifié tard. Une correction change tous les résultats, nous l'avons donc documentée plutôt que de la faire à la hâte.
- **Peut-on faire confiance à l'application ?** Pour la calibration et le contrôle sur les solutions connues, oui, nous l'avons vérifié. Pour les composés co-élués et le dibenz(a,h)anthracène, non, et nous le signalons.
- **Que feriez-vous avec plus de temps ?** Séparer les pics co-élués, valider la méthode sur davantage d'échantillons, et comparer avec d'autres modèles de calibration.
