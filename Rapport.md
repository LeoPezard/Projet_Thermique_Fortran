---
lang: fr
papersize: a4
geometry: margin=1.3cm
fontsize: 15pt

header-includes:
  - \usepackage{float}
  - \usepackage{multirow}
  - \usepackage{graphicx}
  - \usepackage{tabularx}
---

\listoffigures

\listoftables


# Introduction

Après avoir vu en TD le profil de température dans un mur en 1D en régime stationnaire, cette note de calcul a pour objectif de modéliser le comportement thermique d'une cellule plane de batterie lithium-ion par la méthode des volumes finis.

Quatre niveaux de modélisation sont abordés de manière progressive :

- Simulation d'une cellule élémentaire multicouche avec ses cinq couches constitutives, en convection symétrique.
- Détermination d'un milieu homogène équivalent (conductivité et source de chaleur équivalentes) représentant la cellule de façon simplifiée.
- Simulation d'un assemblage de 50 cellules et détermination du coefficient d'échange h minimal pour maintenir la batterie dans des plages de température de sécurité lors d'une charge lente et d'une charge rapide.
- Prise en compte d'une résistance de contact entre chaque cellule impactant ainsi la conduction entre chaque couche.

Le code de calcul est écrit en Fortran 95, basé sur le programme développé en TD (résolution TDMA, volumes finis, conditions aux limites de convection). Le programme a été étendu pour gérer les maillages multicouches et les conductivités harmoniques aux interfaces.

# Description du problème

Les grandeurs connues sont l'épaisseur et les propriétés thermophysiques de chaque
couche, les sources volumiques de chaleur et les conditions aux limites convectives
($h$, $T_{fluide}$). Les grandeurs recherchées sont le profil de température $T(x)$
dans la cellule, la température maximale $T_{max}$, et — pour l'assemblage — le
coefficient d'échange minimal $h_{min}$ permettant de respecter les seuils de sécurité.

##  Géométrie et composition de la cellule

La cellule élémentaire est assimilée à une **plaque plane multicouche 1D**, d'épaisseur totale L = 250 $\mu m$ , composée de cinq couches successives décrites dans le Tableau 1.

```{=latex}
\begin{table}[H]
\centering
\caption{Composition d'une cellule élémentaire Li-Ion et sources volumiques de chaleur}
\begin{tabularx}{\textwidth}{|X|c|c|c|c|}
\hline
\textbf{Couche / Matériau} & \textbf{Épaisseur ($\mu$m)} & \textbf{K (W/m·K)} & \textbf{S lente (W/m³)} & \textbf{S rapide (W/m³)} \\
\hline
Collecteur cathode — Aluminium & 20 & 237,0 & 640 & 6 400 \\
Couche active cathode — NMC/LFP & 100 & 1,58 & 2 720 & 27 200 \\
Séparateur — Polyoléfine microporeuse & 25 & 0,33 & 1 280 & 12 800 \\
Couche active anode — Graphite & 90 & 5,0 & 2 720 & 27 200 \\
Collecteur anode — Cuivre & 15 & 398,0 & 640 & 6 400 \\
\hline
\end{tabularx}
\end{table}
```

##  Phénomènes physiques et hypothèses

La chaleur est générée par effet Joule dans chaque couche. Les hypothèses retenues sont les suivantes :

- Régime permanent (stationnaire) — bilan thermique en équilibre.
- Géométrie 1D plane — pas de variations transversales de température.
- Propriétés thermiques homogènes et indépendantes de la température dans chaque couche.
- Symétrie du refroidissement — convection identique sur les deux faces pour une cellule en cœur d'assemblage.

## Propriétés thermophysiques — justification des conductivités

Les sources volumiques sont fournies dans l'énoncé. Les conductivités thermiques (direction normale aux couches, direction through-plane) ont été recherchées dans la littérature :

```{=latex}
\begin{table}[H]
\centering
\caption{Conductivités thermiques retenues et leur justification bibliographique}
\begin{tabularx}{\textwidth}{|l|c|c|X|}
\hline
\textbf{Matériau} & \textbf{K (W/m·K)} & \textbf{Incertitude} & \textbf{Source / Justification} \\
\hline
Aluminium & 237,0 & ±1\% & Incropera \& DeWitt, Table A.1 [1] — valeur standard bien établie. \\
NMC/LFP cathode & 1,58 & ±30\% & Maleki et al. (1999), J. Electrochem. Soc. [2] — direction through-plane. \\
Séparateur polyoléfine & 0,33 & ±20\% & Shi et al. (2020), J. Power Sources [3] — varie 0,1–0,5 W/mK selon la porosité. \\
Graphite anode & 5,0 & ±40\% & Maleki et al. (1999) [2] — direction through-plane ; varie 1–10 W/mK. \\
Cuivre & 398,0 & ±1\% & Incropera \& DeWitt, Table A.1 [1] — valeur standard. \\
\hline
\end{tabularx}
\end{table}
```

**Note sur l'incertitude :** Les conductivités des électrodes actives (cathode NMC et anode graphite) présentent une variabilité importante selon la formulation, la densification et la direction de mesure. Une incertitude de ±30–40 % sur K entraîne une incertitude similaire sur le gradient de température dans ces couches, mais peu d'effet sur la résistance thermique totale (dominée par les couches épaisses à faible K). L'incertitude principale sur $T_{max}$ est estimée à ±15 % pour la cellule isolée.

##  Conditions aux limites

Pour représenter une cellule en cœur d'un assemblage refroidi des deux côtés (configuration symétrique) :

- Face gauche (Ouest) : convection — flux $\Phi =  h \cdot (T_{face} - T_{fluide})$, avec $T_{fluide} = 40$ °C
- Face droite (Est) : convection — même condition par symétrie.
- Coefficient h = 10 $W/m^2 K$ — refroidissement naturel / convection forcée modérée, typique pour des cellules prismatiques avec plaque froide selon les préconisations NREL 2024 [4].



#  Modèle numérique

##  Méthode des volumes finis

La méthode des volumes finis 1D est appliquée sur un maillage structuré. Pour chaque volume de contrôle i de demi-largeur $dx_i$, le bilan d'énergie s'écrit :

$$a_w \cdot T_{i-1} - a_p \cdot T_i + a_e \cdot T_{i+1} = -b_i$$

avec les conductances nodales calculées par la moyenne harmonique des conductivités des volumes voisins, assurant la continuité du flux aux interfaces entre couches de conductivités différentes :

$$k_e = \frac{2 \cdot K_i \cdot K_{i+1} }{(K_i + K_{i+1})}$$ 

Le terme source intégré sur le volume utilise la méthode des trapèzes (conformément au TD01 modifié) :

$b_i = S_i \cdot dx_i$ (source uniforme par couche)

## Résolution TDMA

Le système tridiagonal est résolu par l'algorithme de Thomas (TDMA), implémenté en
deux passes. La première passe (descente) calcule les facteurs de récurrence $P_i$
et $Q_i$ :

$$P_i = \frac{a_e}{a_p - a_w \cdot P_{i-1}} \qquad Q_i = \frac{b_i + a_w \cdot Q_{i-1}}{a_p - a_w \cdot P_{i-1}}$$

La seconde passe (remontée) reconstruit le profil de température :

$$T_i = P_i \cdot T_{i+1} + Q_i$$

La robustesse est assurée par un test de pivot nul (seuil $10^{-12}$).

La discrétisation des conditions aux limites de convection s'écrit pour le nœud de bord :

$$a_p = a_e + h \qquad b_1 = S_1 \cdot dx_1 + h \cdot T_{\text{fluide}}$$

##  Maillage multicouche

Le nombre de volumes par couche est réparti proportionnellement à l'épaisseur de chaque couche, avec un minimum de 1 volume par couche. Par exemple, pour N = 100 volumes au total, la répartition obtenue est :

```{=latex}
\begin{table}[H]
\centering
\caption{Répartition des volumes finis pour N = 100}
\begin{tabular}{|l|c|c|}
\hline
\textbf{Couche} & \textbf{$N_{vol}$} & \textbf{dx ($\mu$m)} \\
\hline
Collecteur Al cathode   & 8  & 2,5 \\
Couche active cathode   & 40 & 2,5 \\
Séparateur              & 10 & 2,5 \\
Couche active anode     & 36 & 2,5 \\
Collecteur Cu anode     & 6  & 2,5 \\
\hline
\end{tabular}
\end{table}
```

# Objectif 1 — Cellule multicouche : résultats et validation

##  Profil de température de la cellule

Les profils de température obtenus **sans source interne** sont tracés afin de valider l'hypothèse de linéarité de la température dans chaque couche. Le nombre d'échantillons n ne change pas les valeurs des températures mais permet un maillage plus fin du domaine.

Nous obtenons ainsi un profil de température linéaire dans chaque couche de la cellule, avec pour valeurs aux limites suivantes :

```{=latex}
\begin{table}[H]
\centering
\caption{Conditions aux limites — sans source interne}
\begin{tabular}{|l|c|c|c|}
\hline
 & \textbf{h (W/m²K)} & \textbf{T (°C)} & \textbf{q (W)} \\
\hline
Ouest & 10 & 40 & 0 \\
Est   & 10 & 50 & 0 \\
\hline
\end{tabular}
\end{table}
```

```{=latex}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{1_ep_n=100.png}
\caption{Profil de température (n = 100)}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{1_ep_n=1000.png}
\caption{Profil de température (n = 1000)}
\end{figure}
```

En ajoutant **une source interne**, le profil de température à l'intérieur de chaque couche devient parabolique

```{=latex}
\begin{table}[H]
\centering
\caption{Conditions aux limites — avec source interne}
\begin{tabular}{|l|c|c|c|}
\hline
 & \textbf{h (W/m²K)} & \textbf{T (°C)} & \textbf{q (W)} \\
\hline
Ouest & 10 & 40 & 0 \\
Est   & 10 & 40 & 0 \\
\hline
\end{tabular}
\end{table}
```

```{=latex}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{1_ep_source.png}
\caption{Profil de température avec source interne (n = 1000)}
\end{figure}
```

**Analyse :** Le gradient de température dans une cellule individuelle est extrêmement faible en raison de la très faible épaisseur de la cellule (250 $\mu m$) et de sa conductivité relativement élevée. C'est la résistance convective de surface qui contrôle l'écart $T_{max} - T_{fluide}$. Avec $h = 10 W/m^2 K$, même en charge rapide, la cellule reste en dessous des seuils de sécurité grâce à la faible puissance volumique d'une cellule seule.

##  Comparaison solution numérique / solution analytique

**Pour un milieu homogène sans source interne, le profil de température est sous la forme :** 

$$T(x) = \frac{T_{Est}-T_{Ouest}}{e_{milieu}} \cdot x + T_{Ouest}$$

\underline{Vérification de la continuité du flux thermique :}

En régime permanent sans source interne, la loi de Fourier impose un flux constant 
dans tout le domaine :

$$\varphi = -\lambda_j \cdot \frac{\Delta T_j}{e_j} = \text{cste}$$

Le flux traversant la cellule est extrait depuis la face Ouest :

$$\varphi = h \cdot (T(0) - T_{ext,O}) = 10 \cdot (45{,}024 - 40) = 50{,}24 \ W/m^2$$

\underline{Vérification sur la couche NMC (x = 20 µm à x = 120 µm, $\lambda_{NMC} = 1{,}58 \ W/m \cdot K$) :}

La chute de température attendue analytiquement sur cette couche vaut :

$$\Delta T_{NMC}^{analytique} = -\frac{\varphi \cdot e_{NMC}}{\lambda_{NMC}} 
= -\frac{50{,}24 \times 100 \times 10^{-6}}{1{,}58} = -3{,}18 \times 10^{-3} \text{°C}$$

*Le signe négatif de $\Delta T^{\text{analytique}}$  traduit la convention de Fourier (flux positif dans le sens croissant de x) ; la comparaison porte sur les valeurs absolues.*


La chute de température lue sur le profil numérique entre x = 20 µm et x = 120 µm 
vaut :

$$\Delta T_{NMC}^{numérique} = T(120 \ \mu m) - T(20 \ \mu m) = 45{,}028 - 45{,}024 
= 4{,}00 \times 10^{-3} \text{°C}$$

L'écart relatif entre les deux valeurs est de :

$$\epsilon = \frac{|\Delta T^{num} - \Delta T^{ana}|}{|\Delta T^{ana}|} 
= \frac{|4{,}00 - 3{,}18|}{3{,}18} \approx 26 \ \%$$

Cet écart est imputable à la précision de lecture graphique sur un intervalle de 
température très faible (~0,004°C). 

**Pour un milieu homogène avec source uniforme S et convection symétrique** (h, $T_f$), la solution analytique est $T(x) = -\frac{S}{2K}\cdot x \cdot(x-L) + T_{bord}$ avec $T_{bord} = T_f + \frac{S\cdot L}{2 \cdot h}$. La solution numérique concorde avec cette référence, confirmant la bonne implémentation du schéma TDMA et des conditions aux limites convectives.


#  Objectif 2 — Milieu homogène équivalent

## Détermination du milieu équivalent

**Conductivité équivalente $K_{eq}$ :** Les couches thermiques sont en série, donc la résistance totale est la somme des résistances individuelles. La conductivité équivalente est calculée par :

**$R_{tot} = \Sigma \frac{e_j}{K_j}$ ,     $K_{eq} = L_{tot} / R_{tot}$**

**Source équivalente $S_{eq}$ :** Elle est calculée comme la moyenne volumique des sources de chaque couche : $S_{eq} = \Sigma(S_j \cdot e_j) / L_{tot}$, garantissant la conservation de la puissance totale dissipée dans la cellule.


```{=latex}
\begin{table}[H]
\centering
\caption{Paramètres du milieu équivalent homogène}
\begin{tabular}{|l|c|l|}
\hline
\textbf{Grandeur} & \textbf{Valeur} & \textbf{Unité / Remarque} \\
\hline
$L_{tot}$       & 250    & $\mu$m — épaisseur totale d'une cellule \\
$K_{eq}$        & 1,591  & W/m·K — dominé par le séparateur \\
$S_{eq}$ lente  & 2 285  & W/m³ — moyenne volumique \\
$S_{eq}$ rapide & 22 848 & W/m³ — moyenne volumique \\
\hline
\end{tabular}
\end{table}
```


**Remarque :** La conductivité équivalente $K_{eq} = 1,591 W/mK$ est très proche de celle de la couche active cathode (1,58 W/mK), car la résistance thermique est dominée par les couches les plus épaisses à faible conductivité (cathode, séparateur, anode). Les collecteurs métalliques (Al et Cu) ont des conductivités très élevées (>100 W/mK) mais contribuent peu à la résistance thermique totale vu leur faible épaisseur.

## Comparaison multicouche / milieu équivalent

Pour observer le profil de température, sans source interne, ces conditions de simulations ont été utilisées : 

```{=latex}
\begin{table}[H]
\centering
\caption{Conditions aux limites — comparaison sans source interne}
\begin{tabular}{|l|c|c|c|}
\hline
 & \textbf{h (W/m²K)} & \textbf{T (°C)} & \textbf{q (W)} \\
\hline
Ouest & 10 & 40 & 0 \\
Est   & 10 & 50 & 0 \\
\hline
\end{tabular}
\end{table}
```

```{=latex}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{K_eq_1_ep_PAS_S.png}
\caption{Profil avec conductivité équivalente, sans source interne}
\end{figure}
```

Pour observer l'effet de la source interne, la température des fluides aux extrémités est mise égale car l'épaisseur de la cellule ne permettra pas de voir les différences de température à cet ordre de grandeur d'épaisseur : 

```{=latex}
\begin{table}[H]
\centering
\caption{Conditions aux limites — comparaison avec source interne}
\begin{tabular}{|l|c|c|c|}
\hline
 & \textbf{h (W/m²K)} & \textbf{T (°C)} & \textbf{q (W)} \\
\hline
Ouest & 10 & 40 & 0 \\
Est   & 10 & 40 & 0 \\
\hline
\end{tabular}
\end{table}
```

```{=latex}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{K_eq_S_eq_1_ep.png}
\caption{Profil avec conductivité équivalente, avec source interne - charge lente}
\end{figure}
```

Le modèle équivalent reproduit fidèlement le comportement global (température maximale et minimale) avec un écart inférieur à 0,01 °C sur $T_{max}$ dans les deux modes de charge.

**Limites du modèle équivalent :** Le profil de température interne à la cellule est lissé par l'homogénéisation : le modèle équivalent ne peut pas restituer les variations locales de température aux interfaces entre couches (gradient dans le séparateur notamment), ni les effets de concentration de flux aux interfaces de conductivités très différentes. Pour l'assemblage, où l'intérêt est la température macroscopique de l'ensemble, cette approximation est parfaitement acceptable.

### Validation de la conductivité équivalente 

Sans source interne et en régime permanent, le flux est constant dans tout le domaine.
On l'extrait depuis la face Ouest :

$$\varphi = h \cdot (T(0) - T_{ext,O}) = 10 \cdot (45{,}0245 - 40) = 50{,}245 \ W/m^2$$

La conductivité équivalente extraite numériquement vaut :

$$K_{eq}^{num} = \frac{\varphi \cdot L_{tot}}{T(L) - T(0)} 
= \frac{50{,}245 \times 250 \times 10^{-6}}{45{,}0323 - 45{,}0245} 
= \frac{1{,}256 \times 10^{-2}}{7{,}774 \times 10^{-3}} \approx 1{,}616 \ W/m \cdot K$$

Ici, T(L) et T(0) sont récupérées dans le fichier de sortie pour plus de précision.
Comparée à la valeur analytique obtenue par résistances en série :

$$K_{eq}^{ana} = \frac{L_{tot}}{\displaystyle\sum_j \frac{e_j}{\lambda_j}} = 1{,}59 \ W/m \cdot K$$

L'écart relatif est de $\epsilon \approx 1{,}6\%$, imputable à la position des nœuds 
aux centres des demi-mailles aux bords. La validation est concluante.


#  Objectif 3 — Assemblage de 50 cellules - Batterie complète

##  Modèle de l'assemblage

L'assemblage est modélisé comme un milieu homogène continu de même propriétés ($K_{eq}$, $S_{eq}$) sur l'épaisseur totale $L_{asm} = 50 \cdot 250 \mu m = 12,5 mm$. Le refroidissement est symétrique (convection sur les deux faces), avec un fluide à $T_f$ = 40 °C.

##  Détermination du coefficient h minimal

Les simulations sont réalisées pour différentes valeurs de h. Les critères de sécurité retenus sont (d'après l'annexe du sujet et Battery University BU-808) :

- Température interne recommandée en fonctionnement de charge : idéalement 15 à 35 °C
- Limite de conception conservatrice : $T_{int,max} \leq 45$ °C.
- Limite haute absolue en fonctionnement normal : $T_{int,max} \leq 60$ °C
- Seuil critique du matériau : shutdown du séparateur vers 120 à 135 °C

Pour ces différentes simulations, j'ai lancé le code Fortran en modifiant à chaque fois le type de charge (lente/rapide) ainsi que les coefficients h de chaque côté de la cellule. J'ai ensuite extrait la valeur de la température max dans la batterie.

```{=latex}
\begin{table}[H]
\centering
\caption{$T_{max}$ de l'assemblage 50 cellules en fonction de h ($T_{fluide}$ = 40°C)}
\begin{tabularx}{\textwidth}{|X|c|c|c|c|}
\hline
\textbf{h (W/m²K)} & \textbf{$T_{max}$ lente (°C)} & \textbf{$T_{max}$ rapide (°C)} & \textbf{Statut lente} & \textbf{Statut rapide} \\
\hline
5    & 42,88 & 68,69 & OK & $> 60$°C \\
10   & 41,45 & 54,48 & OK & $> 45$°C \\
20   & 40,74 & 47,39 & OK & $> 45$°C \\
50   & 40,31 & 43,13 & OK & OK \\
100  & 40,17 & 41,71 & OK & OK \\
200  & 40,10 & 40,99 & OK & OK \\
500  & 40,06 & 40,56 & OK & OK \\
1000 & 40,04 & 40,42 & OK & OK \\
\hline
\end{tabularx}
\end{table}
```

```{=latex}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{T=f(h).png}
\caption{Température maximale de l'assemblage vs coefficient d'échange h (charge lente et rapide)}
\end{figure}
```

**Résultat clé :** Pour la charge lente, le critère $T \leq 45$ °C est satisfait dès $h = 5 W/m^2 K$ (convection naturelle). Pour la charge rapide, il faut $h \geq 50 W/m^2 K$ pour satisfaire la limite de conception de 45°C - ce qui correspond à une convection forcée modérée (air ou liquide). En pratique industrielle, $h \approx 100-500 W/m^2 K$ est typique pour un refroidissement par plaque froide à eau, ce qui assure une marge de sécurité confortable même en charge rapide.

## Profils de température dans l'assemblage

```{=latex}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{sol_comparaison_plot.png}
\caption{Profils de température dans l'assemblage de 50 cellules — charge lente et charge rapide (h = 5 W/m²K)}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{comparaison_h=50.png}
\caption{Profils de température dans l'assemblage de 50 cellules — charge lente et charge rapide (h = 50 W/m²K)}
\end{figure}
```

Le profil de température est parabolique avec un maximum en cœur d'assemblage, conformément à la solution analytique d'un milieu homogène avec source uniforme et convection symétrique. La position du maximum au centre (x = $\frac{L_{asm}}{2}$) est cohérente avec la symétrie du problème.

# Objectif 4 — Résistance de contact entre cellules

## Modélisation de la résistance de contact

Dans un assemblage réel de batteries lithium-ion, les cellules ne sont pas en contact
parfait les unes avec les autres. L'interface entre deux cellules successives introduit
une **résistance thermique de contact** $R_c$ (en $m^2K/W$), due aux irrégularités de
surface, aux couches d'oxyde et à l'absence de continuité du matériau conducteur.
Cette résistance, souvent négligée en première approximation, peut devenir significative
dans un assemblage dense où les interfaces sont nombreuses. Elle est modélisée ici comme
une couche fictive homogène d'épaisseur $e_c = 1\ \mu m$ et de conductivité :

$$\lambda_c = \frac{e_c}{R_c}$$

Le choix de $e_c = 1\ \mu m$ est arbitraire : seul le rapport $e_c / \lambda_c = R_c$
a une signification physique. Cette épaisseur est retenue car elle est du même ordre
que les couches réelles de la cellule, ce qui assure la cohérence du maillage, et
suffisamment faible pour ne pas perturber l'épaisseur totale de l'assemblage
(50 interfaces $\times$ 1 µm = 50 µm, négligeable devant $L_{asm} = 12{,}5$ mm).
La conductivité associée est donc imposée par :

$$\lambda_c = \frac{e_c}{R_c} = \frac{1 \times 10^{-6}}{1 \times 10^{-3}} = 10^{-3}\ \text{W/m·K}$$

La valeur $R_c = 10^{-3}\ \text{m}^2\text{K/W}$ retenue est un ordre de grandeur
représentatif d'un contact sec entre surfaces métalliques, en accord avec les données
de la littérature [Incropera & DeWitt (Table 3.2, *Thermal contact conductance*)]. La source volumique de la couche fictive est nulle ($S_c = 0$),
l'interface mécanique entre cellules ne constituant pas une zone de génération de chaleur.

Cette couche est donc insérée dans la simulation entre chaque paire de cellules successives dans le maillage multicouche.


## Profil de température sur 2 cellules avec résistance de contact

```{=latex}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{contact_2_cell.png}
\caption{Profil de température sur 2 cellules avec résistance de contact $R_c$}
\end{figure}
```

Le profil de température met en évidence une **discontinuité de pente** aux interfaces
entre cellules. Physiquement, le flux thermique étant continu, la chute de température
à l'interface vaut :

$$\Delta T_{contact} = R_c \cdot \varphi$$

Plus $R_c$ est élevé, plus ce saut de température est marqué. À l'intérieur de chaque
cellule, le profil reste parabolique conformément aux résultats des objectifs précédents.

## Profil de température sur 50 cellules avec résistance de contact
```{=latex}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{50_cell_contact.png}
\caption{Profil de température sur 50 cellules avec résistance de contact $R_c$ sans source interne}
\end{figure}
```

```{=latex}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{Rc_50_Source.png}
\caption{Profil de température sur 50 cellules avec résistance de contact $R_c$ avec source interne}
\end{figure}
```

Sur l'assemblage complet de 50 cellules, les résistances de contact se manifestent sous
forme de **paliers de température** régulièrement espacés, superposés à l'enveloppe
parabolique globale. La température maximale en cœur d'assemblage est ainsi plus élevée
qu'en l'absence de résistances de contact, la résistance thermique totale de l'assemblage
étant augmentée par les 49 interfaces.

## Conductivité équivalente de l'assemblage avec résistances de contact

Il est possible de définir une conductivité équivalente globale de l'assemblage incluant
les résistances de contact, par analogie avec les résistances en série :

$$R_{tot}^{asm} = N_{rep} \cdot \frac{L_{cell}}{K_{eq}} + (N_{rep} - 1) \cdot R_c$$

$$\boxed{K_{eq}^{asm} = \frac{L_{asm}}{R_{tot}^{asm}} = \frac{N_{rep} \cdot L_{cell}}{N_{rep}
\cdot \dfrac{L_{cell}}{K_{eq}} + (N_{rep}-1) \cdot R_c}}$$

Ce milieu homogène équivalent corrigé permet de retrouver la température maximale globale
de l'assemblage sans avoir à résoudre explicitement les interfaces, au prix d'une perte
d'information sur les sauts locaux de température aux interfaces.

**Vérification des cas limites :** pour $R_c \to 0$, on retrouve $K_{eq}^{asm} \to K_{eq}$,
ce qui est cohérent avec l'absence d'interface. Pour $R_c \gg L_{cell}/K_{eq}$, la
conductivité de l'assemblage est dominée par les résistances de contact, conformément
au comportement physique attendu.

```{=latex}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{Rc_Keq_50_lin.png}
\caption{Profil de température sur 50 cellules — comparaison modèle explicite et conductivité équivalente avec $R_c$ avec $T_{fluide} = 40$ °C à l'Ouest et $T_{fluide} = 50$ °C à l'Est}
\end{figure}
```

```{=latex}
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{Rc_Keq_50_Seq.png}
\caption{Profil de température sur 50 cellules — comparaison modèle explicite et conductivité équivalente avec $R_c$ avec $T_{fluide} = 40$ °C de chaque côté}
\end{figure}
```

## Comparaison des températures maximales
```{=latex}
\begin{table}[H]
\centering
\caption{Impact de la résistance de contact sur $T_{max}$ de l'assemblage avec h = 10 $W/M^2K$ par côté}
\begin{tabular}{|l|c|c|}
\hline
\textbf{Configuration} & \textbf{$T_{max}$ lente (°C)} & \textbf{$T_{max}$ rapide (°C)} \\
\hline
Sans résistance de contact         & 41,25 & 54,48 \\
Avec $R_c$ (modèle explicite)      & 45,6 & 57,6 \\
\hline
\end{tabular}
\end{table}
```

La résistance de contact dégrade la conductivité thermique effective de l'assemblage et
élève la température maximale en cœur de batterie. Pour maintenir $T_{max} \leq 45$°C
en charge rapide, le coefficient d'échange $h$ minimal nécessaire est donc **plus élevé**
qu'en l'absence de résistances de contact. Ce résultat souligne l'importance de minimiser
$R_c$ en pratique, par exemple par l'utilisation de pâtes thermiques ou de plaquettes
conductrices aux interfaces entre cellules.


# Simulation bonus — Batterie réaliste de 2,5 cm d'épaisseur (100 cellules) :

À titre de comparaison, une simulation a été réalisée sur un assemblage de 100 cellules ($L_{asm} = 2{,}5$ cm) avec les cas extrêmes $h = 5\ W/m^2K$ et $h = 500\ W/m^2K$ en charge lente et rapide. 

Les résultats confirment les tendances observées : en charge rapide avec $h = 5\ W/m^2K$, $T_{max}$ atteint **97,9°C**, dépassant largement la limite absolue de 60°C et approchant le seuil de shutdown du séparateur (~120°C) — ce régime est donc incompatible avec la sécurité sans refroidissement forcé. 


\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{2,5cm_rapide_h5.png}
\caption{Profil de température — assemblage 2,5 cm, charge rapide, h = 5 W/m²K}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{2,5cm_lente_h5.png}
\caption{Profil de température — assemblage 2,5 cm, charge lente, h = 5 W/m²K}
\end{figure}


Avec $h = 500\ W/m^2K$, $T_{max}$ redescend à **41,6°C** en charge rapide et **40,17°C** en charge lente, confirmant qu'un refroidissement intensif permet de maîtriser thermiquement même un assemblage deux fois plus épais.



\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{2,5cm_rapide_h500.png}
\caption{Profil de température — assemblage 2,5 cm, charge rapide, h = 500 W/m²K}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{2,5cm_lente_h500.png}
\caption{Profil de température — assemblage 2,5 cm, charge lente, h = 500 W/m²K}
\end{figure}


#  Conclusion

##  Synthèse des résultats

Cette étude a permis de simuler avec succès le comportement thermique d'une cellule de batterie lithium-ion par volumes finis 1D. Les résultats principaux sont :

```{=latex}
\begin{table}[H]
\centering
\caption{Synthèse des réponses aux objectifs}
\begin{tabularx}{\textwidth}{|l|X|}
\hline
\textbf{Question posée} & \textbf{Réponse} \\
\hline
Cellule isolée — charge lente  & $T_{max}$ = 40,03°C — très faible $\Delta T$ dans la cellule ($<$ 0,001 K). La résistance convective domine. \\
Cellule isolée — charge rapide & $T_{max}$ = 40,29°C avec h = 10 W/m²K — conforme. Mais l'assemblage de 50 cellules requiert un refroidissement renforcé. \\
Conductivité équivalente $K_{eq}$ & 1,591 W/m·K (résistance en série dominée par le séparateur et les couches actives). \\
Assemblage — charge lente  & $h_{min}$ = 5 W/m²K (convection naturelle suffisante, $T_{max}$ = 42,9°C $<$ 45°C). \\
Assemblage — charge rapide & $h_{min}$ = 50 W/m²K (convection forcée nécessaire, $T_{max}$ = 43,1°C $<$ 45°C). En pratique $h \geq 100$ W/m²K recommandé. \\
Assemblage avec $R_c$ & Résistance de contact $R_c = 10^{-3}\ m^2K/W$ : $T_{max}$
passe de 41,25°C à 45,6°C (lente) et de 54,48°C à 57,6°C (rapide) avec $h = 10\
W/m^2K$. La conductivité effective de l'assemblage est réduite à
$K_{eq}^{asm} < K_{eq}$, nécessitant un $h_{min}$ plus élevé pour rester sous 45°C. \\
\hline
\end{tabularx}
\end{table}
```

##  Limites du modèle et pistes d'amélioration

Plusieurs limitations importantes doivent être soulignées :

- Le modèle est strictement stationnaire : les phases de charge transitoires (démarrage, fin de charge) ne sont pas capturées. Une extension en régime transitoire serait nécessaire pour évaluer les chocs thermiques.
- L'homogénéisation du milieu équivalent efface les détails intra-cellulaires. Pour des analyses de sécurité fine (surchauffe du séparateur), le modèle multicouche est indispensable.
- Les incertitudes sur les conductivités des électrodes actives (±30–40 %) se répercutent directement sur le $\Delta T$ interne à la cellule. Une étude de sensibilité paramétrique serait pertinente.
- La modélisation 1D ignore les effets de bord et les hétérogénéités latérales (inhomogénéités de courant), potentiellement importants pour les grandes cellules.
- La chaleur de polarisation (pertes électrochimiques) et la chaleur entropique ont été négligées. Leur inclusion pourrait modifier les résultats de 5–15 % selon le taux de charge.

**Conclusion générale :** Le modèle développé est suffisamment robuste pour dimensionner le système de refroidissement d'un assemblage de batteries en première approximation. La principale conclusion opérationnelle est que la charge rapide nécessite impérativement une convection forcée avec $h \geq 50 W/m^2 K$ pour maintenir l'assemblage sous 45°C, tandis que la charge lente peut être gérée avec une convection naturelle ou forcée légère.

