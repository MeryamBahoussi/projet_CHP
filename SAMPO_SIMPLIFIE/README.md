
#oui

# SAMPO
Simulation d'Ablation Multi-Phasique avec des écoulements Compressibles

Simulation of Ablation with MultiPhase cOmpressible flows

## Contexte physique
Code d'étude pour l'ablation liquide ou l'ablation multi-espèces.
Permet de calculer des configurations multi-physiques avec des domaines fluide et solide couplés par un processus de changement de phase à la paroi.

### Ablation Liquide 
Dans le cas de l'ablation liquide, le domaine fluide est constitué d'un écoulement diphasique, à savoir un gaz et un liquide.
Le liquide est injecté à la paroi suite à la fusion du domaine solide.
Pour modéliser l'écoulement diphasique, un modèle à interface diffuse est utilisé. Dans la mesure où les deux phases sont considérées comme non miscibles, le modèle à 5 équations d'Allaire et al. ou de Massoni et al. est utilisé.

Pour résoudre numériquement ce modèle, une décomposition d'opérateurs est utilisée. Cela permet de séparer les phénomènes liés à l'acoustique (plus les effets visqueux, thermiques, capillaires) de la partie convective des équations.
Des schémas volumes finis sont utilisés pour les deux parties. Ils s'appuient sur des schémas de type Godunov basés sur des solveurs de Riemann approchés.
Un schéma implicite peut être employé pour la partie acoustique tandis que la partie transport est toujours traitée en explicite.


### Ablation Multi-espèces
Dans le cas de l'ablation multi-espèces, le domaine fluide est constitué d'un mélange d'espèces gazeuses. Des équations de convection-diffusion sont utilisées pour modéliser l'écoulement. Les termes de production d'espèces chimiques ne sont pas pris en compte, ce qui correspond donc au régime de chimie figée dans l'écoulement gazeux.

Pour résoudre numérique ce modèle, un schéma volumes finis direct est utilisé. Une intégration temporelle implicite est réalisée pour s'affranchir des conditions de stabilité des effets dissipatifs. 

### Résolution du domaine solide
Dans tous les cas, l'équation de la chaleur est résolue dans le domaine solide. Les propriétés du matériaux, à savoir sa densité, sa conductivité thermique ainsi que sa capacité calorifique, sont supposées ne pas dépendre de la température. 
Un schéma volumes finis implicite est utilisé dans cette partie.

## Compilation du code
Il faut se placer dans le dossier lancement/ puis exécuter la commande suivante
./lancement.sh -compile
 options possibles
 -clean     [ ]      : pour nettoyer la compilation
 -gcc       [ ]      : utilisation du compilateur gcc
 -debug     [ ]      : mode debuggage

## Exécution du code
Pour utiliser le code, il faut créer l'arborescence suivante
./cas_test/CALCULS

Les fichiers de paramètres pour un cas test *toto* sont à renseigner dans le répertoire cas_test/CALCULS/toto/entrees/

Pour exécuter le cas test *toto* sur *N* processeurs avec le code, il faut lancer la commande suivante
./lancement.sh -d toto -n N
 options possibles
 -d         [defaut] : nom du repertoire du cas test 
 -n         [1]      : nombre de processeurs pour lancer le calcul
 -gcc       [ ]      : utilisation du binaire compile avec le compilateur gcc
 -debug     [ ]      : utilisation du binaire compile en mode debug
 -binaire   [ ]      : utilisation d'un binaire spécifique (par défaut on prend celui de la branche courante GIT)
