#!/bin/bash

function echoHelp
{
    echo "
╔═════════════════════════════════════════════════════
║ Script qui permet d'utiliser le code
║ 4 modes d'utilisation possibles
║ - compilation du code
║ - lancement d'un calcul en interactif
║ - soumission d'un calcul en batch sur une machine distante
║ - lancement de la non regression
║
║══ Pour la compilation (par defaut avec le compilateur ifort)
║ ./lancement.sh -compile
║ options possibles
║ -clean     [ ]      : pour nettoyer la compilation
║ -gcc       [ ]      : utilisation du compilateur gcc
║ -debug     [ ]      : mode debuggage
║ -totalview [ ]      : active aussi le mode debug
║
║══ Pour lancer un cas D sur N procs. Les entrees doivent etre dans le dossier cas_test/CALCULS/D
║ ./lancement.sh -d D -n N
║ options possibles
║ -d         [defaut] : nom du repertoire du cas test 
║ -n         [1]      : nombre de processeurs pour lancer le calcul
║ -omp       [1]      : lancement avec n threads openMP par processus MPI
║ -gcc       [ ]      : utilisation du binaire compile avec le compilateur gcc
║ -debug     [ ]      : utilisation du binaire compile en mode debug
║ -binaire   [ ]      : utilisation d'un binaire spécifique (par defaut on prend celui de la branche courante GIT)
║ -totalview [ ]      : active le mode debug et lance le debugger totalview à l'execution
║ -valgind   [ ]      : pour lancer valgrind
║
║══ Pour soumettre un calcul en batch sur un cluster distant C
║ ./lancement.sh -C [-d][-n][-time] 
║ options possibles
║ -d         [defaut] : nom du repertoire du cas test 
║ -n         [1]      : nombre de processeurs pour lancer le calcul
║ -gcc       [ ]      : compilateur gcc (à la compilation ou l'execution)
║ -debug     [ ]      : utilisation du binaire compile en mode debug
║ -binaire   [ ]      : utilisation d'un binaire spécifique (par defaut on prend celui de la branche courante GIT)
║ -quiet     [ ]      : pour ne pas interagir avec l'utilisateur (utile pour la non reg)
║
║══ Pour lancer la non regression et la comparaison avec la version V sur le scenario S
║ ./lancement.sh -non_reg V S
╚═════════════════════════════════════════════════════
"
}

#####################################
#  Initialisation des options a la valeur par defaut
#####################################
usage=2  # 1:compilation 2:lancement 3:batch 4:non reg 5:aide
Nproc=1
Nomp=1
repertoire=defaut
debug=0
totalview=0
valgrind=0
clean=0
gcc=1
binaire_fourni=0
tlimit_fourni=0
tlimit=18000
mlimit=3000
quiet=0

#   Initialisation des variables locales
local_path=$PWD
code=SAMPO

# machine sur laquelle est lance le script
machine=`echo $HOSTNAME | tr -d '[0-9]'` 
#version de l'os sur lequel on lance le script ${HOSTOSVERS}

#####################################
#  Analyse des arguments
#####################################
if (($# != 0)) ; then
    while [[ "${1:-vide}" != vide ]]
    do
	case $1 in
	    "-d")
		shift
		repertoire=$1;;
	    "-debug")	    
		debug=1;;
	    "-n")
		shift
		Nproc=$1;;
	    "-omp")
		shift
		Nomp=$1;;
	    "-compile")
		usage=1;;
	    "-clean")
		clean=1;;
	    "-totalview")
		totalview=1
		debug=1;;
	    "-valgrind")
		valgrind=1;;
	    "-gcc")
		gcc=1;;
	    "-time")
		tlimit_fourni=1
		shift
		tlimit=$1;;
	    "-quiet")
		quiet=1;;
	    "-non_reg")
		usage=4
		shift
		version=$1
		shift
		target=$1;;
	    "-binaire")
		binaire_fourni=1
		shift
		binaire=$1;;
	    "-h" | "-help")
		usage=5;;
	    *)
		echo "option non prevue "$1;;
	esac
	shift
    done
fi

# positionnement de l'environnement du code
function sourceEnvCode 
{
    if [ "$ENVIRONMENT" != "BATCH" ] ; then
	cd ..
	source ./env_code_$USER.sh
	cd -
	echo ""
	# deja fait dans le fichier batch
    fi
}

# Initialisation des differents chemins en fonction de la machine et du mode de lancement
function definePath()
{
    if [ "$ENVIRONMENT" == "BATCH" ]; then
	echo "on est sur $machine en batch"
	path_cas_test=cas_test
	path_src=$PWD
	path_bin=$PWD
    else # local ou xxx interactif
	path_cas_test=../cas_test/CALCULS
	path_src=../sources
	path_bin=$PWD/../bin
    fi
    path_cas_test=$path_cas_test/$repertoire
}

# recuperation du binaire a utiliser
function getBinaire()
{
    if [[ $binaire_fourni = 0 ]]; then
	optBin="_"
	if (( debug==1 )) ; then
	    optBin=${optBin}g
	fi
	# on fait toujours sans mkl et avec petsc
	# sans petsc ici
	optBin=${optBin}
	if (( gcc==1 )) ; then
	    optBin=${optBin}gcc
	fi
	UNAME_S="`uname -s`"
	if [ "$UNAME_S" == "Linux" ]; then
	    dirBin=LINUX32e_${HOSTOSVERS}$optBin
	elif [ "$UNAME_S" == "Darwin" ]; then
	    dirBin=Darwin$optBin	    
	fi
	TYPEREEL='5'
	namBin=a.$code$optBin${TYPEREEL}_$HOSTARCH
	currentBranch=`git branch | grep \* | cut -d ' ' -f2`
	binaire=$path_bin/$dirBin/$currentBranch/$namBin
    fi
}

function compilation()
{
    echo "
######################################################
#   Compilation du code
######################################################
"       
    cd $path_src
    optMake="-j8 MACHTYPE=$machine"
    if (( debug==1 )) ; then
	optMake="${optMake} debug=yes"
    fi
    if (( gcc == 1 )) ; then
	optMake="${optMake} gcc=yes"
    else
	optMake="${optMake} gcc=no"
    fi
    if (( clean==1 )) ; then
	optMake="${optMake} clean"
    fi
    make $optMake
    echo "Compilation terminee"
}

function nonReg()
{
    echo "
######################################################
#   Lancement de la non regression et de la comparaison
######################################################
"
    cd ../../PrivateTest/Prod/
    ./run.sh $version toto $target
    cd ../Comparaison/
    ./comparaison.sh $version -2 reference $target
}


# Lancement d'un calcul
function lancementCalcul
{
    # test pour savoir si le dossier existe
    if [ ! -d $path_cas_test ]; then
	echo -e "Erreur : le dossier de calcul $path_cas_test n'existe pas ! \n"; exit
    fi
    # test pour savoir si le binaire existe
    if [ ! -f  $binaire ]; then
	echo -e "Erreur : le binaire $binaire n'existe pas ! \n"; exit
    fi
    # Copie du binaire
    cp $binaire $path_cas_test/bin 
    
    #  Copie et archivage des sources et choix de la commande d'execution
    if [ "$ENVIRONMENT" == "BATCH" ] ; then
	echo "on ne copie pas les sources en batch" 
	F90=ccc_mprun
    else
	tar -czf  $path_cas_test/src.tar.gz $path_src/src
	echo "copie des sources terminee" 
	F90=mpirun
    fi

    # Lancement
    cd $path_cas_test
    echo "pwd: $PWD"
    echo -e 'On utilise le binaire ' $binaire '\n'
    
    if (( totalview==1 )) && (( valgrind==1 )) ; then
	echo -e "Erreur : options -totalview et -valgrind incompatibles"
	exit
    fi
    
    if (( valgrind==1 )) ; then
	valgrind --error-limit=no --leak-check=full ./bin 2>&1 | tee ./execution.dat
    elif ((totalview==1)) ; then
	totalview -mpi "Open MPI" -np $Nproc ./bin
    else #lancement normal
	OMP_NUM_THREADS=$Nomp
	$F90 -n $Nproc --oversubscribe ./bin 2>&1 | tee ./execution.dat
	rm ./bin
    fi        
}

case $usage in
    1) # Compilation du code
	sourceEnvCode
	definePath
	compilation	
	;;   
    2) # lancement du code
	sourceEnvCode
	definePath
	getBinaire
	lancementCalcul	
	;;
    3) # Soumission d'un calcul en batch
	sourceEnvCode
	definePath
	getBinaire
	;;
    4) # Lancement de la non regression et de la comparaison
	nonReg ;;
	
    5) # affiche l'aide
	echoHelp ;;
esac
