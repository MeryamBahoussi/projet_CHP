# Makefile adapt� � la construction de programmes Fortran90
CODE=SAMPO
############################################################################

UNAME_S:=$(shell uname -s)
ifeq ($(UNAME_S),Linux)
	UNAME=LINUX32e_${HOSTOSVERS}
	HOSTARCH=x86_64-intel-linux
endif
ifeq ($(UNAME_S),Darwin)
	UNAME=Darwin
	HOSTARCH=x86_64-apple-darwin18
endif
ifeq ($(UNAME_S),CYGWIN_NT-10.0)
	UNAME=Cygwin
	HOSTARCH=x86_64_cygwin
endif

PTMPDIR=/tmp/$(LOGNAME)
TMPDIR=$(PTMPDIR)/COMPIL/tmp

# pour utiliser mkl, petsc et/ou mutation++
mkl=no
petsc=no
mutationpp=no
gcc=yes
cgns=no
CC=gcc
############################################################################
include mkobj.lst
.SUFFIXES : .F .f .F90 .f90 .c .cpp .o

.PHONY : clean all $(OBJDIRS) $(OBJDIR)

DEFRULE=lien

ifeq ($(mpi),yes)
	MPIFLAGS = -lmpi
else
	MPIFLAGS =
endif

ifeq ($(gcc),yes)
	COMPILO=gcc
else
	COMPILO=intel
endif

HELP=ifort -V

ifeq ($(gcc),yes)
	F90 = mpif90
else
	F90 = mpiifort
endif
F90LIN = $(F90)
ifeq ($(mkl),no)
else
	ifeq ($(gcc),yes)
		MKL_INC=-I${MKL_ROOT}/include
	else
                # pour le compilateur Intel pas besoin de sp�cifier le chemin de la MKL
		MKL_INC=
	endif
	MKL_LDFLAGS=
	MKL_LDADD=-lmkl_core -lmkl_intel_lp64 -lmkl_sequential -liomp5
	MKL_FLAGS=-DMKL_ACTIF
endif
ifeq ($(petsc),yes)
	PETSCFLAGS= -DPETSC_ACTIF
	PETSC_LDADD=-lpetsc
	PETSC_INCL=-I${PETSC_INC} -I${PETSC_DIR}/include
endif
ifeq ($(gcc),yes)
        # Librairies optimis�es Intel pour compilation avec gcc
	LDADD +=
endif

ifeq ($(mutationpp),no)
else
	MUTATION_DIR=/usr/local
	MUTATION_INC=-I$(MUTATION_DIR)/include -I$(MUTATION_DIR)/include/mutation++ -I$(MUTATION_DIR)/include/eigen3
	MUTATION_LDFLAGS = -L$(MUTATION_DIR)/lib
	MUTATION_FLAGS = -DMPP_ACTIF
	MUTATION_LDADD = -lmutation++ -lstdc++
endif
ifeq ($(cgns),yes)
	#CGNS_INC=-I$(CGNS_INC)
	CGNS_LDFLAGS = -L$(CGNS_LIB)
	CGNS_FLAGS = -DCGNS_ACTIF
	CGNS_LDADD = -lcgns
endif

CPPINC = ${MUTATION_INC}
F90INC = ${PETSC_INCL} ${MKL_INC} ${CGNS_INC}
LDADD += ${MKL_LDADD} ${PETSC_LDADD} $(MUTATION_LDADD) ${CGNS_LDADD}
LDFLAGS += ${PETSCFLAGS} $(PETSC_LDFLAGS) $(MKL_LDFLAGS) $(MUTATION_LDFLAGS) ${LDADD} ${CGNS_LDFLAGS}
F90FLAGS = $(PETSCFLAGS) $(MUTATION_FLAGS) $(MKL_FLAGS) ${CGNS_FLAGS}
CPPFLAGS = $(MUTATION_FLAGS)

# Pour r�cup�rer le nom de la branche GIT
GIT_CODE_PROD:=$(shell git branch | grep "\*" | sed -e "s%*%%g" | sed -e "s% %%g")
GITHASH_CODE_PROD:=$(shell git log --graph | head -1 | cut -d' ' -f3)

F90FLAGS += -DGIT="'${GIT_CODE_PROD}'" -DGITHASH="'${GITHASH_CODE_PROD}'" $(MPIFLAGS) $(PETSCFLAGS) $(MKL_FLAGS)

#------- MODE DEBUG ? --------
ifeq ($(debug),yes)
	OBJDIR = $(UNAME)_g
	ifeq ($(gcc),yes)
		F90FLAGS += -gdwarf-2 -cpp -ffree-line-length-none -fdollar-ok -Wall -Wextra -Wno-unused-parameter -Wno-unused-dummy-argument -Wno-compare-reals -Wno-tabs
		CPPFLAGS += -gdwarf-2
	else
		F90FLAGS += -gdwarf-2 -DDEBUG_CODE -traceback -fpp #-fp-model precise
		CPPFLAGS += -g
	endif
else
	OBJDIR = $(UNAME)
	ifeq ($(gcc),yes)
                # Mis en forme du code
		F90FLAGS += -cpp -ffree-line-length-none -fdollar-ok               # Optimisations
		F90FLAGS +=
	else
		OPTIM_F90FLAGS = -O2  $(OPT_ARCH) -fpp -fp-model precise
		NOPTIM_F90FLAGS = -O0 -fpp
	endif
endif
# Compatibilite pour initmak.pl a la place de depmod:
TYPEREEL = 5
ifndef TMPDIR
	TMPDIR=/tmp/$(LOGNAME)
endif
ifndef PTMPDIR
	PTMPDIR=$(TMPDIR)
endif
SECTION = $(CODE)/${GIT_CODE_PROD}
# Finalement on compile sur le /tmp
TMPLIB = $(TMPDIR)/$(SECTION)

# Creation de la chaine OUT si non specifiee pour differencier
# les compilations avec ou sans MKL, PETSC et le mode debug ou pas:
ifeq ($(debug),yes)
	DOUT ?=_
	DOUT:=${DOUT}g
	VARIANTE=Debug
else
	VARIANTE=Release
endif

# Changement de logique, on fait toujours la MKL
ifeq ($(mkl),no)
else
	DOUT ?=_
	DOUT:=${DOUT}m
endif
ifeq ($(petsc),yes)
	DOUT ?=_
	DOUT:=${DOUT}p
endif
ifeq ($(gcc),yes)
	DOUT ?=_
	DOUT:=${DOUT}gcc
	VARIANTE+=gcc
else
	VARIANTE+=intel
endif
ifdef DOUT
	OUT ?= $(DOUT)
endif

OBJDIR = $(TMPLIB)/objets$(OUT)$(TYPEREEL)/$(HOSTARCH)
BIBLIO = $(LDFLAGS)
LD=$(F90)
EXEDIR=$(UNAME)$(OUT)
F90FLAGS += -DVARIANTE="'${VARIANTE}'"
LIBUDF=$(BIN)libs

ifeq ($(gcc),yes)
	MODCOMMAND = -I$(OBJDIR)/modules -J$(OBJDIR)/modules
else
	MODCOMMAND = -module $(OBJDIR)/modules
endif

# binaire et librairie du projet, BIN peut etre choisi par l'utilisateur (setenv ou export BIN)
BIN ?= ../bin/$(EXEDIR)/${GIT_CODE_PROD}/
TARGET_BIN = $(BIN)a.$(CODE)$(OUT)$(TYPEREEL)_$(HOSTARCH)
TARGET_LIB = $(OBJDIR)/lib/lib$(CODE).a

# regle par defaut: ne pas installer le produit en meme temps qu'on compile:
all : $(TARGET_LIB) $(TARGET_BIN)
	@echo "done for $@ on $<"

# sources de la librairie
SRCDIRS = src/MAIN src/FLUIDE src/MAILLAGE src/MPI src/SOLVEUR src/SOLIDE src/INTERFACE src/IO

SUFFIXES = .F90 .f90 .F .f .cpp

PROJ_SRC = $(foreach srcdir, $(SRCDIRS), $(foreach suf, $(SUFFIXES), $(wildcard $(srcdir)/*$(suf))))

LIB_SOURCES = $(PROJ_SRC)

OBJDIRS = $(foreach SRCDIR, $(SRCDIRS), $(addprefix $(OBJDIR)/, $(SRCDIR))) \
$(OBJDIR)/modules $(OBJDIR)/lib $(OBJDIR)/bin $(BIN)

# nom des objets de la librairie.
LIB_BASES = $(basename $(LIB_SOURCES))
LIB_OBJECTS = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(LIB_BASES)))

# Liste des fichiers non-optimis�s
NOPTIM_LIST =
#NOPTIM_LIST += src/SOLVEUR/solveur_petsc_global.f90
#NOPTIM_LIST += src/MAILLAGE/m_maillage.f90

#NOPTIM_LIST += src/FLUIDE/acoustique.f90 src/FLUIDE/calcul_fluide.f90 src/FLUIDE/conditions_limites_euler.f90 src/FLUIDE/conditions_limites_implicite.f90 src/FLUIDE/conditions_limites_visqueux.f90 src/FLUIDE/decode.f90 src/FLUIDE/flux_EUCCLHYD.f90 src/FLUIDE/flux_Gallice.f90 src/FLUIDE/flux_HLL.f90 src/FLUIDE/implicite_energie_1d.f90 src/FLUIDE/implicite_energie.f90
#NOPTIM_LIST += src/FLUIDE/flux_Gallice.f90 src/FLUIDE/initialisation_fluide.f90
#NOPTIM_LIST += src/FLUIDE/lagrange_implicite_chalons.f90 src/FLUIDE/lagrange_implicite_direct.f90 src/FLUIDE/m_EOS.f90 src/FLUIDE/m_transport.f90 src/FLUIDE/outils_conditions_limites.f90 src/FLUIDE/outils_pas_acoustique.f90 src/FLUIDE/premier_euler.f90 src/FLUIDE/premier_euler_total.f90 src/FLUIDE/premier_visqueux.f90 src/FLUIDE/second_energie.f90 src/FLUIDE/second_euler.f90 src/FLUIDE/second_forces_ext.f90 src/FLUIDE/second_visqueux.f90
#NOPTIM_LIST += src/INTERFACE/deplacement_paroi.f90  src/INTERFACE/itf_bilans.f90 src/INTERFACE/itf_CL.f90 src/INTERFACE/itf_fusion.f90 src/INTERFACE/itf_initialisation.f90 src/INTERFACE/itf_parametres.f90 src/INTERFACE/itf_sublimation.f90 src/INTERFACE/itf_synchronisation.f90 src/INTERFACE/outils_itf.f90
#NOPTIM_LIST += src/IO/input.f90 src/IO/m_gestion_IO.f90 src/IO/residu.f90 src/IO/sorties.f90
#NOPTIM_LIST += src/MAILLAGE/m_maillage.f90 src/MAILLAGE/m_raffinement_maillage.f90 src/MAILLAGE/outils_maillage.f90
#NOPTIM_LIST += src/MPI/gestion_interbloc_et_patch.f90 src/MPI/m_communication_MPI.f90 src/MPI/m_MPI.f90 src/MPI/synchro_interbloc.f90
#NOPTIM_LIST += src/SOLIDE/conditions_limites_solide.f90 src/SOLIDE/initialisation_solide.f90 src/SOLIDE/premier_thermique.f90 src/SOLIDE/second_thermique.f90 src/SOLIDE/thermique_solide.f90
#NOPTIM_LIST += src/SOLVEUR/mkl_dss.f90 src/SOLVEUR/m_outils_CSR.f90 src/SOLVEUR/m_PETSC.f90 src/SOLVEUR/solveur_1d.f90 src/SOLVEUR/solveur_balayage.f90 src/SOLVEUR/solveur_MKL.f90 src/SOLVEUR/solveur_petsc_global.f90


# Faire le fichier version:
initversion:
	@if test ! -f $(OBJDIR)/Version.o ; then \
	$(MAKE) theversion ; \
	fi
	@if test ! -f $(OBJDIR)/modules/version.mod ; then \
	$(MAKE) theversion ; \
	fi

# regle de construction des repertoires
$(OBJDIRS) $(OBJDIR) :
	@if test ! -d $@ ; then \
	echo "mkdir" $@ ; \
	mkdir -p $@ &&  touch $@/.stamp; \
	fi

infos:
	echo "UNAME=$(UNAME), ARCH=$(ARCH), MACHTYPE=$(MACHTYPE), OUT=$(OUT), F90=$(F90), TARGET_BIN=$(TARGET_BIN)"

# Installation du produit, a mettre en coherence avec le script:
install:
	@echo "Installing $(CODE) on $(prefix) release $(TAG): $(prefix)/a.$(CODE)$(OUT)_$(UNAME)_$(SVN)"
	cp -f --preserve=timestamps $(TARGET_BIN) $(prefix)/a.$(CODE)$(OUT)_$(UNAME)_$(TAG)_$(VERSION)
	chmod a+rX $(prefix)/a.$(CODE)$(OUT)_$(UNAME)_$(TAG)_$(VERSION)
#	mkdir -p $(prefix)
#	cp -pf $(TARGET_BIN) $(prefix)/a.$(CODE)_$(OBJDIR)_$(SVN)

lien: all
	@echo "Linking $(CODE) from $(prefix) release $(TAG): $(prefix)/a.$(CODE)$(OUT)_$(UNAME)_$(SVN)"
	ln -sf $(TARGET_BIN) a.$(CODE)$(OUT)

# r�gles de construction des binaires.
$(TARGET_BIN) : $(BIN)a.$(CODE)$(OUT)$(TYPEREEL)_$(HOSTARCH)

# Faire le fichier de version a la main:
theversion:
	perl ./mkversion.pl -c '$(OBJDIR)' -i 'ifort -V' -h 'code_diphasique_$(HOSTARCH)' -v '$(VERSION)' -d '$(DATE)' -t '$(TAG)'
	@echo "Compilation du fichier Version:"
	$(F90) $(MODCOMMAND) -c $(F90FLAGS) $(OPTIM_F90FLAGS) $(F90INC) -o $(OBJDIR)/Version.o $(OBJDIR)/Version.F90

# r�gle de construction de la librairie: inutile...
$(TARGET_LIB) :	depend $(OBJDIRS) $(LIB_OBJECTS)
	ar -ruv $@ $(LIB_OBJECTS)

$(OBJDIR)/%.o : %.f90
# 	@if echo $(NOPTIM_LIST) | grep -q $< ; then \
# 		echo "Compilation NON OPTIMISEE du fichier: " $< ; \
# 		$(F90) $(MODCOMMAND) -c  $(F90INC) $(F90FLAGS) $(NOPTIM_F90FLAGS) -o $@ $< ; \
# 	else  \
# 		echo "Compilation du fichier: " $< ; \
# 		$(F90) $(MODCOMMAND) -c $(F90INC) $(F90FLAGS) $(OPTIM_F90FLAGS)  -o $@ $< ; \
# 	fi
	@echo "Compilation du fichier: " $< ; \
	$(F90) $(MODCOMMAND) -c $(F90INC) $(F90FLAGS) $(OPTIM_F90FLAGS)  -o $@ $<

$(OBJDIR)/%.o : %.F90
	$(F90) $(MODCOMMAND) -c $(F90FLAGS) $(OPTIM_F90FLAGS) $(F90INC) -o $@ $<

$(OBJDIR)/%.o : %.f
	$(F90) $(MODCOMMAND) -c $(F90FLAGS) $(OPTIM_F90FLAGS) $(F90INC) -o $@ $<

$(OBJDIR)/%.o : %.F
	$(F90) $(MODCOMMAND) -c $(F90FLAGS) $(OPTIM_F90FLAGS) $(F90INC) -o $@ $<

$(OBJDIR)/%.o : %.c
	$(CC) -c $(CFLAGS) -I./includes -o $@ $<

$(OBJDIR)/%.o : %.cpp
	$(CC) -c $(CPPFLAGS) $(CPPINC) -o $@ $<

# r�gles diverses (d�pendences, nettoyage).
ALL_FORTRAN = $(LIB_SOURCES)

clean :
	\rm -f $(TARGET_LIB) \
	$(LIB_OBJECTS) $(OBJDIR)/Version.o \
	$(OBJDIR)/modules/*.mod

#dependance : mkdep.lst mkobj.lst
depend : mkdep.lst mkobj.lst

#  $(ALL_FORTRAN) initmak.pl $(OBJDIR)
mkdep.lst : $(ALL_FORTRAN) initmak.pl $(OBJDIR)
	@echo -e "\n Running LEGACY dependency analysis for $(ALL_FORTRAN)"
	perl ./initmak.pl -d -f -R -L '$(SRCDIRS)'

depfast:
	@echo "Running FAST dependency analysis"
	@perl ./initmak.pl -f -R -L '$(SRCDIRS)'

#inclusion des d�pendances.
include mkdep.lst
