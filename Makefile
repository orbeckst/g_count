# $Id: Makefile,v 1.30 2004/03/18 01:04:28 oliver Exp $
#
#
# Makefile for the compilation of g_count and relatives

# Makefile switches:
# Switch on on commandline by setting them, eg 
#        make DEBUG=1 WITH_MPI=1 g_count 
#
#  DEBUG        compile with debugging, -g and no optimisations
#  WITH_MPI     use mpi compiler instead of gcc
#               (set CC and LD in the 'ifdef WITH_MPI' section)       

# this should be changed to /opt/gromacs/gromacs-3.0 in the future
# synapse
GMX_TOP_DIR        := /sansom/fedpacks
MY_GMX_TOP_DIR     := /sansom/gfio/oliver/Gromacs/Gromacs-3.1.4

# EXEC depends on your machine/OS
GMX_EXEC_PREFIX := $(GMX_TOP_DIR)#
GMX_LIB_DIR     := $(GMX_EXEC_PREFIX)/lib#
#GMX_INCLUDE_DIR := $(GMX_TOP_DIR)/include#
GMX_INCLUDE_DIR := $(GMX_TOP_DIR)/include/gromacs#

# this is only necessary for the creation of etags
# (for compilation it is not important)
GMX_SOURCE_DIR  := /sansom/fedpacks/src/BUILD/gromacs-3.2.1

# install binaries to
BIN_DIR := $(HOME)/bin
#BIN_DIR := $(MY_GMX_TOP_DIR)/i686-pc-linux-gnu/bin

CPPFLAGS += -I$(GMX_INCLUDE_DIR)

ifdef DEBUG
CFLAGS   +=  -DDEBUG -g -Wall  -Wno-unused 
else
CFLAGS   +=  -g -O2 -fomit-frame-pointer -finline-functions \
             -funroll-loops -Wall -Wno-unused 
endif

# static (but g_ri3Dc segfaults at the end)
#LDFLAGS  +=  -lm -Wl,-static -L$(GMX_LIB_DIR) -lmd -lgmx 

# dynamic linking of libgm, require LD_LIBRARY_PATH
#LD_LIBRARY_PATH += /usr/src/gm-1.6.3_Linux/binary/lib/
# for libgm.so on synapse

LDFLAGS  +=  -lm -L$(GMX_LIB_DIR) -lmd -lgmx 

ifdef WITH_MPI
  # on synapse.biop.ox.ac.uk, gmx linked with mpi stuff:
  MPI_D := /usr/local/mpich-gm_GNU
  CC := $(MPI_D)/bin/mpicc
  LD := $(CC)
else
  CC := gcc
  LD := $(CC)
endif

INSTALL := install

AUX_NAMES :=  utilgmx xf count plt
AUX_SRC := $(addsuffix .c, $(AUX_NAMES))
AUX_H   := $(addsuffix .h, $(AUX_NAMES))
AUX_OBJ := $(addsuffix .o, $(AUX_NAMES))

# g_count
G_COUNT     := g_count
G_COUNT_SRC := g_count.c 
G_COUNT_H   := $(AUX_H)
G_COUNT_OBJ := g_count.o 

# g_zcoord
G_ZCOORD     := g_zcoord
G_ZCOORD_SRC := g_zcoord.c 
G_ZCOORD_H   := $(AUX_H)
G_ZCOORD_OBJ := g_zcoord.o 

# g_flux
G_FLUX     := g_flux
G_FLUX_SRC := g_flux.c 
G_FLUX_H   := $(AUX_H)
G_FLUX_OBJ := g_flux.o 


GRID_SRC := grid3D.c xdr_grid.c 
GRID_H   := grid3D.h xdr_grid.h
GRID_OBJ := grid3D.o xdr_grid.o

# g_ri3Dc
G_RI3DC     := g_ri3Dc 
G_RI3DC_SRC := g_ri3Dc.c 
G_RI3DC_H   := $(AUX_H)  $(GRID_H) 
G_RI3DC_OBJ := g_ri3Dc.o $(GRID_OBJ)

# a_ri3Dc
A_RI3DC     := a_ri3Dc 
A_RI3DC_SRC := a_ri3Dc.c 
A_RI3DC_H   := $(AUX_H)  plt.h $(GRID_H)  
A_RI3DC_OBJ := a_ri3Dc.o plt.o $(GRID_OBJ)

DIST_GRIDCOUNT_H   := $(sort  $(A_RI3DC_H) $(G_RI3DC_H) $(GRID_H) $(AUX_H))
DIST_GRIDCOUNT_SRC := $(sort  $(A_RI3DC_SRC) $(G_RI3DC_SRC) $(GRID_SRC) $(AUX_SRC))
DIST_GRIDCOUNT := gridcount.tar.gz 

# template: g_xx
G_XX     := g_xx
G_XX_SRC := g_xx.c 
G_XX_H   := $(AUX_H)
G_XX_OBJ := g_xx.o 

ALL_PROG := $(G_COUNT) $(G_FLUX) $(G_ZCOORD) $(G_RI3DC) $(A_RI3DC)
ALL_SOURCES := $(G_COUNT_SRC) $(G_COUNT_H) \
               $(G_FLUX_SRC) $(G_FLUX_H) \
               $(G_ZCOORD_SRC) $(G_ZCOORD_H) \
	       $(G_RI3DC_SRC) $(G_RI3DC_H) \
	       $(A_RI3DC_SRC) $(A_RI3DC_H) \
	       $(AUX_SRC) $(AUX_H)

define usage
\nIn order to compile programs edit the Makefile for paths to the Gromacs\
\nlibraries. Then do \
\n   make clean; make PROGRAM\
\nwhere PROGRAM can be one of \`$(ALL_PROG)'.\
\nPerhaps you have to edit variables at the top of the Makefile \
\nto make it work. Programs are statically linked so you can take the \
\nbinaries wherever you like (hopefully...).\
\n\
\nIn order to install all programs that you compiled, type\
\n   make BIN_DIR=<your_target_dir> install\
\nor set BIN_DIR at the top of the Makefile.\
\n\
\nInstall targets: \
\n   all          compile \`$(ALL_PROG)' (default)\
\n   install      install all compiled programs in BIN_DIR\
\n   install-grid for \`$(G_RI3DC) $(A_RI3DC)'\
\n   dist-grid    creates a distribution tarball $(DIST_GRIDCOUNT)\
\n   clean        clean object files etc\
\n   distclean    remove every generated file\
\n\
\nMakefile switches:\
\nSwitch on on commandline by setting them, eg \`make WITH_MPI=1 g_count' \
\n  DEBUG        compile with debugging, -g and no optimisations\
\n  WITH_MPI     use mpi compiler instead of gcc\
\n               (set CC and LD in the 'ifdef WITH_MPI' section)\
\n\
\nUse \`make doxygen' to generate documentation in\
\n       /sansom/public_html/localinfo/dev/g_count
endef
# 'emacs font-lock

.PHONY: all help doxygen
all:	$(ALL_PROG)

help:
	@echo -e "$(usage)"

doxygen: doxygen.config
	doxygen doxygen.config

$(AUX_OBJ): $(AUX_SRC) $(AUX_H)



$(G_COUNT): $(G_COUNT_OBJ) $(AUX_OBJ) 
	$(LD) -o $@ $^ $(LDFLAGS)
$(G_COUNT_OBJ): $(G_COUNT_SRC) $(G_COUNT_H)



$(G_ZCOORD): $(G_ZCOORD_OBJ) $(AUX_OBJ) 
	$(LD) -o $@ $^ $(LDFLAGS)
$(G_ZCOORD_OBJ): $(G_ZCOORD_SRC) $(G_ZCOORD_H)


$(G_FLUX): $(G_FLUX_OBJ) $(AUX_OBJ) 
	$(LD) -o $@ $^ $(LDFLAGS)
$(G_FLUX_OBJ): $(G_FLUX_SRC) $(G_FLUX_H)



#.phony: grid
#grid: $(G_RI3DC) $(A_RI3DC)

$(G_RI3DC): $(G_RI3DC_OBJ) $(AUX_OBJ) 
	$(LD) -o $@ $^ $(LDFLAGS)
$(G_RI3DC_OBJ): $(G_RI3DC_SRC) $(G_RI3DC_H)


$(A_RI3DC): $(A_RI3DC_OBJ) $(AUX_OBJ) 
	$(LD) -o $@ $^ $(LDFLAGS)
$(A_RI3DC_OBJ): $(A_RI3DC_SRC) $(A_RI3DC_H)




$(G_XX): $(G_XX_OBJ) $(AUX_OBJ) 
	$(LD) -o $@ $^ $(LDFLAGS)
$(G_XX_OBJ): $(G_XX_SRC) $(G_XX_H)


TAGS: $(ALL_SOURCES)
	etags *.c *.h 
	find $(GMX_SOURCE_DIR) -name '*.[ch]' | xargs etags --append 


.PHONY: clean distclean install install-grid dist-grid
install:
	@echo ">> Installing all compiled programs: "
	for p in $(ALL_PROG); do \
	    if [ -e $$p ]; then  \
	       echo ">>> Check if binary \`$$p' has to be rebuilt..."; \
	       $(MAKE) $$p;      \
	       echo ">>> Installing file \`$$p' ..."; \
	       $(INSTALL) -v -m 755 $$p $(BIN_DIR); \
	    fi; \
	done;

install-grid:  $(G_RI3DC) $(A_RI3DC)
	for p in $^; do \
	    if [ -e $$p ]; then  \
	       echo ">>> Installing file \`$$p' ..."; \
	       $(INSTALL) -v -m 755 $$p $(BIN_DIR); \
	    fi; \
	done;


clean:
	-rm core *.o *.a *~ 

distclean: clean
	-rm $(ALL_PROG)

dist-grid: $(DIST_GRIDCOUNT)

$(DIST_GRIDCOUNT):
	-rm -rf grid
	mkdir grid grid/examples
	cp $(DIST_GRIDCOUNT_H) $(DIST_GRIDCOUNT_SRC)  grid
	cp examples/Makefile.grid examples/slicer.pl grid/examples
	cp README.grid LICENSE grid
	cp Makefile grid
	cd grid; patch < ../Makefile_grid.patch
	tar -zcvf $@ grid

.cvsignore: Makefile
	echo $(ALL_PROG) > $@
	echo "*.log *~ core *.a TAGS tmp" >> $@
	echo "*.xtc *.trr *.tpr *.edr *.ndx" >> $@
	echo "*Test*" >> $@
