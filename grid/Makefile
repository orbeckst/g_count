# $Id: Makefile,v 1.29 2004/03/17 22:04:31 oliver Exp $
#
#
# Makefile to compile the grid counter
# run 'make help' for, um, help

# Makefile switches:
# Switch on on commandline by setting them, eg 
#        make DEBUG=1  g_ri3Dc 
#
#  DEBUG        compile with debugging, -g and no optimisations


GMX_TOP_DIR     := /usr/local/gromacs

# EXEC depends on your machine/OS
GMX_EXEC_PREFIX := $(GMX_TOP_DIR)/i686-pc-linux-gnu#
GMX_LIB_DIR     := $(GMX_EXEC_PREFIX)/lib#
GMX_INCLUDE_DIR := $(GMX_TOP_DIR)/include#

# this is only necessary for the creation of etags
# (for compilation it is not important)
GMX_SOURCE_DIR  := /usr/local/src/gmx

# install binaries to
#BIN_DIR := $(HOME)/bin
BIN_DIR := $(GMX_EXEC_PREFIX)/bin

CPPFLAGS += -I$(GMX_INCLUDE_DIR)

ifdef DEBUG
CFLAGS   +=  -DDEBUG -g -Wall  -Wno-unused 
else
CFLAGS   +=  -g -O2 -fomit-frame-pointer -finline-functions \
             -funroll-loops -Wall -Wno-unused 
endif

LDFLAGS  +=  -lm -L$(GMX_LIB_DIR) -lmd -lgmx 

CC      := gcc
LD      := $(CC)
INSTALL := install

AUX_NAMES :=  utilgmx xf count plt
AUX_SRC := $(addsuffix .c, $(AUX_NAMES))
AUX_H   := $(addsuffix .h, $(AUX_NAMES))
AUX_OBJ := $(addsuffix .o, $(AUX_NAMES))


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

ALL_PROG := $(G_RI3DC) $(A_RI3DC)

define usage
\nIn order to compile programs edit the Makefile for paths to the Gromacs\
\nlibraries. Then do \
\n   make clean; make PROGRAM\
\nwhere PROGRAM can be one of \`$(ALL_PROG)'.\
\nPerhaps you have to edit variables at the top of the Makefile \
\nto make it work. Programs are statically linked so you can take the \
\nbinaries wherever you like (hopefully...).\
\n\
\nIn order to install all programs, type\
\n   make BIN_DIR=<your_target_dir> install\
\nor set BIN_DIR at the top of the Makefile.\
\n\
\nInstall targets: \
\n   all          compile \`$(ALL_PROG)' (default)\
\n   install      install all compiled programs in BIN_DIR\
\n   clean        clean object files etc\
\n   distclean    remove every generated file\
\n\
\nMakefile switches:\
\nSwitch on on commandline by setting them, eg \`make WITH_MPI=1 g_count' \
\n  DEBUG        compile with debugging, -g and no optimisations
endef
# 'emacs font-lock

.PHONY: all help 
all:	$(ALL_PROG)

help:
	@echo -e "$(usage)"

$(AUX_OBJ): $(AUX_SRC) $(AUX_H)


$(G_RI3DC): $(G_RI3DC_OBJ) $(AUX_OBJ) 
	$(LD) -o $@ $^ $(LDFLAGS)
$(G_RI3DC_OBJ): $(G_RI3DC_SRC) $(G_RI3DC_H)


$(A_RI3DC): $(A_RI3DC_OBJ) $(AUX_OBJ) 
	$(LD) -o $@ $^ $(LDFLAGS)
$(A_RI3DC_OBJ): $(A_RI3DC_SRC) $(A_RI3DC_H)


TAGS: $(ALL_SOURCES)
	etags *.c *.h 
	find $(GMX_SOURCE_DIR) -name '*.[ch]' | xargs etags --append 


.PHONY: clean distclean install 

install:  $(G_RI3DC) $(A_RI3DC)
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

.cvsignore: Makefile
	echo $(ALL_PROG) > $@
	echo "*.log *~ core *.a TAGS tmp" >> $@
	echo "*.xtc *.trr *.tpr *.edr *.ndx" >> $@
	echo "*Test*" >> $@
