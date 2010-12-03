# GNU Makefile to compile the g_count and friends
#
#   Copyright (C) 2003-2010 Oliver Beckstein <orbeckst@gmailcom>
#   This program is made available under the terms of the GNU Public License. 
#   See the file COPYING or http://www.gnu.org/copyleft/gpl.html
#
# Edit BIN_DIR and the GMX_* variables to reflect the locations in
# your setup. Run 'make help' for, um, help.

# Makefile switches:
# Switch on on commandline by setting them, eg 
#        make DEBUG=1 WITH_MPI=1 g_count 
#
#  DEBUG        compile with debugging, -g and no optimisations

#====================================================================
# 
# THE FOLLOWING VARIABLES MUST BE CHECKED/SET BY THE USER
# -------------------------------------------------------
#
# See the file INSTALL for instructions.
#
# Include directory for Gromacs header files:
GMX_INCLUDE_DIR = /usr/local/gromacs/include/gromacs

# Set the directories where Gromacs libraries are to be found:
GMX_LIB_DIR     = /usr/local/gromacs/lib

# Install binaries into:
BIN_DIR = /usr/local/gromacs/bin
#
#
#====================================================================

# this is only necessary for the creation of etags
# (for compilation it is not important)
GMX_SOURCE_DIR  := 


CPPFLAGS += -I$(GMX_INCLUDE_DIR)

ifdef DEBUG
CFLAGS   +=  -DDEBUG -g -Wall  -Wno-unused 
else
CFLAGS   +=  -g -O2 -fomit-frame-pointer -finline-functions \
             -funroll-loops -Wall -Wno-unused 
endif


# dynamic linking of libgm, require LD_LIBRARY_PATH
#LD_LIBRARY_PATH += /usr/src/gm-1.6.3_Linux/binary/lib/
# for libgm.so on synapse

LDFLAGS  +=  -lm -L$(GMX_LIB_DIR) -lmd -lgmx

ifeq ($(CC),icc)
   # '-lpthread ' may be required for icc
   LDFLAGS +=  -lpthread
endif

# hmph... seems to be required on Mac OS X, otherwise I get
# "ld: could not find entry point "start" (perhaps missing crt1.o) for inferred architecture i386"
LD := $(CC)


INSTALL := install

AUX_NAMES :=  utilgmx xf count
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


# template: g_xx
G_XX     := g_xx
G_XX_SRC := g_xx.c 
G_XX_H   := $(AUX_H)
G_XX_OBJ := g_xx.o 

ALL_PROG := $(G_COUNT) $(G_FLUX) $(G_ZCOORD)

# removed:                
ALL_SOURCES := $(G_COUNT_SRC) $(G_COUNT_H) \
               $(G_FLUX_SRC) $(G_FLUX_H) \
               $(G_ZCOORD_SRC) $(G_ZCOORD_H) \
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
\n   clean        clean object files etc\
\n   distclean    remove every generated file\
\n\
\nMakefile switches:\
\nSwitch on on commandline by setting them, eg \`make DEBUG=1 g_count' \
\n  DEBUG        compile with debugging, -g and no optimisations\
\n\
\nImportant Makefile parameters:\
\nGMX_LIB_DIR...............$(GMX_LIB_DIR)\
\nGMX_INCLUDE_DIR...........$(GMX_INCLUDE_DIR)\
\n\
\nBIN_DIR...................$(BIN_DIR)\
\n\
\nCC........................$(CC)\
\nLD........................$(LD)\
\nCPPFLAGS..................$(CPPFLAGS)\
\nCFLAGS....................$(CFLAGS)\
\nLDFLAGS...................$(LDFLAGS)
endef
# 'emacs font-lock

define include_check
if test -e $(GMX_INCLUDE_DIR)/xtcio.h; then \
    echo "OK   Gromacs include directory found: $(GMX_INCLUDE_DIR)"; \
else \
    echo "BAD  Gromacs include directory missing. Set GMX_INCLUDE_DIR in Makefile!"; \
fi
endef

_LIBS = $(wildcard $(GMX_LIB_DIR)/libgmx.*)
define lib_check
libs=($(_LIBS)); \
if test "$${#libs[@]}" -gt 0; then \
    echo "OK   Gromacs lib directory found: $(GMX_LIB_DIR)"; \
else \
    echo "BAD  Gromacs lib directory missing. Set GMX_LIB_DIR in Makefile!"; \
fi
endef

.PHONY: all check help
all:	$(ALL_PROG)

help:
	@echo -e "$(usage)"

check:
	@echo "============================================================"
	@echo "Checking if GMX_INCLUDE_DIR and GMX_LIB_DIR are set properly"
	@echo "============================================================"
	@$(call include_check)
	@$(call lib_check)
	@echo "============================================================"
	@echo "If you got any 'BAD' entries please read INSTALL."


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



# $(G_XX): $(G_XX_OBJ) $(AUX_OBJ) 
# 	$(LD) -o $@ $^ $(LDFLAGS)
# $(G_XX_OBJ): $(G_XX_SRC) $(G_XX_H)


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


clean:
	-rm core *.o *.a *~ 

distclean: clean
	-rm $(ALL_PROG) TAGS

.cvsignore: Makefile
	echo $(ALL_PROG) > $@
	echo "*.log *~ core *.a TAGS tmp" >> $@
	echo "*.xtc *.trr *.tpr *.edr *.ndx" >> $@
	echo "*Test*" >> $@

# new numbering scheme: major is Gromacs compatibility
#                       minor are my releases
NAME  := g_count
MAJOR := gmx4.0
MINOR := 6

TAR_NAME := $(NAME)-$(MAJOR)-$(MINOR).tar.bz2
TAR_DIR  := $(NAME)-$(MAJOR)-$(MINOR)

.phony: dist rsync

dist: $(TAR_NAME)
$(TAR_NAME): $(ALL_SOURCES) Makefile README INSTALL LICENCE ChangeLog examples
	rm -rf $(TAR_DIR)
	mkdir $(TAR_DIR)
	cp -r $^ $(TAR_DIR)
	tar --exclude=*~ -jcvf $@ $(TAR_DIR) && rm -rf $(TAR_DIR)

UPLOAD_URI := clathrin:/sansom/public_html/html/sbcb/oliver/download/Gromacs
rsync: $(TAR_NAME)
	rsync -avP $< $(UPLOAD_URI) 
