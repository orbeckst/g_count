## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.
#
# $Id: Makefile.am,v 1.1.1.1 2001/11/12 20:49:01 oliver Exp $
# $Log: Makefile.am,v $
# Revision 1.1.1.1  2001/11/12 20:49:01  oliver
# Transferring g_count into the Biop CVS repository -- this checkin only
# creates the directory structure and adds README and stuff to
# incorporate this into a Gromacs source tree
#
#
# written for compiling Gromacs v3 programmes
# Based on gromacs-3.0.3/src/tools/Makefile
#


AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"

LDADD = ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la
LSQ   = $(top_srcdir)/src/tools/lsq.c

EXTRA_DIST              = README  

# These programs are not compiled or installed by default - you will have to 
# issue "make <program>" and copy the binary to the correct location yourself! 
# Add new entries in Makefile.am!

bin_PROGRAMS         =   g_count g_zcoord 
EXTRA_PROGRAMS       =   g_diff g_msd2

# Note: you don't have to list sources for "prog" if it is the single file prog.c

g_count_SOURCES 	= g_count.c xmgr.c utilgmx.c
g_zcoord_SOURCES        = g_zcoord.c xmgr.c
g_msd2_SOURCES          = g_msd2.c  xmgr.c $(LSQ)


CLEANFILES   = 	*~ \\\#*










