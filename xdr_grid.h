/*
  $Id: xdr_grid.h,v 1.4 2004/03/17 23:26:29 oliver Exp $
  $Log: xdr_grid.h,v $
  Revision 1.4  2004/03/17 23:26:29  oliver
  - realised that I was using the TOKEN-less xdr format although I though
    that I was using the other one, ie xdr_token_grid()
  - completely switched over to TOKEN-less but left some of the TOKEN code
    in to remind myself how to do it
  - tested both format: both work but are incompatible and fail on each others
    data files (even though I removed stringent TOKEN/TOKEN-less checks that
    used to be there from development days)

  Revision 1.3  2002/08/15 23:51:39  oliver
  Changed file format: start directly with the int for VERSION in both formats so that they can be detected gracefully

  Revision 1.2  2002/08/15 14:51:11  oliver
  FINALLY: 3Dgrid output works. I settle for a format with interspersed tokens

  Revision 1.1  2002/08/14 14:29:08  oliver
  Writing grid files as XDR files. The first time the whole thing
  compiles AND does not crash!
  xdr_grid.[ch] contain
     definition and names of token types and the t_XDRgrid structure
     xdr_grid():        my xdr routine for a grid (struct t_XDRgrid)
     xdr_token_grid():  same purpose as xdr_grid but intersperses data 
                        with tokens to aid in debugging/error detection  
     xdr_grid_access(): universal xdr wrapper for all token types in the
                        data file (made possible by some black pointer
  		        magic) -- used by xdr_grid

  If the preprocessor directive GRID_NO_TOKENS is defined, the
  xdr_grid routine is used for reading/writing (this uses the
  XDR_GRID() macro). The version numbers of the token-less format are
  the negative of the base version number.
			
     grid2_serialise():
     grid3_serialise(): turn a 2D or 3D array into a vector

  g_ri3Dc: grid_write() loads an t_XDRgrid structure and calls xdr_grid()

*/

#ifndef _xdr_grid_h
#define _xdr_grid_h

static char *SRCID_xdr_grid_h = "$Id: xdr_grid.h,v 1.4 2004/03/17 23:26:29 oliver Exp $";

#include <stdio.h>
#include <rpc/xdr.h>
#include "fatal.h"
#include "names.h"
#include "smalloc.h"
#include "typedefs.h"
#include "utilgmx.h"

/* Version number of the file format is an integer, which is incremented 
   when changes are made. If the change breaks something then one has to add 
   specific tests on the format in version_check() and return INCOMPATIBLE 
*/
#define GRID_FF_VERSION 1      /* grid xdr file format version */
#define HEADER_MAX 256         /* descriptive header string */
#define NGRID_MAX  512*512*512 /* max size for a grid (so that I can use
				  xdr auto-allocate in xdr_array()*/

/* these must correspond to eGridToken_names[] and eGridType_names[] */
enum eGridToken {
  egtokVERSION,egtokHEADER,egtokGRIDTYPE,egtokDIMENSION,egtokSIZE,egtokDELTA,
  egtokORIGIN,egtokGRID,egtokNR
};
enum eGridType {egtyREGULAR,egtyOTHER,egtyNR};
enum version_ok {INCOMPATIBLE,OLDER,OK,verokNR};

extern char *eGridToken_names[egtokNR+1];
extern char *eGridType_names[egtyNR+1];
extern char *version_ok_names[verokNR+1];

typedef struct {
  int   version;
  char  *header;
  enum eGridType type;
  int   dim;
  int   *size;
  real  *delta;
  real  *origin;
  real  *grid;     /* serialised grid (one vector) because I do not
                      know how to put there a type which handles 2D
                      (**real, 3D ***real etc) */
} t_XDRgrid;

extern enum version_ok version_check(int);
extern bool xdr_grid (XDR *, t_XDRgrid *);
extern bool xdr_grid_access (XDR *,enum eGridToken,void *,int);
extern real *grid3_serialise(real ***,int *);
extern real *grid2_serialise(real  **,int *);
extern real ***grid3_unserialise(real *,int *);
extern real  **grid2_unserialise(real *,int *);

#define GRIDTOKEN(e)   ENUM_NAME((e),egtokNR,eGridToken_names)
#define GRIDTYPE(e)    ENUM_NAME((e),egtyNR,eGridType_names)
#define VERSION_OK(e)  ENUM_NAME((e),verokNR,version_ok_names)


#endif /* _xdr_grid_h */
