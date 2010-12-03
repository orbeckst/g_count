/*
  $Id: xdr_grid.c,v 1.4 2004/03/17 23:26:29 oliver Exp $
  $Log: xdr_grid.c,v $
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
     xdr_grid_access(): universal xdr wrapper for all token types in the
                        data file (made possible by some black pointer
  		      magic) -- used by xdr_grid
     grid2_serialise():
     grid3_serialise(): turn a 2D or 3D array into a vector

  g_ri3Dc: grid_write() loads an t_XDRgrid structure and calls xdr_grid()

*/

#include "xdr_grid.h"


char *eGridToken_names[egtokNR+1] = 
{"VERSION","HEADER","GRIDTYPE","DIMENSION","SIZE","DELTA","ORIGIN","GRID"};
char *eGridType_names[egtyNR+1] = {"regular","other"};
char *version_ok_names[verokNR+1] = {"INCOMPATIBLE","OLDER","OK"};

#ifdef DOUBLE
#define xdr_real xdr_double
#else
#define xdr_real xdr_float 
#endif


enum version_ok version_check (int version) {
  if (abs(version) < abs(GRID_FF_VERSION))  return OLDER;
  return OK;
}
    
  
/* 
   all xdr_* return 1 on success, 0 on error
*/



/*
  Simple grid format without interspersed tokens. This is the default
  format (as in: you need to edit XDR_GRID() in grid3D.h in order to
  use the token-interspersed format). There does not seem to be a big
  advantage in using the token format, although in principle you have
  some error checking capability to detect file corruption.
*/

bool xdr_grid (XDR *xdrs, t_XDRgrid *xg) {
  int ngrid;
  
  ngrid = xg->size[0] * xg->size[1] * xg->size[2];
  return (xdr_int    (xdrs, &xg->version)           &&
	  (version_check(xg->version) > INCOMPATIBLE) &&
	  xdr_string (xdrs, &xg->header, HEADER_MAX) &&
	  xdr_enum   (xdrs,(enum_t *) &xg->type)    &&
	  xdr_int    (xdrs, &xg->dim)               &&
	  xdr_vector (xdrs, (char *) xg->size,   
		      xg->dim, sizeof(int), (xdrproc_t) xdr_int)   &&
	  xdr_vector (xdrs, (char *) xg->delta,  xg->dim, 
		      sizeof(real), (xdrproc_t)xdr_real)   &&
	  xdr_vector (xdrs, (char *) xg->origin, xg->dim, 
		      sizeof(real), (xdrproc_t)xdr_real)   &&
	  xdr_array  (xdrs, (char **) &xg->grid, &ngrid, NGRID_MAX,
		      sizeof(real), (xdrproc_t)xdr_real) );
} 

/* 
   File format with tokens automagically (ie using xdr_grid_access)
   interspersed.  (This is incompatible with the token-less format).
   This is more fancy and allows for additional error checks. It can
   be used by defining 
       #define XDR_GRID(xdr,xg) xdr_token_grid((xdr),(xg)) 
   in grid3D.h. However, during development I accidentally stuck with
   the simple format, which seems to work well enough for a year
   now. Howver, with tokens works, too. I leave in these two routines 
   and all the TOKEN clutter mainly to remind myself how to do it....
 */
bool xdr_token_grid (XDR *xdrs, t_XDRgrid *xg) {
  return (xdr_grid_access(xdrs,egtokVERSION,  &xg->version,    0)    &&
	  (version_check(xg->version) > INCOMPATIBLE)                &&
	  xdr_grid_access(xdrs,egtokHEADER,   &xg->header,     0)    &&
	  xdr_grid_access(xdrs,egtokGRIDTYPE, &xg->type,       0)    &&
	  xdr_grid_access(xdrs,egtokDIMENSION,&xg->dim,        0)    &&
	  xdr_grid_access(xdrs,egtokSIZE,     &xg->size,  xg->dim)   &&
	  xdr_grid_access(xdrs,egtokDELTA,    &xg->delta, xg->dim)   &&
	  xdr_grid_access(xdrs,egtokORIGIN,   &xg->origin,xg->dim)   &&
	  xdr_grid_access(xdrs,egtokGRID,     &xg->grid, 
			  xg->size[0] * xg->size[1] * xg->size[2]) );
}


bool xdr_grid_access (XDR *xdrs,const enum eGridToken current,void *data,int n) {
  /* for some tokens, data points to an array and we need to know the
     size n of it. If *data == NULL memory is allocated here.

     nb:  this routine does READ and WRITE
     nb2: thank gromacs that I have libxrdf.c to see how it is done...
     nb3: stupid xdr routines: CAST EVERYTHING to what they want (see man xdr)
     and http://docs.sun.com/?q=xdr_enum&p=/doc/802-5885/6i9k4u0bp&a=view
  */
  char *sp;   /* needed to write constants thru xdr */
  int  *ip;
  real **rp;
  enum eGridType *enum_tmp=NULL;
  enum_t token;
  
  token = (enum_t)current;
  
  if (current != egtokVERSION) {
    /* file format does not have a VERSION token rather the version is the 
       first entry 
    */
    if (!xdr_enum(xdrs,&token)) return FALSE;
    if (current != (enum eGridToken)token) {
      msg("xdr_grid_access(): Read error: Token expected '%s' != token "
	  "read '%s'\n", GRIDTOKEN(current),GRIDTOKEN((enum eGridToken)token));
      return FALSE;
    };
  };

  switch ((enum eGridToken)token) {
  case egtokVERSION:
  case egtokDIMENSION:
    ip = (int *)data;
    return  xdr_int(xdrs,ip);
  case egtokHEADER:
    sp = (char *)data;
    return xdr_string(xdrs,&sp,HEADER_MAX);
  case egtokGRIDTYPE:
    enum_tmp = (enum eGridType *)data;
    return xdr_enum(xdrs,(enum_t *)enum_tmp);
  case egtokSIZE:
    return xdr_vector(xdrs, (char *)(*(int **)data), (unsigned int)n,
		      sizeof(int), (xdrproc_t)xdr_int);
  case egtokDELTA:
  case egtokORIGIN:
  case egtokGRID:
    rp = (real **)data;
    if (!*rp && xdrs->x_op == XDR_DECODE) {
      snew(*rp,n);
      if (!*rp) return FALSE;
    } 
    else if (xdrs->x_op == XDR_FREE)
      sfree(*rp);

    return xdr_vector(xdrs,(char *)(*rp), (unsigned int)n,
		      sizeof(real), (xdrproc_t)xdr_real);
  default:
    fatal_error(0,"xdr_token_write(): unknown token type %d.\n", token);
    break;
  }
  return FALSE;
};


/* 
   serialise a 3D array of given sizes mx[] 
*/
real *grid3_serialise(real ***g,int *mx) {
  real *v,*tmp;        /* serialised grid */
  int i,j,k;           /* indices of the 3D grid */ 
  
  snew(v,mx[0]*mx[1]*mx[2]);
  if (!v) return NULL;

  tmp = v;
  for(k=0;k<mx[2];k++)
    for(j=0;j<mx[1];j++)
      for(i=0;i<mx[0];i++)
	*(tmp++) = g[i][j][k];
  return v;
}

/* 
   serialise a 2D array of given sizes mx[] 
*/
real *grid2_serialise(real **g,int *mx) {
  real *v,*tmp;        /* serialised grid */
  int i,j;             /* indices of the 2D grid */ 
  
  snew(v,mx[0]*mx[1]);
  if (!v) return NULL;

  tmp = v;
  for(j=0;j<mx[1];j++)
    for(i=0;i<mx[0];i++)
      *(tmp++) = g[i][j];
  return v;
}


real ***grid3_unserialise(real *v,int *mx) {
  real ***g;           /* 3D grid */
  int i,j,k;           /* indices of the 3D grid */ 
  
  g=grid3_alloc(mx[0],mx[1],mx[2]);
  if (!g) return NULL;

  for(k=0;k<mx[2];k++)
    for(j=0;j<mx[1];j++)
      for(i=0;i<mx[0];i++)
	g[i][j][k] = *(v++);
  return g;
}


real **grid2_unserialise(real *v,int *mx) {
  real **g;            /* 2D grid */
  int i,j;             /* indices of the 3D grid */ 
  
  g=grid2_alloc(mx[0],mx[1]);
  if (!g) return NULL;

  for(j=0;j<mx[1];j++)
    for(i=0;i<mx[0];i++)
      g[i][j] = *(v++);
  return g;
}
