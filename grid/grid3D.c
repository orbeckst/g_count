/*
  $Id: grid3D.c,v 1.7 2004/03/17 23:26:29 oliver Exp $
  $Log: grid3D.c,v $
  Revision 1.7  2004/03/17 23:26:29  oliver
  - realised that I was using the TOKEN-less xdr format although I though
    that I was using the other one, ie xdr_token_grid()
  - completely switched over to TOKEN-less but left some of the TOKEN code
    in to remind myself how to do it
  - tested both format: both work but are incompatible and fail on each others
    data files (even though I removed stringent TOKEN/TOKEN-less checks that
    used to be there from development days)

  Revision 1.6  2002/08/20 01:13:50  oliver
  NEW:
  * a_ri3Dc can resize (well, shrink) the grid (only z1, z2, R)
    - needed setup_tgrid() and setup_corners() from g_ri3Dc for this
      dirty hack, so they were moved into grid3D.c
  FIXED:
    calculation of the radius from the grid corners was wrong because I
    used abs() where I should have used fabs(). Didn't have any impact
    so far but for the readjustment of the grid it matters!

  Revision 1.5  2002/08/19 17:20:47  oliver
  IMPORTANT:
    write a true density to the 3D grid file (in nm^-3)
  small things:
  - a few variables (sumP, dV...) are now double because for long
    trajectories the loss in accuracy is substantial
  - declared command line parameters as static
  - moved autoset_cavity() from g_ri3Dc into count.c (.h) so that I can
    use it for g_count as well
  - cleaned up argument list

  Revision 1.4  2002/08/17 17:37:10  oliver
  volume calculations from t_grid packaged as two functions, DeltaV()
  and gridV

  Revision 1.3  2002/08/15 23:51:39  oliver
  Changed file format: start directly with the int for VERSION in both formats so that they can be detected gracefully

  Revision 1.2  2002/08/15 14:51:11  oliver
  FINALLY: 3Dgrid output works. I settle for a format with interspersed tokens

  Revision 1.1  2002/08/14 15:09:59  oliver
  splitting off functions from g_ri3Dc into grid3D.c (and .h) which are common to the grid generator (g_ri3Dc) and analysis program(s) (probably will be only one, a_ri3Dc, codenamed 'dry county')


  collect stuff common to g_ri3Dc and a_ri3Dc

  see comments about the file format in grid3D.h
*/

#include "grid3D.h"

/*
  write the 3D grid to the output file with XDR routines

  load the t_XDRgrid structure; one small problem is that I want the
  data already as one long vector. The grid3_alloc() puts it already
  in one chunk but I want to keep it general so we need some
  grid_serialiser
*/
bool grid_write (FILE *fp,t_tgrid *tg,char *header) {
  XDR         xdrs;
  t_XDRgrid   xg;
  
  xg.version = GRID_FF_VERSION;
  xg.header  = header;
  xg.type    = egtyREGULAR;
  xg.dim     = 3;
  xg.size    = tg->mx;
  xg.delta   = tg->Delta;
  xg.origin  = tg->a;
  
  if(! (xg.grid=grid3_serialise(tg->grid,tg->mx))) 
    fatal_error(0,"grid_write(): Cannot allocate memory for grid "
		"serialisation.\n");

  xdrstdio_create (&xdrs, fp, XDR_ENCODE);
  if (! XDR_GRID(&xdrs,&xg)) {
    msg("grid_write(): XDR write failed.\n");
    return FALSE;
  };
  xdr_destroy(&xdrs);
  sfree(xg.grid);
  return TRUE;
}


/* load the t_tgrid structure from the XDR stream; 
   XDR_GRID() allocates memory for the grid on reading

   NB: We do not check during reading that the dimension of the grid
   is really 3 (and thats the only number that fits into a t_tgrid
   structure), hence this routine is definitely not general ... in
   principle, xdr can auto-malloc and determine the size on the fly
   (using the xdr_array) but this requires a rewrite for t_tgrid, too,
   so that it uses pointers instead of declared array[3] for delta,
   origin, and max */
bool grid_read (FILE *fp,t_tgrid *tg,char *header) {
  XDR         xdrs;
  t_XDRgrid   xg;
  int i;
  
  xg.header  = header;
  xg.size    = tg->mx;
  xg.delta   = tg->Delta;
  xg.origin  = tg->a;
  xg.grid    = NULL;        /* force XDR_GRID() to malloc */

  xdrstdio_create (&xdrs, fp, XDR_DECODE);
  if (! XDR_GRID(&xdrs,&xg)) {  
    msg("grid_read(): XDR read failed.\n");
    return FALSE;
  }; 
  if (abs(xg.version) > abs(GRID_FF_VERSION)) 
    msg("WARNING: grid was WRITTEN with version %d, but is READ "
	"with version %d\n",
	xg.version, GRID_FF_VERSION);
  if (xg.dim != 3)                       /* see CAVEAT above */
    fatal_error(0,"grid_read(): The data file appears to have dim=%d, but this \n"
	          "             routine can only read 3D grids.\n",xg.dim); 

  dmsg("Version (file) %d\n",xg.version);
  dmsg("GridType   %s\n",GRIDTYPE(xg.type));
  dmsg("Dimension  %d\n",xg.dim);
  dmsg("Size       ");
  for(i=0;i<xg.dim;i++) dmsg("%5d ", xg.size[i]);
  dmsg("\n");
  dmsg("Delta      ");
  for(i=0;i<xg.dim;i++) dmsg("%5.3f ", xg.delta[i]);
  dmsg("\n");
  dmsg("Origin     ");
  for(i=0;i<xg.dim;i++) dmsg("%5.3f ", xg.origin[i]);
  dmsg("\n");
  
  if(! (tg->grid=grid3_unserialise(xg.grid,xg.size))) 
    fatal_error(0,"grid_write(): Cannot allocate memory for grid "
		"un-serialisation.\n");

  xdr_destroy(&xdrs);
  return TRUE;
}

gmx_inline double DeltaV (t_tgrid *tg) {
  return (double)tg->Delta[XX] * (double)tg->Delta[YY] * (double)tg->Delta[ZZ];
}

gmx_inline double DeltaA (t_tgrid *tg) {
  return (double)tg->Delta[XX] * (double)tg->Delta[YY]; 
}

gmx_inline double gridV (t_tgrid *tg) {
  return (double)tg->mx[XX] * (double)tg->mx[YY] * (double)tg->mx[ZZ] * DeltaV(tg);
}

gmx_inline double gridA (t_tgrid *tg) {
  return (double)tg->mx[XX] * (double)tg->mx[YY] * DeltaA(tg);
}


bool setup_tgrid (t_tgrid *g, t_cavity *cylinder, rvec Delta) {
  rvec c[8];     /* corners of the grid */
  real du;       /* difference between original a,b and after grid setup */
  int d;         /* coordinate dimensions XX, YY, ZZ */
  int gsize;     /* total mem allocated for grid */

  setup_corners (c, cylinder);
  copy_rvec (c[0], g->a);
  copy_rvec (c[7], g->b);

  copy_rvec (Delta, g->Delta);

  for (d=XX; d<DIM; d++) {
    /* floor() so that we dont get a half empty fringe cell. float
       operations can give a bit funny results here so we see to it
       that representation errors dont screw the
       setup. 100*GMX_REAL_EPS is more empirical than well founded..
    */
    g->mx[d] = (int) floor ((g->b[d] - g->a[d])/Delta[d] + 100*GMX_REAL_EPS);
    /* make corners consistent with grid */
    du = 0.5*( (g->b[d] - g->a[d]) - g->mx[d]*Delta[d]);
    du = ((real)fabs(du) < 100*GMX_REAL_EPS ? 0 : du);
    g->a[d] += du;
    g->b[d] -= du;
  };

  /* readjust cylinder -- we DO change the input parameters if the
     grid spacing is incommensurable with them */
  cylinder->z1 = g->a[ZZ];
  cylinder->z2 = g->b[ZZ];
  cylinder->radius = 0.5 * MIN( fabs(g->b[XX] - g->a[XX]),
				fabs(g->b[YY] - g->a[YY]));
  msg("After grid placement: new values for\n"
       "  radius = %.3f nm\n  z1 = %.3f nm\n  z2 = %.3f nm\n",
       cylinder->radius,cylinder->z1,cylinder->z2);


  gsize = g->mx[XX] * g->mx[YY] * g->mx[ZZ];
  dmsg ("setup_tgrid(): boundaries of the grid: a[3] = (%.3f, %.3f, %.3f)\n"
	"                                       b[3] = (%.3f, %.3f, %.3f)\n",
	g->a[XX], g->a[YY], g->a[ZZ], 
	g->b[XX], g->b[YY], g->b[ZZ]);
  dmsg ("setup_tgrid(): Allocating %d x %d x %d = %d bytes = %.1f MB \n"
	"               for the grid at resolution Delta = (%g, %g, %g) nm.\n",
	g->mx[XX], g->mx[YY], g->mx[ZZ], 
	sizeof(real)*gsize, sizeof(real)*gsize/(1024.0*1024.0), 
	Delta[XX], Delta[YY], Delta[ZZ]);

  g->grid = grid3_alloc(g->mx[XX],g->mx[YY],g->mx[ZZ]);
  return (g->grid != NULL);
};


void setup_corners (rvec c[], t_cavity *cyl) {
  /* we know that there are 8 elements in c[] */
  rvec rx, ry;    /* steps along the grid boundary */

  new_rvec (rx, 2*cyl->radius, 0, 0);
  new_rvec (ry, 0, 2*cyl->radius, 0);

  /* walk along the edges 0 -> 7 
     c[n-1] + r = c[n]
     Actually, we only need c[0] and c[7] for now, but it does not really
     cost anything...
  */
  new_rvec (c[0], cyl->cpoint[XX] - cyl->radius, 
	          cyl->cpoint[YY] - cyl->radius, 
	          cyl->z1);
  rvec_add (c[0], rx, c[1]);
  rvec_add (c[1], ry, c[2]);
  rvec_sub (c[2], rx, c[3]);
  new_rvec (c[4], cyl->cpoint[XX] - cyl->radius, 
	          cyl->cpoint[YY] + cyl->radius, 
	          cyl->z2);
  rvec_sub (c[4], ry, c[5]);
  rvec_add (c[5], rx, c[6]);
  rvec_add (c[6], ry, c[7]);
  return;
};

