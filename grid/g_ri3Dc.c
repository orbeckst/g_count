/*
 * $Id: g_ri3Dc.c,v 1.16 2004/03/17 22:04:31 oliver Exp $
 * $Log: g_ri3Dc.c,v $
 * Revision 1.16  2004/03/17 22:04:31  oliver
 * - Try to use liggmx xvgr routines,
 *   but had to keep xmgropen() which was moved to xf.c (together with all the
 *   other xfarbe functions)
 * - compiles with Gromacs 3.2.1
 * - updated the grid only packaging and checked that it compiles
 *
 * Revision 1.15  2003/10/10 15:19:01  oliver
 * reordered header line
 *
 * Revision 1.14  2003/10/09 09:17:11  oliver
 * reordered status line to be more useful
 *
 * Revision 1.13  2002/08/20 23:09:46  oliver
 * no major things, cleanup & checkin before I go home
 *
 * Revision 1.12  2002/08/20 01:13:50  oliver
 * NEW:
 * * a_ri3Dc can resize (well, shrink) the grid (only z1, z2, R)
 *   - needed setup_tgrid() and setup_corners() from g_ri3Dc for this
 *     dirty hack, so they were moved into grid3D.c
 * FIXED:
 *   calculation of the radius from the grid corners was wrong because I
 *   used abs() where I should have used fabs(). Didn't have any impact
 *   so far but for the readjustment of the grid it matters!
 *
 * Revision 1.11  2002/08/19 17:20:47  oliver
 * IMPORTANT:
 *   write a true density to the 3D grid file (in nm^-3)
 * small things:
 * - a few variables (sumP, dV...) are now double because for long
 *   trajectories the loss in accuracy is substantial
 * - declared command line parameters as static
 * - moved autoset_cavity() from g_ri3Dc into count.c (.h) so that I can
 *   use it for g_count as well
 * - cleaned up argument list
 *
 * Revision 1.10  2002/08/17 17:40:22  oliver
 * * normalisation of the grid is now done in g_ri3Dc for once and all
 *   (see comments). Also write average number of molecules into the header
 * * packaged autosetting the t_cavity structure from the box into
 *   function autoset_cavity() -- should be moved to one of the othe
 *   files (count.c ?) so that I can use it for g_count as well
 *   Returns TRUE if any parameter was autoset
 *
 * Revision 1.9  2002/08/17 01:02:25  oliver
 * improved setting of the grid
 * - getting rid of rounding errors
 * - define limits consistently
 *
 * Revision 1.8  2002/08/16 23:13:03  oliver
 * New:
 * * can use either time step or unitary weights for counting
 *   (select with option)
 *   uses get_timestep() from trjcat.c
 * * defaults for cpoint, radius, z1, z2 are now taken from the
 *   simulation box
 *
 * Fixed:
 * - user supplied subtitle was lost from header
 *
 * Minor mop up
 * - moved file format description to grid3D.h
 * - changed return types to bool or void
 * - removed some unused stuff (bInCircle)
 * - more debug output: file format version and  weight used
 * - reworked header text
 *
 * Revision 1.7  2002/08/15 17:05:36  oliver
 * First official release.
 * - fixed bug with length of xfarbe file header
 *
 * Revision 1.6  2002/08/15 15:27:26  oliver
 * fixed -subtitle; header is now guaranteed to be < HEADER_MAX in the grid file
 * Also renamed MAX_HEADER to HEADER_MAX
 *
 * Revision 1.5  2002/08/15 14:51:11  oliver
 * FINALLY: 3Dgrid output works. I settle for a format with interspersed tokens
 *
 * Revision 1.4  2002/08/14 15:09:59  oliver
 * splitting off functions from g_ri3Dc into grid3D.c (and .h) which are common to the grid generator (g_ri3Dc) and analysis program(s) (probably will be only one, a_ri3Dc, codenamed 'dry county')
 *
 * Revision 1.3  2002/08/14 14:29:08  oliver
 * Writing grid files as XDR files. The first time the whole thing
 * compiles AND does not crash!
 * xdr_grid.[ch] contain
 *    definition and names of token types and the t_XDRgrid structure
 *    xdr_grid():        my xdr routine for a grid (struct t_XDRgrid)
 *    xdr_grid_access(): universal xdr wrapper for all token types in the
 *                       data file (made possible by some black pointer
 * 		      magic) -- used by xdr_grid
 *    grid2_serialise():
 *    grid3_serialise(): turn a 2D or 3D array into a vector
 *
 * g_ri3Dc: grid_write() loads an t_XDRgrid structure and calls xdr_grid()
 *
 * Revision 1.2  2002/08/13 17:36:05  oliver
 * * new program: g_ri3Dc
 *   . count molecules in grid cells
 *   . effective radius pore profile
 *   . volume calculation
 *   . xy and xz density projections (xfarbe format)
 *   No 3D grid output yet (lacking a suitable format)
 * * changed DR_DEFAULT to new value
 * * smaller changes of comments
 *
 * Revision 1.1.2.3  2002/08/12 00:58:55  oliver
 * current working setup
 *
 * Revision 1.1.2.2  2002/08/11 19:39:44  oliver
 * first (kind of) working version
 * - writes output (but no 3D grid--as long as I dont know what to do
 *   with it)
 * - pore profile
 * - xy and xz projections in xfarbe format (nice!)
 *
 * Revision 1.1.2.1  2002/08/10 22:06:28  oliver
 * Start development on g_ri3Dc. Based on g_flux. Have not tried to compile yet
 *
 *
 *
 */
static char *SRCID_g_ri3Dc_c = "$Id: g_ri3Dc.c,v 1.16 2004/03/17 22:04:31 oliver Exp $";

#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
#include "rdgroup.h"
#include "mshift.h"
#include "gstat.h"
#include "names.h"
#include "utilgmx.h"
#include "count.h"
#include "xdr_grid.h"
#include "grid3D.h"

/* g_ri3Dc
   -------

   Purpose: 

   place a rectangular grid over a cylinder (that typically
   encompasses a pore) and count the occupancy of each cell over the
   length of a trajectory.

   Typically application: count water molecules, and determine the
   true solvent accessible volume of a pore.

     V = Sum_{ijk: N(i,j,k)>N_min} deltaV

   ijk       cell (0=<i<M for XX,YY,ZZ)
   N(ijk)    occupancy of cell
   N_min     OCCUPIED_MIN
   delta_V   volume of a cell = delta_x * delta_y * delta_z

     R*(z) = sqrt(V(z)/(pi*delta_z))

   R*(z)     pore profile
   V(z)      volume of a slice centered at height z

   and the total effective radius R*

     R* = sqrt(V/pi*L)
 

   Implementation:

   Corners of the grid are calculated from z1, z2, a point on the pore
   axis (cpoint) and the radius. The grid always includes the full
   circle.

   The corners are numbered from 0 to 7, starting in the lower left corner:

	       4    7  	 The edges of the grid are denoted by 
	      /|   /	 a[d] and b[d] (d=XX,YY,ZZ) and are components of
             5----6	 the corners of the box, eg
	       |  	    a[XX] = c[0][XX],  b[XX] = c[7][XX]
	       3----2	    a[YY] = c[0][YY],  b[YY] = c[7][YY]
		   /   	    a[ZZ] = c[0][ZZ],  b[ZZ] = c[7][ZZ]
       	     0----1

   Typically, the user gives a desired spatial resolution Delta[d]
   which determines the spacing, ie the number of cells mx[d] in
   each dimension d.

     mx[d] = ceil((b[d] - a[d])/Delta[d])

   
   For each molecule's center of mass x (or if counting atoms, just
   its coordinate) the position in the grid is calculated as

     ijk[d] = floor ( (x[d] - a[d])/Delta[d] )

   and if within a cell, the counter 

      occupation[ijk[d]] 

   for this cell is incremented BY THE TIMESTEP (so that we can deal
   with trajectories that have different time steps)

  grid: use XDR routines (similar to the gromacs libxdrf.c but less elaborate)
        advantages:
	- binary format (more efficient than huge ascii files)
	- could implement compression similar to xdr3dfcoord()
	- machine independent
  For the file format see grid3D.h 
*/


void init_t_result (t_result *, t_tgrid *, real **, real **, real **);
static gmx_inline double gridcount (t_tgrid *, rvec, real);
static real get_timestep(char *);

void init_t_result (t_result *res, t_tgrid *tgrid, real **profile, 
		    real **xyproj, real **xzproj) {
  res->volume = 0.0;
  res->reff   = 0.0;
  res->tgrid  = tgrid;
  res->profile= profile;
  res->xyproj = xyproj;
  res->xzproj = xzproj;
  return;
};



/* 
   check if x is in grid and add dt to the cell 

   NB: if float overflow is an issue (very long trajectories) then
   recompile with the -DDOUBLE directive to turn all reals into
   doubles. (This almost always happens for the grid sum sumP so we
   define it as double anyway and return the increment as double, too) 
*/
static gmx_inline double gridcount (t_tgrid *g, rvec x, real dt) {
  int        d;              /* loop over XX, YY, ZZ */
  real       u;              /* x-a, distance of mol from lower left of grid */
  int        ijk[3];         /* index of grid cell */

  for(d=XX; d<DIM; d++) {
    if ((u=x[d] - g->a[d]) < 0 || x[d] - g->b[d] >= 0)   return 0;
    ijk[d]=(int) floor(u/g->Delta[d]);
  }
  /* increase occupancy */
  g->grid[ijk[XX]][ijk[YY]][ijk[ZZ]] += dt;

  return (double) dt;
}

/* from trjcat.c */
static real get_timestep(char *fnm)
{
  /* read first two frames in trajectory 'fnm' to determine timestep */
  
  int        status;
  real       t0, dt;
  t_trxframe fr;
  bool ok;
  
  ok=read_first_frame(&status,fnm,&fr,TRX_NEED_X);
  if(!ok || !fr.bTime)
    fatal_error(0,"\nCouldn't read time from first frame.");
  t0=fr.time;
    
  ok=read_next_frame(status,&fr);
  if(!ok || !fr.bTime) 
    fatal_error(0,"\nCouldn't read time from second frame.");
  dt=fr.time-t0;

  close_trj(status);
  
  return dt;
}

     
int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]g_rid3Dc[TT] places a 3D grid into a simulation box and counts "
    "the number of molecules (typically water) in each cell over a trajectory. "
    "The rectangular grid is large enough to encompass a cylinder of given "
    "radius and height (specify radius R, lower and upper z, and a point on the "
    "pore axis (which is always parallel to the z-axis). If no parameters are "
    "given for the cylinder, a cylinder is fit into the simulation box. "
    "The cylinder (and thus the "
    "grid) is fixed. The spatial resolution can be given either as one "
    "number for "
    "all three dimensions (cartesian coordinates; unit is nm) or separately as "
    "three (NOTE: (1) It is recommended to use at least the same grid spacing "
    "in x and z (which is compatible with a cylindrical symmetry along z)). "
    "(2) The exact dimensions of the grid (radius etc) are re-adjusted to "
    "integral numbers of the grid spacing.)\n"
    "At the end, a density map (t[cell]/T_total) is written to the grid file."
    "[PAR]"
    "The output is the number density, averaged over the trajectory. "
    "Read this file "
    "into a_ri3Dc, the grid analysis program, and produce radial distribution ",
    "functions, density plots, pore profiles etc.\n"
    "The grid data  is written as a binary file (but xdr format, readable ",
    "on any machine with the appropriate version of a_ri3Dc). a_ri3Dc has an option "
    "to write it out as ascii text.\n" 
    "[PAR]",
    "Suggested use for water:\n",
    "create an index file for the water molecules:\n",
    "  [TT]echo -e \"keep 0\\ndel 0\\nr SOL\\nq\\n\" ",
    "| make_ndx -f in.pdb -o sol.ndx[TT]\n",
    "and run [TT]g_ri3Dc -m[TT] on it ([TT]-m[TT] is the default)."
    "[PAR]"
    "Caveats and known limitations:\n"
  };

  static char *bugs[] = {
    "The program guarantees to use the user supplied grid spacing. If the other "
    "dimensions are incommensurable with Delta they are changed to comply.",
    "If you want radial distribution functions (yes, you do!) always use "
    "Delta[XX] == Delta[YY] to keep the cylindrical symmetry", 
    "z-axis is the only allowed axis (and this will probably not "
    "change in the future)",
    "-m behaves different from the standard usage ",
    "within the g_* programs -- it figures out _for itself_ what the "
    "molecules are and does not need MOLECULE numbers but ATOM_IDs.",
    "-m is the DEFAULT behaviour. It works nicely with a SOL index file but "
    "more complicated solvents are untested.",
    "For IONS you haved to use the -nom option!",
    "The XDR file is not compressed.",
    "Even the developer mistypes the name frequently",
  };

  static bool bdtWeight  = TRUE; /* weigh counts by the timestep */ 
  static bool bMolecular = TRUE; /* default is to use molecules    */
  static t_cavity   cavity = {   /* define volume to count mols in */
    {0, 0, 1},            /* axis -- cannot be changed */
    {0, 0, 0},            /* cpoint */
    0,                    /* radius */
    0, 0,                 /* z1 < z2 */
    0                     /* volume - calculate later */
  };
  static rvec Delta = {0.02, 0.02, 0.02}; /* spatial resolution in nm */
  static char buf[HEADER_MAX];        /* additional text for graphs */
  static char *header = buf;          

  t_pargs pa[] = {
    { "-m",      FALSE, etBOOL, {&bMolecular},
      "index contains atoms, but g_ri3Dc counts the molecules"},
    { "-cpoint", FALSE, etRVEC, {&(cavity.cpoint)},
      "Point on the central symmetry axis of the grid"},
    { "-R",      FALSE, etREAL, {&(cavity.radius)},
      "Maximum radius that should be contained in the grid"},
    { "-z1",     FALSE, etREAL, {&(cavity.z1)},
      "Center grid between z1, and ..."},
    { "-z2",     FALSE, etREAL, {&(cavity.z2)},
      "z2 (these boundaries are kept fixed!)"},
    { "-delta",  FALSE, etRVEC, {&(Delta)},
      "Spatial resolution in X, Y, and Z (in nm)"},
    { "-dtweight", FALSE, etBOOL, {&bdtWeight},
      "HIDDENWeigh counts with the time step (yes) or count all equally (no)"},
    { "-subtitle", FALSE, etSTR, {&header},
      "Some text to add to the output graphs"},
  };

  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTPS, NULL, NULL, ffREAD },
    { efNDX, NULL, NULL, ffOPTRD },
    { efDAT, "-grid",    "gridxdr", ffWRITE },
  };

  FILE       *fGrid;         /* 3D grid with occupation numbers */
  t_tgrid    tgrid;          /* all information about the grid */
  t_topology *top;           /* topology                   */
  double     sumP;           /* Sum_xyz P(x,y,z) */
  double     T_tot;          /* total time */   
  real       weight;         /* each count adds one weight to P(x,y,z) */
  rvec       *x,*x_s;        /* coord, with and without pbc*/
  rvec       *xmol_cm=NULL;  /* COM coordinates for molecules */
  matrix     box;            /* box matrix (3x3)           */
  real       t,tlast,dt;     /* time, time of last frame, time step  */
  int        natoms;         /* number of atoms in system  */
  int        status;
  int        i;              /* loopcounters                 */
  char       *grpname;       /* name of the group            */
  int        gnx;            /* number of atoms in group*/
  int        gnmol;          /* number of molecules in group */
  atom_id    *molndx;        /* index of mols in atndx */
  atom_id    *index;         /* atom numbers in index file  */
  t_block    *mols;          /* all molecules in system    */
  char       s_tmp[STRLEN];  /* string for use with sprintf() */

#define NFILE asize(fnm)
#define NPA   asize(pa)

  /* no additional header by default */
  strncpy(header,"",HEADER_MAX);

  CopyRight(stderr,argv[0]);
  fprintf (stderr, "\nVersion: %s\n\n", SRCID_g_ri3Dc_c);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);

  if (bDebugMode()) {
    dfprintf ("%s -- debugging...\n\n", Program());
    dfprintf ("Source version: %s\n", SRCID_g_ri3Dc_c);
    dfprintf ("Version of the grid file format: %d\n", GRID_FF_VERSION);
  };
  dmsg("bdtWeight = %i\n",bdtWeight);


  /* open input files */
  top=read_top(ftp2fn_null(efTPS,NFILE,fnm));
  get_index(&(top->atoms),ftp2fn_null(efNDX,NFILE,fnm),
	    1,&gnx,&index,&grpname);

  /* get info about topology etc */
  mols=&(top->blocks[ebMOLS]);

  /* construct array molndx of mol indices in atndx if -m is set */
  molndx = NULL;
  gnmol  = 0;
  if (bMolecular) {
    msg ("Interpreting indexfile entries as parts of molecules and "
	 "using \ntheir center of mass.\n");
    snew (molndx, mols->nr);
    if ( (gnmol = mols_from_index (index, gnx, mols, molndx, mols->nr)) < 0) {
      fatal_error (1, "Error: could not find  molecules.\n");
    };
    msg ("%-10s%10s%10s\n", "Group", "Molecules", "Atoms");
    msg ("%-10s%10d%10d\n", grpname,  gnmol, gnx);
    msg ("%-10s%10d%10d\n", "System", mols->nr, top->atoms.nr);

    snew (xmol_cm,gnmol);
  };

  dt=get_timestep(ftp2fn(efTRX,NFILE,fnm));

  natoms = read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  snew(x_s,natoms);
  tlast = t - dt;

  /* look at cavity 
     set defaults from the simulation box size
  */
  msg("\n");
  autoset_cavity(&cavity,box,NPA,pa);

  /*  setup grid  */
  if (!setup_tgrid (&tgrid, &cavity, Delta)) 
    fatal_error (-1,"FAILED: setup_tgrid()\n");


/*
  main loop over frames
*/
  sumP = 0;       /* sum_xyz P(x,y,z) */
  T_tot = 0;
  do {
    dt = t - tlast;
    tlast = t;
    T_tot += (double) dt;
    weight = bdtWeight ? dt : 1;     /* weight for a count in P(x,y,z) */
    if (bMolecular) {
      /* prepare COM coordinates for molecules in xmol_cm */
      x2molxcm (top, molndx, gnmol, box, x, x_s, xmol_cm);
      /* loop over all molecules */
      for(i=0; i<gnmol; i++) 
	sumP += gridcount(&tgrid, xmol_cm[i], weight);
    } else {
      /* over all atoms */
      for(i=0; i<gnx; i++) 
	sumP += gridcount(&tgrid, x[index[i]], weight);
    }
  } while(read_next_x(status,&t,natoms,x,box));

  /* normalisation (here: r = (x,y,z) <-> [ijk]
     NB: Integral_V dr f(r) = Sum_r dV f(r)
     
     In general:
     h(r) = 1/T Sum_i dt_i grid[r,t] correctly weighted time average of counts
     <N>_t = Sum_r h(r)              average number of molecules in V over T
     P(r) = h(r)/<N>                 probability for finding a molecule in the 
                                     volume dV at r
     p(r) = P(r)/dV                  probability density
     n(r) = h(r)/dV                  number density
		
     here:
     h(r)  = 1/T grid[r]
     <N>_t = 1/T sumP
     n(r)  = 1/(dV*T) grid[r]
     
     Write the time averaged number density to the grid file. This is
     NOT the probability distribution but it can be easily obtained by
     computing N = Sum_r dV n(r), and p(r) = n(r)/N
  */
  { 
    int i,j,k;
    real dV = DeltaV(&tgrid);
    
    for(i=0;i<tgrid.mx[XX];i++)
      for(k=0;k<tgrid.mx[ZZ];k++)
	for(j=0;j<tgrid.mx[YY];j++)
	  tgrid.grid[i][j][k] /= (dV*T_tot);
  }
  

  strcpy(s_tmp,header);
  snprintf (header, HEADER_MAX,
	    "%s{T}%gps{grp}%s{z1}%.2f{z2}%.2f{Rmax}%.2f"
	    "{<Delta>}%.3f"
            "{<N>}%.2f{XTC}%s"
            "{Delta}(%.3f,%.3f,%.3f){sumP}%g%s",
	    s_tmp, T_tot, grpname, cavity.z1, cavity.z2, cavity.radius, 
	    (tgrid.Delta[XX]+tgrid.Delta[YY]+tgrid.Delta[ZZ])/3.,
            sumP/T_tot, ftp2fn(efTRX,NFILE,fnm), 
            tgrid.Delta[XX],tgrid.Delta[YY],tgrid.Delta[ZZ],
            sumP, (bdtWeight ? "ps" : ""));
  header[HEADER_MAX-1] = '\0';   /* better safe than sorry */
  msg("Writing 3D density grid with header\n%s\n",header);
	  
  fGrid    = ffopen (opt2fn("-grid",    NFILE, fnm), "w");
  if (!grid_write(fGrid,&tgrid,header))
    msg("Writing the 3D grid failed.\n");

  /* clean up a bit */

  sfree(tgrid.grid);

  fprintf(stderr,"\n");
  close_trj(status);
  fclose(fGrid);
  thanx(stdout);

  return 0;
}

