/*
 * $Id: g_flux.c,v 1.6 2004/03/17 22:04:31 oliver Exp $
 */
static char *SRCID_g_flux_c = "$Id: g_flux.c,v 1.6 2004/03/17 22:04:31 oliver Exp $";

#include <math.h>
#include <string.h>
#include <stdarg.h>
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
#include "xvgr.h"
#include "gstat.h"
#include "names.h"
#include "xf.h"
#include "utilgmx.h"
#include "count.h"

/* algorithm for detecting pore crossings ('pX'): "Check on exit"
   For each molecule
   (1) * detect boundary crossing across z=z1 or z2 ('bX')

         min(x[t-1], x[t]) < z < max(x[t-1], x[t]) ? 1 : 0
	 
       * type of the crossing (enter or exit)

	 z1 < x[t] < z2:
                        out -> in  : (O,I)   (molecule ENTERS the pore)
			divergence = 1
         x[t] < z1  or  x[t] > z2: 
                        in  -> out : (I,O)   (molecule EXITS  the pore)
			divergence = -1
       * compute direction u relative to the pore axis:
                                               / +1  up
         u = sgn( Dot(x[t]-x[t-1], axis) ) = {
                                               \ -1  down
   if bX occurred AND divergence == +1 
      store (t, u) 

   (IMPORTANT: only trajectories O->I->O are counted as complete
   transitions through the pore; I->O->I are periodic boundary trajectories.
   Thus, we only need to remember ENTER (influx) states of water molecules)

   (2) detect pX
   if bX occurred AND divergence == -1 
      compare u[t] with the stored u[t-Dt]
      If the direction of the flow is the same (u[t] == u[t-Dt])
         increment counters

      erase (t-DT,u)

  
   Hacks:
   * u==0 used as a signal that (t,u) is not initialised
       
*/

/* residency time stuff -- quick hack */
#define RESIDENCY_BINS 1000    /* max number of bins */
#define RESIDENCY_TMAX 2000.0  /* max time 1000 ps */
                               /* ==> bin width = 2 ps */

#define DR_DEFAULT 0
#define SKIPPED_UPDATE  -1


typedef struct {  
  /* records for every molecule the last ENTER event into the cavity */
  atom_id  id;             /* redundant but clean */
  real     t;              /* enter time in ps */
  int      u;              /* +1: flow parallel to axis, -1: anti parallel
			       0: not initialised */
} t_crossing;

typedef struct {
  /* flux results */
  real tau;              /* (average) residency/crossing time */
  int total;             /* |up| +  |down| */           
  int net;               /*  up  +   down  */           
  int up;                /* flux in +1 direction (>0) */
  int down;              /* flux in -1 direction (<0) */
} t_flux;

typedef struct {
  /* all results in one place */
  real time;             /* time of the current frame in ps */
  real last_update_t;    /* the last time at which update_result() was
                            called */
  t_flux  PHI;            /* flux for this time frame */
  t_flux  CPHI;           /* cumulated PHI; contains total at the end */
  int     ncr;            /* number of entries in crossed */
  atom_id *crossed;       /* all particles that crossed the pore */
  real    **residency;     /* histogram of residency and permeation times 
			     [center of bin] [tau_r] [tau_cross] */
  real    resDelta;       /* bin width of the histogram */
} t_result;

void init_t_result (t_result *, int);
void init_t_flux   (t_flux *);
int update_result (t_result *);
static gmx_inline bool bCrossing(const rvec, const rvec, const t_cavity *);
static gmx_inline int  typeCrossing (const rvec, const t_cavity *);
real pore_crossing (const rvec, const rvec, const real, 
		    t_crossing *, t_cavity *, t_result *);
int store_crossed (const atom_id, t_result *);

void init_t_result (t_result *res, int ncrmax) {
  int i;

  res->time = 0;
  res->last_update_t = -1;     /* must be different from time !! */
  init_t_flux(&(res->PHI));
  init_t_flux(&(res->CPHI));
  res->ncr = 0;
  snew(res->crossed,ncrmax);
  res->residency=grid2_alloc(RESIDENCY_BINS,3);
  res->resDelta = (real)RESIDENCY_TMAX/(real)RESIDENCY_BINS;
  /* load first column with mean value of the bin... laziness rules 
     (note that the array comes zero-initialised !)
   */
  for(i=0;i<RESIDENCY_BINS;i++) 
    res->residency[i][0] = (i+0.5)*res->resDelta;
  return;
};

void init_t_flux (t_flux *f) {
  f->tau   = 0;
  f->total = 0;
  f->net   = 0;
  f->up    = 0;
  f->down  = 0;
  return;
};

int update_result (t_result *result) {
  /* update_result() may be called more than once but changes occur
     only at the very first call. This is not very robust---I
     cannot call update_result repeatedly because I will loose all
     further changes that my occur during this time frame. There is
     nothing to do beyond simply splitting it in two routines: one for
     CPHI and the other for PHI. Or an additional flag to finalise?
  */
  if (result->time == result->last_update_t) return SKIPPED_UPDATE;

  if ((result->PHI.total = abs(result->PHI.up) + abs(result->PHI.down))) {
    result->CPHI.tau   += result->PHI.tau;

    result->PHI.tau    /= (real) result->PHI.total;
    result->PHI.net     = result->PHI.up + result->PHI.down;

    result->CPHI.total += result->PHI.total;
    result->CPHI.net   += result->PHI.net;
    result->CPHI.up    += result->PHI.up;
    result->CPHI.down  += result->PHI.down;
  };
  result->last_update_t = result->time;
  return result->PHI.total;
};



static gmx_inline bool bCrossing(const rvec x, const rvec x_last, 
				 const t_cavity *cavity) {
  real zmin, zmax;
  /* particle crossed a boundary 1 or 2  if
     zmin < z1|2 < zmax
  */
  zmin = MIN(x[ZZ],x_last[ZZ]);
  zmax = MAX(x[ZZ],x_last[ZZ]);
  return (zmin < cavity->z1 && cavity->z1 < zmax) ||
         (zmin < cavity->z2 && cavity->z2 < zmax);
};

static gmx_inline int typeCrossing (const rvec x, const t_cavity *cavity) {
  /* +1: into the cavity   Out->In, 
     -1: out of the cavity In->Out  */
  return (cavity->z1 < x[ZZ] && x[ZZ] < cavity->z2) ? 1 : -1;
};

 

real pore_crossing (const rvec x, const rvec x_last, const real t, 
		   t_crossing *inflx, t_cavity *cavity, t_result *res) {
  int u;           /* +1: flow parallel pore axis, -1: anti parallel   */
  int divergence;  /* +1: flow into the pore,      -1: out of the pore */
  real tau;
  rvec dx;
  rvec *axis;
  int ires;        /* index in the residency time histogram */

  tau = 0;         /* crossing time tau > 0 signals a crossing */

  /*
  if (bDebugMode()) {
    fprintf (debug, "[t=%8.1f mol_id=%5d] "
	            "---- pore_crossing():  DEBUG  ----\n", 
	     t, inflx->id);
    fprintf (debug, 
	     "   x=(%g,%g,%g) x_last=(%g,%g,%g)\n",
	     x[XX],x[YY],x[ZZ],x_last[XX],x_last[YY],x_last[ZZ]);
  };
  */

  if (bCrossing(x, x_last, cavity)) {
    divergence = typeCrossing (x, cavity);

    rvec_sub(x, x_last, dx);
    axis = &(cavity->axis);
    u = (iprod(dx, *axis) > 0 ? 1 : -1);

    if (bDebugMode()) {
      fprintf (debug, "[t=%8.1f mol_id=%5d] ---- boundary crossing ----\n", 
	       t, inflx->id);
      fprintf (debug, 
	       "   x=(%6.4f,%6.4f,%6.4f) x_last=(%6.4f,%6.4f,%6.4f) "
	       "dx=x-x_last=(%6.4f,%6.4f,%6.4f)\n",
	       x[XX],x[YY],x[ZZ],x_last[XX],x_last[YY],x_last[ZZ],
	       dx[XX],dx[YY],dx[ZZ]);
      fprintf (debug, 
	       "   divergence = %d     u = %d    (axis=(%4f,%4f,%4f))\n",
	       divergence, u, (*axis)[XX], (*axis)[YY], (*axis)[ZZ]); 
    };

    if (divergence == 1 ) {
      /* store the ENTER event */
      inflx->t = t;
      inflx->u = u;
    } else {
      /* check for a pore crossing.  If uninitialised (molecule
	 started in cavity), u==0, and thus this test will fail 
         (and the residency time is computed from t=0) */
      tau = (u != 0) ? t - inflx->t : t;
      ires = (int)floor(tau/res->resDelta);  /* index in residency histogram */
      if (ires > RESIDENCY_BINS - 2) {
	ires=RESIDENCY_BINS-1; 
	dmsg("residency histogram: [t=%g] tau=%g > MAX=%g\n", t, tau, RESIDENCY_TMAX); 
	dmsg("residency histogram:        ires=%i  u=%i\n",ires,u);
      } 
      
      if (u == inflx->u) {
	/* complete pore crossing ! */
	t_flux *PHI;
	
	PHI = &(res->PHI);
	PHI->tau  += tau;
	if (u>0) { PHI->up++;}
	else     { PHI->down--;};
	/* CPHI, PHI->tau, PHI->net, PHI->total are updated en bloc in
           main()->update_result() at the end of the frame. */
	list_add_atomid (inflx->id, &(res->ncr), res->crossed);
	/* build histogram of crossing times */
	++(res->residency[ires][2]);
      };
      /* build histogram of residency times (regardless if crossed or
         exited the same way) and crossing times */
      ++(res->residency[ires][1]);

      /* erase (t,u): a valid crossing requires a new ENTER event */
      inflx->t = 0;
      inflx->u = 0;
    };
  };
  return tau;
};
  

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]g_flux[TT] calculates the flux of molecules (or atoms) through "
    "a  cylindrical region as a function of time. It takes an index file "
    "with atomnumbers and chooses the molecules that these atoms belong to "
    "(or optionally (with [TT]-nom[TT]) it simply takes all atoms in the "
    "index) and generates an output file with the "
    "number of molecules/atoms at each time frame.[PAR]",
    "Output:\n"
    "=======-----------------------------------------------------------\n"
    " PHI+,  PHI-              flux parallel/anti-parallel to pore axis\n"
    "                          (PHI+ > 0, PHI- < 0)\n"
    " PHI   =  PHI+  +  PHI-   net flux\n"
    "|PHI|  = |PHI+| + |PHI-|  total flux\n"
    " tau = Sum tau_i/|PHI|    (average) residency/crossing time\n" 
    "\n"
    "                            t\n"
    "Cumulated flux C[PHI(t)] = Sum PHI(t'),\n"
    "                           t'=0\n"
    "                 <tau>(t)= C[tau]/|C[PHI]|\n"
    "------------------------------------------------------------------\n"
    "-dat:  (time t in ps)\n"
    "   t PHI(t) PHI+(t) PHI-(t) |PHI(t)| tau \n"
    "-cdat:\n"
    "   t C[PHI] C[PHI+] C[PHI-] |C[PHI]| <tau>(t)[PAR]",
    "Suggested use for water:\n",
    "create an index file for the water molecules:\n",
    "  [TT]echo -e \"keep 0\\ndel 0\\nr SOL\\nq\\n\" ",
    "| make_ndx -f in.pdb -o sol.ndx[TT]\n",
    "and run [TT]g_flux -m[TT] on it ([TT]-m[TT] is the default)."
  };

  static char *bugs[] = { 
    "g_flux is currently ALPHA development.",
    "-m behaves different from the standard usage "
    "within the g_* programs -- it figures out _for itself_ what the "
    "molecules are and does not need MOLECULE numbers but ATOM_IDs.",
    "-m is the DEFAULT behaviour!",
    "The DEFAULT volume/radius correction is specific for MY setup",
    "Despite appearance it only makes sense to specify a pore axis "
    "approximately parallel to the z-axis because we only really base "
    "the definition of the boundaries on z coordinates",
    "XMGrace files have the wrong suffix (xvg instead of agr)",
  };

  bool bMolecular = TRUE;   /* default is to use molecules    */
  t_cavity   cavity = {   /* define volume to count mols in */
    {0, 0, 1},            /* axis     */
    {0, 0, 0},            /* cpoint   */
    -1,                   /* radius   */
    -1,                   /* z1 < z2 */
    -1,
    0                    /* volume - calculate later */
  };
  real dR = DR_DEFAULT;   /* volume correction */

  t_pargs pa[] = {
    { "-m",      FALSE, etBOOL, {&bMolecular},
      "index contains atoms, but g_flux counts the molecules"},
    { "-axis",   FALSE, etRVEC, {&(cavity.axis)},
      "Vector pointing parallel to the pore axis"},
    { "-cpoint", FALSE, etRVEC, {&(cavity.cpoint)},
      "HIDDENPoint on the pore axis"},
    { "-R",      FALSE, etREAL, {&(cavity.radius)},
      "HIDDENRadius of the cylindrical pore cavity"},
    { "-z1",     FALSE, etREAL, {&(cavity.z1)},
      "Confine waters to have a z coordinate larger than z1, and"},
    { "-z2",     FALSE, etREAL, {&(cavity.z2)},
      "smaller than z2"},
    { "-dR",     FALSE, etREAL, {&dR},
      "HIDDENCorrect water-accessible cavity volume" }
  };

  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTPS, NULL, NULL, ffREAD },
    { efNDX, NULL, NULL, ffOPTRD },
    { efXVG, "-o",    "flux", ffWRITE },
    { efDAT, "-dat",  "flux", ffWRITE },
    { efDAT, "-cdat", "cflux", ffWRITE },
    { efNDX, "-ncr",  "crossed", ffWRITE },
    { efDAT, "-res", "residency", ffWRITE }
  };

  FILE       *out;           /* xmgr file with raw numbers */
  FILE       *fData;         /* current  flux    */
  FILE       *fCData;        /* accumulated flux */
  FILE       *fncr;          /* index of crossing particles */
  FILE       *fres;          /* histogram of residency and crossing times */
  t_topology *top;           /* topology                   */
  rvec       *x,*x_s;        /* coord, with and without pbc*/
  rvec       *x_last;        /* coord of last frame with pbc for
                                molecules (COM) or atoms */
  rvec       *xmol_cm;       /* COM coordinates for molecules */
  matrix     box;            /* box matrix (3x3)           */
  real       t;              /* time  */
  t_simtime  stime;          /* important times in th simulation */
  int        natoms;         /* number of atoms in system  */
  int        status;
  int        i;              /* loopcounters                 */
  int        tot_frames = 0; /* t_tot = tot_frames * delta_t */
  char       *grpname;       /* name of the group            */
  int        gnx;            /* number of atoms in group*/
  int        gnmol;          /* number of molecules in group */
  atom_id    *molndx;        /* index of mols in atndx */
  atom_id    *index;         /* atom numbers in index file  */
  t_block *mols;             /* all molecules in system    */
  t_crossing *influx;        /* array recording all enter events 
				for mols/atoms */
  t_result   result;         /* flux for current t and total flux,
                                particles that crossed the pore */

  char       s_tmp[STRLEN];  /* string for use with sprintf() */

#define NFILE asize(fnm)
#define NPA   asize(pa)

  CopyRight(stderr,argv[0]);
  fprintf (stderr, "\nVersion: %s\n\n", SRCID_g_flux_c);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);

  if (bDebugMode()) {
    dfprintf ("%s -- debugging...\n\n", Program());
  };
  
  /* open input files */

  top=read_top(ftp2fn_null(efTPS,NFILE,fnm));
  get_index(&(top->atoms),ftp2fn_null(efNDX,NFILE,fnm),
	    1,&gnx,&index,&grpname);
  
  /* step size in ps from tpx file 
  dt = dt_tpx (ftp2fn_null(efTPS,NFILE,fnm));
  dmsg ("Read from the topology: stepsize delta_t = %3g ps\n", dt);
  */

  /* get info about topology etc */
  mols=&(top->blocks[ebMOLS]);

  /* construct array molndx of mol indices in atndx if -m is set */
  molndx = NULL;
  gnmol = -1;
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
  };
  

  if (bMolecular) {
    init_t_result(&result,gnmol);
    snew (xmol_cm,gnmol);
    snew (x_last, gnmol);
    snew (influx, gnmol);
    for (i=0; i<gnmol; i++) {
      /*
	influx[i].id  id of the molecule 
	              ??? as seen in the input tpr/xtc file???
	influx[i].t records the time when molecule id entered the pore.
	influx[i].u == 0 is the sign for an uninitialised crossing!
	(valid values for u: +1, -1) 
      */
      influx[i].id = molndx[i];
      influx[i].t  = 0;
      influx[i].u  = 0;
    };
  } else {
    init_t_result(&result,gnx);
    xmol_cm = NULL;
    snew (x_last, gnx);
    snew (influx, gnx);
    for (i=0; i<gnx; i++) {
      influx[i].id = index[i];
      influx[i].t  = 0;
      influx[i].u  = 0;
    };
  };

  natoms = read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  snew(x_s,natoms);

  /* look at cavity 
     set defaults from the simulation box size
  */
  msg("\n");
  autoset_cavity(&cavity,box,NPA,pa);
  cavity.vol = sqr(cavity.radius - dR) * PI * (cavity.z2 - cavity.z1);

  /* open output files */

  sprintf (s_tmp, "R=%4g nm, z\\s1\\N=%4g nm, z\\s2\\N=%4g nm,"
	   " V=%4g nm\\S3\\N, indexgroup %s",
	     cavity.radius, cavity.z1, cavity.z2, cavity.vol,
	     grpname);
  out = xmgropen (opt2fn("-o",NFILE,fnm),
		 "Flux",
		 s_tmp,
		 "Time [ps]",
		 bMolecular ? "number of molecules" : 
		 "number of atoms");
  fData   = ffopen (opt2fn("-dat",  NFILE, fnm), "w");  
  fCData  = ffopen (opt2fn("-cdat", NFILE, fnm), "w");  
  fncr    = ffopen (opt2fn("-ncr",  NFILE, fnm), "w");  

  /* start counting .... */

  /* fill x_last from the first frame  */
  if (bMolecular) {
    x2molxcm (top, molndx, gnmol, box, x, x_s, x_last);
  } else {
    for(i=0; i<gnx; i++) {
      copy_rvec(x[index[i]], x_last[i]);
    };
  };

  /* ... and advance */
  if (! read_next_x(status,&t,natoms,x,box)) {
    msg ("Only one frame in trajectory.\nGAME OVER!\n");
    /* clean up a bit */
    fprintf(stderr,"\n");
    close_trj(status);
    fclose(out);
    fclose(fData);
    fclose(fCData);
    fclose(fncr);
    exit(EXIT_FAILURE);
  };

  /* 
     main loop over frames 
  */
  do {
    /* set initial values of the result struct 
       result.CPHI accumulates
     */
    result.time = t;
    init_t_flux(&(result.PHI));

    if (bMolecular) {
      /* prepare COM coordinates for molecules in xmol_cm */
      x2molxcm (top, molndx, gnmol, box, x, x_s, xmol_cm);

      for (i = 0; i < gnmol; i++) {
	pore_crossing (xmol_cm[i], x_last[i], t, 
		       &(influx[i]), &cavity, &result);
	copy_rvec(xmol_cm[i], x_last[i]);
      };
      /* End loop over all molecules */
    } else {
      for(i=0; i<gnx; i++) {
	pore_crossing (x[index[i]], x_last[i], t, 
		       &(influx[i]), &cavity, &result);
	copy_rvec(x[index[i]], x_last[i]);
      }
      /* End loop over all atoms */
    } 
    
    update_result (&result);

    {
      static char s_PHI[STRLEN], s_CPHI[STRLEN];
      /* output format:
	 t[ps] PHI(t) PHI+(t) PHI-(t) |PHI|  tau 
	 t     C(PHI) C(PHI+) C(PHI-) |CPHI| <tau>
      */
      
      sprintf (s_PHI, "%6d %6d %6d %6d  %6g", 
	       result.PHI.net, result.PHI.up, result.PHI.down, 
	       result.PHI.total, result.PHI.tau);

      sprintf (s_CPHI, "%6d %6d %6d %6d  %6g", 
	       result.CPHI.net, result.CPHI.up, result.CPHI.down, 
	       result.CPHI.total, 	     
	       result.CPHI.total > 0 ? 
	          result.CPHI.tau/((real) result.CPHI.total) : 0); 
      
      fprintf (out,   "%6g %s %s\n", result.time, s_PHI, s_CPHI);
      fprintf (fData, "%6g %s\n",    result.time, s_PHI);
      fprintf (fCData,"%6g %s\n",    result.time, s_CPHI);
    };

    /* tot_frames++; */
  } while(read_next_x(status,&t,natoms,x,box));

  /* post production 
     ---------------
     (1)  sort and write index of all particles that crossed 
  */
  quicksort (result.crossed, 0, result.ncr-1);
  sprintf(s_tmp, "%s_crossed", grpname);
  fwrite_index(fncr, result.crossed, bMolecular ? etxMOL : etxATOM,
	      result.ncr, top, s_tmp);

  /* histogram of residency and crossing times */
  fres = ffopen (opt2fn("-res",  NFILE, fnm), "w");    
  fprintf(fres,"#%7s  %9s  %9s\n","[ps]","N(tau_res)","N(tau_cross)");
  for(i=0;i<RESIDENCY_BINS;i++) 
    fprintf(fres,"%8.2f  %9i  %9i\n",result.residency[i][0],
	    (int)result.residency[i][1],(int)result.residency[i][2]);

  /* clean up a bit */
  fprintf(stderr,"\n");
  close_trj(status);
  fclose(out);
  fclose(fData);
  fclose(fCData);
  fclose(fncr);
  fclose(fres);

  thanx(stdout);
  
  return 0;
}

