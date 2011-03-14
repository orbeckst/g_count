/*
 * $Id: g_count.c,v 1.37 2009/07/02 15:15:12 oliver Exp $
 *
 * This program is based on a Gromacs 2.0 g_* program 
 * see http://www.gromacs.org
 *
 */
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "rdgroup.h"
#include "mshift.h"
#include "xvgr.h"
#include "gstat.h"
#include "princ.h"
#include "names.h"
#include "gmx_fatal.h"
#include "xf.h"
#include "utilgmx.h"
#include "count.h"

/* tau = minimum time spent in cavity to be tracked, assumed to be ps */
#define TAU_DEFAULT   10
/* empirical radius adjustment for cylindrical volume calculation */
#define DR_DEFAULT    0.030

/* "times" are in frames */
typedef struct {      /* record all atoms that spend time in the cavity */
  atom_id   id;       /* number of molecule or atom */
  int       time;     /* time spent so far (measured in (int) frames) */
} t_atmcav;



int track_ndx(atom_id atom, 
	       t_atmcav *tracked, int natoms);
gmx_bool bTracked(t_atmcav *t, int tau_frames);
void do_tracking(FILE *fTrack, FILE *fTDat, t_atmcav *trckd, 
		 t_simtime *t, t_topology *top, char *grpname);
void dump_a_cav(FILE *fp, t_atmcav *trckd, 
		atom_id mxcav, t_simtime *t,
		t_topology *top);


void do_tracking(FILE *fTrack, FILE *fTDat, t_atmcav *trckd, 
		 t_simtime *t, t_topology *top, char *grpname)
{
  t_atmcav *tr;

  int i, j;          /* loop */
  int ncav;          /* number of molecules ever in the cavity */
  int ntrck;         /* number of mol/atoms tracked (ie t > tau) */

  atom_id *atndx = NULL;     /* indices from block.        */
  t_block *mols;             /* all molecules in system    */
  t_atom  *atoms;            /* all atoms */
  t_resinfo *resinfo;        /* resnames and numbers */
		       
  /* initialize, 
     use NULL to ask for the current max index private to track_id 
  */
  ncav  = track_ndx(NO_ATID, NULL, 0);
  ntrck = 0;

  msg ("Writing index of all %d atoms in the cavity.\n", ncav);

  /* needed for translation mols -> atom ids */
  mols  = &(top->mols);
  atndx = mols->index;

  /* translation atom_id --> resnr: 
     resinfo[atoms[iatom].resind].nr
   */
  atoms = top->atoms.atom;
  resinfo = top->atoms.resinfo;
  
  // quicksort (trckd, 0, ncav-1);

  /* ids to new index file */
  fprintf (fTrack, "[ %s_tracked ]\n", grpname);

  for (i = 0; i < ncav; i++) {
    tr = &(trckd[i]);
    if (bTracked(tr, t->tau_frames)) {
      ntrck++;
      /* ADD +1 when WRITING (external) index file !!! */
      fprintf(fTrack, "%5u %s", tr->id + 1, 
	      (ntrck % 15 == 0) ? "\n" : "");
    }
  }
  fprintf (fTrack, "\n\n");

  /* simply dump all entries (same as tau=0) */
  fprintf (fTrack, "[ %s_cavity ]\n", grpname);
  for (i = 0; i < ncav; i++) {
    tr = &(trckd[i]);
    /* ADD +1 when WRITING (external) index file !!! */
    fprintf (fTrack, "%5u %s", tr->id + 1, 
	     ((i+1) % 15 == 0) ? "\n" : "");
  }
  fprintf (fTrack, "\n\n");

  msg ("Number of atoms tracked (tau=%g): %d (%5.2f%% of all atoms "
       "in the cavity)\n", 
       t->tau,
       ntrck, (real) 100*ntrck/ncav);

  if (bDebugMode()) {
    fprintf (debug, "\n# Threshold for counting: tau = %g\n", t->tau);
    fprintf (debug, "# (count if  t_cav > tau ps)\n");
    fprintf (debug, "# Total simulation time: t_tot = %g ps\n", t->t_tot); 
    fprintf (debug, "# Total simulation frames: tot_frames = %d\n", 
	             t->tot_frames); 
    fprintf (debug, "# delta_t = t_tot/(tot_frames-1) = %g ps\n", t->delta_t);
  };

  /* dump all data to fTDat (possibly as input to g_track.pl) */
  /* write statistic: time spent in cavity */
  if (bDebugMode()) {
    fprintf (debug, "# Number of atoms tracked: %d (%5.2f%% of all atoms "
	     "in the cavity)\n", 
	     ntrck, (real) 100*ntrck/ncav);
    fprintf (debug, "# (For the ids internal to %s subtract 1! These are "
	     "ids that appear \n# in the pdb or gro or ndx files.)\n", 
	     Program());
    fprintf (debug, "# Tracking atoms\n# %8s %6s %6s %6s %12s %13s %s\n",
	     "tr->id+1", "a_id", "res_id", "name", "t_cav (ps)", 
	     "t_cav/t_tot", "tracked?");
    dump_a_cav (NULL, trckd, ncav, t, top);
  };
    
  fprintf (fTDat, "[ type=atoms t_tot=%g ]\n", t->t_tot);
  dump_a_cav(fTDat, trckd, ncav, t, top);

  return;    
};


void dump_a_cav (FILE *fp, t_atmcav *trckd, 
		 atom_id mxcav, t_simtime *t,
		 t_topology *top){
  int i;
  gmx_bool bDebugOutput;

  atom_id *atndx;        /* indices from block.        */
  t_block *mols;             /* all molecules in system    */
  t_atom  *atoms;            /* all atoms */		 
  t_resinfo *resinfo;        /* resnames and numbers */

  /* if fp == NULL  given AND debug file opened print debugging table */
  if ((bDebugOutput = !fp)) {
    if (!debug) fp = stderr;
  };

  /* needed for translation mols -> atom ids */
  mols  = &(top->mols);
  atndx = mols->index;
  /* translation atom_id --> resnr 
     resinfo[atoms[iatom].resind].nr
  */
  atoms = top->atoms.atom;
  resinfo = top->atoms.resinfo;

  for (i = 0; i < mxcav; i++) {
    atom_id igmx;
    t_atmcav *tr;
    
    tr = &(trckd[i]);
    igmx = tr->id; 
    
    if (!bDebugOutput) {
    fprintf (fp, "%6u %6u %6s %12.1f %13.6f\n",
	     igmx + 1,
	     resinfo[atoms[igmx].resind].nr + 1,
	     *(top->atoms.atomname[igmx]),
	     tr->time * t->delta_t,
	     (real)tr->time / t->tot_frames);
    } else {
    fprintf (debug, "  %8u %6u %6u %6s %12.1f %13.6f %s\n",
	     tr->id + 1,
	     igmx + 1,
	     resinfo[atoms[igmx].resind].nr + 1,
	     *(top->atoms.atomname[igmx]),
	     tr->time * t->delta_t,
	     (real)tr->time / t->tot_frames,
	     bTracked (tr, t->tau_frames) ?  "#TRACKED#" : "");
    };
  };
  return;
};


/* store an atom_id in tracked[] unless it is already there. always
   return the position of atom in tracked[] 
   
   calling with atom=NO_ATID, tracked[] == NULL returns the current
   number of entries

   calling with a atom already in tracked[] returns its index -- only
   use this if you know that it is already in there, otherwise you
   will put it in */

int track_ndx (atom_id atom, 
	       t_atmcav tracked[], int natoms) {
  int  i;
  static int ntr;   /* remember number of tracked atoms  */

  /* track_id (type, NO_ATID, NULL, ...) returns current total ntr */
  if (!tracked && atom == NO_ATID)   return ntr;
  
  if (ntr > natoms) {
    gmx_fatal(FARGS, "track_ndx(): number of tracked atoms"
		 " ntr=%d exceeds total number of atoms %d.\n",
		 ntr, natoms);
  };

  for (i = 0; tracked[i].id != atom && i < ntr; i++);
  if (i >= ntr) 
    tracked[ntr++].id = atom;

#ifdef DEBUG
  dfprintf ("track_ndx(): atomid = %u, ntr = %d.\n",
	    atom, ntr);
#endif

  return i;
};


gmx_bool bTracked (t_atmcav *atm, int tau_frames) {
  return atm->time > tau_frames;
};


int main(int argc,char *argv[])
{
  const char *desc[] = {
    "[TT]g_countx[TT] counts the numbers of molecules within a cylindrical ",
    "region as a function of time. It takes an index file with ",
    "atomnumbers and generates an output file with the number of atoms "
    "which are present in the cavity a each time frame.\n",
    "The data file contains\n",
    "        t(ps) number concentration(mol/l) [density(g/cm^3)][PAR]",
    "Density and concentration are calculated from the volume of the cavity, ",
    "V = (R*)^2 * pi * (z2 - z1), where R* is the apparent ",
    "accessible radius (corrected with wall atom radius r_W, water "
    "radius r_water, and d(wall-water) from the Lennard-Jones potential:\n "
    "        R* = rho - d(wall-water) + r_water \n"
    "and rho the radius of the pore, measured from the atomic center of a "
    "wall atom."
    "Typically, R = rho - r_wall, so the radius correction is finally\n"
    "        dR = r_wall - d(wall-water) + r_water\n"
    "Thus, the density is calculated for the ",
    "_accessible_ volume.[PAR]"
    "The -track index file contains the atom_ids  of all atoms that",
    "were at least for tau ps in the cavity (entry [ *_tracked ]). ",
    "The entry [ *_cavity ] contains all atoms that were for at ",
    "least one time frame in the cavity.[PAR]",
    "[TT]-tdat[TT] writes a data file containing all atoms in "
    "the cavity and their total residency time. Format:\n",
    "         atom_nr  molecule_nr  name   total_time (ps)  T_sim/total_time"
    "[PAR]",
    "Defaults:\n"
    "If no values are given for cpoint then the center of the initial box is "
    "taken. "
    "If no radius is given, the maximum radius fitting in the box is chosen. "
    "z1 and z2 default to the full box."
    "[PAR]"
    "Suggested use for water:\n",
    "create an index file for the water oxygens:\n",
    "  [TT]echo -e \"keep 0\\ndel 0\\nt OW\\nq\\n\" ",
    "| make_ndx -f in.pdb -o ow.ndx[TT]\n",
    "and run [TT]g_count -nom[TT] on it ([TT]-nom[TT] is the default)."
  };

  const char *bugs[] = { 
    "Counting molecules (hidden option -m) does not work; it used to work "
    "in older versions of g_count (pre-Gromacs 4.0)",
    "You should only supply water oxygens for water. "
    "However, this means that the mass-density is wrong because hydrogens are not taken into account "
    "so you should either just use concentrations/number densities or correct the mass yourself.",
    "The density is calculated as (total mass in cavity)/((R*-r_water)^2*pi*(z2-z1)) -- "
    "so it makes only sense if this approximates the true cavity volume. "
    "R* = R - dR, where dR is the volume/radius correction.",
    "The DEFAULT volume/radius correction dR is specific for methane pseudo "
    "atoms (ffgmx, ffG43a1) and SPC water.",
    "The program is only tested with pore axis parallel to z-axis (0,0,1).",
  };

  gmx_bool bMolecular = FALSE;   /* default is to use atoms    */
  real tau = TAU_DEFAULT;   /* t_cav > tau => mol tracked  */
  t_cavity   cavity = {   /* define volume to count mols in */
    {0, 0, 1},            /* axis     */
    {0, 0, 0},            /* cpoint   */
    -1,                   /* radius   */
    -1,                   /* z1 < z2 */
    -1,
    0                    /* volume - calculate later */
  };
  real dR = DR_DEFAULT;   /* effective radius R* = R - dR [PERSONAL DEFAULT]*/
  int  ngroups = 1;  /* not used >1 */
  t_pargs pa[] = {
    { "-m",      FALSE, etBOOL, {&bMolecular},
      "HIDDENindex contains atoms, but g_count counts the molecules (NOT WORKING?!?)"},
    { "-axis",   FALSE, etRVEC, {&(cavity.axis)},
      "Vector pointing parallel to the pore axis"},
    { "-cpoint", FALSE, etRVEC, {&(cavity.cpoint)},
      "Point on the pore axis"},
    { "-R",      FALSE, etREAL, {&(cavity.radius)},
      "Radius of the cylindrical pore cavity"},
    { "-z1",     FALSE, etREAL, {&(cavity.z1)},
      "Confine waters to have a z coordinate larger than z1, and"},
    { "-z2",     FALSE, etREAL, {&(cavity.z2)},
      "smaller than z2"},
    { "-tau",    FALSE, etREAL, {&tau},
      "Total time (ps) that a particle has to be in the cavity at "
      "least so that it is tracked"},
    { "-dR",     FALSE, etREAL, {&dR},
      "Correct water-accessible cavity volume" },
    { "-ng",       FALSE, etINT, {&ngroups},
      "HIDDENNumber of groups to consider" },    
  };
  FILE       *out;           /* xmgr file with raw numbers */
  FILE       *fData;         /* t n c rho     :generic nxy data file      */
  FILE       *fTDat;         /* atom t_cav    :generic nxy data file      
			        input for g_track.pl */
  FILE       *fDens;         /* xmgr density */
  FILE       *fConc;         /* xmgr concentration */
  FILE       *fTrack;        /* index file of all molecules that were
                                ever in the cavity */
  t_topology top;            /* topology                   */
  int        ePBC;
  char       title[STRLEN];
  rvec       *xtop;
  gmx_bool   bTop;
  rvec       *x,*x_s;        /* coord, with and without pbc*/
  rvec       xcm;            /* center of mass of molecule */
  matrix     box;            /* box matrix (3x3)           */
  real       t,tm;           /* time, total mass molecule  */
  t_simtime  stime;          /* important times in th simulation */
  real       totalmass;      /* mass of all mol./atoms within cavity */
  int        natoms;         /* number of atoms in system  */
  int        nmolecules;     /* number of molecules in system  */
  t_trxstatus *status;
  int        i,j;            /* loopcounters                 */
  int        tot_frames = 0; /* t_tot = tot_frames * delta_t */
  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTPS, NULL, NULL, ffREAD },
    { efXVG, "-o",    "c_numbers", ffWRITE },
    { efXVG, "-dens", "c_density", ffOPTWR },
    { efXVG, "-conc", "c_conc", ffOPTWR },
    { efDAT, "-dat", "cavity", ffWRITE },
    { efDAT, "-tdat", "c_track", ffWRITE },
    { efNDX, NULL, NULL, ffOPTRD },
    { efNDX, "-track", "c_track", ffOPTWR }
  };
  output_env_t oenv;
  int        nmol;           /* number of molecules within the cavity */
  real       conc, dens;     /* concentration and mass density in cavity */
  int        tndx;           /* index in cav_id[][] */
  t_atmcav   *trckd;         /* all data about atoms in the cavity */
  char       s_tmp[STRLEN];  /* string for use with sprintf() */
  char       s_title[STRLEN];/* and another one... */  
  /* from gmx_traj.c */
  const char *indexfn = NULL;
  char       **grpname;
  int        *isize0,*isize;
  atom_id    **index0,**index;
  atom_id    *atndx = NULL;
  t_block    *mols= NULL;
  /* old 3.3.x and modified */
  char       *ggrpname;      /* name of the FIRST group       */
  int        gnx = 0;        /* number of atoms in FIRST group*/
  atom_id    *gindex = NULL; /* index of FIRST group */
  int        gnmol = 0;      /* XXX number of molecules in group */
  atom_id    *molndx = NULL; /* XXX index of mols in atndx */
  t_atom     *atoms;         /* XXX replaces 'a' ???? */

#define NFILE asize(fnm)
#define NPA   asize(pa)

  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs,
		    &oenv);

  if (bDebugMode()) {
    dfprintf ("%s -- debugging...\n\n", Program());
  };
  
  /* open input files */

  /* old:   top=read_top(ftp2fn_null(efTPS,NFILE,fnm)); */
  bTop = read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
  sfree(xtop);

  if (!bTop) {
    gmx_fatal(FARGS, "Need a run input file");
  }

  if (bMolecular) {
    gmx_fatal(FARGS, "Sorry, -m option not working at the moment");
  }
  else {
    indexfn = ftp2fn_null(efNDX,NFILE,fnm);
  }

  if (ngroups != 1) {
    gmx_fatal(FARGS, "Sorry, only a single group currently allowed.");
  }

  snew(grpname,ngroups);
  snew(isize0,ngroups);
  snew(index0,ngroups);

  get_index(&(top.atoms),indexfn,ngroups,isize0,index0,grpname);

  // remnant from old code that also looked at molecules; TODO: cleanup
  isize = isize0;
  index = index0;

  /* ngroups == 1 at moment */
  gnx = isize[0];
  gindex = index[0];
  ggrpname = grpname[0];
  atoms = top.atoms.atom;
  mols = &(top.mols);
  atndx = mols->index;
  
  natoms = read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  snew(x_s,natoms);

  /* look at cavity 
     set defaults from the simulation box size
  */
  msg("\n");
  autoset_cavity(&cavity,box,NPA,pa);

  cavity.vol = sqr(cavity.radius - dR) * PI * (cavity.z2 - cavity.z1);
  msg ("Water accessible volume of cylindrical cavity V=%6g nm^3\n"
       "(solvent accessible volume with correction dR = %g nm)\n", 
       cavity.vol, dR);

  /* check tau */
  if (tau < 0) {
    msg ("Warning:  tau = %4g must be >= 0. "
	 "Using the default tau = %4g.", tau, TAU_DEFAULT);
    tau = TAU_DEFAULT;
  };

  /* open output files */
  sprintf (s_tmp, "R=%4g nm, z\\s1\\N=%4g nm, z\\s2\\N=%4g nm,"
	   " V=%4g nm\\S3\\N, indexgroup %s",
	     cavity.radius, cavity.z1, cavity.z2, cavity.vol,
	     ggrpname);
  sprintf(s_title, "%s in cavity", ggrpname); 
  out = xmgropen (opt2fn("-o",NFILE,fnm),
                  s_title, 
                  s_tmp, 
                  "Time [ps]",
                  bMolecular ? "number of molecules" : "number of atoms");
  fData  = ffopen (opt2fn("-dat",    NFILE, fnm), "w");
  fTDat  = ffopen (opt2fn("-tdat",   NFILE, fnm), "w");
  fTrack = ffopen (opt2fn("-track",  NFILE, fnm), "w");
  fConc  = xmgropen (opt2fn("-conc", NFILE, fnm), 
                     s_title,
		     s_tmp,
		     "Time [ps]",
		     "concentration [mol/l]");
  fDens = xmgropen (opt2fn("-dens", NFILE, fnm), 
                    s_title,
		    s_tmp,
		    "Time [ps]",
		    "density [g/cm\\S3\\N]");
  


  /* store all atom_ids (and molecule numbers) that are ever in the cavity:
  */
  nmolecules = mols->nr; 
  snew(trckd, natoms);
  
  do {
    /* write time in ps */
    fprintf(out, "%10g", t);
    fprintf(fData, "%10g", t);
    fprintf(fConc, "%10g", t);
    fprintf(fDens, "%10g", t);

    nmol = 0;       /* restart counting at each timestep */
    totalmass = 0;  /* mass of all particles in the cavity */

    /* count atoms */
    for(i=0; i<gnx; i++) {
      if (bInCavity(x[gindex[i]], &cavity)) {
	nmol++;
	totalmass += atoms[gindex[i]].m;
	tndx = track_ndx(gindex[i],trckd,natoms);
	(trckd[tndx].time)++; 
      };
    }
    
    /* concentration = number density in mol l^-1 */
    conc = nmol/AvogadroConstant * 1e+24 / cavity.vol;

    /* mass density in g cm^-3 */
    dens = totalmass * AMU * 1e+24 / cavity.vol;

    fprintf (out,  "\t%5d\n", nmol);
    fprintf (fData, "\t%5d\t%6g\t%6g\n", nmol, conc, dens);
    fprintf (fDens, "\t%6g\n", dens);
    fprintf (fConc, "\t%6g\n", conc);

    tot_frames++;
  } while(read_next_x(oenv,status,&t,natoms,x,box));
  close_trj(status);

  /* initialize stime */
  stime.tau        = tau;
  stime.tot_frames = tot_frames;
  stime.t_tot      = t;
  stime.delta_t    = stime.t_tot/(stime.tot_frames-1);
  stime.tau_frames = (int)tau/stime.delta_t;  // compare to t_atmcav.time

  dmsg ("Simulation time: t_tot ==  (tot_frames-1) * delta_t\n"
	"                 %g ps  == %d * %g ps \n",
	stime.t_tot, stime.tot_frames, stime.delta_t); 

  /* [sort and] write tracking index (both if also looking for
     molecules). Only track molecules above tau threshold 
  */
  do_tracking(fTrack, fTDat, trckd, &(stime), &top, ggrpname);

  /* clean up a bit */
  fprintf(stderr,"\n");
  fclose(out);
  fclose(fData);
  fclose(fTrack);
  fclose(fConc);
  fclose(fDens);
  sfree(trckd);

  thanx(stdout);
  
  return 0;
}

