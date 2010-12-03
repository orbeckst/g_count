/*
 * $Id: g_count.c,v 1.35 2009/01/19 16:19:45 oliver Exp $
 *
 * This program is based on a Gromacs 2.0 g_* program 
 * see http://www.gromacs.org
 *
 */
static char *SRCID_g_count_c = "$Id: g_count.c,v 1.35 2009/01/19 16:19:45 oliver Exp $";

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

#define TAU_DEFAULT   10
#define DR_DEFAULT    0.030

typedef struct {      /* record all atoms that spend time in the cavity */
  atom_id   id;       /* number of molecule or atom */
  real      time;     /* time spent so far (measured in (int) frames) */
  real      t_last;   /* last time the particle was seen in the cavity */
  int       crosses;  /* number of enter and exit events */
} t_atmcav;



int track_ndx (enum ndxtype type, atom_id atom, 
	       t_atmcav tracked[], int natoms);
bool bTracked (t_atmcav *t, real t_tot, real tau);
void do_tracking (FILE *fTrack, FILE *fTDat, t_atmcav *trckd[], 
		  bool bMolecular,  t_simtime *t,
		  t_topology *top, char *grpname);
void dump_a_cav (FILE *fp, enum ndxtype type, t_atmcav *trckd, 
		 atom_id mxcav, t_simtime *t,
		 t_topology *top);


void do_tracking (FILE *fTrack, FILE *fTDat, t_atmcav *trckd[], 
		  bool bMolecular,  t_simtime *t, 
		  t_topology *top, char *grpname)
{
  enum ndxtype type;
  t_atmcav *tr;

  int i, j;          /* loop */
  int moleculesize; 
  int ncav[etxNR];   /* number of molecules ever in the cavity */
  int ntrck[etxNR];  /* number of mol/atoms tracked (ie t > tau) */

  atom_id *a, *atndx;        /* indices from block.        */
  t_block *mols;             /* all molecules in system    */
  t_atom  *atoms;            /* all atoms */
		       
  for (type = etxATOM; type < (bMolecular ? etxNR : etxMOL); type++) {
    /* initialize, 
       use NULL to ask for the current max index private to track_id 
    */
    ncav[type]  = track_ndx (type, NO_ATID, NULL, 0);
    ntrck[type] = 0;
  };

  if (bMolecular) {
  msg ("Writing index of all %d atoms and %d molecules in the cavity.\n",
       ncav[etxATOM], ncav[etxMOL]);
  } else {
  msg ("Writing index of all %d atoms in the cavity.\n",
       ncav[etxATOM]);
  };    

  /* needed for translation mols -> atom ids */
  mols  = &(top->blocks[ebMOLS]);
  a     = mols->a;
  atndx = mols->index;
  /* translation atom_id -> resnr 
     top->atoms->atom[atom_id]->resnr
  */
  atoms = top->atoms.atom;

  for (type = etxATOM; type < (bMolecular ? etxNR : etxMOL); type++) {
    /*    quicksort (trckd[type], 0, ncav[type]-1); */

    /* ids to new index file */
    fprintf (fTrack, "[ %s_tracked%s ]\n", grpname,
	     type == etxMOL ? "_molecules" : "");  

    /* if we have molecules then we use the molecule numbers instead
       of the atoms; this makes sure that we track whole molecules 
    */	     
    if (bMolecular && type == etxATOM) {
      for (i = 0; i < ncav[etxMOL]; i++) {
	tr = &(trckd[etxMOL][i]);
	/* for each tracked molecule ... */
	if (bTracked(tr, t->tot_frames, t->tau) ) {
 	  moleculesize = atndx[tr->id + 1] - atndx[tr->id];
	  dfprintf ("tracking trckd[etxMOL][%d].id = %d, moleculesize = %d, "
		"first atom_id in mol = %d\n", 
		i, tr->id, moleculesize, a[atndx[tr->id]]);
	  for (j = 0; j < moleculesize; j++) {
	    /* ... fetch all atoms of this molecule */
	    ntrck[type]++;
	    /* ADD +1 when WRITING (external) index file !!! */
	    fprintf (fTrack, "%5u %s", a[atndx[tr->id]+j] + 1,  
		     (ntrck[type] % 15 == 0) ? "\n" : "");
	  };
	};
      };
    } else {
	/* either etxMOL or etxATOM without molecules */
	for (i = 0; i < ncav[type]; i++) {
	  atom_id igmx;

	  tr = &(trckd[type][i]);
	  igmx = (type == etxMOL ? a[atndx[tr->id]] : tr->id);

	  if (bTracked ( tr, t->tot_frames, t->tau) ) {
	    ntrck[type]++;
	    /* ADD +1 when WRITING (external) index file !!! */
	    fprintf (fTrack, "%5u %s", atoms[igmx].resnr + 1, 
		     (ntrck[type] % 15 == 0) ? "\n" : "");
	  };
	};
    };

    fprintf (fTrack, "\n\n");


    /* simply dump all entries (same as tau=0) */
    fprintf (fTrack, "[ %s_cavity%s ]\n", grpname,
	     type == etxMOL ? "_molecules" : "");  
    for (i = 0; i < ncav[type]; i++) {
      tr = &(trckd[type][i]);
      /* ADD +1 when WRITING (external) index file !!! */
      fprintf (fTrack, "%5u %s", tr->id + 1, 
	       ((i+1) % 15 == 0) ? "\n" : "");
    };
    fprintf (fTrack, "\n\n");

    msg ("Number of %s tracked (tau=%g): %d (%5.2f%% of all %s "
	 "in the cavity)\n", 
	 ENDXTYPE(type), t->tau,
	 ntrck[type], (real) 100*ntrck[type]/ncav[type],
	 ENDXTYPE(type) );
  };


  if (bDebugMode()) {
    fprintf (debug, "\n# Threshold for counting: tau = %g\n", t->tau);
    fprintf (debug, "# (count if  t_cav > tau ps)\n");
    fprintf (debug, "# Total simulation time: t_tot = %g ps\n", t->t_tot); 
    fprintf (debug, "# Total simulation frames: tot_frames = %d\n", 
	             t->tot_frames); 
    fprintf (debug, "# delta_t = t_tot/(tot_frames-1) = %g ps\n", t->delta_t);
  };

  /* dump all data to fTDat (possibly as input to g_track.pl) */
  for (type = etxATOM; type < (bMolecular ? etxNR : etxMOL); type++) {

    /* write statistic: time spent in cavity */
    if (bDebugMode()) {
      fprintf (debug, "# Number of %s tracked: %d (%5.2f%% of all %s "
	       "in the cavity)\n", 
	       ENDXTYPE(type), ntrck[type], (real) 100*ntrck[type]/ncav[type],
	       ENDXTYPE(type) );
      fprintf (debug, "# (For the ids internal to %s subtract 1! These are "
	       "ids that appear \n# in the pdb or gro or ndx files.)\n", 
	       Program());
      fprintf (debug, "# Tracking %s\n# %8s %6s %6s %6s %12s %13s %s\n",
	       ENDXTYPE(type),
	       "tr->id+1", "a_id", "res_id", "name", "t_cav (ps)", 
	       "t_cav/t_tot", "tracked?");
      dump_a_cav (NULL, type, trckd[type], ncav[type], t, top);
    };
    
    fprintf (fTDat, "[ type=%s t_tot=%g ]\n", ENDXTYPE(type), t->t_tot);
    dump_a_cav (fTDat, type, trckd[type], ncav[type], t, top);
  };

  return;    
};


void dump_a_cav (FILE *fp, enum ndxtype type, t_atmcav *trckd, 
		 atom_id mxcav, t_simtime *t,
		 t_topology *top){
  int i;
  bool bDebugOutput;

  atom_id *a, *atndx;        /* indices from block.        */
  t_block *mols;             /* all molecules in system    */
  t_atom  *atoms;            /* all atoms */		 

  /* if fp == NULL  given AND debug file opened print debugging table */
  if ((bDebugOutput = !fp)) {
    if (!debug) fp = stderr;
  };

  /* needed for translation mols -> atom ids */
  mols  = &(top->blocks[ebMOLS]);
  a     = mols->a;
  atndx = mols->index;
  /* translation atom_id -> resnr 
     top->atoms->atom[atom_id]->resnr
  */
  atoms = top->atoms.atom;

  for (i = 0; i < mxcav; i++) {
    atom_id igmx;
    t_atmcav *tr;
    
    tr = &(trckd[i]);
    igmx = (type == etxMOL ? a[atndx[tr->id]] : tr->id); 
    
    if (!bDebugOutput) {
    fprintf (fp, "%6u %6u %6s %12.1f %13.6f\n",
	     igmx + 1,
	     atoms[igmx].resnr + 1,
	     *(top->atoms.atomname[igmx]),
	     tr->time * t->delta_t,
	     tr->time / t->tot_frames);
    } else {
    fprintf (debug, "  %8u %6u %6u %6s %12.1f %13.6f %s\n",
	     tr->id + 1,
	     igmx + 1,
	     atoms[igmx].resnr + 1,
	     *(top->atoms.atomname[igmx]),
	     tr->time * t->delta_t,
	     tr->time / t->tot_frames,
	     bTracked (tr, t->tot_frames, t->tau) ?  "#TRACKED#" : "");
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

int track_ndx (enum ndxtype type, atom_id atom, 
	       t_atmcav tracked[], int natoms) {
  int  i;
  static int ntr[etxNR];   /* remember number of tracked atoms and
                              mols; remember them separately so that
                              we can call the same routine for atoms
                              and mols */

  /* track_id (type, NO_ATID, NULL, ...) returns current total ntr */
  if (!tracked && atom == NO_ATID)   return ntr[type];
  
  if (ntr[type] > natoms) {
    gmx_fatal(FARGS, "track_ndx(): number of tracked %s"
		 " ntr[%s]=%d exceeds total number of %s %d.\n",
		 ENDXTYPE(type), ENDXTYPE(type), ENDXTYPE(type),
		 ntr[type], natoms);
  };

  for (i = 0; tracked[i].id != atom && i < ntr[type]; i++);
  if (i >= ntr[type]) tracked[ntr[type]++].id = atom;

#ifdef DEBUG
  dfprintf ("track_ndx(): atomid = %u, ntr[type=%s] = %d.\n",
	    atom, ENDXTYPE(type), ntr[type]);
#endif

  return i;
};


bool bTracked (t_atmcav *t, real t_tot, real tau) {
  return t->time > tau;
};


int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]g_count[TT] counts the numbers of molecules within a cylindrical ",
    "region as a function of time. It takes an index file with ",
    "atomnumbers and ",
    "chooses the molecules that these atoms belong to ",
    "(or optionally (with [TT]-nom[TT]) it simply takes all atoms in the ",
    "index) and generates an output file with the number of molecules/atoms "
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
    "least one time frame in the cavity.",
    "Entries containing molecule numbers are titled ",
    "[ *_molecules ].[PAR]"
    "[TT]-tdat[TT] writes a data file containing all atoms or molecules in "
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
    "create an index file for the water molecules:\n",
    "  [TT]echo -e \"keep 0\\ndel 0\\nr SOL\\nq\\n\" ",
    "| make_ndx -f in.pdb -o sol.ndx[TT]\n",
    "and run [TT]g_count -m[TT] on it ([TT]-m[TT] is the default)."
    "[PAR]Known problems and CAVEATs:\n"
    "---------------------------",
  };

  bool bMolecular = TRUE;   /* default is to use molecules    */
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

  t_pargs pa[] = {
    { "-m",      FALSE, etBOOL, {&bMolecular},
      "index contains atoms, but g_count counts the molecules"},
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
      "Correct water-accessible cavity volume" }
  };
  FILE       *out;           /* xmgr file with raw numbers */
  FILE       *fData;         /* t n c rho     :generic nxy data file      */
  FILE       *fTDat;         /* atom t_cav    :generic nxy data file      
			        input for g_track.pl */
  FILE       *fDens;         /* xmgr density */
  FILE       *fConc;         /* xmgr concentration */
  FILE       *fTrack;        /* index file of all molecules that were
                                ever in the cavity */
  t_topology *top;           /* topology                   */
  rvec       *x,*x_s;        /* coord, with and without pbc*/
  rvec       xcm;            /* center of mass of molecule */
  matrix     box;            /* box matrix (3x3)           */
  real       t,tm;           /* time, total mass molecule  */
  t_simtime  stime;          /* important times in th simulation */
  real       totalmass;      /* mass of all mol./atoms within cavity */
  int        natoms;         /* number of atoms in system  */
  int        nmolecules;     /* number of molecules in system  */
  int        status;
  int        i,j;            /* loopcounters                 */
  int        tot_frames = 0; /* t_tot = tot_frames * delta_t */
  char       *grpname;       /* name of the group            */
  int        gnx;            /* number of atoms in group*/
  int        gnmol;          /* number of molecules in group */
  int        moleculesize;   /* size of molecule in numbers*/
  atom_id    *molndx;        /* index of mols in atndx */
  atom_id    *index;         /* atom numbers in index file  */
  atom_id    *a, *atndx;     /* indices from block.        */
  t_block *mols;             /* all molecules in system    */
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
  int        nmol;           /* number of molecules within the cavity */
  real       conc, dens;     /* concentration and mass density in cavity */
  int        tndx;           /* index in cav_id[][] */
  t_atmcav   *trckd[etxNR];  /* all data about mols/atoms in the cavity */
  enum ndxtype type;         /* type of index (interpreted as atom or mol) */

  char       s_tmp[STRLEN];  /* string for use with sprintf() */
  char       s_title[STRLEN];/* and another one... */  

  static char *bugs[] = { 
    "-m behaves differently from the standard usage "
    "within the g_* programs -- it figures out _for itself_ what the "
    "molecules are and does not need MOLECULE numbers but ATOM_IDs.",
    "-m is the DEFAULT behaviour!",
    "When counting ions you MUST use -nom !",
    "The density is calculated as (total mass in cavity)/((R-r_water)^2*pi*(z2-z1)) -- "
    "so it makes only sense if this approximates the true cavity volume.",
    "The DEFAULT volume/radius correction is specific for methane pseudo "
    "atoms (ffgmx, ffG43a1) and SPC water.",
    "The program is only tested with pore axis parallel to z-axis (0,0,1).",
  };

#define NFILE asize(fnm)
#define NPA   asize(pa)

  CopyRight(stderr,argv[0]);
  fprintf (stderr, "\nVersion: %s\n\n", SRCID_g_count_c);

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
  a = mols->a;
  atndx = mols->index;

  /* construct array molndx of mol indices in atndx if -m is set */
  molndx = NULL;
  gnmol = -1;
  if (bMolecular) {
    msg ("Interpreting indexfile entries as parts of molecules and "
         "using \ntheir center of mass.\n");
    snew (molndx, mols->nr);
    if ( (gnmol = mols_from_index (index, gnx, mols, molndx, mols->nr)) < 0) {
      gmx_fatal(FARGS, "Error: could not find  molecules.\n");
    };
    msg ("%-10s%10s%10s\n", "Group", "Molecules", "Atoms");      
    msg ("%-10s%10d%10d\n", grpname,  gnmol, gnx);
    msg ("%-10s%10d%10d\n", "System", mols->nr, top->atoms.nr);
  };

  natoms = read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
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
	     grpname);
  sprintf(s_title, "%s in cavity", grpname);
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
  snew (trckd[etxATOM], natoms);
  snew (trckd[etxMOL],  nmolecules);
  
  do {
    /* write time in ps */
    sprintf (s_tmp, "%10g", t);
    fprintf(out, s_tmp);
    fprintf(fData, s_tmp);
    fprintf(fConc, s_tmp);
    fprintf(fDens, s_tmp);

    nmol = 0;       /* restart counting at each timestep */
    totalmass = 0;  /* mass of all particles in the cavity */

    if (bMolecular) {
      rm_pbc(&(top->idef),top->atoms.nr,box,x,x_s); /* remove pbc. */
      
      /* Loop over all molecules. Calculate the center of mass for each
	 molecule. To do so, give an index with all atoms in the molecule
	 to calc_xcm. Index number for atom j of molecule i is 
	 a[atndx[index[i]]+j]. Don't you just love Gromacs? See block.h
	 for the definition of t_block and how to use it. 

	 (yeah - but how to generate the stupid molecules file? It
	 should contain the MOLECULE numbers). Either use a script
	 (~/Gromacs/Scripts/molndx.pl) or use my own
	 molecule index ( mols_from_index(), see above)
 
      */
      for (i = 0; i < gnmol; i++) {
	moleculesize = atndx[molndx[i]+1] - atndx[molndx[i]];
	tm=calc_xcm(x_s, moleculesize, &a[atndx[molndx[i]]], 
		  top->atoms.atom, xcm, FALSE);
	/* We used non-pbc coordinates. Now put cm back in the box */
	for (j = 0; j < DIM; j++) {
	  if (xcm[j] < 0) xcm[j]+= box[j][j];
	  if (xcm[j] > box[j][j]) xcm[j]-=box[j][j];
	}
	/* too much...
	if (bDebugMode) print_ldist (xcm, &cavity, index[i], tm);
	*/
	if ( bInCavity (xcm, &cavity) ) {
	  nmol++;
	  totalmass += tm;

	  tndx = track_ndx (etxMOL, molndx[i], trckd[etxMOL], gnmol);
	  (trckd[etxMOL][tndx].time)++; 
	  for (j = 0; j < moleculesize; j++) {
	    tndx = track_ndx (etxATOM, a[atndx[molndx[i]] + j], 
				       trckd[etxATOM], natoms);
	    (trckd[etxATOM][tndx].time)++; 
	  };
	};
      }
      /* End loop over all molecules */
    } else {
      for(i=0; i<gnx; i++) {
	if ( bInCavity (x[index[i]], &cavity) ) {
	  nmol++;
	  totalmass += top->atoms.atom[index[i]].m;
	  tndx = track_ndx (etxATOM, index[i], 
				     trckd[etxATOM], natoms);
	  (trckd[etxATOM][tndx].time)++; 
	};
      }
      /* End loop over all atoms */
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
  } while(read_next_x(status,&t,natoms,x,box));

  /* initialize stime */
  stime.tau        = tau;
  stime.tot_frames = tot_frames;
  stime.t_tot      = t;
  stime.delta_t    = stime.t_tot/(stime.tot_frames-1);

  dmsg ("Simulation time: t_tot ==  (tot_frames-1) * delta_t\n"
	"                 %g ps  == %d * %g ps \n",
	stime.t_tot, stime.tot_frames, stime.delta_t); 
  
  /* [sort and] write tracking index (both if also looking for
     molecules). Only track molecules above tau threshold 
  */
  do_tracking (fTrack, fTDat, trckd, bMolecular, 
	       &(stime), top, grpname);

  /* clean up a bit */
  fprintf(stderr,"\n");
  close_trj(status);
  fclose(out);
  fclose(fData);
  fclose(fTrack);
  fclose(fConc);
  fclose(fDens);

  thanx(stdout);
  
  return 0;
}

