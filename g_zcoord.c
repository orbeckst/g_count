/*! \file g_zcoord.c
 * \brief  outputs the z coordinates of atoms
 * \author Oliver Beckstein <oliver@bioch.ox.ac.uk>
 * \date   2002
 * $Id: g_zcoord.c,v 1.20 2009/06/24 18:47:02 oliver Exp $
 */
/*! \mainpage g_zcoord development documentation
  *
  * \section intro Introduction
  *
  * This is the introduction.
  *
  * \section install Installation
  *
  * \subsection step1 Step 1: Opening the box
  * 
  * etc...
  */

static char *SRCID_g_zoord_c = "$Id: g_zcoord.c,v 1.20 2009/06/24 18:47:02 oliver Exp $";

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
#include "statutil.h"
#include "rdgroup.h"
#include "mshift.h"
#include "xvgr.h"
#include "gstat.h"
#include "princ.h"
#include "names.h"
#include "matio.h"
#include "xf.h"
#include "count.h"

/* stringification magic 
   use xstr(MACRO) to insert MACRO as a string in the sourcecode
 */
#define xstr(s) str(s)
#define str(s) #s


#define XPM_INITIAL_FRAMES 1000 /*! do dynamic alloc every N frames  */
#define Z_SANITY        10.0    /*! dont believe |z| > this value [nm] */
#define MINZVAL_DEFAULT 0.0
#define MAXZVAL_DEFAULT 6.0
#define XPM_Z_SPACING   0.01   /*! z spacing (bin width) for xpm graph [nm] */

//! explanation for one function - MB doxygen example
/*! 
 * and here a more detailed doxygen description of the
 * extremely interesting dump_data function :))
 */
void dump_data (FILE *f[],  const char *format, ...);


void dump_data (FILE *f[],  const char *format, ...) {
  int i;

  va_list args;
  va_start(args, format);
  
  /* next entry after last one should be NULL. This should even work
     if all entries in f[] are NULL (ie no data output selected)
   */
  while (*f) {
    vfprintf (*(f++), format, args);
  };
  
  va_end(args);
  return;
};


int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_zcoord can output the z(t) coordinates of atoms which are located in a "
    "cylindrical region. "
    "This region is defined by its radius R, the axis, and a point on this axis."
    "It takes an index file with atom, ",
    "and generates an output file with the ",
    "z coordinates for each entry in the index file. Coordinates of particles "
    "that are NOT in the cylinder are set to a special value, which is set "
    "with the [TT]-nan[TT] option.",

    "[PAR]For quicker visualisation the [TT]-mat[tt] option is useful. It creates ",
    "an xpm matrix file. The boundaries of the graph in z are taken from "
    "z1 and z2; the spacing can be set with the hidden option "
    "[TT]-matspacing[TT]. Its default is " xstr(XPM_Z_SPACING) " nm.",

    "[PAR]Due to the huge amount of diskspace that the .dat and .xvg file can "
    "take up these options have been all made optional. YOU NEED TO SET AT LEAST "
    "ONE OF [TT]-dat[TT], [TT]-o[TT], or [TT]-mat[TT] to generate any output."
  };

  t_cavity   cavity = {   /* define volume to count mols in */
    {0, 0, 1},            /* axis */
    {0, 0, 0},            /* cpoint */
    -1,                   /* radius */
    -1,                   /* z1 < z2 */
    -1,
    0                     /* volume - calculate later */
  };

  static real spacing = XPM_Z_SPACING;    /* spacing in nm (!) for conversion to xpm */
  static int  ngroups = 1;      /* not used yet >1 */
  static real invalid = 9999.0; /* mark coordinates of particles not in cylinder */ 
  t_pargs pa[] = {
    { "-axis",   FALSE, etRVEC, {&(cavity.axis)},
      "Vector pointing parallel to the pore axis"},
    { "-cpoint", FALSE, etRVEC, {&(cavity.cpoint)},
      "Point on the pore axis"},
    { "-R",      FALSE, etREAL, {&(cavity.radius)},
      "Radius of the cylindrical pore cavity"},
    { "-z1",     FALSE, etREAL, {&(cavity.z1)},
      "Confine waters to have a z coordinate larger than z1, and"},
    { "-z2",     FALSE, etREAL, {&(cavity.z2)},
      "smaller than z2 (in nm). Also the z-boundaries for the xpm matrix."},
    { "-nan",    FALSE, etREAL, {&invalid},
      "number that is used in output files for coordinates when the particle "
      "is not in the cylinder"},
    { "-matspacing", FALSE, etREAL, {&spacing},
      "HIDDENz spacing (in nm) for the xpm matrix"},
    { "-ng",       FALSE, etINT, {&ngroups},
      "HIDDENNumber of groups to consider" },    
  }; 
  FILE       *out = NULL;    /* xmgr file with coordinates */
  FILE       *fData= NULL;   /* simple data file: t z      */
  t_topology top;            /* topology                   */
  rvec       *xtop;
  int        ePBC;
  char       title[STRLEN];
  bool       bTop;
  rvec       *x,*x_s;        /* coord, with and without pbc*/
  rvec       xcm;            /* center of mass of molecule */
  matrix     box;            /* box matrix (3x3)           */
  real       t,tm;           /* time, total mass molecule  */
  int        natoms;         /* number of atoms in system  */
  int        status;
  int        i,j;            /* loopcounters               */
  int        k_frame;        /* counts frames */

  char       *indexfn;
  char       **grpname;      /* groupnames for all groups */  
  int        *isize;
  atom_id    **index;        /* molecule numbers in index for each group  */
  atom_id    *atndx = NULL;  /* indices from block.        */
  t_block    *mols= NULL;    /* all molecules in system    */  

  char       *ggrpname;      /* name of the FIRST group       */
  int        gnx = 0;        /* number of atoms in FIRST group*/
  atom_id    *gindex = NULL; /* index of FIRST group */

  real maxzval = MINZVAL_DEFAULT; /* max z in nm (!) for conversion to xpm */
  real minzval = MAXZVAL_DEFAULT; /* min z in nm (!) for conversion to xpm */
  int        nbins = 0;           /* number of bins for conv to xpm */
  real       **zbins = NULL; /* bins for conversion to xpm */
  FILE       *xpmOut = NULL; /* xpm matrix file */
  t_rgb      rlo, rhi;       /* low and high colors for matrix */
  int nlevels = 2;
  int max_frames = XPM_INITIAL_FRAMES; 
                            /* current max number of frames in xpm matrix */
  real       *x_axis = NULL;/* x-axis labels */
  real       *y_axis = NULL;/* y-axis labels */
  bool       bXPM = FALSE;  /* produce XPM matrix ?*/
  bool       bXVG = FALSE;  /* no output generated by default */
  bool       bDAT = FALSE; 

  t_filenm fnm[] = {
    { efTRX, "-f",   NULL,     ffREAD  },
    { efTPS, NULL,   NULL,     ffREAD },
    { efXVG, "-o",   "zcoord", ffOPTWR },
    { efXPM, "-mat", "zcoord", ffOPTWR }, /* this is for xpm output */
    { efDAT, "-dat", "zcoord", ffOPTWR },
    { efNDX, NULL,   NULL,     ffOPTRD }
  };

#define NFILE asize(fnm)
#define NPA   asize(pa)
#define NDFILES 2
  FILE       *dfiles[NDFILES+1]; /* output files to dump data to 
				    terminate with NULL */
 
  static char *bugs[] = { 
    "Even if the cylinder axis is not (0,0,1), the top and bottom of the cylinder are"
    "still determined by the z-coordinates z1 and z2.",
    "ATTENTION: no data is generated by default !!!",
    "xpm graph: If no proper (i.e. |z| < " xstr(Z_SANITY) "nm) z-limits are "
    "given they are set to "
    "the arbitrary values z_min=" xstr(MINZVAL_DEFAULT) " and "
    "z_max=" xstr(MAXZVAL_DEFAULT),
  };


  CopyRight(stderr,argv[0]);
  fprintf (stderr, "\nVersion: %s\n\n", SRCID_g_zoord_c);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);
  
  bDAT = (opt2fn_null("-dat", NFILE, fnm)) ? TRUE : FALSE;
  bXVG = (opt2fn_null("-o",   NFILE, fnm)) ? TRUE : FALSE;
  bXPM = (opt2fn_null("-mat", NFILE, fnm)) ? TRUE : FALSE;

  if (! (bXPM  || bXVG || bDAT)) {
    /* check if this is a pointless run */
    msg("No output selected. I stop here.\n");
    return (1);
  }

  bTop = read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
  sfree(xtop);

  if (!bTop) {
    gmx_fatal(FARGS, "Need a run input file");
  }

  indexfn = ftp2fn_null(efNDX,NFILE,fnm);

  if (ngroups != 1) {
    gmx_fatal(FARGS, "Sorry, only a single group currently allowed.");
  }

  snew(grpname,ngroups);
  snew(isize,ngroups);
  snew(index,ngroups);
  get_index(&(top.atoms),indexfn,ngroups,isize,index,grpname);
  
  /* ngroups == 1 at moment */
  gnx = isize[0];
  gindex = index[0];
  ggrpname = grpname[0];
  mols = &(top.mols);
  atndx = mols->index;
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  snew(x_s,natoms);


  /* look at cavity */
  if (cavity.z1 >= cavity.z2) {
    gmx_fatal(FARGS, "z-boundaries have worng values: They must be z1 < z2.");
  }
  
  /* look at cavity 
     set defaults from the simulation box size
  */
  msg("\n");
  autoset_cavity(&cavity,box,NPA,pa);
  cavity.vol = sqr(cavity.radius) * PI * (cavity.z2 - cavity.z1); // not needed, just for completeness

  if (bXPM) {
    /* for xpm output: set limits to z1 and z2 */
    /* mild sanity check */
    minzval = cavity.z1 < -Z_SANITY  ? MINZVAL_DEFAULT : cavity.z1;
    maxzval = cavity.z2 >  Z_SANITY  ? MAXZVAL_DEFAULT : cavity.z2;

    /* xpm graphics: initialize bins  */
    nbins = (int)ceill((maxzval - minzval) / spacing);
    snew(zbins,max_frames);
    snew(x_axis,max_frames);

    msg("Limits for xpm matrix: %g nm < z < %g nm\n", minzval, maxzval);
  }

  /* initialize dfiles: use it to dump data to different files simultanously */
  for (i = 0; i < NDFILES+1; i++) dfiles[i] = NULL;
  i = 0;       /* use to count how many output files are open */
  if (bXVG) {
    out    = xmgropen(opt2fn("-o",NFILE,fnm),
		      "z coordinate",
		      "Coordinate",
		      "Time (ps)","z (nm)");
    dfiles[i++] = out;
  }
  if (bDAT) {
    fData  = ffopen (opt2fn("-dat", NFILE, fnm), "w");
    dfiles[i++] = fData;
  }


  k_frame = 0;     /*! counts the frames */
  do {
    dump_data (dfiles,  "%10g", t);

    if (bXPM) {
      // new timeframe
      // check if we need more memory
      if (k_frame >= max_frames) {
        max_frames += XPM_INITIAL_FRAMES;
        if (!(srenew(zbins,max_frames)) ||
            !(srenew(x_axis,max_frames)) ) {
          msg("Out of memory for xpm graph: I abort and try to process what "
              "I have got so far (frame = %g).\n", k_frame);
          break;
        }
      }
      snew(zbins[k_frame],nbins);
      x_axis[k_frame]=t/1000;  // in ns

      // initialize bins to zero
      for (i=0; i<nbins; i++) 
        zbins[k_frame][i] = 0.0;
    }

    /* We are not interested in molecules, just in the z-coordinates 
       of atoms */
    for(i=0; i<gnx; i++) {
      if ( bInCavity (x[gindex[i]], &cavity) ) {  // only look at cylinder
	if (bXPM) {
	  if ( (x[gindex[i]][ZZ] >= minzval) && (x[gindex[i]][ZZ] <= maxzval))
	    zbins[k_frame][(int)floorl((x[gindex[i]][ZZ] - minzval) / spacing)]=1.0;
	}
      } 
      else {
	/* particle OUTSIDE cylinder: mark coordinates as 'invalid'
	   but keep them in the output file */
	for(j=0; j<3; j++) {
	  x[gindex[i]][j] = invalid;
	}
      }
      dump_data (dfiles,"\t%10g", x[gindex[i]][ZZ]);
    }

    /* some printout for debugging/development - MB 
    for (i=0; i<nbins; i++) 
    for (i=0; i<10; i++) 
      printf("%1d",(int)zbins[k_frame][i]);
    printf("\n");
    */

    dump_data (dfiles, "\n");
    k_frame++;
  } while(read_next_x(status,&t,natoms,x,box));

  if (bXPM) {
    /* colors for matrix */
    rlo.r = 1.0, rlo.g = 1.0, rlo.b = 1.0; // no water color
    rhi.r = 0.4, rhi.g = 0.4, rhi.b = 1.0; // water color
  
    /* create labels */
    snew(y_axis,nbins+1);
    for(i=0; i<=nbins; i++) {                   // n+1 labels with MAT_SPATIAL_Y
      y_axis[i]=10 * (minzval + (i * spacing)); // Yes, I WANT Angstroms !!!
    }
    /* write xpm matrix */
    xpmOut = fopen(opt2fn("-mat", NFILE, fnm),"w");
    write_xpm(xpmOut,
	      MAT_SPATIAL_Y,
	      "g_zcoord matrix",     // title
	      "existence",           // label for the continuous legend
	      "timeframe / ns",      // label for the x-axis
	      "z coordinate / A",    // label for the y-axis - where the F**k is the
                          	      // bloody Angstrom
	      k_frame, nbins,        // size of the matrix
	      x_axis,            // the x-ticklabels
	      // would be nice to have the times here
	      y_axis,            // the y-ticklabels
	      // would be nice to have the zcoordinate label in nm or A
	      zbins,             // element x,y is matrix[x][y]
	      0.0,               // output lower than lo is set to lo
	      1.0,               // output higher than hi is set to hi
	      rlo,               // rgb value for level lo
	      rhi,               // rgb value for level hi
	      &nlevels);         // number of color levels for the output
    fclose(xpmOut);
  };

  /* clean up a bit */
  fprintf(stderr,"\n");
  close_trj(status);

  if (bDAT) fclose(fData);
  if (bXVG) fclose(out);
  
  if (bXPM) {
    // free memory from xpm matrix data
    for (i=0; i<k_frame; i++)
      sfree(zbins[i]);
    sfree(zbins);
    sfree(x_axis);
  }


  thanx(stdout);
  
  return 0;
}

