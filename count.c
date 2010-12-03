/* $Id: count.c,v 1.7 2009/06/15 15:49:57 oliver Exp $
   
   specialised functions which are used in in g_count and g_flux

*/

#include "count.h"

char *endxtype_names[etxNR+1] = {
  "atoms", "molecules", NULL };
char *eDensUnit_names[eduNR+1] = {
  "nm^-3", "n(SPC,T=300K,P=1bar)", "mol/l", "Angstrom^-3", NULL };
/* conversion divisors: n' = n/DensUnit, see comments in count.h  */
real DensUnit[eduNR] = {
  1.000000,  32.3227, 0.602214, 1000.0};  


bool autoset_cavity (t_cavity *cavity,matrix box,int npa,t_pargs pa[]) {
  bool bAutoSet = FALSE;

  if(!opt2parg_bSet("-cpoint",npa,pa)) {
    new_rvec(cavity->cpoint,0.5*box[XX][XX],0.5*box[YY][YY],0.5*box[ZZ][ZZ]);
    bAutoSet = TRUE;
    msg("auto-set: cpoint = (%.3f,%.3f,%.3f)\n",
	cavity->cpoint[XX],cavity->cpoint[YY],cavity->cpoint[ZZ]);
  }
  
  if(!opt2parg_bSet("-R",npa,pa)) {
    cavity->radius=MIN(
		MIN(cavity->cpoint[XX],fabs(cavity->cpoint[XX]-box[XX][XX])),
		MIN(cavity->cpoint[YY],fabs(cavity->cpoint[YY]-box[YY][YY])));
    bAutoSet = TRUE;
    msg("auto-set: radius R=%.4f nm\n",cavity->radius);
  }

  if(!opt2parg_bSet("-z1",npa,pa)) {
    cavity->z1 = 0;
    bAutoSet = TRUE;
    msg("auto-set: z1 = %.3f nm\n",cavity->z1); 
  }
  if(!opt2parg_bSet("-z2",npa,pa)) {
    cavity->z2 = box[ZZ][ZZ];
    bAutoSet = TRUE;
    msg("auto-set: z2 = %.3f nm\n",cavity->z2); 
  }

  if (cavity->z1 > cavity->z2) {
    real tmp = cavity->z1;
    cavity->z1 = cavity->z2;
    cavity->z2 = tmp;
    bAutoSet = TRUE;
    msg ("Warning: your input z1 > z2 made no sense. I swapped them.\n");
  };
  cavity->vol = PI * cavity->radius * cavity->radius * (cavity->z2 - cavity->z1);

  return bAutoSet;
}


void x2molxcm (t_topology *top, int ePBC, atom_id *molndx, int gnmol, matrix box, 
	       rvec *x, rvec *x_s, rvec *xmol_cm) {
  /* put molecular center of mass coordinates into xmol_cm (must have
     at least gnmol entries allocated!) */
  t_block   *mols;             /* all molecules in system    */
  atom_id   *a, *atndx;        /* indices from block.        */
  rvec      xcm;               /* center of mass */
  int       moleculesize;           
  int       i,j;

  /* remove pbc. */
  rm_pbc(&(top->idef),ePBC,top->atoms.nr,box,x,x_s); 

  /* Loop over all molecules. Calculate the center of mass for each
     molecule. To do so, give an index with all atoms in the molecule
     to calc_xcm. Index number for atom j of molecule i is 
     a[atndx[index[i]]+j]. Don't you just love Gromacs? See block.h
     for the definition of t_block and how to use it. 
     
     (yeah - but how to generate the stupid molecules file? It
     should contain the MOLECULE numbers). Either use a script
     (~/Gromacs/Scripts/molndx.pl) or use my own
     molecule index ( mols_from_index() )
  */

  mols=&(top->mols);
  atndx = mols->index;
  // a = mols->a;

  for (i = 0; i < gnmol; i++) {
    moleculesize = atndx[molndx[i]+1] - atndx[molndx[i]];
    // calc_xcm(x_s, moleculesize, &a[atndx[molndx[i]]], 
    calc_xcm(x_s, moleculesize, &(atndx[molndx[i]]),       // XXX BROKEN
	     top->atoms.atom, xcm, FALSE);
  
    /* We used non-pbc coordinates. Now put cm back in the box */
    for (j = 0; j < DIM; j++) {
      if (xcm[j] < 0) xcm[j]+= box[j][j];
      if (xcm[j] > box[j][j]) xcm[j]-=box[j][j];
      xmol_cm[i][j] = xcm[j];
    };
  };
  return;
};


real npbc2com (t_topology *top, atom_id *molndx, matrix box, rvec *x_s, 
	       atom_id i, rvec xcm) { 
  /* calculate the center of mass of molecule i from its non-pbc
     coordinates 
     returns: total mass of the molecule
  */

  t_block   *mols;             /* all molecules in system    */
  atom_id   *a, *atndx;        /* indices from block.        */
  real      tm;                /* molecule mass */
  int       moleculesize;           
  int       j;


  mols=&(top->mols);
  // a = mols->a;
  atndx = mols->index;

  moleculesize = atndx[molndx[i]+1] - atndx[molndx[i]];
  tm=calc_xcm(x_s, moleculesize, &(atndx[molndx[i]]),           // XXX BROKEN
	      top->atoms.atom, xcm, FALSE);
  
  /* We used non-pbc coordinates. Now put cm back in the box */
  for (j = 0; j < DIM; j++) {
    if (xcm[j] < 0) xcm[j]+= box[j][j];
    if (xcm[j] > box[j][j]) xcm[j]-=box[j][j];
  }
  return tm;
};



bool bInCavity (const rvec x, const t_cavity *cavity) {

  if (x[ZZ] > cavity->z1 && x[ZZ] < cavity->z2) {
    return (ldist (x, cavity->axis, cavity->cpoint) < cavity->radius);
  }
  return FALSE;
};

int mols_from_index (atom_id *index, int gnx, t_block *mols, 
		     atom_id *molndx, int max_mol) {
  /* make a list of all molecule indices which contain atoms from the index */
  
  int i, k, l;
  int nr_mol;    /* number of molecules in molndx */
  
  nr_mol = 0;

  /* XXX: COMPLETELY BROKEN */

/*   for (i = 0; i < gnx; i++) { */
/*     /\* check every single atom in the index file against all entries */
/*        in the molecules block (BUG: only find corresponding molecule */
/*        if we hit the first entry in a[], ie the one that the */
/*        mol->index points to. Should really sweep thru the atoms for each */
/*        molecule(size) as well)  */
/*     *\/ */
/*     k = 0; */
/*     while ( mols->a[mols->index[k]] != index[i] && k < mols->nr ) */
/*       { k++; }; */
/*     if (k < mols->nr) { */
/*       /\* found one: add found molindex k to molndx if not already there *\/ */
/*       l = nr_mol; */
/*       while ( molndx[--l] != k && l >= 0); */
/*       if (l < 0) { */
/* 	assert (nr_mol < max_mol); */
/* 	molndx[nr_mol++] = k; */
/* #ifdef DEBUG */
/* 	dfprintf ("mols_from_index(): molndx[%5d] = %5d\n", nr_mol-1, k); */
/* #endif */
/*       }; */
/*     }; */
      
/*   }; */

  return nr_mol;
}

void print_ldist (const rvec x, const t_cavity *cavity, const atom_id idx, const real mass) {
  char *s;
  bool incav;

  incav = bInCavity (x, cavity);

  s = incav ? "IN CAVITY" : "";
  fprintf (stderr, "atom[%5d]=[%6f, %6f, %6f]:"
	           "dist = %4f m=%4f u, bCav [%d] %s\n",
	   idx, x[XX], x[YY], x[ZZ],
	   ldist (x, cavity->axis, cavity->cpoint), mass, incav, s);
};

#define LNDXMX 15    /* entrie per line in an index file */

void fwrite_index (FILE *fp, atom_id *list, enum ndxtype list_type,
		 int nmx, t_topology *top, char *grpname){
  /* write an index group to file fp. Interprete entries of list as
     MOLECULES (type == etxMOL or atoms (etxATOM) 

     (1) list contains MOLECULES (bMolecular == 1)
         * write residue numbers to index group
	 * lookup atoms for each molecule
	   -> write atom numbers 
     (2) list contains ATOMS
         only write atoms to the index group

     Routine derived from do_tracking() (what a mess...).
 
*/

  int i,j;
  bool bDebugOutput;

  atom_id *atndx;        /* indices from block.        */
  t_block *mols;             /* all molecules in system    */
  t_atom  *atom;             /* all atoms */		 
  int     moleculesize;
  enum ndxtype type;

  if (!fp) {
    msg ("fwrite_index() in file %s, line %g: attempt to write without "
	  "opening a file.\n", __FILE__, __LINE__);
    return;
  };
  assert (list_type == etxMOL || list_type == etxATOM);

  msg ("Writing index of all %d %s that crossed the pore.\n",
       nmx, ENDXTYPE(list_type));


  /* needed for translation mols -> atom ids */
  mols  = &(top->mols);
  atndx = mols->index;
  /* translation atom_id -> resnr 
     top->atoms->atom[atom_id]->resnr
  */
  atom = top->atoms.atom;

  for (type = etxATOM; type <= list_type; type++) {
    /* ids to new index file */
    fprintf (fp, "[ %s%s ]\n", grpname,
	     type == etxMOL ? "_residues" : "");  

    /* if we have molecules then we use the molecule numbers to get
       the atoms */
    if (list_type == etxMOL && type == etxATOM) {
      int counter;   /* number of entries */
      for (i = counter = 0; i < nmx; i++) {
	moleculesize = atndx[list[i] + 1] - atndx[list[i]];
	for (j = 0; j < moleculesize; j++, counter++) {
	  /* ... fetch all atoms of this molecule */
	  /* ADD +1 when WRITING (external) index file !!! */
	  fprintf (fp, "%5u %s", atom[atndx[list[i]]+j].atomnumber + 1,   // XXX possibly BROKEN
		   NEWLINE(LNDXMX,counter));
	};
      };
    } else {
	/* either etxMOL or etxATOM without molecules */
	for (i = 0; i < nmx; i++) {
	  atom_id igmx;
	  igmx = (type == etxMOL ? atom[atndx[list[i]]].atomnumber : list[i]); // XXX poss. BROKEN for etxMOL
	  /* ADD +1 when WRITING (external) index file !!! */
	  fprintf (fp, "%5u %s", atom[igmx].resnr + 1, NEWLINE(LNDXMX,i));
	};
    };
    fprintf (fp, "\n\n");
  }; /* end for loop over type */

  return;    
};



void dump_atomlist (FILE *fp, enum ndxtype type, atom_id *list, 
		 int nmx, t_topology *top){
  /* write an index group to file fp. Interprete entries of list as
     MOLECULES (type == etxMOL or atoms (etxATOMS) */

  int i, j;          /* loop */
  int moleculesize; 
  bool bDebugOutput;

  atom_id *a, *atndx;        /* indices from block.        */
  t_block *mols;             /* all molecules in system    */
  t_atom  *atoms;            /* all atoms */		 

  /* if fp == NULL  given AND debug file opened print debugging table */
  if ((bDebugOutput = !fp)) {
    if (!debug) fp = stderr;
  };

  /* needed for translation mols -> atom ids */
  mols  = &(top->mols);
  atndx = mols->index;
  /* translation atom_id -> resnr 
     top->atoms->atom[atom_id]->resnr
  */
  atoms = top->atoms.atom;

  for (i = 0; i < nmx; i++) {
    atom_id igmx;

    igmx = (type == etxMOL ? atndx[list[i]] : list[i]);  // XXX BROKEN for etxMOL
    
    if (!bDebugOutput) {
    fprintf (fp, "%6u %6u %6s\n",
	     igmx + 1,
	     atoms[igmx].resnr + 1,
	     *(top->atoms.atomname[igmx]));
    } else {
    fprintf (debug, "  %8u %6u %6u %6s\n",
	     list[i] + 1,
	     igmx + 1,
	     atoms[igmx].resnr + 1,
	     *(top->atoms.atomname[igmx]));
    };
  };

  return;
};
