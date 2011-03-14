/* $Id: count.c,v 1.8 2009/07/02 16:05:48 oliver Exp $
   
   specialised functions which are used in the gridcounters and also 
   in g_count and g_flux

   Copyright (C) 2003-2007 Oliver Beckstein <orbeckst@jhmi.edu>
   This program is made available under the terms of the GNU Public License. 
   See the file LICENSE or http://www.gnu.org/copyleft/gpl.html

*/

#include "count.h"

char *endxtype_names[etxNR+1] = {
  "atoms", "molecules", NULL };
char *eDensUnit_names[eduNR+1] = {
  "nm^-3", "n(SPC,T=300K,P=1bar)", "mol/l", "Angstrom^-3", NULL };
/* conversion divisors: n' = n/DensUnit, see comments in count.h  */
real DensUnit[eduNR] = {
  1.000000,  32.3227, 0.602214, 1000.0};  


gmx_bool autoset_cavity (t_cavity *cavity,matrix box,int npa,t_pargs pa[]) {
  gmx_bool bAutoSet = FALSE;

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



gmx_bool bInCavity (const rvec x, const t_cavity *cavity) {

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
  gmx_bool incav;

  incav = bInCavity (x, cavity);

  s = incav ? "IN CAVITY" : "";
  fprintf (stderr, "atom[%5d]=[%6f, %6f, %6f]:"
	           "dist = %4f m=%4f u, bCav [%d] %s\n",
	   idx, x[XX], x[YY], x[ZZ],
	   ldist (x, cavity->axis, cavity->cpoint), mass, incav, s);
};

#define LNDXMX 15    /* entrie per line in an index file */

void fwrite_index (FILE *fp, atom_id *list, 
		 int nmx, t_topology *top, char *grpname){
  /* write an index group to file fp */
  int i;

  if (!fp) {
    msg ("fwrite_index() in file %s, line %g: attempt to write without "
	  "opening a file.\n", __FILE__, __LINE__);
    return;
  };

  msg ("Writing index of all %d atoms that crossed the pore.\n", nmx);
  for (i = 0; i < nmx; i++) {
    /* ADD +1 when WRITING (external) index file !!! */
    fprintf (fp, "%5u %s", list[i] + 1, NEWLINE(LNDXMX,i));
  };
  fprintf (fp, "\n\n");
};


