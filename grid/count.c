/* $Id: count.c,v 1.5 2002/08/20 01:13:50 oliver Exp $
   
   specialised functions which are used in in g_count and g_flux

   $Log: count.c,v $
   Revision 1.5  2002/08/20 01:13:50  oliver
   NEW:
   * a_ri3Dc can resize (well, shrink) the grid (only z1, z2, R)
     - needed setup_tgrid() and setup_corners() from g_ri3Dc for this
       dirty hack, so they were moved into grid3D.c
   FIXED:
     calculation of the radius from the grid corners was wrong because I
     used abs() where I should have used fabs(). Didn't have any impact
     so far but for the readjustment of the grid it matters!

   Revision 1.4  2002/08/19 17:16:36  oliver
   - define units to measure a density in; currently available: unity
     (nm^-3 -- canonical gromacs units), SPC (bulk SPC water at 300K,
     1bar), molar (in mol/l), Angstroem (in A^-3)
   - moved autoset_cavity() from g_ri3Dc into count.c (.h) so that I can
     use it for g_count as well

   Revision 1.3  2002/08/13 17:36:05  oliver
   * new program: g_ri3Dc
     . count molecules in grid cells
     . effective radius pore profile
     . volume calculation
     . xy and xz density projections (xfarbe format)
     No 3D grid output yet (lacking a suitable format)
   * changed DR_DEFAULT to new value
   * smaller changes of comments

   Revision 1.1.2.5  2002/08/11 14:45:16  oliver
   rename write_index() to fwrite_index() to avoid collision with a
   gromacs function in index.c

   Revision 1.1.2.4  2002/06/04 21:48:01  oliver
   current working version of g_flux

   Revision 1.1.2.3  2002/05/29 18:43:23  oliver
   * g_flux writes index of all mols/atoms that crossed the pore
   * new function list_add_atomid (const atom_id, int *, atom_id *) in utilgmx.c:
     add an atom_id to list if its not allready in there. Sort list later with
     quicksort(list, 0, last)
   * write_index() in count.c: clean up of do_tracking() from g_count.c;
     writes the contents of list into an index (either as residue numbers
     ('molecules') or atom ids
   In g_flux:
   * store net flux in t_flux
   * new: update_result(); updates t_result from internal values
   * added the list of particles that crossed to t_result
   * inline a few functions
   * HIDE -dR option
   * prune unused variables
   * only allocate for gnx atoms (not natoms)

   Revision 1.1.2.2  2002/05/27 23:52:04  oliver
   first working version of g_flux (working means: compiles, does not core dump, results look vaguely sensible)

   Revision 1.1.2.1  2002/05/27 12:21:25  oliver
   Moved common functions and definition for g_count, g_sldist, g_flux into count.[ch]. g_count compiles

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


void x2molxcm (t_topology *top, atom_id *molndx, int gnmol, matrix box, 
	      rvec *x, rvec *x_s, rvec *xmol_cm) {
  /* put molecular center of mass coordinates into xmol_cm (must have
     at least gnmol entries allocated!) */
  t_block   *mols;             /* all molecules in system    */
  atom_id   *a, *atndx;        /* indices from block.        */
  rvec      xcm;               /* center of mass */
  int       moleculesize;           
  int       i,j;

  /* remove pbc. */
  rm_pbc(&(top->idef),top->atoms.nr,box,x,x_s); 

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

  mols=&(top->blocks[ebMOLS]);
  a = mols->a;
  atndx = mols->index;

  for (i = 0; i < gnmol; i++) {
    moleculesize = atndx[molndx[i]+1] - atndx[molndx[i]];
    calc_xcm(x_s, moleculesize, &a[atndx[molndx[i]]], 
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


  mols=&(top->blocks[ebMOLS]);
  a = mols->a;
  atndx = mols->index;

  moleculesize = atndx[molndx[i]+1] - atndx[molndx[i]];
  tm=calc_xcm(x_s, moleculesize, &a[atndx[molndx[i]]], 
	      top->atoms.atom, xcm, FALSE);
  
  /* We used non-pbc coordinates. Now put cm back in the box */
  for (j = 0; j < DIM; j++) {
    if (xcm[j] < 0) xcm[j]+= box[j][j];
    if (xcm[j] > box[j][j]) xcm[j]-=box[j][j];
  }
  return tm;
};



bool bInCavity (rvec x, t_cavity *cavity) {

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

  for (i = 0; i < gnx; i++) {
    /* check every single atom in the index file against all entries
       in the molecules block (BUG: only find corresponding molecule
       if we hit the first entry in a[], ie the one that the
       mol->index points to. Should really sweep thru the atoms for each
       molecule(size) as well) 
    */
    k = 0;
    while ( mols->a[mols->index[k]] != index[i] && k < mols->nr )
      { k++; };
    if (k < mols->nr) {
      /* found one: add found molindex k to molndx if not already there */
      l = nr_mol;
      while ( molndx[--l] != k && l >= 0);
      if (l < 0) {
	assert (nr_mol < max_mol);
	molndx[nr_mol++] = k;
#ifdef DEBUG
	dfprintf ("mols_from_index(): molndx[%5d] = %5d\n", nr_mol-1, k);
#endif
      };
    };
      
  };

  return nr_mol;
}

void print_ldist (rvec x, t_cavity *cavity, atom_id idx, real mass) {
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

  atom_id *a, *atndx;        /* indices from block.        */
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
  mols  = &(top->blocks[ebMOLS]);
  a     = mols->a;
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
	  fprintf (fp, "%5u %s", a[atndx[list[i]]+j] + 1, 
		   NEWLINE(LNDXMX,counter));
	};
      };
    } else {
	/* either etxMOL or etxATOM without molecules */
	for (i = 0; i < nmx; i++) {
	  atom_id igmx;
	  igmx = (type == etxMOL ? a[atndx[list[i]]] : list[i]);
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
  mols  = &(top->blocks[ebMOLS]);
  a     = mols->a;
  atndx = mols->index;
  /* translation atom_id -> resnr 
     top->atoms->atom[atom_id]->resnr
  */
  atoms = top->atoms.atom;

  for (i = 0; i < nmx; i++) {
    atom_id igmx;

    igmx = (type == etxMOL ? a[atndx[list[i]]] : list[i]); 
    
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
