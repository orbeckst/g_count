/* $Id: utilgmx.c,v 1.5 2008/01/14 15:14:05 oliver Exp $
   everyday functions
   $Log: utilgmx.c,v $
   Revision 1.5  2008/01/14 15:14:05  oliver
   cleaned up and compiles with Gromacs 3.3.1; breaks backward compatibility with Gromacs 3.2.1

   Revision 1.4  2002/08/13 17:36:05  oliver
   * new program: g_ri3Dc
     . count molecules in grid cells
     . effective radius pore profile
     . volume calculation
     . xy and xz density projections (xfarbe format)
     No 3D grid output yet (lacking a suitable format)
   * changed DR_DEFAULT to new value
   * smaller changes of comments

   Revision 1.2.2.2  2002/08/11 19:35:23  oliver
   grid allocation functions for real variables in  2D and 3D
   (from Numerical Recipes & simplified + use snew())

   Revision 1.2.2.1  2002/05/29 18:43:23  oliver
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

   Revision 1.2  2001/02/19 14:26:22  oliver
   moved dt_tpx to utilgmx

   Revision 1.1  2001/02/18 16:38:13  oliver
   Initial revision

*/

static char *SRCID_utilgmx_c = "$Id: utilgmx.c,v 1.5 2008/01/14 15:14:05 oliver Exp $";

#include "utilgmx.h"
#include "fatal.h"
#include "futil.h"
#include "maths.h"
#include "smalloc.h"

real ***grid3_alloc(int nx, int ny, int nz) {
  /* allocate a 3D array, size nx*ny*nz which is referenced through pointers.
     By using snew() the memory is already 0-initialised.
     range is g[0..nx-1,0..ny-1,0..nz-1].
     x: rows
     y: columns
     z: depth
     This is a trimmed down version of "Numerical Recipes in C", p944 , f3tensor() 

     This memory/pointer stuff is mind boggling: Draw it on
     paper. Basically, g[] is a vector, pointing to the rows of an
     array g[][], nx*ny (allocated as one block). These entries again
     point to rows in the 'real' chunk of memory nx*ny*nz, g[][][].

     Because we have allocate the memory in blocks we can initialise
     all pointers in a fairly tricky (to me, anyway) fashion.

     The original routine allows for another memory range of the
     indices as well.  

  */
  int  i,j;
  real ***g = NULL;

  /* allocate a row of 'ptr to ptr to rows' (= index by x values), length nx */
  snew(g,nx);
  if (!g) gmx_fatal(FARGS,"grid3_alloc(): Cannot allocate memory for nx=%d\n",nx);
  
  /* allocate nx 'ptrs to rows' (length ny each) and set ptrs to them */
  snew(g[0],nx*ny);
  if (!g) gmx_fatal(FARGS,"grid3_alloc(): Cannot allocate memory for nx*ny=%d*%d\n",
		      nx,ny);
    
  /* allocate nx*ny 'rows' (length nz each; these memory locations hold the reals) */
  snew(g[0][0],nx*ny*nz);
  if (!g) gmx_fatal(FARGS,"grid3_alloc(): Cannot allocate memory for data "
		      "(nx*ny*nz=%d*%d%d)\n",nx,ny,nz);

  /* initialise all pointers (starting from the already initialised
     initial locations), "leaping" through the array 
  */
  for(j=1; j<ny; j++) g[0][j]=g[0][j-1]+nz;
  for(i=1; i<nx; i++) {
    g[i]=g[i-1]+ny;
    g[i][0]=g[i-1][0]+ny*nz;
    for(j=1; j<ny; j++) g[i][j]=g[i][j-1]+nz;
  }

  /* return ptr to an array of ptr to rows */
  return g;
}

real **grid2_alloc(int nx, int ny) {
  /* allocate a 3D array, size nx*ny*nz which is referenced through pointers.
     By using snew() the memory is already 0-initialised.
     range is g[0..nx-1,0..ny-1].
     x: rows
     y: columns
     This is a trimmed down version of "Numerical Recipes in C", p944 , f3tensor() 
  */
  int  i,j;
  real **g = NULL;

  /* allocate a row of 'ptr to ptr to rows' (= index by x values), length nx */
  snew(g,nx);
  if (!g) gmx_fatal(FARGS,"grid2_alloc(): Cannot allocate memory for nx=%d\n",nx);
  
  /* allocate nx 'ptrs to rows' (length ny each) and set ptrs to them */
  snew(g[0],nx*ny);
  if (!g) gmx_fatal(FARGS,"grid2_alloc(): Cannot allocate memory for nx*ny=%d*%d\n",
		      nx,ny);

  for(i=1;i<nx;i++) g[i]=g[i-1]+ny;
  
  /* return ptr to an array of ptr to rows */
  return g;
}


void quicksort(atom_id v[], int left, int right)
{
  int i,last;

  if (left >= right)                    /* Do nothing if array contains */
    return;                             /* fewer than two elements      */
  aswap(v,left,(left+right)/2);         /* Move partition element       */
  last=left;                            /* to v[0]                      */
  for(i=left+1; (i<=right); i++)        /* partition                    */
    if (v[i] < v[left]) {
      last++;
      aswap(v,last,i);                  /* watch out for macro trick    */
    }
  aswap(v,left,last);                   /* restore partition element    */
  quicksort(v,left,last-1);
  quicksort(v,last+1,right);
}


void msg (const char *format, ...)
{
  va_list args;
  va_start(args, format);

  vfprintf (stderr, format, args);
  if (bDebugMode()) {
    vfprintf(debug, format, args);
  };

  va_end(args);
  return;
};

void dmsg (const char *format, ...)
{
  va_list args;
  va_start(args, format);

  if (bDebugMode()) {
    vfprintf(stderr, format, args);
    vfprintf(debug,  format, args);
  };

  va_end(args);
  return;
};
  

void dprintf (const char *format, ...)
{
  va_list args;
  va_start(args, format);
  
  if (bDebugMode()) {
    vfprintf(stderr, format, args);
  };
  va_end(args);
  return;
};

void dfprintf (const char *format, ...)
{
  /* write to debug file (see fatal.h) */
  va_list args;
  va_start(args, format);
  
  if (bDebugMode()) {
    vfprintf(debug, format, args);
  };
  va_end(args);
  return;
};


real ldist (rvec x, rvec p, rvec c) {
  /* calculate the (perpendicular) distance of point x from the line
     through c pointing along p 
  */
  real lambda;
  rvec d, u;

  /* u = x - c */
  rvec_sub (x, c, u);
  
  /* d = [(u.p)/p^2]p  - u */
  lambda = iprod (u, p) / iprod (p, p);
  svmul (lambda, p, d);
  rvec_dec (d, u);

  return norm (d);
};

real dt_tpx (char *fn) {
  int step, natoms;
  real t, lambda;
  t_inputrec  ir;  

  /* discard almost all info */
  read_tpx(fn,&step,&t,&lambda, &ir, NULL, &natoms, NULL, NULL, NULL, NULL);

  return ir.delta_t * ir.nstxtcout;
};

int list_add_atomid (const atom_id id, int *nid, atom_id *list) {
  /* append id to list and increase nid if not already in there 
     return the new number of ids, nid, if any changes occurred, 0 otherwise
     
     ATTENTION: nid is modified in the caller!!
   */
  int j;
  for (j=0; j < *nid  &&  id != list[j]; j++);
  if (j >= *nid) {
    list[(*nid)++] = id;
    return *nid;
  };
  return 0;
};
