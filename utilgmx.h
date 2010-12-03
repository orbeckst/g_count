/* $Id: utilgmx.h,v 1.6 2003/10/09 09:15:27 oliver Exp $
   $Log: utilgmx.h,v $
   Revision 1.6  2003/10/09 09:15:27  oliver
   use exoisting def of MIN/MAX if exists

   Revision 1.5  2002/08/13 17:36:05  oliver
   * new program: g_ri3Dc
     . count molecules in grid cells
     . effective radius pore profile
     . volume calculation
     . xy and xz density projections (xfarbe format)
     No 3D grid output yet (lacking a suitable format)
   * changed DR_DEFAULT to new value
   * smaller changes of comments

   Revision 1.3.2.4  2002/08/11 19:35:23  oliver
   grid allocation functions for real variables in  2D and 3D
   (from Numerical Recipes & simplified + use snew())

   Revision 1.3.2.3  2002/06/04 21:48:01  oliver
   current working version of g_flux

   Revision 1.3.2.2  2002/05/29 18:43:23  oliver
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

   Revision 1.3.2.1  2002/05/27 23:52:04  oliver
   first working version of g_flux (working means: compiles, does not core dump, results look vaguely sensible)

   Revision 1.3  2001/02/19 14:31:52  oliver
   added "tpxio.h" for dt_tpx()

   Revision 1.2  2001/02/19 14:26:40  oliver
   moved dt_tpx to utilgmx

   Revision 1.1  2001/02/18 16:38:24  oliver
   Initial revision

*/

#ifndef _utilgmx_h
#define _utilgmx_h

static char *SRCID_utilgmx_h = "$Id: utilgmx.h,v 1.6 2003/10/09 09:15:27 oliver Exp $";

#include <stdarg.h>
#include "typedefs.h"
#include "vec.h"
#include "tpxio.h"

#define PI 3.14159265358979323844
#define AMU 1.66054e-27                 /* 1u = AMU kg */
#define AvogadroConstant 6.02214e+23    /* mol^-1 */

#define aswap(v,i,j) {  \
  atom_id temp;         \
                        \
  temp=v[i];            \
  v[i]=v[j];            \
  v[j]=temp;            \
}

#ifndef MIN
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif
#define NEWLINE(period, count) ((count)+1) % (period) ? "" : "\n"

extern real ***grid3_alloc(int nx, int ny, int nz);
extern real **grid2_alloc(int nx, int ny);
extern void msg     (const char *format, ...);
extern void dmsg    (const char *format, ...);
extern void dprintf (const char *format, ...);
extern void dfprintf (const char *format, ...);
extern void quicksort (atom_id v[], int left, int right);
extern real ldist (rvec x, rvec p, rvec c);
extern real dt_tpx (char *fn);
extern int  list_add_atomid (const atom_id, int *, atom_id *);

static gmx_inline void new_rvec (rvec r, real x, real y, real z)
{
  r[XX] = x;
  r[YY] = y;
  r[ZZ] = z;
};

#endif /* _utilgmx_h */
