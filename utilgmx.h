/* $Id: utilgmx.h,v 1.8 2009/07/22 15:31:14 oliver Exp $ */

#ifndef _utilgmx_h
#define _utilgmx_h

static char *SRCID_utilgmx_h = "$Id: utilgmx.h,v 1.8 2009/07/22 15:31:14 oliver Exp $";

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
extern real ldist (const rvec x, const rvec p, const rvec c);
extern real dt_tpx (char *fn);
extern int  list_add_atomid (const atom_id, int *, atom_id *);

static inline void new_rvec (rvec r, real x, real y, real z)
{
  r[XX] = x;
  r[YY] = y;
  r[ZZ] = z;
};

#endif /* _utilgmx_h */
