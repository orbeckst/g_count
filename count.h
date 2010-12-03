/* $Id: count.h,v 1.4 2002/08/19 17:16:36 oliver Exp $
   
   specialised functions which are used in in g_count and g_flux

   $Log: count.h,v $
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

   Revision 1.1.2.6  2002/08/11 14:45:16  oliver
   rename write_index() to fwrite_index() to avoid collision with a
   gromacs function in index.c

   Revision 1.1.2.5  2002/06/04 21:48:01  oliver
   current working version of g_flux

   Revision 1.1.2.4  2002/05/29 18:43:23  oliver
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

   Revision 1.1.2.3  2002/05/28 11:37:09  oliver
   put DR_DEFAULT back into g_count (and remove from count.h)

   Revision 1.1.2.2  2002/05/27 23:52:04  oliver
   first working version of g_flux (working means: compiles, does not core dump, results look vaguely sensible)

   Revision 1.1.2.1  2002/05/27 12:21:25  oliver
   Moved common functions and definition for g_count, g_sldist, g_flux into count.[ch]. g_count compiles

*/

#ifndef _count_h
#define _count_h

static char *SRCID_count_h = "$Id: count.h,v 1.4 2002/08/19 17:16:36 oliver Exp $";

#include <assert.h>
#include "typedefs.h"
#include "names.h"
#include "vec.h"
#include "princ.h"
#include "rmpbc.h"
#include "readinp.h"
#include "utilgmx.h"

/* types of index: atoms or molecules */
enum ndxtype    {etxATOM, etxMOL, etxNR};  
/* units to mesure a number density n (given in nm^-3) in:
   n' = n/DensUnit

   SPC: 
    rho0 = 0.9669;      # bulk SPC water at 300K (Gromacs 3, ffgmx)
                        # new: from max of P(N), labbook II, p159:
                        #            0.966894±0.000462 g cm-³
                        #      from N/<V>, II, p158
                        #            0.966458±0.00458  g cm-³
    # volume of one water molecule at T=300K, P=1bar (Labbook II p159)
    # in nm³
    v_water   = 0.030938;
    1/v_water = 32.3227

   MOLAR:
    1 nm^-3 = 1/N_Avogadro * (10^-8 dm)^-3

   ANG: in Angstroem^-3
    1 nm^-3 = 1 (10 A)^-3 
 */
enum eDensUnit {eduUNITY, eduSPC, eduMOLAR, eduANG, eduNR};
extern char *endxtype_names[etxNR+1];
extern char *eDensUnit_names[eduNR+1];
extern real DensUnit[eduNR];


#define ENDXTYPE(e)      ENUM_NAME(e,etxNR,endxtype_names)
#define EDENSUNITTYPE(e) ENUM_NAME(e,eduNR,eDensUnit_names)


typedef struct {
  rvec    axis;       /* vector parallel to pore axis */
  rvec    cpoint;     /* any point on the axis */
  real    radius;     /* max radius of the cylindrical pore */
  real    z1, z2;     /* confined between z1 and z2 */
  real    vol;        /* volume V=r^2 pi (z2-z1) */
} t_cavity;


typedef struct {      /* times in the simulation */
  int   tot_frames;   /* number of frames read */
  real  t_tot;        /* simulation time (ps) */
  real  delta_t;      /* time step (ps) */
  real  tau;          /* threshold for tracking a particle: 
			 t_cav > tau */
} t_simtime;


extern bool autoset_cavity (t_cavity *,matrix,int,t_pargs []);
extern bool bInCavity (rvec x, t_cavity *cavity);
extern void x2molxcm (t_topology *top, atom_id *molndx, int gnmol, matrix box, 
		      rvec *x, rvec *x_s, rvec *xmol_cm);
extern real npbc2com (t_topology *, atom_id *, matrix, rvec *, 
		      atom_id, rvec);
extern int mols_from_index (atom_id *index, int gnx, t_block *mols, 
		     atom_id *molndx, int max_mol);
extern void print_ldist (rvec x, t_cavity *cavity, atom_id idx, real mass);
extern void fwrite_index (FILE *, atom_id *, enum ndxtype,  int, t_topology *, char *);
extern void dump_atomlist (FILE *, enum ndxtype, atom_id *, int, t_topology *);

#endif   /* _count_h */
