/* $Id: count.h,v 1.7 2009/07/02 16:05:48 oliver Exp $
   
   specialised functions which are used in in g_count and g_flux
*/

#ifndef _count_h
#define _count_h

static char *SRCID_count_h = "$Id: count.h,v 1.7 2009/07/02 16:05:48 oliver Exp $";

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
extern bool bInCavity (const rvec x, const t_cavity *cavity);
extern void x2molxcm (t_topology *top, int ePBC, atom_id *molndx, int gnmol, matrix box, 
		      rvec *x, rvec *x_s, rvec *xmol_cm);
extern real npbc2com (t_topology *, atom_id *, matrix, rvec *, 
		      atom_id, rvec);
extern int mols_from_index (atom_id *index, int gnx, t_block *mols, 
		     atom_id *molndx, int max_mol);
extern void print_ldist (const rvec x, const t_cavity *cavity, const atom_id idx, const real mass);
extern void fwrite_index (FILE *, atom_id *, int, t_topology *, char *);
extern void dump_atomlist (FILE *, enum ndxtype, atom_id *, int, t_topology *);

#endif   /* _count_h */
