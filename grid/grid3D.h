/*
  $Id: grid3D.h,v 1.9 2004/03/17 23:26:29 oliver Exp $
  $Log: grid3D.h,v $
  Revision 1.9  2004/03/17 23:26:29  oliver
  - realised that I was using the TOKEN-less xdr format although I though
    that I was using the other one, ie xdr_token_grid()
  - completely switched over to TOKEN-less but left some of the TOKEN code
    in to remind myself how to do it
  - tested both format: both work but are incompatible and fail on each others
    data files (even though I removed stringent TOKEN/TOKEN-less checks that
    used to be there from development days)

  Revision 1.8  2002/08/20 01:13:50  oliver
  NEW:
  * a_ri3Dc can resize (well, shrink) the grid (only z1, z2, R)
    - needed setup_tgrid() and setup_corners() from g_ri3Dc for this
      dirty hack, so they were moved into grid3D.c
  FIXED:
    calculation of the radius from the grid corners was wrong because I
    used abs() where I should have used fabs(). Didn't have any impact
    so far but for the readjustment of the grid it matters!

  Revision 1.7  2002/08/19 17:20:47  oliver
  IMPORTANT:
    write a true density to the 3D grid file (in nm^-3)
  small things:
  - a few variables (sumP, dV...) are now double because for long
    trajectories the loss in accuracy is substantial
  - declared command line parameters as static
  - moved autoset_cavity() from g_ri3Dc into count.c (.h) so that I can
    use it for g_count as well
  - cleaned up argument list

  Revision 1.6  2002/08/17 17:37:10  oliver
  volume calculations from t_grid packaged as two functions, DeltaV()
  and gridV

  Revision 1.5  2002/08/16 23:18:28  oliver
  corrected file format description (no VERSION token)

  Revision 1.4  2002/08/15 23:51:39  oliver
  Changed file format: start directly with the int for VERSION in both formats so that they can be detected gracefully

  Revision 1.3  2002/08/15 17:05:36  oliver
  First official release.
  - fixed bug with length of xfarbe file header

  Revision 1.2  2002/08/15 14:51:11  oliver
  FINALLY: 3Dgrid output works. I settle for a format with interspersed tokens

  Revision 1.1  2002/08/14 15:09:59  oliver
  splitting off functions from g_ri3Dc into grid3D.c (and .h) which are common to the grid generator (g_ri3Dc) and analysis program(s) (probably will be only one, a_ri3Dc, codenamed 'dry county')


  collect stuff common to g_ri3Dc and a_ri3Dc



  grid: use XDR routines (similar to the gromacs libxdrf.c but less elaborate)
        advantages:
	- binary format (more efficient than huge ascii files)
	- could implement compression similar to xdr3dfcoord()
	- machine independent

  Format: (The tokens like HEADER are actually of type enum
  eGridToken; exception: there is no VERSION token; the file
  starts with an int.
  
  (VERSION)   int            # version info (to catch the inevitable 
                               format change in the future)
  HEADER    str[HEADER_MAX]  # descriptive txt
  GRIDTYPE  enum             # currently we only have regular grids
  DIMENSION int              # dim of array
  SIZE      int int ...      # number of cells per dim  \
  DELTA     real real ...    # grid width per dim        } allows reconstruction 
  ORIGIN    real real ...    # lower left corner of grid/
  GRID      real ..................................
                 (0,0,0) (1,0,0) (2,0,0) ... (N1,0,0)       \
                 (0,1,0) (1,1,0) (2,1,0) ... (N1,1,0)       \
		 ...                                        \
                 (0,N2,0) (1,N2,0) (2,N2,0) ... (N1,N2,0)   \
                 (0,0,1) (1,0,1) (2,0,1) ... (N1,0,1)       \
                 (0,1,1) (1,1,1) (2,1,1) ... (N1,1,1)       \
		 ...                                        \
                 (0,N2,0) (1,N2,0) (2,N2,0) ... (N1,N2,0)   \
		 ... 
		 ...
		 .... (N1,N2,N3)

		 (ordered by z-slices)

*/


#ifndef _grid3D_h
#define _grid3D_h

#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "fatal.h"
#include "count.h"
#include "xdr_grid.h"

#define OCCUPIED_MIN 0       /* occupancy must be greater than this to
                                  be counted as occupied for volume
                                  calculations */
#define XDR_GRID(xdr,xg) xdr_grid((xdr),(xg))
/* alternative (and incompatible) format, which works, too. 
   See xdr_grid.c for more comments.
#define XDR_GRID(xdr,xg) xdr_token_grid((xdr),(xg)) 
*/

typedef struct {
  /* everything about the grid */
  real ***grid;    /* MX x MY x MZ */
  int  mx[3];      /* MX, MY, MZ */
  real a[3], b[3]; /* bottom left and top right corner of the grid in nm */ 
  rvec Delta;      /* spatial resolution in nm */
} t_tgrid;

typedef struct {
  /* all results in one place */
  real      volume;      /* calculated volume */
  real      reff;        /* effective radius */
  t_tgrid *tgrid;        /* full grid structure */
  real   **profile;      /* (z, R*(z)) */
  real   **xyproj;       /* projection of the occupancy in the xy plane */
  real   **xzproj;       /* projection of the occupancy in the xz plane */
} t_result;

extern bool grid_write (FILE *,t_tgrid *,char *);
extern bool grid_read  (FILE *,t_tgrid *,char *);
extern gmx_inline double DeltaV (t_tgrid *);
extern gmx_inline double DeltaA (t_tgrid *);
extern gmx_inline double gridV  (t_tgrid *);
extern gmx_inline double gridA  (t_tgrid *);
extern bool setup_tgrid (t_tgrid *, t_cavity *, rvec);
extern void setup_corners (rvec [], t_cavity *);

#endif  /* _grid3D_h */
