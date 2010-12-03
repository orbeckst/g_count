/* $Id: plt.c,v 1.3 2004/02/10 01:50:23 oliver Exp $
   $Log: plt.c,v $
   Revision 1.3  2004/02/10 01:50:23  oliver
   added average density output

   Revision 1.2  2003/11/25 17:44:19  oliver
   - plt output (and ascii dump) have now proper units
   - moved changed write_plt into density_write_plt and moved ascii dump into
     plt.c: density_write_ascii

   Revision 1.1  2003/11/24 19:21:18  oliver
   output for gOpenMol binary density seems to work


   writing density in plt format gOpenMol
   http://www.csc.fi/gopenmol/developers/plt_format.phtml

   Stuff ripped from Christoph Freudenberger's g_sdf.c
   g_sdf.c,v 1.25 2003/08/05 11:00:00 cfreuden

   (c) 2004 Oliver Beckstein <oliver@biop.ox.ac.uk>
   This program is made available under the terms of the GNU Public License.

*/
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "futil.h"
#include "statutil.h"
#include "grid3D.h"

static char *SRCID_plt_c = "$Id: plt.c,v 1.3 2004/02/10 01:50:23 oliver Exp $";

static void i_write(FILE *output, int value)
{
  fwrite(&value,sizeof(int),1L,output);
}

static void f_write(FILE *output,float value)
{
  fwrite(&value,sizeof(float),1L,output);
}

/* write plt binary density file for gOpenMol
   Normalise each grid point by dividing by unit.
*/
void density_write_plt (char *fnm, t_tgrid *tgrid, real unit)
{
  FILE *fp;
  int  i,j,k;
  real norm;

  norm = 1./unit;
  fp = ffopen (fnm, "wb");

  /* rank */
  i_write(fp,3);

  /* Type of surface */
  i_write(fp,42);

  /* Zdim, Ydim, Xdim */
  for (i=ZZ; i>=XX; i--)
    i_write(fp,tgrid->mx[i]);

  /* [Z,Y,X][min,max] (box corners in Angstroem)
     see grid3D.h for definition of t_grid (note: distances in nm)
  */
  for (i=ZZ; i>=XX; i--) {
    f_write(fp,10. * tgrid->a[i]);
    f_write(fp,10. * tgrid->b[i]);
  }

  for(k=0;k<tgrid->mx[ZZ];k++) {
    for(j=0;j<tgrid->mx[YY];j++) {
      for(i=0;i<tgrid->mx[XX];i++) {
        f_write(fp,norm * tgrid->grid[i][j][k]);
      }
    }
  }
  ffclose(fp);
};


void density_write_ascii (char *fnm, t_tgrid *tgrid, real unit)
{
  FILE *fp;
  int  i,j,k;
  real norm;

  norm = 1./unit;
  fp = ffopen (fnm, "w");

  for(k=0;k<tgrid->mx[ZZ];k++) {
    for(j=0;j<tgrid->mx[YY];j++) {
      for(i=0;i<tgrid->mx[XX];i++) {
        fprintf(fp,"%.6f ",norm * tgrid->grid[i][j][k]);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

