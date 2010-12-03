/* $Id: plt.h,v 1.3 2004/02/10 01:50:23 oliver Exp $
   $Log: plt.h,v $
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

#ifndef _plt_h
#define _plt_h


void density_write_plt   (char *, t_tgrid *, real);
void density_write_ascii (char *, t_tgrid *, real);

#endif /* _plt_h */
