/*
 * $Id: xf.c,v 1.2 2004/03/18 01:04:28 oliver Exp $
 * 
 * please refer to 
 * http://www.fhi-berlin.mpg.de/grz/pub/xfarbe/xfarbe-2.6/html/index.html 
 * for details about xfarbe and the file format
 */
static char *SRCID_xf_c = "$Id: xf.c,v 1.2 2004/03/18 01:04:28 oliver Exp $";

#include <string.h>
#include <ctype.h>
#include "sysstuff.h"
#include "string2.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "smalloc.h"
#include "utilgmx.h"
#include "xf.h"


/* minimum parameter file; colours are defined in xf.h */
bool xf_write_XFarbe (int ncolours) {
  FILE *XFarbe;
  int i;
  
  if (!(XFarbe=fopen(XFARBE,"w"))){
    msg("error: Cannot open file " XFARBE "for writing.\n");
    return FALSE;
  }
  if (ncolours > XFARBE_NCOLS) {
    msg("warning: The requested number of levels %d is not available. You will \n"
	"         only get %d instead.\n",ncolours,XFARBE_NCOLS);
    ncolours = XFARBE_NCOLS;
  }

  msg("Writing resource file './" XFARBE "' for xfarbe with %d colours.\n",ncolours);
  dmsg("Use it with one of\n"
      "bash$ export XAPPLRESDIR=.\n"
      "tcsh> setenv XAPPLRESDIR .\n");
  fpxfi("ncol",ncolours);
  fpxfi("mode",0);
  fpxf ("autolev","off");
  fpxf ("xyratio","on");
  fpxf ("legend", "on");
  fpxf ("legend_format","%4.2f");
  fpxf ("x-annotation", "on");   /* annotation on but increment 0 */  
  fpxf ("y-annotation", "on");   /* allows data probing in xfarbe */
  fpxfi("x-annotation-incr",0);  /* _with_ real x and y values    */
  fpxfi("y-annotation-incr",0);  /* describe in */
  /* http://www.fhi-berlin.mpg.de/grz/pub/xfarbe/xfarbe-2.6/html/index.html#Data%20Probing */ 
  for(i=0;i<ncolours;i++) 
    fprintf(XFarbe,"xfarbe.col%d:   #%.6X\n",i,xfarbe_colours[i]);

  fclose(XFarbe);
  return TRUE;
}


FILE *xf_header (const char *fname,const char *header,
                    const int nx,const int ny) {      
  /* Open data file in xfarbe format from the grid tgrid.  See 
     http://www.fhi-berlin.mpg.de/gnz/pub/xfarbe/xfarbe-2.6/html/#Format%20of%20Data%20Files 
     for the file format definition
  */
  FILE *fOut;
  char xfheading[XFARBE_MAX_HEADER];

  if (!(fOut = ffopen (fname, "w"))) return NULL;
  if (strlen(header)>0) {
    strncpy(xfheading,header,XFARBE_MAX_HEADER-1);
    xfheading[XFARBE_MAX_HEADER-1] = '\0';
    fprintf(fOut,"%s",xfheading);
  };
  fprintf(fOut,"\n");
  fprintf(fOut,"%d %d\n",nx,ny);
  return fOut;
};


/* append int+1 levels (upto value max) to the file */
void xf_append_levels (FILE *fData, int ncols, real max) {
  int i;
  real delta;

  delta = max/(ncols-1);
  fprintf(fData,"%d\n",ncols-1);        /* number of levels */
  fprintf(fData, "%d  %f\n",0,delta/2); /* first level filters out noise */
  for(i=1;i<ncols-1;i++)                /* equally spaced */
    fprintf(fData, "%d  %f\n",i,i*delta);
  fprintf(fData, "%d\n", ncols-1);      /* everything else */
  return;
}

void xf_append_axes_annotation (FILE *fData, real x1, real x2, 
			     real y1, real y2) {
  fprintf(fData,"%f  %f\n",x1, x2);
  fprintf(fData,"%f  %f\n",y1, y2);
  return;
}



/* create xmgr file with subtitle
   (xvgropen() does not have subtitles)
*/
FILE *xmgropen(char *fn,char *title,char *subtitle,
	       char *xaxis,char *yaxis)
{
  FILE *xmgr;
  time_t t;
  
  xmgr=(FILE *)ffopen(fn,"w");
  fprintf(xmgr,"# Grace project file");
  fprintf(xmgr,"# This file was created by %s\n",Program());
  time(&t);
  fprintf(xmgr,"# All this happened at: %s",ctime(&t));
  fprintf(xmgr,"#\n");
  fprintf(xmgr,"@version 50102\n"); 
  fprintf(xmgr,"@with g0\n"); 
  fprintf(xmgr,"@    title \"%s\"\n",title);
  fprintf(xmgr,"@    subtitle \"%s\"\n",subtitle);
  fprintf(xmgr,"@    xaxis  label \"%s\"\n",xaxis);
  fprintf(xmgr,"@    yaxis  label \"%s\"\n",yaxis);
  fprintf(xmgr,"@target G0.S0\n"); 
  fprintf(xmgr,"@type nxy\n"); 
  
  return xmgr;
}

