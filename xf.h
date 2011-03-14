/*
 * $Id: xf.h,v 1.2 2004/03/18 01:04:28 oliver Exp $
 *
 * please refer to 
 * http://www.fhi-berlin.mpg.de/grz/pub/xfarbe/xfarbe-2.6/html/index.html 
 * for details about xfarbe and the file format
 */

#ifndef _xf_h
#define _xf_h

#ifdef CPLUSPLUS
extern "C" {
#endif

/***************************************************
 *            xfarbe routines (2D density plots)
 ***************************************************/
extern gmx_bool xf_write_XFarbe (int);
/* write a parameter file with 16 colours from xfrabe_colours */
extern FILE *xf_header (const char *,const char *,const int,const int);
/* write open file, write heading and dimensions X and Y; returns file
   handle for open file 
*/  
extern void xf_append_levels (FILE *, int, real);
/* append level definitions to a XFarbe data file 
   ncols - 1 levels equally spaced up to max,
   above max, all is coloured in the last colour, number ncols-1
 */
extern void xf_append_axes_annotation (FILE *, real x1, real x2, 
					 real y1, real y2);
/* label x axis from x1 to x2, &idem y axis */ 

/*
  write customised XFarbe file 
  use with 'export XAPPLRESDIR=.' which activates the level definitions at 
  the end of the data file

  extract colours with:
     gawk '{gsub("\#","0x"); printf ("%s, ", $2)}' cols
  from the xfarbe v2.5 default XFarbe file

  colours: 0xRRGGBB (2 bytes per red, green, blue)

  http://www.ndirect.co.uk/~thomas.green/javascripts/hexAndRGB.html
*/

static u_int xfarbe_colours[] = 
      {0x000099, 0x1111DD, 0x5555FF, 0x6677DD, 0x8899AA, 0xAABB66, 
       0xCCDD33, 0xEEEE00, 0xFFFF00, 0xFFDD00, 0xFFBB00, 0xFF9900,  
       0xFF6600, 0xEE1100, 0xBB0000, 0x770000};

/* it is not easy finding good colours---the xfarbe default is already
   pretty decent...  */
/*
  static u_int xfarbe_colours[] =
       {0x1400FF, 0x0010FF, 0x0066FF, 0x00A3FF, 0x00E0BB, 
       0x00F0AA, 0x33FF88, 0x00FF00, 0x14FF00, 0x52F000, 
       0x8FEE00, 0xCCDD22, 0xFFC000, 0xFFB800, 0xFF7A00,
       0xFF3D00, 0xFF0000, 0xFF007A, 0xFF00B8, 0xFF00F5 };
  */

#define XFARBE_MAX_HEADER 79 /* if the header line is longer, xfarbe
                                does not read the file properly */
#define XFARBE_NCOLS asize(xfarbe_colours)  
#define XFARBE "XFarbe"      /* name of the parameter file */
#define fpxf(resource,str) fprintf(XFarbe,"xfarbe.%s: %s\n",(resource),(str))
#define fpxfi(resource,i)  fprintf(XFarbe,"xfarbe.%s: %d\n",(resource),(i))

/*************************************************************
 * xmgrace routine 
 * replacement for xvgropen()
 */
extern FILE *xmgropen(const char *fn,const char *title,const char *subtitle,
                      const char *xaxis,const char *yaxis);
/* Open a file, and write a title, subtitle, and axis-labels in Xmgr format */


#ifdef CPLUSPLUS
}
#endif

#endif	/* _xf_h */
