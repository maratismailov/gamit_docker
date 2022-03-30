#include "x1000.h"
 
void
DrawMark( ix, iy, type )

int ix, iy ;   /* Screen coordinates of mark to be drawn */
int type   ;   /* Type of mark (There are 5 types compatible
                  the NCAR Graphics PolyMarkers 
                  1    .
                  2    +
                  3    *
                  4    o
                  5    X
                  All others will default to type 1 */

{

/* Local variables */
XSegment mark[10] ; /* Each of the line segments used to form
                        the mark */
int  nsegs         ; /* Number of segments to use in the mark */
int  mx, my        ; /* Start point when rectangle points are drawn */
int  mwidth, mheight ; /* Mark width and height */

switch ( type )
  {
     case -2:        /*Plus sign, form the mark structure */
         mark[0].x1 = ix - usr_charsize/2 ; mark[0].y1 = iy ;
         mark[0].x2 = ix + usr_charsize/2 ; mark[0].y2 = iy ;
         mark[1].x1 = ix ; mark[1].y1 = iy + usr_charsize/2 ;
         mark[1].x2 = ix ; mark[1].y2 = iy - usr_charsize/2 ;
         XDrawSegments(display, win, gc, mark, 2);
         break ;
      case -3:      /*Asterix shape  */
         mark[0].x1 = ix - usr_charsize/3 ; mark[0].y1 = iy ;
         mark[0].x2 = ix + usr_charsize/3 ; mark[0].y2 = iy ;
         mark[1].x1 = ix ; mark[1].y1 = iy + usr_charsize/3 ;
         mark[1].x2 = ix ; mark[1].y2 = iy - usr_charsize/3 ;
         mark[2].x1 = ix - usr_charsize/3 ; mark[2].y1 = iy - usr_charsize/3;
         mark[2].x2 = ix + usr_charsize/3 ; mark[2].y2 = iy + usr_charsize/3;
         mark[3].x1 = ix - usr_charsize/3 ; mark[3].y1 = iy + usr_charsize/3;
         mark[3].x2 = ix + usr_charsize/3 ; mark[3].y2 = iy - usr_charsize/3;
         XDrawSegments(display, win, gc, mark, 4);
         break ;
      case -4:      /* Circle */
         mwidth = (3*usr_charsize)/2 ;
         mheight =(3* usr_charsize)/2 ;
         XDrawArc(display, win, gc, ix-mwidth/2, iy-mwidth/2, mwidth, mheight, 0, 360*64);
         break ;
      case -5:      /* X symbol */
         mark[0].x1 = ix - usr_charsize/2 ; 
         mark[0].y1 = iy + usr_charsize/2 ;
         mark[0].x2 = ix + usr_charsize/2 ; 
         mark[0].y2 = iy - usr_charsize/2 ;
         mark[1].x1 = ix - usr_charsize/2 ; 
         mark[1].y1 = iy - usr_charsize/2 ;
         mark[1].x2 = ix + usr_charsize/2 ; 
         mark[1].y2 = iy + usr_charsize/2 ;
         XDrawSegments(display, win, gc, mark, 2);
         break ;
      default:
         mwidth = usr_charsize/8  + 1 ; mheight = usr_charsize/8 + 1;
         mx = ix - mwidth/2 ; my = iy - mheight/2 ;
         XDrawRectangle(display, win, gc, mx, my, mwidth, mheight );
  }
}
          
