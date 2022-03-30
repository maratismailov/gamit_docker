#include "x1000.h"
  
void
jr2dr( dxp, dyp )

float *dxp, *dyp;   /* Draw a line from current position to world 
                     relative coordinates dxp, dyp */

{
/* Local variable definitions */
int ix, iy;   /* Pixel coordinates of the new world coordinates
                   for the end of the line */
int  cwpx(), cwpy() ; /* Functions which return pixel coordinates using
                        all of the current world and virtual settings */
float cpux(), cpuy(); /* Return world coordinates for a set of pixel
                         coordinates */
float wcx, wcy ;  /* World corrdinates of current cusor positions */ 
void check_window() ; /* MOD TAH 210523: Added declaration */

/* Check status of windows */
/* MOD TAH 210523: Removed void; this is function call */
check_window();

/* Get the world coordinates of the current location of the
   cursor */
wcx = cpux( Curx ) ;
wcy = cpuy( Cury ) ;

/* Get the new location of the cursor in pixel coordinate*/
ix = cwpx( wcx + *dxp ) ;
iy = cwpy( wcy + *dyp ) ;

/* Now draw the line using current attributes */
XDrawLine(display, win, gc, ix, iy, Curx, Cury );

/* Save the current postion */
Curx = ix ; Cury = iy ;

}

