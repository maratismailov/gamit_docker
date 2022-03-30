#include "x1000.h"
 
void
j2drw_( xp, yp )

float *xp, *yp;   /* Draw a line from current position to world 
                     coordinates xp, yp */

{
/* Local variable definitions */
int ix, iy;   /* Pixel coordinates of the new world coordinates
                   for the end of the line */
int cwpx(), cwpy() ; /* Functions which return pixel coordinates using
                        all of the current world and virtual settings */
void check_window() ; /* MOD TAH 210523: Added declaration */

/* Check status of windows */
/* MOD TAH 210523: Removed void; this is function call */
check_window();

/* Get the current location of the cursor in pixel coordinate*/
ix = cwpx( *xp ) ;
iy = cwpy( *yp ) ;

/* Now draw the line using current attributes */
XDrawLine(display, win, gc, ix, iy, Curx, Cury );
/* XFlush(display); Skip and see waht happens */

/* Save the current postion */
Curx = ix ; Cury = iy ;

}

