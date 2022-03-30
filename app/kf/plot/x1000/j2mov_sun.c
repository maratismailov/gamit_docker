#include "x1000.h"

void
j2mov_( xp, yp )

float *xp, *yp;   /* Move the cursoreto world 
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

/* Save the current postion */
Curx = ix ; Cury = iy ;

}

