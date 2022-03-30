#include "x1000.h"

void
j2mrk( xp, yp, type )

float *xp, *yp;   /* mrke the cursoreto world 
                     coordinates xp, yp */
int *type ;       /* Point type.  If the value is negative
                     then a polymarker is put out, if 
                     positive, then the type is assumed to
                     to an ASCII characeter code and this
                     if written out. */

{
/* Local variable definitions */
int ix, iy;   /* Pixel coordinates of the new world coordinates
                   for the end of the line */
int cwpx(), cwpy() ; /* Functions which return pixel coordinates using
                        all of the current world and virtual settings */
void DrawMark(), DrawChar(); /* Functions which return type of point */

/* Check status of windows */

void check_window();

/* Get the current location of the cursor in pixel coordinate*/
ix = cwpx( *xp ) ;
iy = cwpy( *yp ) ;

/* Save the current postion */
Curx = ix ; Cury = iy ;

/* Now see what type of point this is. */
if( *type <= 0 ) DrawMark( ix, iy, *type ) ;
else             DrawChar( ix, iy, *type ) ;

}

