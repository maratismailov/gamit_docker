#include "x1000.h"

/* Function to return the pixel coordinates of the world point
   xp.  If clipping in turned on then the edge of the current window
   is returned */

int cwpx ( xp ) 

float xp ;        /* World X-coordinate to be converted to pixels */

{
/* Local variables */
float dxp_bl ; /* Offset of world point from bottom left cornor */
float dpdxw  ; /* Ratio of pixels to x world coordinates */
float px_l, px_r ; /* Pixel coordinates of left and right edges of
                      of the view port */

int px ;       /* Pixel coordinates of point */

/* Get world coordinate relative to left side */
dxp_bl = xp - xwl ;

/* Get scale of pixels to world coordinates */
dpdxw = (float) Window_width*(xvr-xvl)/(xwr-xwl) ;

/* Get the pixel coordinates of the left and right edge of viewport */
px_l = xvl * (float) Window_width ;
px_r = xvr * (float) Window_width ;

/* Now get the pixel coorinates of the point */
px = (int) (px_l + dpdxw*dxp_bl + 0.5) ;

return ( px ) ;
}
