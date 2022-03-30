#include "x1000.h"

/* Function to return the world coordinates of the pixel point
   ix.  If clipping in turned on then the edge of the current window
   is returned */

float cpux ( ix ) 

int ix ;        /* Pixel coordinates to be converted to world 
                   coordinates  */

{
/* Local variables */
float dix_bl ; /* Offset of pixel point from bottom left cornor */
float dpdxw  ; /* Ratio of pixels to x world coordinates */
float px_l   ; /* Pixel coordinates of left edge of
                  of the view port */

float px ;       /* World coordinates of point */

/* Get scale of pixels to world coordinates */
dpdxw = (float) Window_width*(xvr-xvl)/(xwr-xwl) ;

/* Get the pixel coordinates of the left and right edge of viewport */
px_l = xvl * (float) Window_width ;

/* Get pixel position relative to left edge of view port. */
dix_bl = (float) ix - px_l ;

/* Now get the world coorinates of the point */
px = xwl + dix_bl/dpdxw ;

return ( px ) ;
}
