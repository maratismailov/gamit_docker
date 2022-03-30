#include "x1000.h"

/* Function to return the world coordinates of the piyel point
   iy.  If clipping in turned on then the edge of the current window
   is returned */

float cpuy ( iy ) 

int iy ;        /* Piyel coordinates to be converted to world 
                   coordinates  */

{
/* Local variables */
float diy_bl ; /* Offset of piyel point from bottom left cornor */
float dpdyw  ; /* Ratio of piyels to y world coordinates */
float py_b   ; /* Piyel coordinates of bottom edge of
                  of the view port */

float py ;       /* World coordinates of point */

/* Get scale of piyels to world coordinates */
dpdyw = (float) Window_height*(xvt-xvb)/(xwt-xwb) ;

/* Get the piyel coordinates of the left and right edge of viewport */
py_b = (1.0 - xvb) * (float) Window_height ;

/* Get piyel position relative to left edge of view port. */
diy_bl = py_b - (float) iy ;

/* Now get the world coorinates of the point */
py = xwb + diy_bl/dpdyw ;

return ( py ) ;
}
