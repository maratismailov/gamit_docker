#include "x1000.h"

/* Function to return the pixel coordinates of the world point
   yp.  If clipping in turned on then the edge of the current window
   is returned */

int cwpy ( yp ) 

float yp ;        /* World Y-coordinate to be converted to pixels */

{
/* Local variables */
float dyp_bl ; /* Offset of world point from bottom left cornor */
float dpdyw  ; /* Ratio of pixels to y world coordinates */
float py_b, py_t ; /* Pixel coordinates of bottom and top edges of
                      of the view port */

int py ;       /* Pixel coordinates of point */

/* Get world coordinate relative to bottom side */
dyp_bl = yp - xwb ;

/* Get scale of pixels to world coordinates */
dpdyw = (float) Window_height*(xvt-xvb)/(xwt-xwb) ;

/* Get the pixel coordinates of the bottom and top edge of viewport */
py_b = (1.0-xvb) * (float) Window_height ;
py_t = (1.0-xvt) * (float) Window_height ;

/* Now get the pixel coorinates of the point */
py = (int) (py_b - dpdyw*dyp_bl + 0.5) ;

return ( py ) ;
}
