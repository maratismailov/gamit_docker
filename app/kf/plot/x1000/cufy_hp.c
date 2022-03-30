#include "x1000.h"

float cufy( yp )
 
float *yp ;  /* User y coorindate */

{
int cwpy() ; /* Convert world to pixel */
int iy ;     /* Pixel coordinayte of point */
float cu ;   /* Ratio of pixel to window */

iy = cwpy( *yp ) ;

cu = -(float) iy / (float) Window_height ;

return (cu ) ;

}
