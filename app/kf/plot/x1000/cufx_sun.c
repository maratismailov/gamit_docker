#include "x1000.h"

float cufx_( xp )

float *xp ;  /* User x coorindate */

{
int cwpx() ; /* Convert world to pixel */
int ix ;     /* Pixel coordinayte of point */
float cu ;   /* Ratio of pixel to window */

ix = cwpx( *xp ) ;

cu = (float) ix / (float) Window_width ;

return (cu ) ;

}
