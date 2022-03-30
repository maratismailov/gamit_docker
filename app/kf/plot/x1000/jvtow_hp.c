#include "x1000.h"

void
jvtow( vx, vy, wx, wy, wz )

float *vx, *vy ; /* Virtual coordinates of point */
float *wx, *wy, *wz ; /* World coordinates of point */

{

int ix, iy ; /* Pixel coordinates (which we then convert
                to world coordinates) */
float cpux() , cpuy() ; /* Functions to convert pixels to
                world coordinates) */

ix = *vx * (float) Window_width ;
iy = (1.0 - (*vy)) * (float) Window_height ;

*wx = cpux(ix) ;
*wy = cpuy(iy) ;
*wz = 0.0 ;
}

