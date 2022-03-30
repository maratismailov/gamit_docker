#include "x1000.h"
 
void
jjust( wj, hj )

/* Routine to save the justification for width and height */

float *wj, *hj ; /* Width and height justification (each value
                    0.0-1.0 for left to right and bottom to top
                    justification */
{
jcenter = *wj;
jheight = *hj;
}

