#include "x1000.h"
   
void
jview( vl, vr, vb, vt )

float *vl, *vr, *vb, *vt ; /* Virtual coordinates of the edges of 
                              the current viewport.  These are
                              simply saved in the global variable
                              list */

{
/* Save the values */
xvl = *vl ; xvr = *vr ; xvb = *vb ; xvt = *vt ;
}

