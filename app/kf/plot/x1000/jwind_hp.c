#include "x1000.h"
 
void
jwind( wl, wr, wb, wt )

float *wl, *wr, *wb, *wt ; /* World coordinates of the edges of 
                              the current viewport.  These are
                              simply saved in the global variable
                              list */

{
/* Save the values */
xwl = *wl ; xwr = *wr ; xwb = *wb ; xwt = *wt ;
}

