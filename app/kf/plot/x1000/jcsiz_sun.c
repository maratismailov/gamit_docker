#include "x1000.h"

void
jcsiz_( width_mm, height_mm, gap_mm )

float *width_mm ;  /* Width oof characters in mm (only one used) */
float *height_mm, *gap_mm ;  /* Height and gap between characters 
                                Neither of which is used */

{
/* Convert mm to pixels approximately, and save in the user
   character size.  This value will be used by error bars, and
   tick sizes.  The size of the font being used will be used for
   spacing text */

usr_charsize = (int) *width_mm / 0.50;

}

