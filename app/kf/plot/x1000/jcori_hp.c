#include "x1000.h"
 
void
jcori( ox,oy,oz, px,py,pz )

/* Routine to save the orientation of text.  Saved as the tangent
   of the baseline orientation */

float *ox,*oy,*oz;  /* Orientation vecotr of baseline of text,
                     z component ignored */
float *px,*py,*pz;  /* Orientation of vector perpendicular to text,
                     (all components ignored, assummed to be 
                      rotation of baseline */

{
/* Compute the tangent of the orientation of the baseline
   vector.  Check that not too large for vertical lines (in
   this case save a large number) */
float mox ; /* Modified version of ox so that we donot divide 
               by zero */
mox = *ox ; 
if( mox == 0.0 ) mox = *oy / 1.e6 ;
jor_tan = (*oy)/(mox) ;

} 
