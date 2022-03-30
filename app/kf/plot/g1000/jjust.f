CTITLE JJUST
 
      subroutine jjust( wj, hj )
 
*     Sets the justification of strings.  There is large quantification
*     problem here.  All we do is look at xj orientation and set
*     JCENTER to closest value we can.
 
      include 'g1000.h'
 
* PASSED VARIABLES
 
*   wj, hj  - Requested justification in width and height.
 
 
 
      real*4 wj, hj
 
*     Set closest value and hope for the best
 
      jcenter = nint(wj-0.5)
      gheight = hj - 0.5
 
      return
      end
 
