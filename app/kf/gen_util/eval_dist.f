CTITLE EVAL_DIST
 
      subroutine eval_dist( X1, X2, dist)

      implicit none 
 
*     This routine computes the distabce between two points
*     whose cartesian coordinates are saved in X1 and X2
 
* PASSED VARIABLES
 
*   X1(3)       - Cartesian coordinates of site 1
*   X2(3)       - Cartesian coordinates of site 2
*   dist        - The distance between them (cord length)
 
      real*8 X1(3), X2(3), dist
 
***** Use simple formula
 
      dist = sqrt( (X2(1)-X1(1))**2 + 
     .             (X2(2)-X1(2))**2 + 
     .             (X2(3)-X1(3))**2   )
 
****  Thats all
      return
      end
 
