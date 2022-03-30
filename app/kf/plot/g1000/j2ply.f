CTITLE J2PLY
 
      subroutine j2ply( num, xp, yp )
 
*     Draw a line between the points given in xp and yp arrays
*     The specific routine used depends on the dash pattern which
*     was set in jlstl.  Direct call to curve or curved.
 
      include 'g1000.h'
 
* PASSED VARIABLES
 
*   num     - NUmber of points to plot
 
 
      integer*4 num
 
*   xp(num), yp(num)    - User coordinates of points
 
 
      real*4 xp(num), yp(num)
 
****  Check dash type to see which routine we should use
 
      if( jdash.gt.1 ) then
          call curved( xp, yp, num)
      else
          call curve( xp, yp, num)
      end if
 
****  Thats all
      return
      end
 
