CTITLE J2PLY
 
      subroutine j2ply( num, xp, yp )
 
*     Fortran implemenation of j2ply.  Emulated by a
*     sequence of calls to j2drw.  This may cause some problems
*     with dashed lines but avoids internal storgage for the
*     arrays in the x1000 package.
 
* PASSED VARIABLES
 
*   num     - Number of points to be plotted
 
      integer*4 num
 
*   xp(num) - X coordinates (world coordinates)
*   yp(num) - Y coordinates (world coordinates)
 
      real*4 xp(num), yp(num)
 
* LOCAL VARIABLES
 
*   i       - Loop counter
 
      integer*4 i
 
***** First move to first point
      call j2mov( xp(1), yp(1) )
 
*     Now draw the line
      do i = 2, num
          call j2drw( xp(i), yp(i))
      end do
 
****  Thats all
      return
      end
 
