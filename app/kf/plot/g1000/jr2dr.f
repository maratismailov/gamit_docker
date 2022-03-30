CTITLE JR2DR
 
      subroutine jr2dr( dxp, dyp )
 
*     Draws line from current position to current position
*     plus dxp, dyp.  Use Line or lined
*     once we get the current position of the cursor.
 
      include 'g1000.h'
 
* PASSED VARIABLES
 
*   dxp, dyp        - USer coordinates of relative draw
 
 
      real*4 dxp, dyp
 
* LOCAL VARIABLES
 
*   ix,iy       - Current pen position in plotter coordinates
 
 
      integer*4 ix,iy
 
*   xs, ys      - Pen poistion in user coordinates
*   xp, yp      - Final position
*   cpux, cpuy  - Routine to convert plotter coordinates to user
*               - cooridinates.
 
 
 
      real*4 xs, ys, xp, yp, cpux, cpuy
 
****  Get the current pen position, in plotter coordinates
 
      call mxmy(ix,iy)
 
*     Convert plotter coordinates to user coordinates
 
      xs = cpux(ix)
      ys = cpuy(iy)
 
*     Get final position
      xp = xs + dxp
      yp = ys + dyp
 
****  Depending of dash type use line or lined.
 
      if( jdash.gt.1 ) then
          call lined( xs, ys, xp, yp )
      else
          call line( xs, ys, xp, yp )
      end if
 
***** Thats all
      return
      end
 
