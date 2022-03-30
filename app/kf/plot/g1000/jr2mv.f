CTITLE JR2MV
 
      subroutine jr2mv( dxp, dyp )
 
*     Moves the pen to user coordinates xp, yp.  Uses the SPPS
*     routine FRSTPT.
 
* PASSED VARIABLES
 
*   dxp, dyp        - User coorinates change in user coordinates
 
 
      real*4 dxp, dyp
 
* LOCAL VARIABLES
 
*   ix,iy       - Current pen position in plotter coordinates
 
 
      integer*4 ix,iy
 
*   xp, yp      - Final position
*   xs, ys      - Start coordinates (curtent pen position)
*   cpux, cpuy  - Routine to convert plotter coordinates to user
*               - cooridinates.
 
 
 
      real*4 xp, yp, xs, ys, cpux, cpuy
 
****  Get the current pen position, in plotter coordinates
 
      call mxmy(ix,iy)
 
*     Convert plotter coordinates to user coordinates
 
      xs = cpux(ix)
      ys = cpuy(iy)
 
*     Get final position
      xp = xs + dxp
      yp = ys + dyp
 
*     Call move to absolute position
 
      call frstpt( xp, yp )
 
*     Thats all
      return
      end
 
