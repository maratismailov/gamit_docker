CTITLE J2DRW
 
      subroutine j2drw( xp, yp )
 
*     Draws line from current position to xp, yp.  Use Line or lined
*     once we get the current position of the cursor.
 
      include 'g1000.h'
 
* PASSED VARIABLES
 
*   xp, yp      - USer coordinates to draw to.
 
 
      real*4 xp, yp
 
* LOCAL VARIABLES
 
*   ix,iy       - Current pen position in plotter coordinates
 
 
      integer*4 ix,iy
 
*   xs, ys      - Pen poistion in user coordinates
*   cpux, cpuy  - Routine to convert plotter coordinates to user
*               - cooridinates.
 
 
 
      real*4 xs, ys, cpux, cpuy
 
****  Get the current pen position, in plotter coordinates
 
      call mxmy(ix,iy)
 
*     Convert plotter coordinates to user coordinates
 
      xs = cpux(ix)
      ys = cpuy(iy)
 
****  Depending of dash type use line or lined.
 
      if( jdash.gt.1 ) then
          call lined( xs, ys, xp, yp )
      else
          call line( xs, ys, xp, yp )
      end if
 
***** Thats all
      return
      end
 
