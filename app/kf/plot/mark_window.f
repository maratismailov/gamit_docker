CTITLE MARK_WINDOW
 
      subroutine mark_window

      implicit none
 
 
*     Routine to mark the window used for polynomial fitting
*     This is done by drawing a box around the window.
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   i,j     - Loop counters
 
      integer*4 i,j

*   xp(5), yp(5)    - Four cornors of the box
 
      real*4 xp(5), yp(5)
 
***** Set up the box to be draw
 
      xp(1) = poly_window(1)
      yp(1) = poly_window(3)
 
      xp(2) = poly_window(2)
      yp(2) = poly_window(3)
 
      xp(3) = poly_window(2)
      yp(3) = poly_window(4)
 
      xp(4) = poly_window(1)
      yp(4) = poly_window(4)
 
      xp(5) = xp(1)
      yp(5) = yp(1)
 
***** Set the line_type
      if( line_type(1).ne.0 ) then
          call jlstl(line_type)
      else
          call jlstl(1)
      end if
 
***** Draw box and flush
      call j2ply( 5, xp, yp)
      call jmcur
 
***** Thats all
      return
      end
 
