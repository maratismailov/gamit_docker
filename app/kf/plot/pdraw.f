CTITLE PDRAW
 
      subroutine pdraw
 
 
*     Routine to draw a line representing the polynomial estimated from
*     the data.  The line is drawn only over the extent of of the window.
*     If a longer line is desired then, the window can be changed before
*     the line is drawn.
 
      include 'plot_param.h'
      include 'plot_com.h'

      integer*4 max_plyt

      parameter ( max_plyt = 200  )
 
*   i,j,k       - Loop counters
 
      integer*4 i,j
 
*   xp(max_plyt), yp(max_plyt)  - The lines to be drawn (max_plyt is
*               - given in the parameter file (must not exceed scratch
*               - common area).
 
      real*4 xp(max_plyt), yp(max_plyt)
 
*   start       - Start x value
*   step        - Step size to be used for line
*   value       - Value of polynomial at current x value
 
      real*8 start, step, value
 
***** Set up the arrays to be plotted
 
      start = poly_window(1)
      step  = (poly_window(2)-poly_window(1))/(max_plyt-1)
 
      do i = 1, max_plyt
          xp(i) = start + step*(i-1)
          value = 0
          do j = 0, poly_order
              if( j.eq.0 .and. xp(i).eq.0 ) then
                  value = poly_vec(1)
              else
                  value = value + poly_vec(j)*xp(i)**j
              end if
          end do
*                                 ! Account for any rotation of plot
          yp(i) = sign_y*value
          xp(i) = sign_x*xp(i)
      end do
 
***** Now set line type.  If none set use 1
      if( line_type(1).gt.0 ) then
          call jlstl(line_type)
      else
          call jlstl(1)
      end if
 
***** Draw line
      call gsclip(1)
      call j2ply( max_plyt, xp, yp )
*                     ! Flush buffer so that user sees line
      call gsclip(0)
      call jmcur
 
***** Thats all
      return
      end
 
