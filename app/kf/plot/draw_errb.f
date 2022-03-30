CTITLE    ..................................................................
 
      subroutine draw_errb(x_data, y_data)
c
c     Subroutine to draw error bars.  Only one type of
c     error bar is supported at this time. This type is a straight
c     line which jumps over the character plotted.
* MOD TAH 910312: added errb_scale to allow rescaling of the error bars
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
c
*                         ! the common block
      include 'plot_com.h'
c
c Variables
c ---------
c x_data  -- the x data and its sigma (may be zero)
c y_data  -- the y data and its sigma (may be zero)
c
      real*4 x_data(1), y_data(2)
 
c
c
c Local variables
c ---------------
c xp, yp -- arrays used to hold the ema data in main memory
c gap_x, gap_y -- the size of the gap in the error bars so that
c     they dont pass through the plotted point
c
      real*4 xp(2), yp(2), gap_x, gap_y
 
      integer*4 i
 
c
c.... Get the ema data into main memory.  Note there is no need to
c     change the sign of these data because this sign change will be
c     handled in the draw tic routine.
      do i = 1,2
         xp(i) = x_data(i)
         yp(i) = y_data(i)

*        Multiply error bars by scaling factor. (Factotrs input with
*        errbar command.
         if( i.eq.2 ) then
            xp(i) = xp(i)*errb_scale(1)
            yp(i) = yp(i)*errb_scale(2)
         end if
      end do
c
c.... Make sure sigmas are positive
      xp(2) = abs(xp(2))
      yp(2) = abs(yp(2))
c
c.... Get the size of the gap in the error bar (if point plotted)
      if( point_type.ne.0 ) then
*                                                     ! points are 2mm*2mm
         call conv_mmtow(1.0, 0.0, x_size_mm, gap_x)
         call conv_mmtow(0.0, 1.0, y_size_mm, gap_y)
*               ! set gap to zero
      else
         gap_x = 0.0
         gap_y = 0.0
      end if
 
*     Increase the gap a bit
      gap_x = gap_x * 1.2
      gap_y = gap_y * 1.2
c
c.... Now draw the x and y error bars
*                              ! sigma for x data so draw
      if( xp(2).ne.0. ) then
         call draw_tic(xp(1)+gap_x,yp(1), xp(2)-gap_x,0.0, 0.0, 0.0)
         call draw_tic(xp(1)-gap_x,yp(1), -(xp(2)-gap_x),0.0, 0.0, 0.0)
      end if
c
*                              ! sigma for y data so draw
      if( yp(2).ne.0. ) then
         call draw_tic(xp(1),yp(1)+gap_y, 0.0,yp(2)-gap_y, 0.0, 0.0)
         call draw_tic(xp(1),yp(1)-gap_y, 0.0,-(yp(2)-gap_y), 0.0, 0.0)
      end if
c
c.... If the errbar type 2 is selected then put additional tick markes
c     at the top and bottom of the error bar
      if( errb_type.eq.2 ) then
*                                 ! put tick at left and right
         if( xp(2).ne.0 ) then
            call draw_tic(xp(1)+xp(2),yp(1)-gap_y,0.0,2*gap_y, 0.,0.)
            call draw_tic(xp(1)-xp(2),yp(1)-gap_y,0.0,2*gap_y, 0.,0.)
         end if
c
*                                 ! put tick at top and bottom
         if( yp(2).ne.0 ) then
            call draw_tic(xp(1)-gap_x,yp(1)+yp(2),2*gap_x,0.0, 0.,0.)
            call draw_tic(xp(1)-gap_x,yp(1)-yp(2),2*gap_x,0.0, 0.,0.)
         end if
c
*                                 ! type two error bars
      end if
c
      return
      end
 
