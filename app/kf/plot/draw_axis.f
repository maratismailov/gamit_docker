CTITLE    .................................................................
 
      subroutine draw_axis(st_x, st_y, en_x, en_y)
c
c     Routine to draw and axis from st_x,st_y to en_x,en_y
c
c Include files
c -------------
*                        ! the parameter file
      include 'plot_param.h'
c
*                        ! the common block
      include 'plot_com.h'
c
c Variables
c ---------
c st_x,st_y the start coordinates in world cooridinates
c en_x,en_y the end cooridinates in world coordiantes
c
      real*4 st_x, st_y, en_x, en_y
 
c
c.... Move to start
      call j2mov(st_x*sign_x, st_y*sign_y)
c
      call j2drw(en_x*sign_x, en_y*sign_y)
c
      return
      end
 
