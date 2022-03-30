CTITLE    ...............................................................
 
      subroutine draw_tic(st_x,st_y, of_x,of_y, fin_x,fin_y)
c
c     Routine to draw a tic mark at absolute position st_x, st_y
c     to relative position of_x,of_y and then to move the cursor
c     relative position fin_x,fin_y
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
c st_x, st_y -- the start position (absolute)
c of_x, of_y -- the movement for tic (relative)
c fin_x, fin_y -- the position to leave the cursor in
c
      real*4 st_x, st_y, of_x, of_y, fin_x, fin_y
 
c
c.... Move the cusor to start position
      call j2mov(st_x*sign_x, st_y*sign_y)
      call s2mov(st_x*sign_x, st_y*sign_y)
c
c.... draw tic relative
      call jr2dr(of_x*sign_x, of_y*sign_y)
      call sr2mv(of_x*sign_x, of_y*sign_y)
c
c.... Now move to finial position
      call jr2mv(fin_x*sign_x, fin_y*sign_y)
      call sr2mv(fin_x*sign_x, fin_y*sign_y)
c
      return
      end
 
