CTITLE    .................................................................
 
      subroutine axis_setup( ch_x, ch_y, size_mm, tic_len)
c
c     Routine to set for drawing an axis.  The line type is set
c     and the view and scales are checked.  The length of the
c     ticmarks in world coordinates is also computed.
c
c Variables
c ---------
c ch_x,ch_y -- size of characters in the x and y directions for
c     the tic mark (one should be set to zero)
c size_mm -- size of page in mm
c tic_len -- length of tic marks in world coordinates
 
      real*4 ch_x,ch_y, size_mm, tic_len
c
c
c.... Set scales and view
      call set_view
      call set_scale
c
c.... set line type
C     call jlstl(1)
c
c.... Get the size to tic mark, should be half the size of charcter
      call conv_mmtow(ch_x/2,ch_y/2,size_mm,tic_len)
c
c
      return
      end
 
