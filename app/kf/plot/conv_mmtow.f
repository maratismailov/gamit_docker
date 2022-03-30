CTITLE    ....................................................................
 
      subroutine conv_mmtow( xs, ys, size, length)
c
c     Routine to convert xs or ys in mm to a change in world coordiantes
c     Either xs or ys should be zero.
c
c Variables
c ---------
c xs, ys -- the vaules in mm to converted to world coodinates.
c     One of these values should be zero.
c     ===================================
c size  -- the size of device in mm in the direction to be converted
c length -- the length of either xs or ys in world coordinates
c
      real*4 xs, ys, size, length
 
c
c local variables
c ---------------
c xw,yw,zw -- one set of world coordiantes
c xw2,yw2,zw2 -- a second set of world coordinates
c lenv -- the length in virtual coordinates
c
      real*4 xw,yw,zw, xw2,yw2,zw2, lenv
 
c
c scracth common
c --------------
c
      common xw,yw,zw, xw2,yw2, zw2, lenv
 
c
c.... Get of origin of virtual coordinates in world coordinates
      call jvtow(.1,.1, xw,yw,zw)
c
c.... Now compute change for either xs, ys
*                          ! tic mark on x axis
      if( xs.eq.0 ) then
         lenv = ys/size
         call jvtow(.1, .1+lenv, xw2,yw2,zw2)
         length = yw2 - yw
*                            ! tic markt on y axis
      else
         lenv = xs/size
         call jvtow(.1+lenv,.1,  xw2, yw2,zw2)
         length = xw2 - xw
      end if
c
      return
      end
 
