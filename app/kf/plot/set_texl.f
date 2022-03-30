      subroutine set_texl(ofx,ofy)
c
c
c     Routine to set up the size and justification of low quality
c     text characters
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
c
*                         ! the plot_com commom
      include 'plot_com.h'
c
c Variables
c ---------
c ofx,ofy -- the offsets in world coordinates to the position of the
c     point to make the charcater centered.
c
      real*4 ofx,ofy
 
c
c.... The prodecures depend of the device
      call set_charsz(charsz_x, charsz_y )
c
c.... set character offsets to zero and then adjust values if need be
      ofx = 0.0
      ofy = 0.0
      call jjust( 0.33333, 0.25)
c
      if( device_id.eq.1 .or. device_id.eq.2 ) then
*                                 ! HP 2648, HP2623
c
          call jjust( 0.5, 0.4)
c
         return
      end if
c
      return
      end
 
