CTITLE GET_EMA_8
 
      subroutine get_ema_8(bak_array, value)

      implicit none 
c
c     Routine to return real*8 value for bak_array
c
c Variables
c ---------
c bak_array -- the record from bak_array (only the first element is
c     real*8
c value -- the real*8 value to be returned
c
      real*8 bak_array(1), value
c
c.... Just save the value
      value = bak_array(1)
c
      return
      end
 
