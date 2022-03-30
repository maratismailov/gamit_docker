CTITLE GET_EMA_4
 
      subroutine get_ema_4(bak_array, value, field)

      implicit none 
c
c     Routine to return real*4 values from bak array
c
c Variables
c ---------
c bak_array -- the bak array (may be only part of it)
c value  -- the value read from the bak_array
c field   -- the field information
c
      real*4 bak_array(1)
 
c
      real*8 value
 
c
      integer*4 field
 
c
c
c.... Get the value
      value = bak_array(field)
c
      return
      end
 
 
