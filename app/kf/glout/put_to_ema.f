CTITLE PUT_TO_EMA
 
      subroutine put_to_ema( ema_value, value, iel )

      implicit none 
 
*     Routine to move value to EMA
 
*   iel     - Element number to be moved
 
      integer*4 iel
 
*   ema_value(1)    - Ema arrya
*   value           - Main memory value
 
 
      real*8 ema_value(1), value
 
 
      if( iel.gt.0 ) then
      end if
 
      return
      end
 
