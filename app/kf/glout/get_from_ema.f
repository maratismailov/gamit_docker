CTITLE GET_FROM EMA
 
      subroutine get_from_ema ( value, vector, iel )

      implicit none  
 
*     routine to get the iel'th element from ema VECTOR and save
*     in value.  Real*8 values only.
 
*   iel     - Elememt from vector to be gotten
 
      integer*4 iel
 
*   value   - Main memory value to be returned
*   vector(1)   - EMA vector from which value will be extracted
 
      real*8 value, vector(1)
 
 
      if( iel.ge.1 ) then
          value = vector(iel)
      end if
 
***** Thats all
      return
      end
 
