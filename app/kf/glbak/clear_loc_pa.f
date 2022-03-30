CTITLE CLEAR_LOC_PARN
 
      subroutine clear_loc_parn( loc_parn, num )

      implicit none  
 
*     Routine to clear the local parmeter number array
 
*   i           - Loop counter
*   loc_parn(1) - Array which will contain parameter numbers
*   num         - Number of values to clear
 
      integer*4 i, loc_parn(1), num
 
***** Loop and clear values
 
      do i = 1, num
          loc_parn(i) = 0
      end do
 
***** Thats all
      return
      end
 
