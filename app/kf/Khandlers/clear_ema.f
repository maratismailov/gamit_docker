CTITLE CLEAR_EMA
 
      subroutine clear_ema ( start, num )
 

      implicit none
 
*     Routine to clear i*4 num words of ema from start to start+num
 
*   start(1j)   - first word to clear
 
      integer*4 start(1)
 
*   i           - Loop counter
*   num         - Number of words to clear
 
      integer*4 i, num
 
 
***** Clear with a do loop
 
      do i = 1, num
          start(i) = 0
      end do
 
***** Thats all
      return
      end
 
