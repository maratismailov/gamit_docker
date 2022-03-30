ctitle
 
      Subroutine clear_close( sum_close, var_close, num_close,
     .                        use_close )

      implicit none 
 
 
*--------------------------------------------------------------------
*     Routine to clear the values used in computing the closures
*     around triplets of baselines.
*
*     T,Herring                   OCT 8, 1986
*--------------------------------------------------------------------
 
*   i           - Loop counter
*   num_close   - number of values in closure computation
*   use_close(3)- baseline numbers of the data in the closure
 
      integer*4 i, num_close, use_close(3)
 
*   sum_close   - Sum of the delays or rates around figure
*   var_close   - variance of the sum
 
      real*8 sum_close, var_close
 
***** SET the values, all to zero
 
      num_close = 0
      do i = 1,3
          use_close(i) = 0
      end do
 
      sum_close = 0.d0
      var_close = 0.d0
 
***** Thats all
      return
      end
 
*-------------------------------------------------------------------
