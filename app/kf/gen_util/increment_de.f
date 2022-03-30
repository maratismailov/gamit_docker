ctitle
 
      subroutine increment_deletes( kunw, delete_count )
 

      implicit none 
 
*     Routine to increment the delete_count based on KUNW flag.
*     KUNW may have multiple bits turned on, indicating all of
*     problems associated with a given observation.  The
*     delete_count is only incremented for the most signficant
*     bit turned on (starting from bit 1.)
 
*   flag        - Start as a copy of kunw and is divided by
*               - two iteratively until a set bit is found.
 
*   kunw        - The Kalman filter unweighted flag, compiled
*               - after all calibrations have been made.
*   delete_count(1) - the counter which will keep track of
*               - the main reasons data are unweighted.
*               - (See OBS_VALUE.FTNI for corresondence between
*               - values in delete_count and meanings)
 
*   set_bit     - Indicate number of bit which is set (first
*               - bit, starting from Bit 1)
 
      integer*4 flag, kunw, delete_count(1), set_bit
 
*   good        - Indicates that set bit has not yet been found
 
      logical good
 
***** If kunw is zero, then just return, we dont need to count
*     it
 
*                                 ! Do not need to count
      if( kunw.eq.0 ) return
 
*                         ! copy kunw so we can work on it
      flag    = kunw
*                         ! set good at first
      good    = .true.
*                         ! start with bit one
      set_bit = 1
 
*                         ! Loop until set bit is found
      do while ( good )
 
*                                     ! bit one is set, we found most
          if( mod(flag,2).ne.0 ) then
*                                     ! significant bit.
              delete_count(set_bit) = delete_count(set_bit) + 1
              good = .false.
*                                     ! divide flag by two and try again
          else
              flag = flag/2
*                                     ! increment delete count index
              set_bit = set_bit + 1
          end if
*                         ! looping to find set bit.
      end do
 
***** Thats all
      return
      end
 
