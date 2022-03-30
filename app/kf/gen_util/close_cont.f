ctitle
 
      subroutine close_cont( s1,s2,l, save_site, nstep, residual,
     .            variance, use_close, num_close, nepoch_obs,
     .            sum_close, var_close )

      implicit none 
 
*-------------------------------------------------------------------
*     Routine to check the list of the residuals and see if sites
*     S1 and S2 are present.  If they are the residual is added
*     to the closure sum (with the appropriate sign), the variance
*     is incrememted, the number of values in the closure (num_close)
*     is incremented and the baseline number is saved.
*
*     T. Herring                  OCT 8, 1986
*------------------------------------------------------------------
 
 
*   i       - A loop counter for looking for the sites
*   l       - The index number for delay or rates
*   nepoch_obs  - number of data at this epoch
*   nstep       - nunber of data types per baseline
*   num_close   - Number of baselines in closure
*   S1, S2  - the site numbers to be checked
*   save_site(2,1)  - the list of baselines at this epoch
*   use_close(3)    - baseline numbers in this closure
 
      integer*4 i, l, nepoch_obs, nstep, num_close, S1, S2,
     .    save_site(2,1), use_close(3)
 
*   residual(1) - the residuals at this epoch
*   sign    - Sign to be used in the closure
*   sum_close   - sum of the residuals in the triplet
*   var_close   - variance of the closure
*   variance(1) - the variances of the residuals
 
      real*8 residual(1), sign, sum_close, var_close, variance(1)
 
*   finished    - Indicates that baseline found
 
      logical finished
 
***** START, See if we can find this baseline (S1/S2)
 
      i = l - nstep
      finished = .false.
*                                 ! loop until we find baseline or run
      do while ( .not. finished )
*                              ! out of data
 
          i = i + nstep
          sign = 0.d0
          if ( S1.eq.save_site(1,i) .and. S2.eq.save_site(2,i) ) then
*                          ! Forward baseline
              sign = 1.d0
          end if
 
          if ( S1.eq.save_site(2,i) .and. S2.eq.save_site(1,i) ) then
*                          ! Reverse baseline
              sign = -1.d0
          end if
 
*         See if we found baseline
*                                  ! Found, sum in values
          if( sign.ne.0.d0 ) then
              finished = .true.
              num_close = num_close + 1
              use_close(num_close) = i
              sum_close = sum_close + sign*residual(i)
              var_close = var_close + variance(i)
          end if
 
*         See if we have run out of data to check
          if( i+nstep.gt.nepoch_obs ) then
              finished = .true.
          end if
 
      end do
 
***** Thats all
      return
      end
 
c....................................................................
