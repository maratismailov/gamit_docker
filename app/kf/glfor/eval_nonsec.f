CTITLE EVAL_NONSEC

      subroutine eval_nonsec(is, epoch, num_nonsec, param_nonsec,
     .                       apr_val_nonsec, dsol, tran_est )

      implicit none 

*     Routine to evaluate the value of the non-secular contributions
*     to the apriori station position
* NOD TAH 031024: Passed tran_est into the routine.  If Bit 16 is
*     set, then log coefficients have been estimated anf therefore
*     should not be added to the station coordinates.  

      include '../includes/const_param.h'

* PASSED VARIABLES
* is  -- Site number
* num_nonsec -- Number of non-secular terms
* param_nonsec(2,num_nonsec) -- Parameters for non-secular terms
* tran_est  -- If Bit 16 is set then log terms are estimated and
*     therefore should not be added to coordinates.  Used in site_OC
*     so that we do not doubly correct the coordinates for the log term

      integer*4 is, num_nonsec, param_nonsec(2,num_nonsec), 
     .          tran_est

* epoch    -- Time at which postion is required
* apr_val_nonsec(8,num_nonsec)
* dsol(3)  -- The change in XYZ at this time

      real*8 epoch, apr_val_nonsec(8,num_nonsec), dsol(3)

* LOCAL VARIABLES
* i,j -- Loop variables

      integer*4 i, j

* dt  -- Change in time in years between solution epoch and offset
*        rate term
* arg -- Argument for the periodic, exponential or logarithm terms

      real*8 dt, arg 

* fatal -- String with any fatal error message written into it

      character*64 fatal

* kbit -- Check if bit is set
      logical kbit

****  Loop over the non-secular terms seeing which ones apply to this
*     site.
      do j = 1,3
         dsol(j) = 0.d0
      end do

      do i = 1, num_nonsec
         if( is.eq.param_nonsec(1,i) ) then
*           OK, this term applies.  Based on type see what we should
*           do.
            if ( param_nonsec(2,i).eq.1 ) then
*              offset and rate change.  See if we are past the time
*              of the change
               if( epoch.ge.apr_val_nonsec(1,i) ) then
*                  We are past the change point.  Compute the change
*                  in position
                   dt = (epoch-apr_val_nonsec(1,i))/365.25d0
                   do j = 1,3
                      dsol(j) = dsol(j) + apr_val_nonsec(2*j+1,i) +
     .                                    apr_val_nonsec(2*j+2,i)*dt
                   end do
               end if
            else if( param_nonsec(2,i).eq.2 ) then
*              Periodic term.  Get the argument and evaluate cose and sine
               arg = 2*pi*(epoch-apr_val_nonsec(1,i))/
     .                           apr_val_nonsec(2,i)
               do j = 1,3
                  dsol(j) = dsol(j) + cos(arg)*apr_val_nonsec(2*j+1,i)+
     .                                sin(arg)*apr_val_nonsec(2*j+2,i)
               end do
            else if( param_nonsec(2,i).eq.3 ) then
*              Exponential decay term
               if( epoch.ge. apr_val_nonsec(1,i)) then
                  arg = -(epoch-apr_val_nonsec(1,i))/apr_val_nonsec(2,i)
                  do j = 1,3
                     dsol(j) = dsol(j) + 
     .                        (1-exp(arg))*apr_val_nonsec(j+2,i)
                  end do
               end if
            else if( param_nonsec(2,i).eq.4 ) then
*              Logarithmic  term
               if( epoch.ge.apr_val_nonsec(1,i) .and.
     .             .not.kbit(tran_est,16) ) then
* MOD TAH 030124: Changed argument to match enfit definition
* MOD TAH 030610: Changed definition again to match tsview and enfit
* MOD TAH 031024: Added check on tran_est to see if log terms have 
*                 been estimated.
                  arg = 1+(epoch-apr_val_nonsec(1,i))/
     .                     apr_val_nonsec(2,i)
                  do j = 1,3
                     dsol(j) = dsol(j) + log(arg)*apr_val_nonsec(j+2,i)
                  end do
               end if
            else 
*              Strange term that should not exit.
               write(fatal,220) is, i, param_nonsec(2,i)
 220           format('Invalid non-secular term for site # ',i4,
     .                ' term ',i5,' Saved type ',i5)
               call report_stat('FATAL','GLOBK','eval_nonsec',' ',
     .                           fatal,0)
            end if
         end if
      end do

***** Thats all
      return
      end





