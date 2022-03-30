CTITLE GLB_APR_CODES
 
      subroutine glb_apr_codes( parn, dim1, dim2, num, type, used, all )

      implicit none  
 
*     Routine to build up the apriori code numbers for the global
*     solution apriori values.  (See GLB_HEADER_DEF.FTNI for code
*     numbers.  NOTE:  The type passed should equal the first
*     type for this parameter type.  The compensation for counters
*     looping from 1 to num is done in this routine.
*     For the coefficients bit 16 is set for out of phase, and
*     not set for in phase terms
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
 
*   dim1        - First dimension for parn array
*   dim2        - Second dimension dimension for parn array
*   i,j         - Loop counters
*   num         - Number of values from parn to use
*   parn(dim1,dim2,1) - Parameter numbers for this type of parameter.
*               - (If the value is zero, then this parameter was
*               - not estimated and it will be recorded here)
*   type        - Code type for these parameters
*   used(*)     - Bit mapped array which indicates if element
*                 used.  Only applied if num is negative of
*                 number to be checked.
*   k           - Actual counter for number when some elememts
*                 will be suppressed because they are not used.
*   all  -- Set to 1 if all values are to be output even if not estimated
*           Otherwise value only output if estimated (apriori velocity
*           will be output even if not estimated, atmosphere only if
*           estimated).
 
      integer*4 dim1, dim2, i,j, num, parn(dim1,dim2,1), type,
     .          used(*), k, all

      logical kbit 

***** Loop over the parameters and see which ones are not estimated
      k = 0 
      do i = 1, abs(num)
          if( num.gt.0 .or. (num.lt.0 .and. kbit(used,i))) then 
              k = k + 1
          end if
          do j = 1, dim1
 
*                                           ! Not estimated, so save
*             Save all apriori values
              if( num.gt.0 .or. (num.lt.0 .and. kbit(used,i))) then 

*                 See if we need to check if estimated
                  if( (all.eq.0 .and. parn(j,1,i).ne.0) .or.
     .                 all.eq.1 ) then  
                     num_apr_codes = num_apr_codes + 1
*                    Use k instead of do loop index i.
                     apr_codes(num_apr_codes) = type+j-1 + 256*k
                     org_apr_codes(num_apr_codes) = type+j-1 + 256*i
                  endif
              end if
          end do
      end do
 
***** Thats all
      return
      end
 
