CTITLE GLB_GLB_CODES
 
      subroutine glb_glb_codes( parn, dim1, dim2, num, type, used )

      implicit none  
 
*     Routine to build up the code numbers for the global
*     solution estimated values.  (See GLB_HEADER_DEF.FTNI for code
*     numbers.  NOTE:  The type passed should equal the first
*     type for this parameter type.  The compensation for counters
*     looping from 1 to num is done in this routine.
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
 
*   dim1        - First dimension for parn array
*   dim2        - Second dimension dimension for parn array
*   i,j         - Loop counters
*   num         - Number of values from parn to use
*   parn(dim1,dim2,1) - Parameter numbers for this type of parameter.
*               - (If the value is non-zero, then this parameter was
*               - estimated and it will be recorded here)
*   range       - Range of do loop if num is passed as zero
*   type        - Code type for these parameters
*   used(*)     - Bit set if element really used. 
*   k           - Actual number when some values are suppressed.
 
      integer*4 dim1, dim2, i,j, num, parn(dim1,dim2,1), range, type,
     .          used(*), k

      logical kbit
 
***** Loop over the parameters and see which ones are not estimated
 
*                             ! Local verus Global tides to be treated
      if( num.eq.0 ) then
          range = 1
      else
          range = abs(num)
      end if

      k = 0
 
      do i = 1, range
          if( num.gt.0 .or. (num.lt.0 .and. kbit(used,i))) then
              k = k + 1
          end if
          do j = 1, dim1
 
*                                           ! Not estimated, so save
              if( parn(j,1,i).ne.0 ) then
 
*                                         ! Do not put upper byte in
                  if( num.eq.0 ) then
                      glb_codes(parn(j,1,i)) = type+j-1
                  else if( num.gt.0 .or.
     .                    (num.lt.0 .and. kbit(used,i))) then
*                     use k instead of do loop index i
                      glb_codes(parn(j,1,i)) = type+j-1 + 256*k
                  end if
              end if
          end do
      end do
 
***** Thats all
      return
      end
 
