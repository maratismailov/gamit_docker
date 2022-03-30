CTITLE GLB_COEFF_CODES
 
      subroutine glb_coeff_codes( parn, dim1, num, type )

      implicit none  
 
*     Routine to build up the code numbers for the global
*     coefficient type parameters.  (See GLB_HEADER_DEF.FTNI for code
*     numbers.  NOTE:  The type passed should equal the first
*     type for this parameter type.  The compensation for counters
*     looping from 1 to num is done in this routine.
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
 
*   dim1        - First dimension for parn array
*   i,j         - Loop counters
*   num         - Number of values from parn to use
*   parn(dim1,1) - Parameter numbers for this type of parameter.
*               - (If the value is zero, then this parameter was
*               - not estimated and it will be recorded here)
*   type        - Code type for these parameters
 
      integer*4 dim1, i,j, num, parn(dim1,1), type
 
***** Loop over coefficients and see which ones were estimated.
 
      do i = 1, num
          do j = 1, dim1
 
*                                         ! Not estimated, so save
              if( parn(j,i).ne.0 ) then
 
                  glb_codes(parn(j,i)) = type + 256*i
*                                                         ! Set in
                  call sbit(glb_codes(parn(j,i)),16,j-1)
*                                         ! or out phase bit
              end if
          end do
      end do
 
 
***** Thats all
      return
      end
 
