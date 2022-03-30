
CTITLE glb_svs_CODES
 
      subroutine glb_svs_codes( parn, dim1, dim2, num, type, apr )

      implicit none  
 
*     Routine to build up the code numbers for the satellite
*     orbit parameters in the global
*     solution estimated values.  (See GLB_HEADER_DEF.FTNI for code
*     numbers.  NOTE:  The type passed should equal the first
*     type for this parameter type.  The compensation for counters
*     looping from 1 to num is done in this routine.

*     This routine uses the extended from of the codes (i.e., the
*     higher level bits (about the I*2 portion) are set.
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
 
*   dim1        - First dimension for parn array
*   dim2        - Second dimension dimension for parn array
*   i,j         - Loop counters
*   num         - Number of values from parn to use
*   parn(dim1,dim2,1) - Parameter numbers for this type of parameter.
*               - (If the value is non-zero, then this parameter was
*               - estimated and it will be recorded here)
*   type        - Code type for these parameters
*   apr  -- Set to 1 if apriori is to be saved, otherwize code for
*           estimated parameter is saved
 
 
      integer*4 dim1, dim2, i,j, num, parn(dim1,dim2), type, apr
 
***** Loop over the parameters and see which ones are not estimated
 
      do i = 1, num    
          do j = 1, dim1
              if( apr.eq.1 ) then   ! Save the apriori 
                  num_apr_codes = num_apr_codes + 1
                  apr_codes(num_apr_codes) = type + j*256 + i*65536
                  org_apr_codes(num_apr_codes) = type + j*256 + i*65536
              else
*                 Save the parameter estimate codes                  
                  if( parn(j,i).ne.0 ) then
                      glb_codes(parn(j,i)) = type + j*256 + i*65536
                  end if
              endif
          end do
      end do

***** Thats all
      return
      end
 
