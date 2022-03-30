CTITLE GLB_PAR_COUNT
      subroutine glb_par_count( apr, mar, np, nm, parn)
 
      implicit none 

 
*     Routine to increment the number of global parameters and
*     markov parameters, and to save the apriori variances and the
*     markov step variances.
* MOD TAH 190528: Added fearue where mar process noise may be negative
*     to decouple rate from value.  (Used to separate UT1 from LOD integrated).
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   nm      - number of markov parameters
*   np      - number of global parameters
*   parn    - parameter number for this pararmeter
 
      integer*4 nm, np, parn
 
*   apr     - apriori sigma for this parameter
*   mar     - Markov sigma for this parameter.  Either of these
*           - two values being non-zero will cause a parameter
*           - to be estimated.
 
      real*4 apr, mar
 
 
***** See if apr or mar are non-zero
 
*                                         ! Yes, we will estimate this
      if( apr.ne.0 .or. mar.ne.0 ) then
*                                         ! parameter
          np = np + 1
          call check_glb_max('global parameters', np, max_glb_parn)
 
*         Check for an origin fix. i.e. we give negative sigma to fix
*         quantity.  Now set sigma to zero
          if( apr.gt.0 ) then
              cov_apr(np) = apr**2
          else
              cov_apr(np) = 0.d0
          end if
 
          parn        = np
 
*         Now see if markov as well
          if( mar.ne.0 ) then
              nm = nm + 1
              call check_glb_max('markov parameters', nm, max_glb_mar)
              cov_mar(nm) = abs(mar)
              ind_mar(nm) = np
          end if
      end if
 
***** Thats all
      return
      end
 
