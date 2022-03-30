CTITLE SOURCE_OC
 
      subroutine source_OC ( comp, source, sol_obs )

      implicit none  
 
*     Routine to compute source o minus c.  Initially only accounts
*     for position and velocity.  Later should include any affects
*     of changes in nutation series if these values are not estimated.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

      include '../includes/glb_hdr_def.h'
 
*   comp        - Component (1=RA, 2=Dec)
*   source      - global source number
 
      integer*4 comp, source
 
*   sol_obs     - Parameter to be updated
 
      real*8 sol_obs

*   deriv(2)    - Derivative of RA or Dec with respect to nutation angles
*   dnut(2)     - Differences in nutation angles between SOLVK value and
*                 current solution.

      real*8 deriv(2), dnut(2)
 
 
***** Remove the C part of the estimates
      if( comp.le.2 ) then
          sol_obs = sol_obs - apr_val_source(comp,1,source) -
     .            apr_val_source(comp,2,source)*(gepoch_expt-
     .                                source_epoch(source))/365.25d0
          if( .not.nut_ang_est ) then
              call nut_ang_part(source, deriv, 
     .              apr_val_source(1,1,source), comp, gepoch_expt)
              dnut(1) = vnut_ang_apr(1) - cnut_ang_apr(1)
              dnut(2) = vnut_ang_apr(2) - cnut_ang_apr(2)
              sol_obs = sol_obs - dnut(1)*deriv(1) - dnut(2)*deriv(2)
          end if
      else
          sol_obs = sol_obs - apr_val_source(comp-2,2,source)
      end if
 
***** Thats all
      return
      end
 
