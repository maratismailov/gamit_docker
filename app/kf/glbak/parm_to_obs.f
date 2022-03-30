CTITLE PARM_TO_OBS
 
 
      real*8 function parm_to_obs( parn, sol_parm)

      implicit none  
 
*     Function to copy the PARN'th element of sol_parm to parm_to_obs
*     if PARN is greater than zero.  Otherwise parm_to_obs is set
*     to zero.
 
 
*   parn        - Parameter number to be copied
 
      integer*4 parn
 
*   sol_parm(1) - Vector containing the corrections to the
*               - parmeters
 
      real*8 sol_parm(1)
 
 
****  If parn is greater than zero then copy
      if( parn.gt.0 ) then
          parm_to_obs = sol_parm(parn)
      else
          parm_to_obs = 0.d0
      end if
 
****  Thats all
      return
      end
 
