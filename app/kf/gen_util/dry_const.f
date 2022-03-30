CTITLE DRY_CONST
 
      subroutine dry_const( in, dry_zen )
 
      implicit none 

 
*     Routine to return the values of user dry zenith delay term.
*     These values are stored in the VALUES block of the data file
*
*                                         11:10 PM SUN., 15 Feb., 1987
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Index to site in baseline
 
      integer*4 in
 
*   dry_zen(2)  - Dry zenith delay and rate (ps and fs/sec)
 
      real*8 dry_zen(2)
 
***** Save the values
 
      dry_zen(1) = user_dry_con(site(in))
      dry_zen(2) = 0.d0
 
***** Thats all
      return
      end
 
