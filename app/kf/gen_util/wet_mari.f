CTITLE WET_MARI
 
      subroutine wet_mari( in, wet_map, dry_zen, wet_zen)

      implicit none 
 
 
*     Routine to compute the wet marini mapping function.
*                                  1:06 PM  MON., 16  FEB., 1987
*
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Site index either 1 or 2
 
      integer*4 in
 
*   dry_zen(2)  - Dry zenith delay term (part of A in Marini formula)
*   wet_map(2)  - Delay and rate mapping function. (ps/s and (fs/s)/ps)
*   wet_zen(2)  - Wet zenith delay term (part of A in Marini formula)
 
 
      real*8 dry_zen(2), wet_map(2), wet_zen(2)
 
***** For the moment use, DRY_MARI mapping function
 
      call dry_mari( in, wet_map, dry_zen, wet_zen )
 
***** Thats all
      return
      end
 
