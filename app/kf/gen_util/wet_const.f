CTITLE WET_CONST
 
      subroutine wet_const( in, wet_zen )
 
      implicit none 

 
*     Routine to return the values of user wet zenith delay term.
*     These values are stored in the VALUES block of the data file
*
*                                         11:10 PM SUN., 15 Feb., 1987
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Index to site in baseline
 
      integer*4 in
 
*   wet_zen(2)  - wet zenith delay and rate (ps and fs/sec)
 
      real*8 wet_zen(2)
 
***** Save the values
 
      wet_zen(1) = user_wet_con(site(in))
      wet_zen(2) = 0.d0
 
***** Thats all
      return
      end
 
