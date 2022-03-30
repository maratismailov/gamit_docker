CTITLE WET_CSC
 
      subroutine wet_csc( in, wet_map )
 
      implicit none 

 
*     Routine to compute the Cosecant wet mapping function.
*                                 10:25 AM  MON., 23  FEB., 1987
 
      include '../includes/kalman_param.h'
      include '../includes/obs_data.h'
 
*   in      - Site index in baseline (1 or 2 )
 
      integer*4 in
 
*   wet_map(2)  - Delay and rate mapping function
*   sine        - Sine of elevation angle
 
      real*8 wet_map(2), sine
 
***** START, compute the values we need
 
      sine = sin( elev(in,1) )
 
      wet_map(1) = 1.d0/sine
      wet_map(2) = -1.d0/sine**2 * cos(elev(in,1)) *
     .             elev(in,2)*1000.d0
 
***** Thats all
      return
      end
 
