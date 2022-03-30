CTITLE WET_CHAO
 
      subroutine wet_chao( in, wet_map )
 
      implicit none 

 
*     Routine to compute the CHAO wet mapping function.
*                                  2:55 PM  MON., 16  FEB., 1987
 
      include '../includes/kalman_param.h'
      include '../includes/obs_data.h'
 
*   in      - Site index in baseline (1 or 2 )
 
      integer*4 in
 
*   A,B     - Constants in CHAO formula
*   beta    - Part of mapping funtion
*   wet_map(2)  - Delay and rate mapping function
*   gamma   - Part of mapping function
 
      real*8 A,B, beta, wet_map(2), gamma
 
***** START, compute the values we need
 
      A = 0.00035d0
      B = 0.017  d0
 
      gamma = A   /( tan(elev(in,1)) + B)
      beta  = 1.d0/( sin(elev(in,1)) + gamma )
 
      wet_map(1) = beta
      wet_map(2) = -beta**2* ( cos(elev(in,1)) -
     .             A/( tan(elev(in,1)) + B)**2/ cos(elev(in,1))**2) *
     .             elev(in,2)*1000.d0
 
***** Thats all
      return
      end
 
