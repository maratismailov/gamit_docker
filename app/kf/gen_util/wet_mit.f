
CTITLE wet_mit
 
      subroutine wet_mit( in, wet_map )

      implicit none 
 
 
*     Routine to compute the MIT  wet mapping function., with
*     constant coefficients
*                              May 22, 1990.
 
      include '../includes/kalman_param.h'
      include '../includes/obs_data.h'
 
*   in      - Site index in baseline (1 or 2 )
 
      integer*4 in
 
*   A,B,C   - Constants in MIT formula
*   beta    - Part of mapping funtion
*   wet_map(2)  - Delay and rate mapping function
*   gamma   - Part of mapping function
 
      real*8 A,B,C, beta, wet_map(2), gamma, sine, cose, topcon 
 
***** START, compute the values we need
*     Coefficients based on Raytrace at Albany and Fairbanks.
*     Preliminary values.
 
      A = 0.0005741d0
      B = 0.0015470d0
      C = 0.0488185d0
 
      beta  = B/( sin(elev(in,1)) + C )
      gamma = A/( sin(elev(in,1)) + beta)
      sine  = sin( elev(in,1) )
      cose  = cos( elev(in,1) )

      topcon = (1.d0 + A/(1.d0 + B/(1.d0 + C)))
 
      wet_map(1) = topcon / ( sine + gamma )

      wet_map(2) = -topcon / ( sine + gamma )**2 *
     .            ( cose - A/( sine + beta)**2 * cose *
     .            ( 1.d0 - B/( sine + C )**2 ) ) *
*                                              !Units (fs/s)/s
     .            elev(in,2) * 1000.d0
 
 
***** Thats all
      return
      end
 
 
