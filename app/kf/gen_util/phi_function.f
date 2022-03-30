CTITLE PHI_FUNCTION
 
      subroutine phi_function(latitude,ellip_hgt,f)

      implicit none 
 
*     J.L. Davis                   3:21 PM  MON., 13  JULY, 1987
*
*     Routine to compute Saastamoinen's function for the decrease
*     in gravity due to height/latitude
 
 
*       ellip_hgt           - Ellipsoidal height, meters
*   ,   latitude            - Latitude, radians
*   ,   f                   - Function value
 
      real*4 ellip_hgt, latitude, f
 
      f = 1 - 0.00266 *cos(2 * latitude)
     .                     - 0.00028D-03 * ellip_hgt
 
      END
 
