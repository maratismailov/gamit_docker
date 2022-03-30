CTITLE
 
      subroutine met_seasonal( T_C, P, rh, tbias, epoch, lat, hgt)
 
      implicit none 

*     Routine to return an estimate of temperature, pressure, and
*     relative humidity based on time, latitude, and height
 
      include '../includes/const_param.h'
 
* PASSED VARIABLES
 
*   t_c     - Temperature in C
*   P       - Pressure in mbar
*   rh      - Relative humidity in 0-1
*   lat     - Latitude of site in radians
*   hgt     - Height of site in m.
*   tbias   - Bias in surface temperature (Should not be added to
*             t_c returned from this routine).
*   t_k     - Temperature Kelvin
 
      real*4 t_c, P, rh, lat, hgt, tbias, t_k
 
*   epoch   - JD of current measurement.
*   dt      - Argument for seasonal term
 
      real*8 epoch, dt
 
*****
*     Get and estimate of temperature first.  Get seasonal argument
      dt = mod ((epoch-dj2000)/365.25,1.d0)*2*pi

      t_c   = ( -20.5 + 48.4*cos(lat)  - 3.1*hgt/1000.d0 ) +
     .         (-14.3 + 3.3*hgt/1000.d0)*sin(lat)*cos(dt)  +
     .         ( -4.7 + 1.1*hgt/1000.d0)*sin(lat)*sin(dt)  
 
      tbias = (  -3.6 - 0.3*cos(lat)  + 2.3*hgt/1000.d0 ) +
     .         ( -1.7 - 1.0*hgt/1000.d0)*sin(lat)*cos(dt)  +
     .         ( -0.2 + 0.8*hgt/1000.d0)*sin(lat)*sin(dt)  
 
*     Now based on standard lapse rate compute the pressure
      t_k = t_c + 273.15
      P = 1013.25d0 * (t_k/(t_k + 6.5d-3*hgt))**5.26d0
 
*     Set the relative humidity
      rh = 0.5
 
****  Thata all
      return
      end
 
