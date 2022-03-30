CTITLE WET_PRESSURE
 
      subroutine wet_pressure(temp_C,rel_hum,wet_pp)

      implicit none 
 
*     J.L. Davis                   3:02 PM  MON., 13  JULY, 1987
*
*     Calculates the wet partial pressure in mbars iven the
*     temperature and relative humidity
 
 
*       rel_hum         - Relative humidity (0-1)
*   ,   temp_C          - Temperature in degrees C
*   ,   temp_K          - Temperature in Kelvins
*   ,   temp_0          - 0 degrees C in Kelvins
*   ,   wet_pp          - Partial pressure of H2O vapor, mbars
 
      real*4 rel_hum, temp_C, temp_K, temp_0, wet_pp
 
      data temp_0 / 273.15 /
 
***** Calculate the absolute temperature
      temp_K = temp_C + temp_0
 
***** Calculate the saturation pressure in mbars
* MOD TAH 880120: multiplied the saturated wet pressure by the rel_hum to
*     to get actual wet pressure.
 
      wet_pp = 6.11 * (temp_K / temp_0) ** (-5.3)
     .              *  exp(25.2 * temp_C / temp_K) * rel_hum
 
      END
 
