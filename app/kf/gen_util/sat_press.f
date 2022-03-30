CTITLE 'SAT_PRESS'
 
 
      real*8 function sat_press(temp_C)
 
      implicit none 

*     J.L. Davis                   3:42 PM  SAT., 20  JUNE, 1987
*
*     Calculate the saturation pressure of water vapor in mbar for
*     a given temperature.
 
*       temp_C          - Temperature in deg C
*   ,   temp_K          - Temperature in Kelvins
*   ,   temp_0          - Zero deg C in Kelvins
 
      real*8 temp_C, temp_K, temp_0
 
      data temp_0 / 273.15D0 /
 
***** Calculate the absolute temperature
      temp_K = temp_C + temp_0
 
***** Calculate the saturation pressure
      sat_press = 6.11 * (temp_K / temp_0) ** (-5.3)
     .                     * exp(25.2 * temp_C / temp_K)
 
      end
 
