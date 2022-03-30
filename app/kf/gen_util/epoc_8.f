CTITLE EPOC_8
 
      subroutine epoc_8( imon, iday, iyr, ihr, imin, jd)

      implicit none 
 
*     Rooutine to return the date corresponding to Julian date 
*     This version uses JD_to_YMDHMS.
 
* Variables
*   imon    - Month
*   iday    - Day of month
*   iyr     - Year (either 19xx or xx and assumes 1900)
*   ihr     - Hour,
*   imin    - Min
 
      integer*4 imon, iday, iyr, ihr, imin
 
*   jd      - Julian date and fracttion corresponding to the
*           - above
 
      real*8 jd
 
* Local Variables
 
*   date(5) - 5 elments of year, month, day, hour, min
 
      integer*4 date(5)
 
*   sectag  - Seconds tag. Set to zero
 
      real*8 sectag
 
****  Setup date and use JD_to_YMDHMS
      call jd_to_ymdhms( jd, date, sectag )
      iyr = date(1) 
      imon = date(2)
      iday = date(3)
      ihr  = date(4)
      imin = date(5)
 
****  Thats all
      return
      end
 
 
