CTITLE FJLDY_8
 
 
      real*8 function fjldy_8( imth, iday, iyr )

      implicit none 
 
*     Routine to compute the julian date.  Here we simply use
*     the more up to date YMDHMS_to_JD routine.
 
* VARIABLES
 
*   iyr     - Year (may be 19xx or just xx and assumed to be 1900)
*   imth    - Month
*   iday    - Day of month
 
      integer*4 iyr, imth, iday
 
* LOCAL VARIABLES
 
*   date(5) - Year, month, day, hour, min
      integer*4 date(5)
 
*   sectag  - Seconds tag. (NOT used here)
*   jd      - Estimated julian date
 
      real*8 sectag, jd
 
****  Copy the values over and call ymdhms_to_jd
      date(1) = iyr
      date(2) = imth
      date(3) = iday
      date(4) = 0
      date(5) = 0
      sectag = 0.d0
 
      call ymdhms_to_jd( date, sectag, jd )
 
      fjldy_8 = jd
 
****  Thats all
      return
      end
 
