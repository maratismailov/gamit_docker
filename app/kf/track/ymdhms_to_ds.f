CTITLE YMDHMS_TO_DS
 
      subroutine ymdhms_to_ds(ymdhm, sec, doy, sod)

      implicit none
 
*     Rouitne to convert year, day-of-year and
*     second-of-day to  yy mm dd hh mm ss
 
*   ymdhm(5)     - Full yy mm dd hh mm of date
 
      integer*4 ymdhm(5)
 
*   sec      - ss: Seconds of date 

*   sod         - second-of-day
*   doy(2)         - year, day-of-year
 
      integer*4  doy(2)
 
      real*8 sod, sec
      
* LOCAL VARIABLES
 
*   date(5)     - Full yy mm dd hh mm
 
      integer*4 date(5)
 
*   sectag      - Seconds tag for jd
*   jd          - date convert to juliane date
*   jan1_jd     - JD on January 1
 
      real*8 sectag, jd, jan1_jd

*     Get sod
      sod = ymdhm(4)*3600 + ymdhm(5)*60 + sec

*     seperate ymdhm to ymd and hms
 
      date(1) = ymdhm(1)
      date(2) = ymdhm(2)
      date(3) = ymdhm(3)
      date(4) = 0
      date(5) = 0
      sectag  = 0.d0

*     Get year
      doy(1) = ymdhm(1)  
 
*     Get julian date
      call ymdhms_to_jd( date, sectag, jd)
 
*     Now do Jan 1
      date(2) = 1
      date(3) = 1
      call ymdhms_to_jd( date, sectag, Jan1_jd)
 
*     Save the day of year
      doy(2)  = jd - Jan1_jd + 1
 
****  Thats all
      return
      end
