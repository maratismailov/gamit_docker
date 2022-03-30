CTITLE DS_TO_YMDHMS
 
      subroutine ds_to_ymdhms(doy, sod, ymdhm, sec )

      implicit none
 
*     Rouitne to convert year, day-of-year and
*     second-of-day to  yy mm dd hh mm ss
 
*   ymdhm(5)     - Full yy mm dd hh mm of date
 
      integer*4 ymdhm(5)
 
*   sec      - ss: Seconds of date 

*   sod         - second-of-day
*   doy(2)         - year, day-of-year
 
      integer*4  doy(2),i 
 
      real*8 sod, sec
      
* LOCAL VARIABLES
 
*   date(5)     - Full yy mm dd hh mm
 
      integer*4 date(5)
 
*   sectag      - Seconds tag for jd
*   jd          - date convert to juliane date
*   jan1_jd     - JD on January 1
 
      real*8 sectag, jd, day1_jd, day0_jd

*     Get year
      ymdhm(1) = doy(1) 
 
*     Get month
      call yds_to_jd( doy(1),doy(2),sod, jd)
      
      do i=1,12
        date(1) = doy(1)
        date(2) = i
        date(3) = 0
        date(4) = 0
        date(5) = 0
        sectag = 0.0
        call ymdhms_to_jd( date, sectag, day1_jd)
        if ( day1_jd.gt.jd) then
          ymdhm(2) = i - 1
          ymdhm(3) = int(jd - day0_jd) + 1
        end if
        day0_jd = day1_jd
      end do

      if (ymdhm(2).eq.0) then
          ymdhm(2) = 12
          ymdhm(3) = int(jd - day0_jd) + 1
      end if
      
      
      ymdhm(4) = int(sod/3600)
      ymdhm(5) = int ( (sod - ymdhm(4)*3600)/60)
      
      sec = sod - ymdhm(4)*3600 - ymdhm(5)* 60
      
         

 
****  Thats all

       
       return
      end
