      real*8 function decyrs( iyr,idoy,sod )
                   
c     Get decimal years from year, day-of-year, sec-of-day
c     taking into account leap year differences (works same as 
c     /kf/utils function 'doy'.

c     R. King 10 May 2003

      implicit none                  

      integer*4 iyr,idoy,imonth,iday,julday,jd

      real*8 sod,xjd      
 
      call monday( idoy,imonth,iday,iyr )
      jd = julday( imonth,iday,iyr ) 
c     convert PEP JD to true JD for TAH routines
      xjd = dble(jd) + sod/86400.d0 -0.5d0
      call jd_to_decyrs( xjd,decyrs )
c      print *,'DECYRS iyr idoy sod imonth iday jd decyrs '
c     .               ,iyr,idoy,sod,imonth,iday,jd,decyrs

      return
      end



