      Subroutine pjdhms(aa,iyr,idoy,imonth,iday,ihr,imin)      

c     Given PEP Julian day PJD,
c     returns year,month,day,hr,min  
c     Used initially (exclusively? by ARC for eclipse printout

      integer*4 jd,iyr,idoy,ihr,imin,iday,imonth

      real*8 aa,bb,hr,xmin

      bb = aa
c      bb = aa - 0.5d0
      jd = idint(bb)
      call dayjul(jd,iyr,idoy)
      hr =  24.d0 * (bb - jd)
      xmin = idint(hr)
      xmin = 60.d0 * ( hr - xmin)
      ihr = idint(hr)
      imin = idint(xmin)
      call monday(idoy,imonth,iday,iyr)
c      write(*,600) iyr,idoy,idint(hr),idint(xmin)
c 600  format("year|day|hour|min>>> ",4i8)

      return
      end

