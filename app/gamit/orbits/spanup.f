      Subroutine SPANUP ( it0, t0, span )

c     Update a time in the (clumsy, USA) form of mon/day/year  hr/min/sec,
c     used in some of ARC, by adding 'span' in seconds

c     Rewritten by R. King 23 January 2001 to use library routines (motivation
c     is to avoid a subtle bug with the g77 compiler)

      implicit none

      integer*4 it0(3),year,month,day,doy,jd,julday,hr,min

      real*8 t0(3),span,sod,sec

      year = it0(3)
      month = it0(1)
      day = it0(2) 
      sod = t0(1)*3600.d0 + t0(2)*60.d0 + t0(3)
        
      jd = julday(month,day,year)
      call timinc(jd,sod,span)
      call dayjul( jd,year,doy)
      call monday(doy,month,day,year)  
      call ds2hms(year,doy,sod,hr,min,sec)
              
      it0(1) = month
      it0(2) = day
      it0(3) = year  
    
      t0(1) = dble(hr)
      t0(2) = dble(min)
      t0(3) = sec

      return
      end
