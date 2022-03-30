      subroutine timeup(itb,tb,ihup,minup,secup)
C
C Add hour,min,sec to particular date
C    itb   : month,day,year
C    tb    : hour,min,sec
C    ihup  : hours to update
C    minup : minutes to update
C    secup : seconds to update
C
      integer jdb,jd,ih,min,julday,isec,ihup,minup
     .        ,idoy,iy,itb(3)
      real*8 fjdb,fract,sec,secup,tb(3)
C
       jdb=julday(itb(1),itb(2),itb(3))
       fjdb=dble(jdb)
       fjdb=fjdb+(tb(1)+tb(2)/60.d0+tb(3)/3600.d0)/24.d0
C      add time
       fjdb=fjdb+(dble(ihup)+dble(minup)/60.d0+secup/3600.d0)/24.d0
       jd=int(fjdb)
       call dayjul(jd,iy,idoy)
       itb(3)=iy
       call monday(idoy,itb(1),itb(2),itb(3))
       fract=(fjdb-dble(jd))*86400.d0
       call ds2hms(iy,idoy,fract,ih,min,sec)
C      write(6,*) 'ih,min,sec',ih,min,sec
       tb(1)=dble(ih)
       tb(2)=dble(min)
       sec=sec+0.1d0
       isec=int(sec)
       tb(3)=dble(isec)

      return
      end


