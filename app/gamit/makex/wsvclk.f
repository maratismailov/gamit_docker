      subroutine wsvclk( usvclk,nprn,iwkn,sow,bc_clock )
c
c     Write a J file record from a the SV clock information from
c     a RINEX Navigation message file

c     R.W. King January 1990, from code in MAKEJ
c

      character*53     afmt

      integer*4        usvclk,itflag,iwkn,iyr,idoy,ihr,min,iprn,ihnsec

      real*8           sow,sec,utcoff,bc_ephem(16), bc_clock(6)
     .               , xeaf0,xeaf1,xeaf2
c

      data afmt/'(i4,1x,i3,2i3,1x,f7.4,3x,i3,1x,f11.4,2x,i2,2x,3d16.8)'/
c

      xeaf0= bc_clock(2)
      xeaf1= bc_clock(3)
      xeaf2= bc_clock(4)
      itflag= 1
      call timcon(itflag,iwkn,sow,iyr,idoy,ihr,min,sec,utcoff)


c        Write a record of the J-file

       write(usvclk,afmt) iyr,idoy,ihr,min,sec,iwkn,sow,nprn
     .                  , xeaf0,xeaf1,xeaf2


cd    stop
      end
