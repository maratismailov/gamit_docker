      subroutine wrxnav ( lu,nprn,iewkn,bc_ephem,bc_clock,subfr1 )
c
c     Write a RINEX 2 navigation file from input FICA data
c     R. W. King  February 1990  - from MAKEX routine WORBIT

      integer*4 lu,nprn,iewkn,iyear,imonth,iday,idoy,ihr,min
     .        , itflag,iyr2
      real*8   bc_ephem(16),bc_clock(6),subfr1(8)
     .       , sec,utcoff
     .       , xetoc,xeaf0,xeaf1,xeaf2,aode,xecrs,xedn,xem0,xecuc,xeecc
     .       , xecus,xeart,xetoe,xecic,xeom0,xecis,xei0,xecrc,xew,xeomd
     .       , xeidt,cflgl2,weekno,pflgl2,svaccr,svhlth,tgd,aodc
     .       , trans_sow,spare


c     note especially the units on bc_ephem(4) (13) and (14)
c     they are in meters, not kilometers!


c       Clock parameters

c             fica_f(13)
      xetoc = bc_clock(1)
c             fica_f(16)
      xeaf0 = bc_clock(2)
c             fica_f(15)
      xeaf1 = bc_clock(3)
c             fica_f(14)
      xeaf2 = bc_clock(4)


c       Broadcast ephemeris quantities

c          Orbit line 1

c             fica_f(26)  (--aode actually in sub-frame 2)
      aode  = subfr1(7 )
c             fica_f(27)
      xecrs = bc_ephem(14)
c             fica_f(28)
      xedn  = bc_ephem( 3)
c             fica_f(29)
      xem0  = bc_ephem( 2)

c          Orbit line 2

c             fica_f(30)
      xecuc = bc_ephem(11)
c             fica_f(31)
      xeecc = bc_ephem( 5)
c             fica_f(32)
      xecus = bc_ephem(12)
c             fica_f(33)
      xeart = bc_ephem( 4)

c          Orbit line 3

c             fica_f(34)
      xetoe = bc_ephem( 1)
c             fica_f(46)
      xecic = bc_ephem(15)
c             fica_f(47)
      xeom0 = bc_ephem( 8)
c             fica_f(48)
      xecis = bc_ephem(16)

c          Orbit line 4
c             fica_f(49)
      xei0  = bc_ephem( 6)
c             fica_f(50)
      xecrc = bc_ephem(13)
c             fica_f(51)
      xew   = bc_ephem(10)
c             fica_f(52)
      xeomd = bc_ephem( 9)

c          Orbit line 5

c             fica_f(54)
      xeidt = bc_ephem( 7)
c             fica_f(7)
      cflgl2= subfr1(1)
c             fica_f(6)
      weekno= iewkn
c             fica_f(11)
      pflgl2= subfr1(5)

c          Orbit line 6

c             fica_f(8)
      svaccr= subfr1(2)
c             fica_f(9)
      svhlth= subfr1(3)
c             fica_f(12)
      tgd   = subfr1(6)
c             fica_f(10)
      aodc  = subfr1(4)    

c          Orbit line 7 (RINEX 2 only)

c             fica_f(3)                 
      trans_sow = subfr1(8)
      spare = 0.d0 

                  

c        Convert clock epoch from week, sec-of-week to calender day, hr, min, sec
c        (still GPS time, not UTC)

      itflag = +4     
      call timcon ( itflag,iewkn,xetoc,iyear,idoy,ihr,min,sec,utcoff ) 
      call monday ( idoy,imonth,iday,iyear )


c       Write the epoch line with SV clock values
                
c     RINEX file wants a 2-digit year
      iyr2 = mod(iyear,100)
      write(lu,10) nprn,iyr2,imonth,iday,ihr,min,sec
     .            , xeaf0,xeaf1,xeaf2
   10 format(i2,5i3,f5.1,3d19.12)     


c        Write broadcast orbit lines  1 - 6

      write(lu,20) aode,xecrs,xedn,xem0
      write(lu,20) xecuc,xeecc,xecus,xeart
      write(lu,20) xetoe,xecic,xeom0,xecis
      write(lu,20) xei0,xecrc,xew,xeomd
      write(lu,20) xeidt,cflgl2, weekno,pflgl2
      write(lu,20) svaccr,svhlth,tgd,aodc 
      write(lu,20) trans_sow,spare,spare,spare

   20 format(3x,4d19.12)


c      bc_ephem( 1)  = fica_f(34)
c      bc_ephem( 2)  = fica_f(29)
c      bc_ephem( 3)  = fica_f(28)
c      bc_ephem( 4)  = fica_f(33)
c      bc_ephem( 5)  = fica_f(31)
c      bc_ephem( 6)  = fica_f(49)
c      bc_ephem( 7)  = fica_f(54)
c      bc_ephem( 8)  = fica_f(47)
c      bc_ephem( 9)  = fica_f(52)
c      bc_ephem(10)  = fica_f(51)
c      bc_ephem(11)  = fica_f(30)
c      bc_ephem(12)  = fica_f(32)
c      bc_ephem(13)  = fica_f(50)
c      bc_ephem(14)  = fica_f(27)
c      bc_ephem(15)  = fica_f(46)
c      bc_ephem(16)  = fica_f(48)
c
c      xetoe  =bc_ephem(1)
c      xem0   =bc_ephem(2)
c      xedn   =bc_ephem(3)
c      xeart  =bc_ephem(4)
c      xeecc  =bc_ephem(5)
c      xei0   =bc_ephem(6)
c      xeidt  =bc_ephem(7)
c      xeom0  =bc_ephem(8)
c      xeomd  =bc_ephem(9)
c      xew    =bc_ephem(10)
c      xecuc  =bc_ephem(11)
c      xecus  =bc_ephem(12)
c      xecrc  =bc_ephem(13)
c      xecrs  =bc_ephem(14)
c      xecic  =bc_ephem(15)
c      xecis  =bc_ephem(16)

      return
      end
