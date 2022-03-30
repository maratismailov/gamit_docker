      subroutine blk1201
     .   ( debug,iflag,offsl1,offsl2,arcvr,asite,asampl
     .   , ff,fi,fc,nf,ni,nc )

c     Read a FICA BLK 1201
c     R. King 19 Dec 1990


      logical debug

c     FICA read flag (0=header, 1=data, 2=ephemeris)
      integer iflag

C     L1,L2 phase center offset from mark- Up, N, E
      real*8       offsl1(3), offsl2(3)

C     site name and receiver serial number
      character*16 asite, arcvr

c     receiver collection interval
      real*8 asampl

      include '../includes/makex.h'

c    Provisional MIT Definition of reduced Trimble Data Block  King 90-12-19
c
c   FIC Block Name         : Station information
c   FIC Block Number       : 1201
c   Original Block Number  : None
c   Original Block Source  : CIGNET file
c
c   Number of Floating Point Items : 7
c   Number of Integer Items        : 0
c   Number of Character Items      : 8

c   Item Numbers           Item Description                        Units
c   _______________        _________________________________       ______________
c
c   Floating Point Items
c
c   1                      Minimac software version                none
c   2                      L1 phase center above mark              meters
c   3                      L2 phase center above mark              meters
c   4                      antenna base above mark                 meters
c   5                      East antenna offset                     meters
c   6                      North antenna offset                    meters
c   7                      Data collection interval                seconds

c   Integer Items
c
c   none
c
c   Character Items
c
c   1-4                    32 character station name
c   5-8                    Operator name
c


c     FICA ARRAYS
      real*8 ff(maxfic)
      integer*4 fi(maxfic)
      character*8 fc(maxfic)
C     number of elements in FICA arrays
      integer*4 nf,ni,nc

      integer*4 i,icount
      character*8 blank8

      save icount

      data icount/0/,blank8/'        '/

      icount = icount + 1
                            
c*    Dummy statements for unused calling arguments     
      if( debug) then
        print *,'nc ni nf ',nc,ni,nf
        print *,'fi',(fi(i),i=1,nc)
       endif


      do i=1,3
      offsl1(1) = ff(2)
      offsl1(2) = ff(6)
      offsl1(3) = ff(5)
      offsl2(1) = ff(3)
      offsl2(2) = ff(6)
      offsl2(3) = ff(5)
      enddo

      arcvr  = fc(1)//blank8
      asite  = fc(3)//fc(4)
      asampl = ff(7)

      if (icount .gt. 2) then
         write (uinfor,114) asite
         write (uinfor,116) arcvr

         write (uscren,114) asite
         write (uscren,116) arcvr


 114     format (1x,'Operator input site code            ',a16)
 116     format (1x,'Operator input receiver serial num  ',a16)
      endif

      iflag = 0

      return
      end
