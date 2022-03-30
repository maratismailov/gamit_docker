      subroutine blk1001
     &   ( debug,iflag,offsl1,offsl2,arcvr,asite,asampl
     .   , ff,fi,fc,nf,ni,nc )

c     Read a FICA BLK 1001
c     R. King 20 Feb 1989

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

c    MIT Definition of MACROMETER II Block 1001        89-02-20

c   FIC Block Name         : MACROMETER II  Input Data
c   FIC Block Number       : 1001
c   Original Block Number  : None
c   Original Block Source  : GAMIT X-File
c
c   Number of Floating Point Items :  9
c   Number of Integer Items        : 13
c   Number of Character Items      : 11

c   Item Numbers           Item Description                        Units
c   _______________        _________________________________       ______________
c
c   Floating Point Items

c   1,2,3                  Antenna L1 phase center wrt mark         meters
c                              Up, North, East
c   4,5,6                  Antenna L2 phase center wrt mark         meters
c   7                      L1 phase multiplication factor
c   8                      L2 phase multiplication factor
c   9                      Seconds of first epoch                   sec
c
c   Integer Items
c   1                      Number of satellites observed
c   2-7                    Satellite PRN numbers
c   8-11                   First epoch                              yr,doy,hr,min
c   12                     Collection interval                      sec
c   13                     Number of epochs
c
c   Character Items
c   1                      Receiver serial number
c   2                      Receiver software version
c   3-4                    Site name
c   6-7                    Latitude of site                        deg,min,sec
c   8-9                    Longitude (W) of site                   deg,min,sec
c   10-11                  Height of site                          meters



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
                   

c*    Dummy statements for unused calling arguments     
      if( debug) then
        print *,'nc ',nc 
        print *,'ni ',ni
        print *,'nf ',nf
       endif

      icount = icount + 1

      do i=1,3
      offsl1(i) = ff(i)
      offsl2(i) = ff(i+3)
      enddo

      arcvr  = fc(1)//blank8
      asite  = fc(3)//fc(4)
      asampl = dble( fi(12) )

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
