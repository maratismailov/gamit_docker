      subroutine blk1101
     .   ( debug,iflag,offsl1,offsl2,sample_int,arcvr,asite
     .   , ff,fi,fc,nf,ni,nc )

c     Read a FICA BLK 1101
c     R. King 19 Dec 1990

      logical debug

c     FICA read flag (0=header, 1=data, 2=ephemeris)
      integer iflag

C     L1,L2 phase center offset from mark- Up, N, E
      real*8       offsl1(3), offsl2(3)

c     16-character station and receiver names returned (to match similar blocks
c     for other receivers, but currently, only asite is set, from fc(1-2).
      character*16 asite,arcvr
c   

c     receiver sampling interval
      real*8 sample_int

      include '../includes/makex.h'

c    Provisional MIT Definition of reduced FICA MINIMAC Data Block  Feigl 89-05-01
c
c   FIC Block Name         : Station information
c   FIC Block Number       : 1101
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
c**      character*8 blank8

      save icount

      data icount/0/
c**     data blank8/'        '/     

c*    Dummy statements for unused calling arguments     
      if( debug) then
        print *,'nc ni nf ',nc,ni,nf   
        print *,'fi ',(fi(i),i=1,ni)

        print *,'fc ',(fc(i),i=1,nc)
       endif

                           
c     what is icount for?  rwk 2 Mar 2000
      icount = icount + 1

      do i=1,3
      offsl1(1) = ff(2)
      offsl1(2) = ff(6)
      offsl1(3) = ff(5)
      offsl2(1) = ff(3)
      offsl2(2) = ff(6)
      offsl2(3) = ff(5)
      enddo
            
c ** I don't know where these arose, but they don't match the definition 
c ** above or the FICA files we have from CIGNET ARGO files
c      arcvr  = fc(1)//blank8
c      asite  = fc(3)//fc(4)
      asite = fc(1)//fc(2)
      arcvr = ' '
      sample_int = ff(7)
       
c      if (icount .gt. 2) then
         write (uinfor,114) asite
         write (uinfor,116) arcvr
c         write (uscren,114) asite
c         write (uscren,116) arcvr
c
c
 114     format (1x,'Operator input site code            ',a16)
 116     format (1x,'Operator input receiver serial num  ',a16)
c      endif

      iflag = 0

      return
      end
