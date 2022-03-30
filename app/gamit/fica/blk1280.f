      subroutine blk1280
     &   ( debug,iflag,svid,tmode,igpswk,gpssec
     &   , dofl1,dofl2,prgl1,prgl2,denrat,iqvec,week,week_set
     &   , ff,fi,fc,nf,ni,nc )

c     Read a FICA BLK 1280 for the Trimble (from CIGNET)
c     R. King   19 Dec 1990
c
c
c--------------------------------------------------------------------------------
c
c    Provisional MIT Definition of reduced FICA Trimble Data Block  King 90-12-19 89-05-01
c
c   FIC Block Name         : Tracking data
c   FIC Block Number       : 1280
c   Original Block Number  : None
c   Original Block Source  : CIGNET file
c
c   Number of Floating Point Items : Variable (minimum 4)
c   Number of Integer Items        : Variable (minimum 6)
c   Number of Character Items      : 0

c   Item Numbers           Item Description                        Units
c   _______________        _________________________________       ______________
c
c   Floating Point Items
c
c   1                      Epoch of data                           GPS sec of week (GPST)
c   2                      Epoch of data                           second of epoch (GPST)
c   3-4                    blank, for readability                  0.0
c      For each satellite channel
c   5 (9, 13, 17, ...)     L1 carrier Doppler phase                full cycles
c                          (positive sign for increasing Doppler shift)
c                          (positive sign for decreasing Phaset)
c   6 (10, 14, 18, ...)    L2 carrier Doppler phase                half cycles
c                          (positive sign for increasing Doppler shift)
c                          (positive sign for decreasing Phaset)
c   7 (11, 15, 19, ...)    L1 pseudorange                          kilometers
c   8 (12, 16, 20, ...)    blank                                   -9999.0
c
c   Integer Items
c
c   1                      Number of channels (nsat) to follow
c   2                      GPS week
c   3                      Year
c   4                      Day of year
c   5                      Hours
c   6                      Minutes
c       For each satellite channel
c   7 (10, 13, 16, ...)    SV id (PRN) of tracker
c   8 (11, 14, 17, ...)    Carrier signal-to-noise ratio, L1       dB-Hz
c   9 (12, 15, 18, ...)    Carrier signal-to-noise ratio, L2       dB-Hz



      include '../includes/makex.h'

c     passed values
      logical     debug,
C     .true. once we have a week num
     &            week_set
C     1 for data, 2 for ephemeris
      integer*4   iflag,
C     sat id num PRN
     &            svid(maxchn),
C     tracking mode
     &            tmode(maxchn),
C     gps week number for epoch
     &            igpswk(maxchn),
C     quality vector L1,L2
     &            iqvec(maxchn,2)
C     gps sec of week for epoch
      real*8      gpssec(maxchn),
C     L1 doppler phase (DATA!)
     &            dofl1(maxchn),
C     L2 doppler phase (DATA!)
     &            dofl2(maxchn),
C     L1 pseudorange
     &            prgl1(maxchn),
C     L2 pseudorange
     &            prgl2(maxchn),
C     signal to noise ratio L1,L2
     &            denrat(maxchn,2)

c     other values
c


c     the FICA arrays
      real*8 ff(maxfic)
      integer*4 fi(maxfic)
      character*8 fc(maxfic)
C     number of elements in array
      integer*4 nf,ni,nc
C     week number
      integer*4 week

      integer*4 itrak,i,j
                   
c*    Dummy statements for unused calling arguments     
      if( debug) then
        print *,'nc ni nf ',nc,ni,nf
        print *,'fc',(fc(i),i=1,nc)
       endif




c     set quality codes, tracking mode,  and snr to
c         default (no data) values
            do i=1,maxchn
               svid(i)     = 0
C              not defined for Trimble
               iqvec(i,1)  = 65535
C              not defined for Trimble
               iqvec(i,2)  = 65535
C              not defined for Trimble
               tmode(i)    = 58
               denrat(i,1) = 0.d0
               denrat(i,2) = 0.d0
            enddo


c     determine the number of valid channels
      itrak= fi(1)

c     GPST week number and second-of-week
      igpswk(1) = fi(2)
      gpssec(1) = ff(1)
      week      = igpswk(1)
      week_set  = .true.


c     loop over valid channels
      do 100 i=1,itrak
         j = i-1

c        prn numbers of tracked SVs
         svid(i) = fi(7+3*j)

c        signal-to-noise ratio
         denrat(i,1) = dble(fi(8+3*j))
         denrat(i,2) = dble(fi(9+3*j))

c        phase
c        X-files follow the "pseudorange" convention for phase sign:
c        positive phase for increasing pseudorange
c        positive phase for decreasing Doppler shift
c        L1 phase
         dofl1(i) = ff(5+4*j)
c        L2 phase in full cycles
         dofl2(i) = ff(6+4*j)/2.0d0

c        pseudorange
c        L1 pseudorange in meters
         prgl1(i) = ff(7+4*j) * 1000.0d0
c        L2 pseudorange is not defined
         prgl2(i) = 0.0d0

c        set quality codes to zero to avoid losing data
         iqvec(i,1) = 0
         iqvec(i,2) = 0
 100  continue


CD     if (debug) print *,'BLK1280: week ',week,' sec of week ',gpssec(1)
CD     if (debug) print *,'BLK1280: PRNs           ',(svid(i),i=1,4)
CD     if (debug) print 203,' BLK1280: L1 phases      ',(dofl1(i),i=1,4)
CD     if (debug) print 203,' BLK1280: L2 phases      ',(dofl2(i),i=1,4)
CD     if (debug) print 203,' BLK1280: L1 amplitudes ',(denrat(i,1),i=1,4)
CD     if (debug) print 203,' BLK1280: L2 amplitudes ',(denrat(i,2),i=1,4)

CD 203  format (a,4(1x,1pe22.14))


      iflag = 1


      return
      end
