      subroutine blk670
     &   ( debug,iflag,svid,tmode,igpswk,gpssec
     &   , dofl1,dofl2,prgl1,prgl2,denrat,iqvec,week,week_set
     &   , ff,fi,fc,nf,ni,nc )

c     Read a FICA BLK 670 for the MACROMETER II
c     R. King 5 October 88
c
c     MACROMETER times are UTC; output is GPST
c     MACROMETER phases are double frequency; output is not
c
c--------------------------------------------------------------------------------
c
c    Provisional MIT Definition of reduced FICA MACROMETER II Data Block  King 88-10-04
c
c   FIC Block Name         : Tracking data (reduced)
c   FIC Block Number       : 670
c   Original Block Number  : None
c   Original Block Source  : GESAR/ASDAP
c
c   Number of Floating Point Items : Variable (minimum 4)
c   Number of Integer Items        : Variable (minimum 6)
c   Number of Character Items      : 0
c
c
c   Item Numbers           Item Description                        Units
c   _______________        _________________________________       ______________
c
c   Floating Point Items
c   1                        Seconds of epoch                      UTC sec
c   2-4                      Dummy, for human readability
c      For each satellite channel
c   5 (7,  9, 11, 13, 15)    2 * (L1 carrier Doppler phase)        cycles
c   6 (8, 10, 12, 14, 16)    2 * (L2 carrier Doppler phase)        cycles
c
c   Integer Items
c   1                        Number of channels (nsat) to follow
c   2                        Year
c   3                        Day of year
c   4                        Hours
c   5                        Minutes
c   6                        Epoch number
c       For each satellite channel
c   7 (10, 13, 16, 19, 22)   SV id (PRN) of tracker
c   8 (11, 14, 17, 20  23)   Carrier amplitude L1
c   9 (12, 15, 18) 21, 24)   Carrier amplitude L2

c----------------------------------------------------------------------------------



      include '../includes/makex.h'

c     passed values
C     .true. to print things
      logical     debug,
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

      real*8  sec,utcoff


c     the FICA arrays
      real*8 ff(maxfic)
      integer*4 fi(maxfic)
      character*8 fc(maxfic)
C     number of elements in array
      integer*4 nf,ni,nc
C     week number
      integer*4 week

      integer*4 ihr,min,iyr,idoy,itrak,i


c     for trick equivalencing of quality vectors.
c     warning, may be machine dependent
      integer*2 iqv2(2)
      integer*4 iqv4

      equivalence (iqv4,iqv2(1))
               

*    Dummy statements for unused calling arguments     
      if( debug) then
        print *,'nc ni nf ',nc,ni,nf
        print *,'fc',(fc(i),i=1,nc) 
        print *,'prgl1 prgl2 ',prgl1,prgl2
       endif



c     set quality codes, tracking mode,  and snr to
c         default (no data) values
            do i=1,maxchn
               svid(i)     = 0
C              not defined for MACROMETER
               iqvec(i,1)  = 65535
C              not defined for MACROMETER
               iqvec(i,2)  = 65535
C              not defined for MACROMETER
               tmode(i)    = 58
               denrat(i,1) = 0.d0
               denrat(i,2) = 0.d0
            enddo


c     determine the number of valid channels
        itrak= fi(1)

c     set the time tags
        iyr       = fi(2) 
        idoy      = fi(3)
        ihr       = fi(4)
        min       = fi(5)
        sec       = ff(1)

CD     print 10,debug,iyr,idoy,ihr,min,sec
CD 10    format(' In BLK670: debug,iyr,idoy,ihr,min,sec',L4,4I4,F8.4 )
c       convert from UTC yy/ddd hh:mm:ss to GPST wkn/sow
        iflag = -2
        call timcon
     .    (iflag,igpswk(1),gpssec(1),iyr,idoy,ihr,min,sec,utcoff)
        week      = igpswk(1)
        week_set  = .true.


c     loop over valid channels

      do i=1,itrak

c       prn numbers of tracked SVs
            svid(i) = fi(3*i+4)

c       signal-to-noise ratio
            denrat(i,1) = dble(fi(3*i+5))
            denrat(i,2) = dble(fi(3*i+6))

c       phase
c           L1 doppler phase
            dofl1(i) = 0.5d0 * ff(2*i+3)
c           L2 doppler phase
            dofl2(i) = 0.5d0 * ff(2*i+4)

c       set quality codes to zero to avoid losing data
            iqvec(i,1) = 0
            iqvec(i,2) = 0

      enddo


CD     if (debug) print *,'BLK670: week ',week,' sec of week ',gpssec(1)
CD     if (debug) print *,'BLK670: PRNs           ',(svid(i),i=1,4)
CD     if (debug) print 203,' BLK670: L1 phases      ',(dofl1(i),i=1,4)
CD     if (debug) print 203,' BLK670: L2 phases      ',(dofl2(i),i=1,4)
CD     if (debug) print 203,' BLK670: L1 amplitudes ',(denrat(i,1),i=1,4)
CD     if (debug) print 203,' BLK670: L2 amplitudes ',(denrat(i,2),i=1,4)

CD 203  format (a,4(1x,1pe22.14))


         iflag = 1


      return
      end
