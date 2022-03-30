      subroutine blk6
     &   ( debug,iflag,svid,tmode,igpswk,gpssec
     &   , dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     &   , iwkn_file,week_set,iwkn_start,sow_start
     &   , ff,fi,fc,nf,ni,nc )

c     Read a FICA BLK 6 (GESAR) or BLK 55 (CORE) for the TI 4100
c     R. King from K. Feigl code in (old) DOFICA  3 October 88
c
c   All times in GPS time
c
c     FICA BLOCK 6
c     map of the floating point items:
c          ff(i)      is
c          1              pseudorange of FTF validity   sec
c          2              FTF bias offset               sec
c          3              user epoch time of pseudorng  GPS sec of week
c          4-11           L1,L2 carrier signal to noise db-Hz
c                         tracker, frequency
c          12-15          L1 pseudorange                kilometers
c          16-19          L2 pseudorange                kilometers
c          20-23          L1 carrier Doppler phase at   cycles
c                         code FTF
c          24-27          L2 carrier Doppler phase at   cycles
c                         code FTF
c          28             recvr internal temperature    deg C
c     map of the integer items:
c          fi(i)      is:                           units
c          1-4            SV PRN of each tracker
c          5-8            tracker mode
c          9-16           L1,L2 quality vector
c                         (tracker,frequency)
C
c     FICA BLOCK 55 (SAME AS 6 FOR THE GOOD STUFF)
C
c     map of the floating point items:
c          ff(i)      is
c          1              pseudorange of FTF validity   sec
c          2              FTF bias offset               sec
c          3              user epoch time of pseudorng  GPS sec of week
c          4-11           L1,L2 carrier signal to noise db-Hz
c                         tracker, frequency
c          12-15          L1 pseudorange                kilometers
c          16-19          L2 pseudorange                kilometers
c          20-23          L1 carrier Doppler phase at   cycles
c                         code FTF
c          24-27          L2 carrier Doppler phase at   cycles
c                         code FTF
c          28-31          L1 carrier velocity at code   Hz
c                         FTF
c          32-35          L2 carrier velocity at code   Hz
c                         FTF
c          36-39          averaged line of site accel   m/s**2
c          40             recvr internal temperature    deg C
c     map of the integer items:
c          fi(i)      is:                           units
c          1-4            SV PRN of each tracker
c          5-8            tracker mode
c          9-12           L1 quality vector
c          13-16          L2 quality vector


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
     &            iqvec(maxchn,2),
C     week number for start of scenario
     &            iwkn_start,
C     week number from file name
     &            iwkn_file

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
     &            denrat(maxchn,2),
C     second of week for sceno. start
     &            sow_start

C     seconds since start of scenario
      real*8       dt
C     second of week for data
     .,            sow_ttag
C     function for differencing times
     .,            secdif
C     unbiased L1 phase
     .,            biasph_l1
C     unbiased L2 phase
     .,            biasph_l2


c     the FICA arrays
      real*8 ff(maxfic)
      integer*4 fi(maxfic)
      character*8 fc(maxfic)
C     number of elements in FICA arrays
      integer*4 nf,ni,nc
C     pointers into FICA arrays
      integer*4 ptri,ptrf

      integer*4 i,j
C     week number for data time tag
      integer*4 iwkn_ttag


c     for trick equivalencing of quality vectors.
c     warning, may be machine dependent
      integer*2 iqv2(2)
      integer*4 iqv4

      equivalence (iqv4,iqv2(1))
                  

*    Dummy statements for unused calling arguments     
      if( debug) then
        print *,'nc ni nf ',nc,ni,nf
        print *,'fc',(fc(i),i=1,nc)
       endif



CD     print *,'BLK6 nf,ni,nc',nf,ni,nc
c      if(debug) print *,'BLK6 week_set iwkn_file ',week_set,iwkn_file
      if (.not.week_set ) then
         if (debug) print *,'BLK6: week not set, returning.'
         return
      else
c        set the time tags
         sow_ttag  = ff(3)
         iwkn_ttag = iwkn_file  

c        Compute the time interval for the "frequency plan correction"
c        seconds from the start of the scenario
         dt = secdif (iwkn_ttag,sow_ttag,iwkn_start,sow_start)

         if (debug) print *, 'BLK6: dt (seconds) = ', dt

         do 20 i=1,4
            gpssec(i) = sow_ttag
            igpswk(i) = iwkn_ttag
c
c           tracker mode
            tmode(i) = fi(4+i)
c
c           satellite id
            svid(i) = fi(i)
c
c           L1 pseudorange (meters)
            prgl1(i) = ff(11+i) * 1000.0d0
c
c           L2 pseudorange (meters)
            prgl2(i) = ff(15+i) * 1000.0d0
c
c           L1 phase biased according to "frequency plan"
            biasph_l1 = ff(19+i)
c           add 6000 to remove this bias
            dofl1(i) = biasph_l1 + 6000.d0 * dt
c
c           L2 phase biased according to "frequency plan"
            biasph_l2 = ff(23+i)
c           subtract 7600 to remove this bias
            dofl2(i) = biasph_l2 - 7600.d0 * dt

            if (debug) then
               print *, 'BLK6:  biased phases L1, L2 ',
     .         biasph_l1, biasph_l2
            endif
   20    continue

c        signal to noise ratios and quality vectors
            ptrf = 4
            ptri = 9
C           loop over L1 and L2
            do 50  j=1,2
C              loop over 4 trackers
               do 40  i=1,4
                   denrat(i,j) = ff(ptrf)
c                  only the least significant 2 bytes
c                  can mean anything, so chop of most
c                  significant bytes via tricky equivalence
                   iqv4  = fi(ptri)
                   iqvec(i,j) = iqv2(2)
                   ptrf = ptrf + 1
                   ptri = ptri + 1
   40       continue
   50    continue


      if (debug) print *,'BLK6: week, sec of week ',igpswk(1),gpssec(1)
      if (debug) print *,'BLK6: PRNs           ',(svid(i),i=1,4)
      if (debug) print *,'BLK6: tracker modes  ',(tmode(i),i=1,4)
      if (debug) print 203,' BLK6: L1 unbiased ph ',(dofl1(i),i=1,4)
      if (debug) print 203,' BLK6: L2 unbiased ph ',(dofl2(i),i=1,4)
      if (debug) print 203,' BLK6: L1 pseudoranges',(prgl1(i),i=1,4)
      if (debug) print 203,' BLK6: L2 pseudoranges',(prgl2(i),i=1,4)
 203  format (a,4(1x,1pe22.14))


         iflag = 1

      endif


      return
      end
