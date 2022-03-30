      subroutine blk80
     &   ( debug,iflag,svid,tmode,igpswk,gpssec
     &   , dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     &   , iwkn_file,week_set,iwkn_start,sow_start
     &   , ff,fi,fc,nf,ni,nc)

c     Read a FICA BLK 80 (GESAR or CORE) for the TI 4100
c     R. King 3 October 88
c     K. Feigl Feb 89
c
c   All times in GPS time
c
c--------------------------------------------------------------------------------
c
c    ARL approved Definition of reduced FICA TI4100 Data Block  Feigl 88/02/18
c
c   FIC Block Name         : Tracking data (reduced size)
c   FIC Block Number       : 80
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
c   1                      Epoch of data                           GPS sec of week
c   2-4                    Dummy, for human readability
c      For each satellite channel
c   5 (9, 13, 17)          L1 carrier Doppler phase (biased)       cycles
c   6 (10, 14, 18)         L2 carrier Doppler phase (biased        cycles
c   7 (11, 15, 19)         L1 pseudorange                          kilometers
c   8 (12, 16, 20)         L2 pseudorange                          kilometers
c
c   Integer Items
c   1                      Number of channels (nsat) to follow
c   2                      GPS week
c   3-6                    Dummy, for human readability
c       For each satellite channel
c   7 (13, 19, 25)         SV id (PRN) of tracker
c   8 (14, 20, 26)         Tracker mode
c   9 (15, 21, 27)         Quality vector for L1
c  10 (16, 22, 28)         Quality vector for L2
c  11 (17, 23, 29)         Carrier signal-to-noise ratio, L1       dB-Hz
c  12 (18, 24, 30)         Carrier signal-to-noise ratio, L2       dB-Hz

c----------------------------------------------------------------------------------



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

c     other values
c
c     FICA arrays
c     the FICA arrays
      integer*4 nf,ni,nc
      real*8 ff(maxfic)
      integer*4 fi(maxfic)
      character*8 fc(maxfic)

c     other values
C     seconds since start of sceno
      real*8       dt
C     second of week for data
     .,            sow_ttag
C     function for differencing times
     .,            secdif
C     unbiased L1 phase
     .,            biasph_l1
C     unbiased L2 phase
     .,            biasph_l2

      integer*4   iwkn_ttag
C     number of sats in block
     .,            itrak
     .,            i


c     for trick equivalencing of quality vectors.
c     warning, may be machine dependent
      integer*2 iqv2(2)
      integer*4 iqv4

      equivalence (iqv4,iqv2(1))
          
      logical week_warning
      integer*4 len,rcpar
      character*16 prog_name   
      character*256 message              
           
      data week_warning/.false./

c     get calling program name for report_stat
      len = rcpar(0,prog_name)
   
*    Dummy statements for unused calling arguments     
      if( debug) then
        print *,'nc ni nf ',nc,ni,nf
        print *,'fc',(fc(i),i=1,nc)
       endif



c     set quality codes, tracking mode, and snr to
c     default (no data) values
      do i=1,maxchn
         svid(i)     = 0
         iqvec(i,1)  = 65535
         iqvec(i,2)  = 65535
         tmode(i)    = 58
         denrat(i,1) = 0.d0
         denrat(i,2) = 0.d0
      enddo

c     determine the number of valid channels
      itrak= fi(1)

c     set the time tags
      iwkn_ttag = fi(2)
      sow_ttag  = ff(1)
      week_set  = .true.

      if (iwkn_ttag .ne. iwkn_file) then
         write (message,10) iwkn_ttag,iwkn_file
         write (uinfor,10) iwkn_ttag,iwkn_file
 10      format (1X,'BLK80 WARNING Week number from Block 80 is ',i5,
     .   'But from file name it is', i5) 
         if( .not.week_warning ) then 
          call report_stat('WARNING',prog_name,'makex/blk70',message,0)
          week_warning = .true.
         endif
      endif

c     are these arrays because of PHASER ?
      do i=1,maxchn
         gpssec(i) = sow_ttag
         igpswk(i) = iwkn_ttag
      enddo

c     Compute the time interval for the "frequency plan correction"
c     seconds from the start of the scenario
      dt = secdif (iwkn_ttag,sow_ttag,iwkn_start,sow_start)

      if (debug) print *, 'In BLK80: dt (seconds) = ', dt

c     loop over valid channels

      do i=1,itrak
c         prn numbers of tracked SVs
          svid(i) = fi(6*i+1)

c         tracker mode
          tmode(i) = fi(6*i+2)

c         quality vectors
c         only the least significant 2 bytes
c         can mean anything, so chop of most
c         significant bytes via tricky equivalence
          iqv4  = fi(6*i+3)
          iqvec(i,1) = iqv2(2)
          iqv4  = fi(6*i+4)
          iqvec(i,2) = iqv2(2)

c         signal-to-noise ratio
          denrat(i,1) = dble(fi(6*i+5))
          denrat(i,2) = dble(fi(6*i+6))

c         L1 phase biased according to "frequency plan"
          biasph_l1 = ff(4*i+1)
c         add 6000 to remove this bias
          dofl1(i) = biasph_l1 + 6000.d0 * dt

c         L2 phase biased according to "frequency plan"
          biasph_l2 = ff(4*i+2)
c         subtract 7600 to remove this bias
          dofl2(i) = biasph_l2 - 7600.d0 * dt

          if (debug) then
             print *, 'In BLK80: biased phases L1, L2 ',
     .       biasph_l1, biasph_l2
          endif

c         pseudorange
c         L1 pseudorange (meters)
          prgl1(i) = ff(4*i+3) * 1000.0d0
c
c         L2 pseudorange (meters)
          prgl2(i) = ff(4*i+4) * 1000.0d0
      enddo

CD     print *,'BLK80: week ',igpswk(1),' sec of week ',gpssec(1)
CD     print *,'BLK80: PRNs           ',(svid(i),i=1,4)
CD     print *,'BLK80: tracker modes  ',(tmode(i),i=1,4)
CD     print 203,' BLK80: L1 unbiased ph ',(dofl1(i),i=1,4)
CD     print 203,' BLK80: L2 unbiased ph ',(dofl2(i),i=1,4)
CD     print 203,' BLK80: L1 pseudoranges',(prgl1(i),i=1,4)
CD     print 203,' BLK80: L2 pseudoranges',(prgl2(i),i=1,4)
CD 203  format (a,4(1x,1pe22.14))


      iflag = 1


      return
      end
