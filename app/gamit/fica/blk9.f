      subroutine blk9
     &   ( debug,iflag,nprn,bc_ephem,bc_clock,subfr1
     &   , week,week_set,ff,fi,fc,nf,ni,nc )
c
c     R. King from K. Feigl code in (old) DOFICA  3 October 88
c

C     FICA BLOCK 9
c     map of the floating point items:
c          ff(i)          is
c          (sub frame 1)
c          1              TLM word (preamble)                dimless
c          2              TLM word (message)                 dimless
c          3              HOW word (time)                    GPS sec of week
c          4              Data ID                            dimless
c          5              Sub-frame ID (1)                   dimless
c          6              full week number (10 bit)          GPS weeks
c          7              C/A and/or P flag, L2 flag         dimless
c          8              SV accuracy                        dimless
c          9              SV health                          dimless
c          10             Age of Data clock                  sec
c          11             L2 P data flag                     dimless
c          12             group delay differential           sec
c          13             clock epoch                        GPS sec of week
c          14             clock drift rate                   sec/(sec**2)
c          15             clock drift                        sec/sec
c          16             clock bias                         sec
c          17             (not used)
c          18             (not used)
c          19             Tracker                            dimless
c          20             SV PRN                             dimless
c          (subframe 2)
c          21             TLM word preamble                  dimless
c          22             TLM word message                   dimless
c          23             HOW word (time)                    dimless
c          24             Data id                            dimless
c          25             Sub-frame id should be 2           dimless
c          26             Age of data (ephemeris)            GPS sec of week
c          27             Radial sine correction (CRS)       meters
c          28             Correction to mean motion          radians/sec
c          29             Mean anomaly at epoch              radians
c          30             in-track cosine amp. (CUC)         radians
c          31             eccentricicty                      dimless
c          32             in-track sine amp (CUS)            radians
c          33             square root of sem-major axis      (meters)**1/2
c          34             time of (ephemeris) epoch          GPS sec of week
c          35             fit interval flag                  dimless
c          36             unused (set to 0.0d0)
c          37             unused (set to 0.0d0)
c          38             unused (set to 0.0d0)
c          39             unused (set to 0.0d0)
c          40             unused (set to 0.0d0)
c          (subframe 3)
c          41             TLM word preamble                  dimless
c          42             TLM word message                   dimless
c          43             HOW word (time)                    dimless
c          44             Data id                            dimless
c          45             Sub-frame id (should be 3)         dimless
c          46             inclintation cosine correction CIC rads
c          47             right ascension of ascending node  rads
c          48             incl. sine correction (CIS)        rads
c          49             inclination                        rads
c          50             radial cosine adj                  rads
c          51             argument of perigee                rads
c          52             right ascension of ascending node  rads/sec
c          53             age of data (ephemeris)            GPS sec of week
c          54             inclination time derivative        rads/sec
c          55             unused (set to 0.0d0)
c          56             unused (set to 0.0d0)
c          57             unused (set to 0.0d0)
c          58             unused (set to 0.0d0)
c          59             unused (set to 0.0d0)
c          60             unused (set to 0.0d0)
c

c CORRESPONDANCE table Prescott Phaser, King Phaser, Feigl Fica

c mt, trakid,    satid,     af0,       af1,       af2,       xtoc
c mt, trakid,    NPRN,      XEAF0,     XEAF1,     XEAF2,     XETOC
c     ff(19),    ff(20),    ff(16),    ff(15),    ff(14),    ff(13)

c satid,         xaode,     crs,       deltan,    m0
c NPRN1,         XEADE,     XECRS,     XEDN,      XEM0
c                ff(26),    ff(27),    ff(28),    ff(29)
c bc_ephem:                    14         3          2

c cuc,           e,         cus,       asqrt
c XECUC,         XEECC,     XECUS,     XEART
c ff(30),        ff(31),    ff(32),    ff(33)
c bc_ephem(11)      5          12         4

c xtoe,          cic,       crc,       cis,       idot
c XETOE,         XECIC,     XECRC,     XECIS,     XEIDT
c ff(34),        ff(46),    ff(50),    ff(48),    ff(54)
c bc_ephem(1)      15        13         16          7

c i0,           omega0,     somega,    omedot
c XEI0,         XEOM0,      XEW,       XEOMD
c ff(49),       ff(47),     ff(51),    ff(52)
c bc_ephem(6)      8           10         9

      include '../includes/makex.h'

c     passed values
      logical     debug,
C     .true. once we have a week num
     &            week_set
C     1 for data, 2 for ephemeris
      integer*4   iflag,
C     for broadcast ephem
     &            nprn
C     broadcast ephemeris params
      real*8      bc_ephem(16),
C     broadcast clock params
     &            bc_clock(6),
C     additional quantities from sub-frame 1
     .            subfr1(8)


c     other values
c
c     ephemeris parameters
      real*8           XEAF0, XEAF1, XEAF2, XETOC
     1                 ,XEADE, XECRS, XEDN, XEM0
     2                 ,XECUC, XEECC, XECUS, XEART
     3                 ,XETOE, XECIC, XECRC, XECIS, XEIDT
     4                 ,XEI0, XEOM0, XEW, XEOMD
     5                 ,cflgl2, svaccr, svhlth, aodc, pflgl2, tgd, aode
     6                 ,howsow

c     the FICA arrays
      real*8 ff(maxfic)
      integer*4 fi(maxfic)
      character*8 fc(maxfic)
C     number of elements in array
      integer*4 nf,ni,nc
C     week number, from block 9
      integer*4 week,
C     week number buffer
     .          iweek_buff,
C     tracker number
     .          trakid     
     .         , i
                    

*    Dummy statements for unused calling arguments     
      if( debug) then
        print *,'nc ni nf ',nc,ni,nf     
        print *,'fi',(fi(i),i=1,ni)
        print *,'fc',(fc(i),i=1,nc)
       endif



c        if we already have a week, don't change unless it
c        it is reasonable

         iweek_buff = int(ff(6))
         if (week_set) then
            if (iweek_buff .eq. week .or. iweek_buff .eq. week +1) then
C              reasonable week, OK, update
               week = iweek_buff
            else
C              bad week, abort to read another record
               return
            endif
         else
            week = iweek_buff
            week_set = .true.
         endif

        if(debug) then
          print 10, ff(6),iweek_buff,week,week_set
10        format(' In BLK 9:  ff(6),iweek_buff,week,week_set='
     1        ,  D16.8,2I5,L4 )
          print 11, ff(20),ff(34),ff(15),ff(16),ff(17)
11        format(' prn,xetoe,clock terms:',5d12.5)
          endif

c        tracker id ok
           trakid = int(ff(19))
c        sat number ok
           nprn = int(ff(20))
c        t oc  -
           xetoc = ff(13)
c        af2   clock drift rate ok
           xeaf2 = ff(14)
c        af1   clock drift ok
           xeaf1 = ff(15)
c        af0   clock bias  ok
           xeaf0 = ff(16)
c        AODE  seconds ok
           xeade = ff(26)
c        C rs  meters ok
           xecrs = ff(27)
c        Delta n ok
           xedn = ff(28)
c        M0 ok
           xem0 = ff(29)
c        C uc ok
           xecuc =  ff(30)
c        e ok
           xeecc = ff(31)
c        C us ok
           xecus = ff(32)
c        Square root A ok
           xeart = ff(33)
c        toe ok
           xetoe = ff(34)
c        C ic ok
           xecic = ff(46)
c        Omega 0 ok
           xeom0 = ff(47)
c        C is  ok
           xecis = ff(48)
c        i0 ok
           xei0 = ff(49)
c        C rc ok
           xecrc = ff(50)
c        Small omega ok
           xew = ff(51)
c        Omega dot ok
           xeomd = ff(52)
c        I dot ok
           xeidt = ff(54)

c        ----these quanties passed to FICA_TO_RINEX but not to MAKEX
c        L2 code flag
           cflgl2 = ff(7)
c        SV accuracy
           svaccr = ff(8)
c        SV health
           svhlth = ff(9)
c        Clock age of data
           aodc   = ff(10)
c        L2 P data flag
           pflgl2 = ff(11)
c        Group delay differential
           tgd    = ff(12)
c        Ephemeris age of data
           aode   = ff(26)
c        Time of frame - HOW word
           howsow = ff(3)


c        store these guys in arrays for easy passing
C        broadcast ephemeris parameters
         bc_ephem(1)= xetoe
         bc_ephem(2)= xem0
         bc_ephem(3)= xedn
         bc_ephem(4)= xeart
         bc_ephem(5)= xeecc
         bc_ephem(6)= xei0
         bc_ephem(7)= xeidt
         bc_ephem(8)= xeom0
         bc_ephem(9)= xeomd
         bc_ephem(10)=xew
         bc_ephem(11)=xecuc
         bc_ephem(12)=xecus
         bc_ephem(13)=xecrc
         bc_ephem(14)=xecrs
         bc_ephem(15)=xecic
         bc_ephem(16)=xecis

c        broadcast clock parameters
         bc_clock(1)=xetoc
         bc_clock(2)=xeaf0
         bc_clock(3)=xeaf1
         bc_clock(4)=xeaf2
C        xeadc age of clock data.  parameter not used.
         bc_clock(5)=0.d0
C        xetdg Ionospheric correction of clock data. Not used.
         bc_clock(6)=0.d0

c        additional broadcast and clock parameters from sub-frame 1
         subfr1(1) = cflgl2
         subfr1(2) = svaccr
         subfr1(3) = svhlth
         subfr1(4) = aodc
         subfr1(5) = pflgl2
         subfr1(6) = tgd
         subfr1(7) = aode
         subfr1(8) = howsow

      iflag = 2
         
      return
      end
