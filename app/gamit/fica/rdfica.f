      subroutine read_ficaf
     &   (debug,iflag,fend,ferr,nprn,svid,tmode,igpswk,
     &   gpssec,dofl1,dofl2,prgl1,prgl2,denrat,iqvec,
     &   bc_ephem,bc_clock)

c
c     read a FICA format data file open on logical unit lu
c
c     Note that not all the returned values will be valid, depending
c     on the value of iflag. For example, bc_ephem and bc_clock
c     will be bogus if iflag = 1, and the phase data will be
c     meaningless if iflag = 2. Iflag is returned to the main program
c     according to what kind of record was read in.
c
c     assume all times in GPS time
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
C
c               TI-NAVIGATOR/ARL:UT FIC Block Definition  02-03-88
c
c    FIC Block Name         : Receiver-Ranging Measurements Data Block
c    FIC Block Number       : 401
c    Original Block Number  : 1
c    Original Block Source  : TI-NAVIGATOR Operating Software
c
c    Number of Floating Point Items : 111
c    Number of Integer Items        : 44
c    Number of Character Items      : 0
c
c
c    Item Numbers           Item Description                        Units
c    _______________        _________________________________       ______________
c
c    Floating Point Items: ff(i)
c
c    i
c    1                      Current FTF                             Seconds
c    2                      System FTF count                        Seconds
c    3 *                    User Epoch Time                         Seconds
c    4                      Current GPS Week                        GPS Weeks
c    5                      FTF of code phase                       Seconds
c    6 *                    FTF of carrier start phase              Seconds
c    7                      FTF Offset                              Seconds
c    8                      FTF interval of carrier phase           Seconds
c    9                      Internal temperature                    Deg. C
c   10                      Wideband AGC L1                         dB
c   11                      Wideband AGC L2                         dB
c        (Tracker 1)
c    12                     Zcount for L1                           Seconds
c    13                     X1 for L1                               Pchips
c    14                     P-Phase for L1                          Pchips
c    15                     Zcount for L2                           Seconds
c    16                     X1 for L2                               Pchips
c    17                     P-Phase for L2                          Pchips
c    18                     Code Frequency for L1                   Hz
c    19                     Code Frequency for L2                   Hz
c    20 *                   Carrier Phase at start FTF(L1)          Cycles
c    21 *                   Carrier Phase at start FTF(L2)          Cycles
c    22                     Carrier Phase at start FTF(L2)          Cycles
c    23                     Carrier Phase at start FTF(L2)          Cycles
c    24                     Carrier Velocity at code FTF(L1)        Hz
c    25                     Carrier Velocity at code FTF(L2)        Hz
c    26                     Carrier Phase at end FTF(L1)            Cycles
c    27                     Carrier Phase at end FTF(L2)            Cycles
c    28                     Avg. Line-of-Sight Acceleration         m/s**2
c    29                     C/N0  L1                                0.5 dB/Hz
c    30                     C/N0  L2                                0.5 dB/Hz
c    31                     Code loop DLL sum bandwidth             Hz
c    32                     Code loop difference bandwidth          Hz
c    33                     Code loop PLL sum bandwidth             Hz
c    34                     Code loop PLL diff. bandwidth           Hz
c    35                     Narrowband AGC L1                       dB
c    36                     Narrowband AGC L2                       dB
c           (Tracker 2)
c    37                     Zcount for L1                           Seconds
c    38                     X1 for L1                               Pchips
c    39                     P-Phase for L1                          Pchips
c    40                     Zcount for L2                           Seconds
c    41                     X1 for L2                               Pchips
c    42                     P-Phase for L2                          Pchips
c    43                     Code Frequency for L1                   Hz
c    44                     Code Frequency for L2                   Hz
c    45 *                   Carrier Phase at start FTF(L1)          Cycles
c    46 *                   Carrier Phase at start FTF(L2)          Cycles
c    47                     Carrier Phase at start FTF(L2)          Cycles
c    48                     Carrier Phase at start FTF(L2)          Cycles
c    49                     Carrier Velocity at code FTF(L1)        Hz
c    50                     Carrier Velocity at code FTF(L2)        Hz
c    51                     Carrier Phase at end FTF(L1)            Cycles
c    52                     Carrier Phase at end FTF(L2)            Cycles
c    53                     Avg. Line-of-Sight Acceleration         m/s**2
c    54                     C/N0  L1                                0.5 dB/Hz
c    55                     C/N0  L2                                0.5 dB/Hz
c    56                     Code loop DLL sum bandwidth             Hz
c    57                     Code loop difference bandwidth          Hz
c    58                     Code loop PLL sum bandwidth             Hz
c    59                     Code loop PLL diff. bandwidth           Hz
c    60                     Narrowband AGC L1                       dB
c    61                     Narrowband AGC L2                       dB
c           (Tracker 3)
c    62                     Zcount for L1                           Seconds
c    63                     X1 for L1                               Pchips
c    64                     P-Phase for L1                          Pchips
c    65                     Zcount for L2                           Seconds
c    66                     X1 for L2                               Pchips
c    67                     P-Phase for L2                          Pchips
c    68                     Code Frequency for L1                   Hz
c    69                     Code Frequency for L2                   Hz
c    70 *                   Carrier Phase at start FTF(L1)          Cycles
c    71 *                   Carrier Phase at start FTF(L2)          Cycles
c    72                     Carrier Phase at start FTF(L2)          Cycles
c    73                     Carrier Phase at start FTF(L2)          Cycles
c    74                     Carrier Velocity at code FTF(L1)        Hz
c    75                     Carrier Velocity at code FTF(L2)        Hz
c    76                     Carrier Phase at end FTF(L1)            Cycles
c    77                     Carrier Phase at end FTF(L2)            Cycles
c    78                     Avg. Line-of-Sight Acceleration         m/s**2
c    79                     C/N0  L1                                0.5 dB/Hz
c    80                     C/N0  L2                                0.5 dB/Hz
c    81                     Code loop DLL sum bandwidth             Hz
c    82                     Code loop difference bandwidth          Hz
c    83                     Code loop PLL sum bandwidth             Hz
c    84                     Code loop PLL diff. bandwidth           Hz
c    85                     Narrowband AGC L1                       dB
c    86                     Narrowband AGC L2                       dB
c           (Tracker 4)
c    87                     Zcount for L1                           Seconds
c    88                     X1 for L1                               Pchips
c    89                     P-Phase for L1                          Pchips
c    90                     Zcount for L2                           Seconds
c    91                     X1 for L2                               Pchips
c    92                     P-Phase for L2                          Pchips
c    93                     Code Frequency for L1                   Hz
c    94                     Code Frequency for L2                   Hz
c    95 *                   Carrier Phase at start FTF(L1)          Cycles
c    96 *                   Carrier Phase at start FTF(L2)          Cycles
c    97                     Carrier Phase at start FTF(L2)          Cycles
c    98                     Carrier Phase at start FTF(L2)          Cycles
c    99                     Carrier Velocity at code FTF(L1)        Hz
c    100                    Carrier Velocity at code FTF(L2)        Hz
c    101                    Carrier Phase at end FTF(L1)            Cycles
c    102                    Carrier Phase at end FTF(L2)            Cycles
c    103                    Avg. Line-of-Sight Acceleration         m/s**2
c    104                    C/N0  L1                                0.5 dB/Hz
c    105                    C/N0  L2                                0.5 dB/Hz
c    106                    Code loop DLL sum bandwidth             Hz
c    107                    Code loop difference bandwidth          Hz
c    108                    Code loop PLL sum bandwidth             Hz
c    109                    Code loop PLL diff. bandwidth           Hz
c    110                    Narrowband AGC L1                       dB
c    111                    Narrowband AGC L2                       dB
c
c
c
c    Integer Items:
c
c         (Tracker 1)
c    1                      PRN                                     Dimensionless
c    2                      Code Indicator                          0=C/A,1=P
c    3 *                    Tracking Mode                           Dimensionless
c                            (See pg. A-14 in TI-4100 Owner's
c                            Manual for description)
c    4                      Antenna number                          Dimensionless
c    5                      Search pass count                       Dimensionless
c    6-9                    Unused
c    10  *                  Quality Vector L1                       Dimensionless
c                            (See pgs. A13-A14 in TI-4100 Owner's
c                            Manual for description)
c    11  *                  Quality Vector L2                       Dimensionless
c                            (See pgs. A13-A14 in TI-4100 Owner's
c                            Manual for description)
c         (Tracker 2)
c    12                     PRN                                     Dimensionless
c    13                     Code Indicator                          0=C/A,1=P
c    14  *                  Tracking Mode                           Dimensionless
c                            (See pg. A-14 in TI-4100 Owner's
c                            Manual for description)
c    15                     Antenna number                          Dimensionless
c    16                     Search pass count                       Dimensionless
c    17-20                  Unused
c    21  *                  Quality Vector L1                       Dimensionless
c                            (See pgs. A13-A14 in TI-4100 Owner's
c                            Manual for description)
c    22  *                  Quality Vector L2                       Dimensionless
c                            (See pgs. A13-A14 in TI-4100 Owner's
c                             Manual for description)
c         (Tracker 3)
c    23                     PRN                                     Dimensionless
c    24                     Code Indicator                          0=C/A,1=P
c    25                     Tracking Mode                           Dimensionless
c                            (See pg. A-14 in TI-4100 Owner's
c                            Manual for description)
c    26                     Antenna number                          Dimensionless
c    27                     Search pass count                       Dimensionless
c    28-31                  Unused
c    32 *                   Quality Vector L1                       Dimensionless
c                            (See pgs. A13-A14 in TI-4100 Owner's
c                             Manual for description)
c    33 *                   Quality Vector L2                       Dimensionless
c                            (See pgs. A13-A14 in TI-4100 Owner's
c                             Manual for description)
c         (Tracker 4)
c    34                     PRN                                     Dimensionless
c    35                     Code Indicator                          0=C/A,1=P
c    36 *                   Tracking Mode                           Dimensionless
c                            (See pg. A-14 in TI-4100 Owner's
c                            Manual for description)
c    37                     Antenna number                          Dimensionless
c    38                     Search pass count                       Dimensionless
c    39-42                  Unused
c    43 *                   Quality Vector L1                       Dimensionless
c                            (See pgs. A13-A14 in TI-4100 Owner's
c                             Manual for description)
c    44 *                   Quality Vector L2                       Dimensionless
c                            (See pgs. A13-A14 in TI-4100 Owner's
c                             Manual for description)

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
c          34             time of epoch                      GPS sec of week
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

c CORRESPONDANCE table Prescott Phaser, King Phaser, Feigl Fica

c mt, trakid,    satid,     af0,       af1,       af2,       xtoc
c mt, trakid,    NPRN,      XEAF0,     XEAF1,     XEAF2,     XETOC
c     ff(19),    ff(20),    ff(16),    ff(15),    ff(14),ff(13)

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

      integer*4   rcvr_chan
      parameter   (rcvr_chan=4)

c     passed values
C     .true. at end of file
      logical     fend,
C     .true. to print things
     &            debug,
C     .true. on error
     &            ferr,
C     .true. once we have a week num
     &            week_set
C     logical unit of file (assumed open)
      integer*4   lu,
C     1 for data, 2 for ephemeris
     &            iflag,
C     for broadcast ephem
     &            nprn,
C     sat id num PRN
     &            svid(rcvr_chan),
C     tracking mode
     &            tmode(rcvr_chan),
C     gps week number for epoch
     &            igpswk(rcvr_chan),
C     quality vector L1,L2
     &            iqvec(rcvr_chan,2)
C     gps sec of week for epoch
      real*8      gpssec(rcvr_chan),
C     L1 doppler phase (DATA!)
     &            dofl1(rcvr_chan),
C     L2 doppler phase (DATA!)
     &            dofl2(rcvr_chan),
C     L1 pseudorange
     &            prgl1(rcvr_chan),
C     L2 pseudorange
     &            prgl2(rcvr_chan),
C     signal to noise ratio L1,L2
     &            denrat(rcvr_chan,2),
C     broadcast ephemeris params
     &            bc_ephem(16),
C     broadcast clock params
     &            bc_clock(6)

c     other values
c
c     ephemeris parameters
      real*8           XEAF0, XEAF1, XEAF2, XETOC
     1                 ,XEADE, XECRS, XEDN, XEM0
     2                 ,XECUC, XEECC, XECUS, XEART
     3                 ,XETOE, XECIC, XECRC, XECIS, XEIDT
     4                 ,XEI0, XEOM0, XEW, XEOMD

      real*8       dsecs,tmcor,dummy,ftf_pseudo,ftf_bias,dt
      real*8       l1dop,l2dop

C     record type 1,2 or 9
      integer*4    mt,mt1,mt2,trakid
C     error count
      integer*4    errct
      integer*4    i,j,iday,nprn1,idummy,
C     day of gps week
     &             idogpsw
      integer*4    ios,
C     pointers into FICA arrays
     &             ptri,ptrf
      character*80 aline
      character*2  blkid

c
      real*8 ff(500)
C     offset from GPS time to UTC time
      real*8      utc_offset
      integer*4 fi(500)
      character*8 fc(500)
C     number of elements in array
      integer*4 num_f,num_i,num_c
      integer*4 blk_id
C     week number, must be obtained
      integer*4 week
C     from almanac (block 62)


c     for trick equivalencing of quality vectors.
c     warning, may be machine dependent
      integer*2 iqv2(2)
      integer*4 iqv4

c     functions
      real*8 get_pseudo_range

      INCLUDE 'makex.ins'
      equivalence (iqv4,iqv2(1))
      save week,week_set

C     speed of light in m per sec (WGS-72)
      data veloc/299792458.d0/

      lu = u_ficaf
      iflag = -1
      errct = 0
      ferr = .false.
      fend = .false.


 100  continue
      call read_fica(lu,blk_id,'BLK ',ff,fi,fc,
     &   num_f,num_i,num_c,ios)

C     end of file
      if (ios .eq. 16#FFFFFFFF) then
         fend = .true.
         return
      else if (ios .ne. 0 .and. errct .lt. 50) then
         write (6,105)
         if (q_batch) write (u_infor,105)
 105     format (x,'READ_FICA: 1 Bad line in FICA file: ')
         call handle_error(ios)
         errct = errct+1
         goto 100
      else if (errct.gt.50) then
         write (6,110)
         if (q_batch) write (u_infor,110)
 110     format (x,'READ_FICAF: REPEATED ERROR IN READ_FICAF: ')
         call handle_error(ios)
         ferr = .true.
         iflag = -1
         return
      endif

c     phase data are in block 6 or 55
      if (week_set .and. (blk_id .eq. 6 .or. blk_id .eq. 55)) then
c        signal to noise ratios and quality vectors
            ptrf = 4
            ptri = 9
            do 50  j=1,2
C     loop over 4 trackers
               do 40  i=1,4
                   denrat(i,j) = ff(ptrf)
c                  only the least significant 2 bytes
c                  can mean anything, so chop of most
c                  significant bytes via tricky equivalence
                   iqv4  = fi(ptri)
                   iqvec (i,j) = iqv2(2)
                   ptrf = ptrf + 1
                   ptri = ptri + 1
   40       continue
   50    continue

c        Get a valid time tag for the "frequency plan correction"
c        Ordinary FICA files have the time tag expressed as both:
c            1) FTF times   ff(1) + ff(2)
c            2) sec of week ff(3)
c        However, FICA files made by NGS2FIC have only (2), the time
c        tag in seconds of week.
c        To preserve compatibility, use the old way, unless
c        the values don't exist, in which case, get dt another way

c        old way
c        pseudorange FTF of validity
         ftf_pseudo = ff(1)

c        FTF bias offset
         ftf_bias   = ff(2)

c        more PHASER magic
         dt = ftf_pseudo + ftf_bias
c
c        new way for NGS2FIC'd files
         if (ftf_pseudo .eq. 0.0d0 .or.
     .       ftf_bias   .eq. 0.0d0 .or.
     .       dt         .le. 1.0e-20) then
            dt = ff(3)
         endif

         do 20 i=1,4
c           second of week, but week number isn't available!
c           if (ff(3) .lt. gpssec(1)) week = week + 1
            if (ff(3) .lt. gpssec(1)) week_set = .false.
            gpssec(i) = ff(3)
            igpswk(i) = week
c
c           tracker mode
            tmode(i) = fi(4+i)
c
c           sattelite id
            svid(i) = fi(i)
c
c           L1 pseudorange (meters)
            prgl1(i) = ff(11+i) * 1000.0d0
c
c           L2 pseudorange (meters)
            prgl2(i) = ff(15+i) * 1000.0d0
c
c           L1 doppler phase
            l1dop = ff(19+i)
c           + 6000 added for consistency with phaser
c           according to the "frequency plan"
            dofl1(i) = l1dop + 6000.d0 * dt
c
c           L2 doppler phase
            l2dop = ff(23+i)
c           according to the "frequency plan"
c           - 7600 added for consistency with phaser
            dofl2(i) = l2dop - 7600.d0 * dt
   20    continue
         iflag = 1

c     TI ROM phase data are in block 401
      else if (blk_id .eq. 401) then
         ptrf = 11
         ptri = 0
         do 140  i=1,4
c           week number
            igpswk(i) = int(ff(4))
            if (igpswk(i) .gt. 0 .and. igpswk(i) .lt. 900) then
               week_set = .true.
            else
               week_set = .false.
            endif

c           us_ftf   = longinteger_to_i4 ( user_buffer( 5) )
c           us_tov   = longinteger_to_i4 ( user_buffer( 7) )
c           code_ftf = longinteger_to_i4 ( user_buffer(10) )

c           code_tov = us_tov + (code_ftf-us_ftf)

            gpssec(i) = ff(3) + (ff(5) - ff(6))
            if (gpssec(i) .lt. ff(3)) then
               print *,'Panic: gpssec < user time',gpssec(i),ff(3)
               week_set = .false.
            endif

c           phase, pseudorange, C/N ratios and quality vectors
C           loop over L1 and L2
            do 150  j=1,2
c              carrier/noise ratios
               denrat(i,j) = ff(ptrf+17+j)

c              quality vector
c              only the least significant 2 bytes
c              can mean anything, so chop of most
c              significant bytes via tricky equivalence
               iqv4  = fi(ptri+9+j)
               iqvec (i,j) = iqv2(2)

 150        continue
c           L1 pseudorange
C           code tov
            prgl1(i) = get_pseudo_range(gpssec(i),
C                                       Z _count
     .                                  ff(ptrf+1),
C                                       X1 count
     .                                  ff(ptrf+2),
C                                       P phase
     .                                  ff(ptrf+3))
c           L2 pseudorange
C           code tov
            prgl2(i) = get_pseudo_range(gpssec(i),
C                                       Z _count
     .                                  ff(ptrf+4),
C                                       X1 count
     .                                  ff(ptrf+5),
C                                       P phase
     .                                  ff(ptrf+6))

c           L1 doppler phase
            dofl1(i) = ff(ptrf + 9)

c           L2 doppler phase
            dofl2(i) = ff(ptrf + 10)

c           sattelite id (PRN)
            svid(i) = fi(ptri+1)

c           tracker mode
            tmode(i) = fi(ptri+3)

            ptrf = ptrf + 25
            ptri = ptri + 11
 140     continue

         iflag = 1
      else if (blk_id .eq. 9) then
         week = int(ff(6))
         week_set = .true.
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
         bc_clock(5)=0.d0
C     xetdg Ionospheric correction of clock data. Not used.
         bc_clock(6)=0.d0
         iflag = 2
C     go read another record
      else
         goto 100
      endif

      if (debug) print *,'RDFICA: PRNs           ',(svid(i),i=1,4)
      if (debug) print *,'RDFICA: tracker modes  ',(tmode(i),i=1,4)
      if (debug) print *,'RDFICA: L1 phases      ',(dofl1(i),i=1,4)
      if (debug) print *,'RDFICA: L2 phases      ',(dofl2(i),i=1,4)
      if (debug) print *,'RDFICA: L1 pseudoranges',(prgl1(i),i=1,4)
      if (debug) print *,'RDFICA: L2 pseudoranges',(prgl2(i),i=1,4)

      return
      end

      real*8 function get_pseudo_range
     .       (code_tov,z_count,x1_count,p_phase)
cc---------------------------------------------------------------------------
c
c      subroutine get_pseudo_range( i2_data, pseudo_range, code_tov )
c
*     A Note on TI time counters:
*     The TI 4100 keeps count of the number of 20 millisec time frames
*     since it was turned on.  These counts are called fundamental time
*     frames (FTF).  There are three separate types of FTF counts.
*     USTOV -- the gps time since the start of the GPS week.  This
*              is a large value (order of days).
*     USFTF -- the FTF corresponding to USTOV.  This is a small value
*              (order of a 10 minutes to severval hours).
*     current FTF -- the FTF of the current operation.  This value is
*              of the same order as USFTF.


cc***************************************************************************
cc     Routine to extract and reconstruct the pseudo range from
cc     the data given in the TI 4100 data.
cc
cc     The I2_data contains:
cc     I2_data     Variable
cc     (words)
cc     1  (i*4)    Z-count units 1.5 secs
cc     3  (i*4)    X1 count (Pchips = 1/10.23 MHz)
cc     5  (i*2)    P_phase  (2.**-16 Pchips)
cc
cc     The Z-count, X1-count, and P_phase give the GPS time encoded
cc     on the satellite signal (measured from start of week).  The
cc     pseudo range is obtained by differencing this value from the
cc     current GPS time as recorded by the station clock (code_TOV)
cc
cc     Becuase of the large numbers involved the Z-count is differenced
cc     from the code_tov before the X1-count and the P_phase are added.
cc
cc     T.Herring                   12:32 PM  THU., 22  MAY , 1986
cc                                 last modified <860827.1639>
cc***************************************************************************
c
c      implicit none
c
c$include gps_parameters.FTNI
c
c      integer*2
c     .    i2_data(1)  ! the pseudo range data buffer
c
c      integer*4
c     .    code_tov    ! the GPS time of code measurement (units FTF's
c                      ! since start of GPS week.
c     .,   longinteger_to_i4   ! function to convert longinteger to I*4
c     .,   P_phase     ! the Pcode phase (2.**-16 Pchips)
c     .,   X1_count    ! the X1_count (Pchips)
c     .,   Z_count     ! the Z count (1.5 seconds)
c
c      real*8
c     .    msb_range       ! the difference of z-count and the code_tov
c                          ! (these are the large number parts of the
c                          ! epochs).
c     .,   pseudo_range    ! the reconstructed pseudo range (m)
c
c***** START, get the components of the pseudo range
c
c      Z_count  = longinteger_to_i4( i2_data(1) )
c      X1_count = longinteger_to_i4( i2_data(3) )
c
c***** Get p_phase and remove the two-complemnt
c      P_phase  = i2_data(5)
c      if( P_phase.lt.0 ) P_phase = 65536.d0 + P_phase
c
c***** Now reconstruct value
c
c      msb_range    =  (code_tov - Z_count*1.5d0*50.d0) * 20.d-3  ! seconds
c
c      pseudo_range =  ( msb_range - X1_count/Pchip_rate
c     .                            - P_phase/65536.d0/Pchip_rate  )
c     .                             * speed_of_light
c
c***** Thats all
c
c      return
c      end
c

      real*8
C     the difference of z-count and the code_tov
     .    msb_range
C     (these are the large number parts of the

C     epochs).

C     the reconstructed pseudo range (meters)
     .,   pseudo_range
C     the GPS time of code measurement (units FTF's
     .,   code_tov
C     since start of GPS week.

C     the Pcode phase (2.**-16 Pchips)
     .,   P_phase
C     reference frequency
     .,   Pchip_rate
C     the X1_count (Pchips)
     .,   X1_count
C     the Z count (1.5 seconds)
     .,   Z_count
C     speed of light in meters per sec (WGS-72)
     .,   speed_of_light

c     might want to check on these constants

C     speed of light in meters per sec (WGS-72)
      data velo/299792458.d0/
C     reference frequency
      data pchip_rate/10.23d06/



C     seconds
      msb_range    =  (code_tov - Z_count*1.5d0*50.d0) * 20.d-3

      pseudo_range =  (msb_range - X1_count/Pchip_rate
     .                           - P_phase/65536.d0/Pchip_rate  )
     .                           * speed_of_light

      get_pseudo_range = pseudo_range
      return
      end

