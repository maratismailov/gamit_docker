      subroutine blk401
     &   ( debug,iflag,svid,tmode,igpswk,gpssec
     &   , dofl1,dofl2,prgl1,prgl2,denrat,iqvec,week,week_set
     &   , ff,fi,fc,nf,ni,nc )

c
c   ??  assume all times in GPS time
c
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
c    5 *                    FTF of code phase                       Seconds
c    6 *                    FTF of carrier start phase              Seconds
c    7                      FTF Offset                              Seconds
c    8                      FTF interval of carrier phase           Seconds
c    9                      Internal temperature                    Deg. C
c   10                      Wideband AGC L1                         dB
c   11                      Wideband AGC L2                         dB
c        (Tracker 1)
c    12    1                Zcount for L1                           Seconds
c    13    2                X1 for L1                               Pchips
c    14    3                P-Phase for L1                          Pchips
c    15    4                Zcount for L2                           Seconds
c    16    5                X1 for L2                               Pchips
c    17    6                P-Phase for L2                          Pchips
c    18    7                Code Frequency for L1                   Hz
c    19    8                Code Frequency for L2                   Hz
c    20 *  9                Carrier Phase at start FTF(L1)          Cycles
c    21 * 10                Carrier Phase at start FTF(L2)          Cycles
c    22   11                Carrier Phase at start FTF(L2)          Cycles
c    23   12                Carrier Phase at start FTF(L2)          Cycles
c    24   13                Carrier Velocity at code FTF(L1)        Hz
c    25   14                Carrier Velocity at code FTF(L2)        Hz
c    26   15                Carrier Phase at end FTF(L1)            Cycles
c    27   16                Carrier Phase at end FTF(L2)            Cycles
c    28   17                Avg. Line-of-Sight Acceleration         m/s**2
c    29   18                C/N0  L1                                0.5 dB/Hz
c    30   19                C/N0  L2                                0.5 dB/Hz
c    31   20                Code loop DLL sum bandwidth             Hz
c    32   21                Code loop difference bandwidth          Hz
c    33   22                Code loop PLL sum bandwidth             Hz
c    34   23                Code loop PLL diff. bandwidth           Hz
c    35   24                Narrowband AGC L1                       dB
c    36   25                Narrowband AGC L2                       dB
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


      include '../includes/makex.h'
        
c     variables for error reporting
      integer*4 len,rcpar
      character*80 prog_name
      character*256 message

c     passed values
      logical     debug,
C     .true. once we have a week num
     &            week_set
C     logical unit of file (assumed open)
      integer*4   lu,
C     1 for data, 2 for ephemeris
     &            iflag,
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
     &            denrat(maxchn,2),
C     gps time of record transmission (sec) ?
     .            codsec
c

C     error count
      integer*4    i,j,iprn
C     pointers into FICA arrays
      integer       ptri,ptrf

      real*8
C     the difference of z-count and the codetov
     .    msbrange
C     (these are the large number parts of the

C     epochs).

C     the pcode phase (pchips)
     .,   pphase
C     the x1count (pchips)
     .,   x1count
C     the z count (seconds)
     .,   zcount
C     total delay (seconds)
     .,   tau

      logical isbad


c     the FICA arrays
      real*8 ff(maxfic)
      integer*4 fi(maxfic)
      character*8 fc(maxfic)
      integer*4 nf,ni,nc
C     week number, from block 9
      integer*4 week


c     There are three time tags in this business
c     1) start
c     2) code
c     3) end

c     all these are dimensioned for L1, L2
      real*8
C     the average acceleration (cycles/s**2)
     .    accel(2)
C     the time interval for shifting the phase
     .,   dtime
C     data.


C     phase at code FTF from end FTF (cycles)
     .,   phsc2(2)
C     phase at code FTF from start FTF (cycles)
     .,   phsc1(2)


C     the carries rate at code FTF (cycles/s)
     .,   ratec(2)
C     the carrier rate at end FTF (cycles/s)
     .,   rate2(2)
C     the carrier rate at start FTF (cycles/s)
     .,   rate1(2)
C     rate of change in range
     .,   delrng(2)
C     Pseudo Range at C (meters)
c *    .,   pratc(2)   --not used

C     d(ion)/dt in m**-1/ s**-1
     .,   dotion
C     d(range)/dt in m/s
     .,   dotrng



c     for trick equivalencing of quality vectors.
c     warning, may be machine dependent
      integer*2 iqv2(2)
      integer*4 iqv4

      equivalence (iqv4,iqv2(1))

c     Official WGS-72 speed of light in meters per second
      real*8 vlight
      parameter (vlight = 299792458.d0)

c     Reference freqency
      real*8 fpchip
      parameter (fpchip = 10.23d06)

c     L1 and L2 frequencies in Hz
      real*8 frq1,frq2
      parameter (frq1  = 154.0d0 * fpchip)
      parameter (frq2  = 120.0d0 * fpchip)

      real*8 wavel1,wavel2
      parameter (wavel1 = vlight/frq1)
      parameter (wavel2 = vlight/frq2)
                
c     wide lane factor
C     facwl = [f(L1)-f(L2)]/[f(L1)+f(L2)] = 17/137
      real*8 facwl
      parameter (facwl = 17.0d0/137.0d0)

c     gear ratio g = 60/77
      real*8 gear
      parameter (gear  = 60.0d0/77.0d0)
        
c     get the calling program name
      len = rcpar(0,prog_name)

      lu = uficaf
           

*    Dummy statements for unused calling arguments     
      if( debug) then
        print *,'nc ni nf ',nc,ni,nf
        print *,'fc',(fc(i),i=1,nc)     
        print *,'week ',week
       endif




c     TI ROM phase data are in block 401
      ptrf = 11
      ptri = 0
      do 140  i=1,4
         isbad = .false.
c        week number
         igpswk(i) = int(ff(4))
         if (igpswk(i) .gt. 0 .and. igpswk(i) .lt. 900) then
            week_set = .true.
         else
            week_set = .false.
         endif

C        TIME TAGS
c        Code time (0.00 seconds)
         codsec = ff(3) + (ff(5) - ff(2))

c        this should end in .08 and be the time at the end of
         gpssec(i) = ff(3) + (ff(5) - ff(2)) + ff(8)/2.d0

         if (gpssec(i) .lt. ff(3)) then
            write(message,'(a,2f12.2)') 'gpssec < user time '
     .                                  , gpssec(i),ff(3)
            write(uinfor,'(a)') message 
           call report_stat('WARNING',prog_name,'makex/blk401',' '
     .                     ,message,0)
            week_set = .false.
         endif

c        difference in seconds between code FTF and start FTF
         dtime = ff(8)/2.0d0
CD        print *,'BLK401: Delta time ',dtime

c        PHASE DATA
c        We must move the interpolate the phase to the code epoch.
c        To do this, get rate and acceleration terms.
c        Carrier rate at code FTF (unbiased).
         ratec(1) = ff(ptrf+13) + 6000.0d0
         ratec(2) = ff(ptrf+14) - 7600.0d0

c        acceleration in cycles per second
         accel(1) = ff(ptrf+17) / wavel1
         accel(2) = ff(ptrf+17) / wavel2

c        seconds between start and code, code and end
         dtime = ff(8)/2.0d0

c        rate at start FTF
         rate1(1) = ratec(1) - accel(1) * dtime
         rate1(2) = ratec(2) - accel(2) * dtime

c        rate at end FTF
         rate2(1) = ratec(1) + accel(1) * dtime
         rate2(2) = ratec(2) + accel(2) * dtime

c        calculate phase at code FTF based on phase at start FTF
         phsc1(1) = ff(ptrf+9) + 6000.0d0 * (gpssec(i) - 2.0d0 * dtime)
     .            + rate1(1) * dtime
     .            + accel(1) * dtime**2/2.0d0
         phsc1(2) = ff(ptrf+10) - 7600.0d0 * (gpssec(i) - 2.0d0 * dtime)
     .            + rate1(2) * dtime
     .            + accel(2) * dtime**2/2.0d0

c        calculate phase at code FTF based on phase at end FTF
         phsc2(1) = ff(ptrf+15) + 6000.0d0 * gpssec(i)
     .            - rate2(1) * dtime
     .            + accel(1) * dtime**2/2.0d0
         phsc2(2) = ff(ptrf+16) - 7600.0d0 * gpssec(i)
     .            - rate2(2) * dtime
     .            + accel(2) * dtime**2/2.0d0

c        compare the two
         if (dabs(phsc2(1) - phsc1(1)) .gt. 5.0d-2 .and.
     .       dabs(phsc2(2) - phsc1(2)) .gt. 5.0d-2 ) then
            isbad = .true.  
            iprn = i*11 - 10
            write(message,'(a,i2,a,f12.2)' ) 
     .           'Phase problem for PRN ',fi(iprn),' at ',gpssec(i) 
            call report_stat('WARNING',prog_name,'makex/blk401',' '
     .                   ,message,0)  
            write(uinfor,'(a,f12.2,a,i3,a,f6.3,4f20.4)') 
     .           'Phase problem at ',gpssec(i),'  PRN ',fi(iprn)
     .          ,' Values: ',dtime,phsc1(1),phsc2(1),phsc1(2),phsc2(2)
         endif

c        put the rate back in according to the "frequency plan"
c        This phase should be increasing when the pseudorange is.

         dofl1(i) = ff(ptrf+15) + 6000.0d0 * gpssec(i)
         dofl2(i) = ff(ptrf+16) - 7600.0d0 * gpssec(i)

c        C/N ratios and quality vectors
         do 150  j=1,2
c           carrier/noise ratios
            denrat(i,j) = ff(ptrf+17+j)

c           quality vector
c           only the least significant 2 bytes
c           can mean anything, so chop of most
c           significant bytes via tricky equivalence
c           If there is a problem decoding the phase,
c           flag it as a machine problem.
            iqv4  = fi(ptri+9+j)
            if (.not. isbad) then
               iqvec (i,j) = iqv2(2)
            else
               iqvec (i,j) = -1
            endif
 150     continue

c        Using the rate of change in the phase, estimate the rate of change
c        in the pseudorange.  Because of the ionoshpere, the relationship
c        is NOT just del(range) = - d(phase) * wavelength
c        We are also ignoring the acceleration term.
c        This calculation assumes the Doppler convention for signs,
c        i.e. wavelength * d(phase)/dt = - d(range)/dt
c
c        Estimate the rate of change of the true range (m/s).
         dotrng = (ratec(1)/wavel1 - ratec(2)/wavel2)
     .          / (1.0d0/wavel2**2 - 1.0d0/wavel1**2)

c        Estimate the rate of change of the ionospheric phase delay (1/m * 1/s)
         dotion = (wavel1*ratec(1) - wavel2*ratec(2))
     .          / (wavel1**2 - wavel2**2)

c        Finally, estimate the rates of change in the 2 pseudodranges
         delrng(1) = dotrng + wavel1**2 * dotion
         delrng(2) = dotrng + wavel2**2 * dotion


c        L1 PSEUDORANGE
         zcount  = ff(ptrf+1)
         x1count = ff(ptrf+2)
         pphase  = ff(ptrf+3)
         msbrange  = codsec - zcount

c        calculate total delay in seconds
         tau =   msbrange - x1count/fpchip - pphase/fpchip

c        reconstruct pseudorange in meters
         prgl1(i) = tau * vlight + dtime * delrng(1)

         if (i .eq. 1) then
CD           print *,'BLK401: P1 at Code ',pratc(1)
CD           print *,'BLK401: P1 delrng  ',dtime * delrng(1)
         endif


c        L2 PSEUDORANGE
         zcount  = ff(ptrf+4)
         x1count = ff(ptrf+5)
         pphase  = ff(ptrf+6)
         msbrange    =  codsec - zcount

c        calculate total delay in seconds
         tau =   msbrange - x1count/fpchip - pphase/fpchip

c        L2 PSEUDORANGE
         zcount  = ff(ptrf+4)
         x1count = ff(ptrf+5)
         pphase  = ff(ptrf+6)
         msbrange  = codsec - zcount

c        calculate total delay in seconds
         tau =   msbrange - x1count/fpchip - pphase/fpchip

c        reconstruct pseudorange in meters
         prgl2(i) = tau * vlight + dtime * delrng(2)

         if (i .eq. 1) then
CD           print *,'BLK401: P2 at Code ',pratc(2)
CD           print *,'BLK401: P2 delrng  ',dtime * delrng(2)
         endif


c        tracker mode
         tmode(i) = fi(ptri+3)

c        sattelite id (PRN)
         svid(i) = fi(ptri+1)

c        print some quantities of interest:
         if (i .eq. 1) then
CD           print *,'BLK401: L1 rate (e)',(ff(ptrf+15)-ff(ptrf+ 9))
CD    .                                   /(2.0*dtime)
CD           print *,'BLK401: L1 rate (m)',ratec(1)
CD           print *,'BLK401: L2 rate (e)',(ff(ptrf+16)-ff(ptrf+10))
CD    .                                   /(2.0*dtime)
CD           print *,'BLK401: L2 rate (m)',ratec(2)
CD           print *,'BLK401: L1 P1      ',dofl1(i),prgl1(i)
CD           print *,'BLK401: L2 P2      ',dofl2(i),prgl2(i)
CD           print *,'BLK401: P diff (c) ',pratc(2) - pratc(1)
CD           print *,'BLK401: P diff (2) ',prgl2(i) - prgl1(i)

c           widelane number
CD           widel = dofl2(i) - dofl1(i)
CD    .            + facwl * (prgl1(i)/wavel1 + prgl2(i)/wavel2)
CD           print *,'BLK401: Wide lane  ',widel
         endif


c        update pointers to begining of items for next tracker
         ptrf = ptrf + 25
         ptri = ptri + 11
 140  continue

      iflag = 1

      return
      end





