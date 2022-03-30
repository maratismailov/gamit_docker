      subroutine settim( debug,rcvrsw,swver,iwkn0,sow0,slop)

c     Calculate the start time using the receiver type and version
c     R. W. King   20 February 1989

      implicit none

      character*3 rcvrsw
      character*512 message
      integer*4 iwkn0,jd,julday,igpsdow
      real*4 swver
      real*8 sow0,utcoff,taiutc,slop
      logical debug
c*     number of channels in the receiver -- no longer used, 
c*    just check for nprn > maxchn   RWK 151015
c*      integer*4  nchan

c     case conversion functions:
      character*3 lowerc

      include '../includes/makex.h'

c     Set default slop in setting epoch search window
          
c rwk 040214: This changed from 1 ms to 100 ms for safety
c      slop = 0.001d0                                     
      slop = 0.1d0

      if(debug) then  
        print *,'SETTIM: rcvrsw swver ',rcvrsw,swver
        print *,'SETTIM: (1) iwkn0 = ', iwkn0
        print *,'SETTIM: (1) sow0  = ', sow0
        print *,'SETTIM: (1) slop   = ', slop
        endif

c     Set parameters for the TI 4100
c     ******************************

c     The TI 4100 interprets the input start time as GPS time, so
c     we need only to subtract the small (<0.3 s) offset at which
c     TI samples its data. For GESAR versions 1.0 and later and most
c     versons of CORE, this is 0.92 seconds;  for GESAR V 25.2, it
c     is 0.26 s; for certain weird case of core 4.1 (here designated
c     CORE 4.11), the offset is 0.
c     e.g., the requested start of 7 hrs 58 min on day 353 of 1986
c     has a TI receiver time of 28680.08 s GPS seconds of day
c     but 28675.08 s UTC.
c       CORE 4.12 is an arbitrary version number for the incorrectly
c     implemented frequency plan of the Yellowknife wk 533 data.
c     The time-tags are assumed to be the same as 4.1 -- mhm 900705

      if (lowerc(rcvrsw).eq.lowerc('COR')) then
cc         nchan = 4
         if (nint(10*swver) .eq. 252) then
C           subtract 0.26 seconds
            call secsum(iwkn0,sow0,-0.26d0,iwkn0,sow0)
         else if (nint(10*swver) .eq. 48) then
C           add 0.08 seconds
            call secsum(iwkn0,sow0,+0.08d0,iwkn0,sow0)
         else if (nint(100*swver).eq. 411) then
C           subtract 1.0 seconds
            call secsum(iwkn0,sow0,-1.0d0,iwkn0,sow0)
         else if (nint(100*swver).eq. 413) then
C           subtract 1.0 seconds and allow 59.0 and 59.08
            slop = 0.1d0
            call secsum(iwkn0,sow0,-1.0d0,iwkn0,sow0)
         else if (nint(10*swver).eq.41 .or. nint(10*swver).eq.57) then
c           subtract 0.92 seconds
            call secsum(iwkn0,sow0,-0.92d0,iwkn0,sow0)
         else
            write (uinfor,100) rcvrsw,swver
 100        format(1x,'SETTIM WARNING: unknown CORE version ',a3,f10.4,
     .   /,1x,'Current known versions: v. 25.2,4.1,4.11,4.12,4.81,5.7'
     .  ,/,1x,'Assuming that data are sampled like GESAR v. 1.0')
            write (message,105) rcvrsw,swver
 105        format('Unknown CORE version ',a3,f10.4,
     .      ' Current known versions: v. 25.2,4.1,4.11,4.12,4.81,5.7',
     .      ' Assuming that data are sampled like GESAR v. 1.0')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
            call secsum(iwkn0,sow0,-0.92d0,iwkn0,sow0)
         endif

      else if (lowerc(rcvrsw).eq.lowerc('GES')) then
cc         nchan = 4
         if(nint(10*swver) .eq. 10 .or.
     .      nint(10*swver) .eq. 11 .or.
     .      nint(10*swver) .eq. 12 .or.
     .      nint(10*swver) .eq. 13 .or.
     .      nint(10*swver) .eq. 14 .or.
     .      nint(10*swver) .eq. 15 .or.
     .      nint(10*swver) .eq. 16 .or.
     .      nint(10*swver) .eq. 19) then
            call secsum(iwkn0,sow0,-0.92d0,iwkn0,sow0)
         else
            write (uinfor,110) rcvrsw,swver
 110        format(1x,'SETTIM WARNING: unknown GESAR version ',a3,f10.4,
     .         /,1x,'Current known versions: v. 25.2, 1.0-1.5',/,
     .           1x,'Assuming that data are sampled like GESAR v. 1.0')
            write (message,115) rcvrsw,swver
 115        format('Unknown GESAR version ',a3,f10.4,
     .           ' Current known versions: v. 25.2, 1.0-1.5',
     .           ' Assuming that data are sampled like GESAR v. 1.0')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
            call secsum(iwkn0,sow0,-0.92d0,iwkn0,sow0)
         end if

      else if (lowerc(rcvrsw).eq.lowerc('ROM')) then
cc          nchan = 4
c         version 1.11 artificial designation for TISTRX translation
         if (nint(100*swver).eq. 111) then
C           subtract 1.0 seconds
            call secsum(iwkn0,sow0,-1.0d0,iwkn0,sow0)
         else
C           ROM is just like GESAR subtract 0.92 seconds to end up at 59.08 GPST
            call secsum(iwkn0,sow0,-0.92d0,iwkn0,sow0)
            write (uinfor,120)
 120        format
     . (1x,'SETTIM: all ROM except 1.11 assumed sampled like GESAR')
            write (message,125)
 125        format('All ROM except 1.11 assumed sampled like GESAR')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
         endif     

c     Set Parameters for LEICA's
c     **************************
C I have almost no idea of Leica firmware histories, therefore:
c Calling Leica Series 200 single-frequency receivers  v.2.00 -- YB  7/11/99
c Calling Leica Series 200 dual-frequency receivers    v.2.50 -- YB  7/11/99
c   revised to include 2.10 (actual firmware)                    RWK 12/8/09
c Calling Leica Series 300 dual-frequency receivers    v.3.50 -- YB  7/11/99
c Calling Leica Series 9400 single-frequency receivers v.4.00 -- YB  7/11/99
c Calling Leica Series 9500 dual-frequency receivers   v.4.50 -- YB  7/11/99 
c Calling Leica Series 500 single-frequency receivers  v.5.00 -- YB  7/11/99
c Calling Leica Series 500 dual-frequency receivers    v.5.50 -- YB  7/11/99
c Calling Leica MC1000 machine control receivers       v.6.00 -- YB  7/11/99
c Calling Leica CRS1000 control station receivers      v.6.50 -- YB  7/11/99 
c Calling Leica GRX1200 (LC1200) firmware CC00         v 7.00 -- RWK 8/23/04   
c Calling Leica GRX1200 (LC1200) fimware 2.14/2.121    v.7.01 -- RWK  3/7/06
c Calling leica GX1230 (LC1230) firmware 1.35          v.7.10 -- RWK 11/9/04 
c Calling Leica GRX1200PRO (LC12PR) and GRX1200GGPRO 
c      (LC14PR) firmware 5.62                          v.7.20 -- RWK  7/2/09
c The GRX1200GGPRO also has 7.50-7.80 naturally        v.7.50 -- RWK 3/30/10
c      and 2.14/2.121                                  v.7.30 -- RWK  7/8/10
c The GMX902 (labeled MC1000 in the header) has 8.03   v 8.03 ---RWK  3/7/07
c The Lecia GR 25 has 120 channels; firmware 3.30      v 9.03 ---RWK 10/9/15
       else if (lowerc(rcvrsw).eq.lowerc('LEI')) then  
         if (nint(100*swver).lt.210) then
cc            nchan = 6 
            slop = 1.d0            
         else if(nint(100*swver).ge.250.and.nint(100*swver).lt.500) then
cc            nchan = 9
            slop = 1.d0
         else if(nint(100*swver).ge.500) then
cc            nchan = maxchn                   
c More stable oscillators          
c           change to be more conservative
c            slop = 0.01d0                
            slop = 0.1d0
         else
            write (uinfor,150) rcvrsw,swver
 150        format(1x,'SETTIM WARNING: unknown LEICA version ',
     . a3,f10.4,/,1x,
     .'See list in settim.f'
     . ,/,1x,'Assume sampling like LEICA 6.0')
            write (message,155) rcvrsw,swver
 155        format('Unknown LEICA version ',a3,f10.4,
     . ' Current known versions:  2.0-8.03',
     . ' Assuming that data are sampled like LEICA 6.5')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
cc            nchan = 12 
c**         change this to be conservative
c            slop = 0.01d0
            slop = 0.1d0

         endif


c     Set Parameters for Macrometers and Mini-Macs
c     ********************************************

      else if (lowerc(rcvrsw).eq.lowerc('MAC')) then
cc         nchan = 6
c        The MACROMETER II interprets the start time as UTC, so we need to
c        add the GPST-UTC difference as well as the -1.525 s offset from
c        the even (UTC) minute.
c        calculate GPST - UTC
         igpsdow = int(sow0/86400.d0)
         jd  = julday(1,5,80) + 7*iwkn0 + igpsdow
c        get the UTC offset according to JD in GPS time
         utcoff = taiutc(jd) - 19.d0
         call secsum( iwkn0,sow0,utcoff-1.525d0, iwkn0,sow0 )

      else if (lowerc(rcvrsw).eq.lowerc('MIN')) then
cc         nchan = 8
c        Minimac data comes from NGS in GPST
c        However, there is a millisecond timing problem in the
c        acquisition software MMAT v. 1.49 such that the data
c        are really acquired at xx.001 seconds
c        Input RINEX files (through Werner Gurtner's translator)
c         have had the proper time conversions for all Mini-Mac versions

         if (nint(100.0d0 * swver) .eq. 149) then 
c**   rwk 010815:  I think this is wrong; since the RINEX files have the correct
c**                time tag (xxx.001), we should be searching around that tag
c**           if(.not.qrinex)
c**    .       call secsum( iwkn0,sow0,+1.0d-3, iwkn0,sow0 )  
c** change to this:
             call secsum( iwkn0,sow0,+1.0d-3, iwkn0,sow0 )

         else if (nint(100.0d0 * swver) .eq. 150) then
c           Special GOTEX Tsukuba version for Peter Morgan
            if(.not.qrinex)
     .       call secsum( iwkn0,sow0,+5.001d00, iwkn0,sow0 )
         else if (nint(100.0d0 * swver) .ge. 161 .and.
     .            nint(100.0d0 * swver) .le. 164) then
            continue
         else if (nint(100.0d0 * swver) .eq. 159 .or.
     .            nint(100.0d0 * swver) .eq. 189 ) then
c            RINEX files from Germany, sampling is 59.001
c            1.89 is mistagged 59.000 (but corrected in MAKEX immediately
c            after call to RRINEX); 1.59 is ok after XTORX
             call secsum( iwkn0,sow0,-0.999d0, iwkn0,sow0 )
        else if (nint(100.0d0 * swver) .eq.196 ) then
c            1.96 is artificial for 1.95 but with a 3.0 s offset from UTC,
c            making a 5.0 s offset from GPST (data from China 1993)
             call secsum( iwkn0,sow0,5.0d0,iwkn0,sow0 )
          else
            write (uinfor,200) rcvrsw,swver
 200     format(1x,'SETTIM WARNING: unknown MINI-MAC version ',a3,f10.4,
     .      /,1x,'Current known versions: 149,150,159,161-164,189,196',/
     .         ,1x,'Assuming that data are sampled like version 164')
            write (message,205) rcvrsw,swver
 205     format('Unknown MINI-MAC version ',a3,f10.4,
     .         ' Current known versions: 149,150,161-164,189,196',
     .         ' Assuming that data are sampled like version 164')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
         endif

c     Set Parameters for the Trimble
c     ******************************
                                             
c        Info from Brian Frohring of Trimble Navigation, e-mail to rwk 28 Mar 97
c
c        Receiver   First firmware version   Last firmware version
c        --------   ----------------------   ---------------------
c       4000SLD                                 3.31A 27-May-93
c       4000SST                                 4.83  02-Sep-93
c       4000SSE      5.50 30-Jul-92             6.01 ??
c       4000SSi      7.00    
c       40000MSGR    7.00    
c       7400MSI      2.21?   
c       ALLOY        5.33 

c       For the BD-750 (MSAG), treat as similar to 7400, as per e-mail of 20 March 2002 to rwk
c       from Warren Gallaher of UNAVCO (see http://www.trimble.com/bd750_ts.asp?Nav=Collection-4891

c       For the "Total station 4700", (NAV) firmware versions found by Matt van Domselaar   

c       For the 5700 and NetRS, use the 1.xx as for the 4700, as per email from Paul Jamason 2 May 2002 
c       For the MS750, also use 1.xx, as per email from Martin Lidberg 27 Jan 2003    
c       For the NETRS the early firmware is 0.7-0, so map this into 0.70  
c       The 4600LS single frequency receiver seems to use 2.50, which might conflict with the
c          firmware for the 7500MSI, Bd-750, and MS-750; map this artificially into 2.99 
c       The later NETR[S,3,4,5,8,9] receivers use 4.xx, so allow this to maintain.  This means
c          treating early 4.xx firmware (ST/SST) the same as the NETRS with respect to 
c          channdel. We still reserve the (artificial) 4.00 for the 'serial P-code' (ST) receiver 
c          (see lib/settyp.f)
c
c      In general, make the following assumptions:
c        3.xx     -> SL/SLD
c        4.xx     -> ST/SST
c        4.00     -> serial P-code ST
c        [567].xx -> SE/SSE
c        7.xx     -> SSi, MSGR
c        1.xx     -> 4700,  5700, NETRS, BD-750, MS-750                                  
c        2.xx     -> 7400MSI, BD-750, R7
c        0.xx     -> NETRS           
c        4.xx     -> NETRS, NETR3, NETR4, NETR5, NETR8, NETR9, TRM852, TRM4-3 
c        5.xx     -> ALLOY 

      else if (lowerc(rcvrsw).eq.lowerc('TRM')) then
cc         nchan = 12
c         write(6,*) 'swver',swver    
         if( nint(100*swver).eq.299 ) then
c           artificial code for L1-only TR4600LS, override obervable types in settyp.f
            slop = 0.3d0 
c           I'm guessing on this; maybe only 9 ---rwk 070607
cc            nchan = 12
         elseif (nint(100*swver).ge.312 .and. 
     .         nint(100*swver).le.758 ) then
            call secsum(iwkn0,sow0,-0.00d0,iwkn0,sow0)
            slop = 0.300d0
c           need to increase to take into account poorly tuned oscillators
c           with long tracks
            slop = 1.d0                  
cc            nchan = 24
cc rwk commented out 101119
cc         elseif (nint(100*swver).gt.70.and.nint(100*swver).lt.300) then
cc           for 4700, 5700, 7400, BD-750 (MSAG), or NETRS  
cc           --use this also, tentatively, for the MS-750, though it may have only 9 channels
cc            nchan = 12
cc            slop = 0.300d0
         else
            write (uinfor,300) rcvrsw,swver
 300        format(1x,'SETTIM WARNING: unknown TRIMBLE version ',
     . a3,f10.4,/,1x,
     .'Current known versions: 0.7 1.xx, 2.xx, 3.12-3.25,3.28,4.10-7.29'
     . ,/,1x,'Assuming that data are sampled like TRIMBLE 4.64')
            write (message,305) rcvrsw,swver
 305        format('Unknown TRIMBLE version ',a3,f10.4,
     .     ' Current known versions: 1.x,2.x,3.x,4.10-7.29'
     .     ,' Assuming that data are sampled like TRIMBLE 4.64')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
c           receiver resets whenever it thinks it has drifted more than
c           256 milliseconds off, so allow 300 ms
            slop = 0.300d0
         endif
                

c     Set parameters for the Ashtech
c     *******************************

      else if (lowerc(rcvrsw).eq.lowerc('ASH') ) then 
c
c---------NOTE added by rwk 021008---------------------------------  

c    The comments below developed over time, according to 
c    our understanding of the Ashtech firmware.  It appears that
c    some of these designations are misleading, so we now recommend
c    using only the following:
c
c      L-XII, M-XII, LM-XII, LM-XII0C, all codeless, with firmware 6Axx, 7Axx 
c        except for 7AP6:  ASHL12  1.0
c                                                               
c      L-XII3, code-tracking (P12) with firmware 7A26 (or 7A??)         :  ASHM12  2.0
c      L-XII3, code-tracking (P12) with firmware 6Cxx, 6Gxx, 6Ixx, 6Mxx :  ASHM12  6.0  
c
c      GG-XXIV single-frequency GPS/Glonass with firmware GM00          :  ASHGGS  3.0 
c
c      P-XII3 (P12) with firmware 7Bxx:  ASHP12 7.20 
c
c      Z-12:  ASHZ12 8.0  (or more specific designations listed below) 
c
c -----------------------------------------------------------------------------
c
c I am calling Ashtech codeless data v 1.00 -- YB
c I am calling the P-code Ashtechs v, 2.00 -- YB  
c    rwk 021008:  It's not clear what this firmware was meant to correspond
c      to since not until the beta version of the P-12 (firmware 6G, 6M)
c      did we appear to have a true P-code capability without biases. 
c I am calling the single-frequency GG24 firmware GM00 3.0 -- rwk 2002/12/20   
c I am calling the beta version P-code Ashtech v 6.00 -- YB   
c  (this version corrects for interchannel biases & can track
c   up to 7-8 satellites simultaneously, first used in SIO1 to SIO2
c   tests in June 1992, at Vandenberg PGGA, Palos Verdes May 1992,
c   and Indonesia 1992)
c I am calling version 7B P-code Ashtech v.7.20 -- YB
c  (this version will automatically switch to codeless tracking under AS,
c   and can track up to 12 satellite channels - received firmware on July 7, 1993))
c Line added by G. Ferhat to make 7A = 7.10.  RWK  940412) 
c    rwk 021008: The '7Axx' firmware, except for 7AP6, was used with codeless
c      receivers (ASHL12), and should be designated 1.0, not 7.1.
c      I don't know about the 7C and 8M, listed below. 
c I am calling version 7C P-code Ashtech v.7.30 -- YB 10/15/99
c I am calling version 8M P-code Ashtech v.7.80 -- YB 10/15/99
c I am calling Z-12 firmware 1B10  v.7.90 -- YB 10/15/99
c I am calling Z-12 firmware 1B12  v.7.92 -- YB 10/15/99
c I am calling Z-12 firmware 1C00  v.8.00 -- YB
c I am calling Z-12 firmware 1C01  v.8.01 -- YB
c I am calling Z-12 firmware 1C63  v.8.02 -- YB 10/15/99
c I am calling Z-12 firmware 1C75  v.8.03 -- YB 10/15/99
c I am calling Z-12 firmware 1C80  v.8.04 -- YB 10/15/99
c I am calling Z-12 firmware 1C85  v.8.05 -- YB 10/15/99
c I am calling Z-12 firmware 1D00  v.8.10 -- YB
c I am calling Z-12 firmware 1E00  v.8.20 -- YB  
c I am calling Z-12 firmware 1E24  v.8.22 -- YB 
c I am calling Z-12 firmware 1E76  v.8.23 -- YB  
c I am calling Z-12 firmware 1E8D  v.8.24 -- rwk
c I am calling Z-12 firmware 1E95  v.8.25 -- YB
c I am calling Z-12 firmware 1E81  v.8.26 -- YB 10/15/99
c I am calling Z-12 firmware 1E82  v.8.27 -- YB 10/15/99
c I am calling Z-12 firmware 1E86  v.8.28 -- YB 10/15/99
c I am calling Z-12 firmware 1E97  v.8.29 -- YB 10/15/99
c I am calling Z-12 firmware 1F00  v.8.30 -- YB  
c I am calling Z-12 firmware 1F39  v.8.32 -- YB  
c I am calling Z-12 firmware 1F50  v.8.35 -- YB
c I am calling Z-12 firmware 1F60  v.8.36 -- YB
c I am calling Z-12 firmware 1G00  v.8.40 -- rwk
c I am calling Z-12 firmware 1H00  v.8.50 -- YB 10/15/99
c I am calling Z-12 firmware 1I00  v.8.60 -- rwk
c I am calling Z-12 firmware 1J00  v.8.70 -- YB
c I am calling Z-12 firmware 2J00  v.8.72 -- YB
c I am calling Z-12 firmware 3J00  v.8.73 -- YB
c I am calling Z-12 firmware 4J00  v.8.74 -- YB
c I am calling Z-12 firmware 1K00  v.8.80 -- YB 10/15/99
c I am calling Z-12 firmware 1L00  v.8.85 -- YB 10/15/99   
c I am calling Z-12 firmware 1SWE1D0 v.8.86 -rwk 1/16/04    
c I am calling Z-12 firmware 1L01  v.8.87 ---rwk06/23/04
c I am calling Z-12 firmware 1Y04  v.8.90 -- YB 10/15/99
c I am calling Z-12 firmware 1Y05  v.8.91 -- YB 10/15/99
c I am calling Z-12 firmware 1Y06  v.8.92 -- YB 10/15/99  
c I am calling Z-12 firmware 1Y07  v.8.93 -- rwk 9/16/03  
c I am calling Z-12 firmware CA00  v.8.95 -- YB 10/15/99
c I am calling Z-12 firmware CB00  v.9.00 -- YB
c I am calling Z-12 firmware CC00  v.9.10 -- YB
c I am calling Z-12 firmware CJ00  v.9.15 -- rwk 2016/04/21 
c I am calling Z-12 firmware CD00  v.9.20 -- YB
c I am calling Z-12 firmware RC00  v.9.40 -- YB 10/15/99
c I am calling Z-12 firmware RC05  v.9.45 -- YB 10/15/99
c I am calling Z-12 firmware RD00  v.9.50 -- YB 10/15/99
c I am calling UZ-12 (Z-FX) firmware UB00 v.9.70 -- YB
c I am calling UZ-12 (Z-FX) firmware UC00 v.9.71 -- YB 10/15/99
c I am calling UZ-12 (Z-FX) firmware UD00 v.9.72 -- YB 10/15/99
c I am calling UZ-12 (Z-FX) firmware UE00 v.9.73 -- YB 10/15/99 
c I am calling the single-frequency GG24 firmware GM00 3.0 -- rwk 2002/12/20   
c I am calling UZ-12 (Z-FX) firmware UF00 v 9.74 -- rwk 2003/09/16   
c I am calling UZ-12 (Z-FX) firmware UG00 v 9.75 -- rwk 2003/09/16   
c I am calling Z-18  firmware V5100 (alias 0051)      v.9.77 -- rwk 2003/09/16  
c I am calling Z-18  firmware V5100 (alias 0057)      v.9.78 -- rwk 2003/09/16   
c I am calling Z-18  firmware V6000 (alias 0060)      v.9.79 -- rwk 2003/09/16   
c I am calling Z-18  firmware V6400 (alias 0064) v.9.80 -- rwk 2001/01/09 
c I am calling Z-18  firmware V6500 (alias 0065) v.9.81 -- rwk 2002/08/18   
c I am calling Z-XII3T firmware IL01-1D04-MCF-12MX v 9.90 -- rwk 2002/09/18   
c I am calling UZ-12 (Z-FX) firmware ZB00 v 9.91 -- rwk 2002/12/04   
c I am calling UZ-12 (Z-FX) firmware ZC00 v 9.92 -- rwk 2003/09/16   
c I am calling UZ-12 (Z-FX) firmware CJ00 v.9.93 -- rwk 2003/09/16     
c I am calling UZ-12 (Z-FX) firmware CN00 v.9.94 -- rwk 2003/09/16        
c I am calling UZ-12 (Z-FX) firmware ZE00 v.9.95 -- rwk 2003/09/16    
c I am calling UZ-12 (Z-FX) firmware ZE21 v.9.96 -- rwk 2003/09/16        
c I am callling the Z-Max firmware Mxxx   v. 9.96---rwk 2007/01/18
c I am callling Z-X firmware ZB00         v.9.97 -- rwk 2010/07/08 
c I am calling the Z-Max firmware 0A01    v. 9.98 --rwk 2013/6/5
c I am calling the PF500 firmware s763Gx24 v.10.01 --rwk 2013/6/5  
c   From C. Vigny:
c       M014 released June 2 2004
c       MB00 released November 2 2004
c       MC00 released May 18 2005
c       MD00 released June 6 2006

         if (     nint(100*swver) .eq. 100
     .       .or. nint(100*swver) .eq. 200  
     .       .or. nint(100*swver) .eq. 300
     .       .or. nint(100*swver) .eq. 600
     .       .or. nint(100*swver) .eq. 710
     .       .or. nint(100*swver) .eq. 720
     .       .or. nint(100*swver) .eq. 730
     .       .or. nint(100*swver) .eq. 780
c  Consolidate all Z-12 firmware and later versions YB 10/15/99)
     .       .or. (nint(100*swver).ge.790.and.nint(100*swver).le.974) )
     .      then
cc             nchan = 12
             slop = 0.300d0
         elseif(nint(100*swver).eq.980 .or. nint(100*swver).eq.981) then  
cc             nchan = 18
             slop = 0.300d0  
         elseif ( nint(100*swver) .eq.990 ) then  
cc             nchan = 12   
c**          it may be possible to tighten this since the clock is calibrated to UTC
             slop = 0.300d0     
         elseif ( nint(100*swver) .ge.996 .and.
     .            nint(100*swver) .le.1001 ) then  
cc             nchan = 24   
c**          it may be possible to tighten this since the clock is calibrated to UTC
             slop = 0.300d0  
         else
            write (uinfor,310) rcvrsw,swver
 310        format(1x,
     . 'SETTIM WARNING: unknown ASHTECH version ',
     . a3,f10.4,/,1x,'Current known versions:v.1.00,2.00,3.00,6.00'
     .    , '7.10,7.20,7.30,7.80,7.90,7.92'
     .    , '8.00,8.01-8.05,8.10,8.20-8.29,8.30-8.36,8.40,8.50,8.60'
     .    , '8.70,8.72-8.74,8.80,8.85,8.90-8.92,8.95'
     .    , '9.00,9.10,9.20,9.40,9.45,9.50,9.70-9.73,9.80,9.81,9.90'
     .    , '9.91,9.92,9.93,9.94,9.95,9.96,9.97,9.98,10.01'
     .    ,' Assuming that data are sampled like Ashtech 9.20')
            write (message,311) rcvrsw,swver
 311        format(1x,
     . 'SETTIM WARNING: unknown ASHTECH',
     . a3,f10.4,1x,'Assuming that data are sampled like Ashtech 9.20')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
            slop = 0.300d0
cc            nchan = 12
         endif 
    
                 

c     Set parameters for the Sercel
c     *******************************

      else if (lowerc(rcvrsw).eq.lowerc('SRT')) then
cc         nchan = 5
         slop = 0.3d0
         call secsum(iwkn0,sow0,-0.00d0,iwkn0,sow0)
      else if (lowerc(rcvrsw).eq.lowerc('SRN')) then
cc         nchan = 5
         slop = 0.3d0
         call secsum(iwkn0,sow0,-0.00d0,iwkn0,sow0)  

c     Set parameters for the Topcons
c     *******************************      

c     Most of these are relabeled versions of receivers from other manufacturers;
c     for all of the ones we've identified, there is now a translation in 
c     MAKEXP from 'TOP' to the original name.  The newest 'Hiper' geodetic receivers
c     are only issued by Topcon, however, so there is no translation to JPS.
c     rwk 030318; list updated rwk 070607 
c            
c TOPPR1        TOPCON GP-R1         C |  L1 receiver(Ashtech)      
c TOPRID      * TOPCON GP-R1D        C |  L1/L2 receiver(Ashtech L-XII)   
c TOPRIP      * TOPCON GP-R1DP       C |  L1/L2 receiver(Ashtech P-XII3)    
c TOPR1Y      * TOPCON GP-R1DY       N |  L1/L2 receiver(Ashtech Z-XII3)       
c TOPPS1        TOPCON GP-S1         C |  L1 receiver(Ashtech)                 
c TOPTRB      * TOPCON TURBO-SII     N |  L1/L2 receiver(AOA RASCAL-8))           
c TOPSX1        TOPCON GP-SX1        C |  L1 receiver(Trimble)         
c TOPSSI      * TOPCON TT4000SSI     P |  L1/L2 receiver(Trimble)    
c TOPDX1      * TOPCON GP-DX1        P |  L1/L2 receiver(Trimble 4700)         
c TOPRIS      * TOPCON GP-R1SD         |  Deprecated code (Apr 2001).  Do not use. 
c TOPLEG      * TOPCON LEGACY GD     N |  L1/L2 receiver (Javad Legacy, no IGS code)  
c TPSLEG      * TPS LEGACY           N |   same as TOPCON LEGACY GD  
c JPSLGE      * JPS E_GGD            N |  Legacy E" 160mm Eurocard-based GPS/GLONASS dual requency receiver
c TPSLGE      * TPS E_GGD            N |    same as JPS E_GGD 
c TPSREG        TPS REGENCY          N |  GPS/GLONASS dual- or single-frequency receiver w/ internal (Regant) chokering antenna, same as JPS REGENCY
c TOPODY        TOPCON ODYSSEY       N |  GPS/GLONASS dual- or single-frequency receiver, same as JPS ODYSSEY (no IGS code)   
c TPSODE      * TPS ODYSSEY_E        N |  GPS/GLONASS L1/L2 with integrated CE device
c JPSEUR      * JPS EUROCARD         N |  GPS and/or GLONASS dual- or single-freq receivers including Eurocard-based Legacy
c TOPHIP      * TOPCON HIPER         N |  L1 receiver (Javad, no IGS code)
c TOPHGD      * TOPCON HIPER-GD      N |  L1/L2 receiver (Javad, no IGS code) 
c TPSHGD      * TPS HIPER-GD         N |   same as TOPHGD / TPS HIPER-GD (Javad)
c TOPHGG      * TOPCON HIPER-GGD     N |  L1/L2 GPS/Glonass receiver (Javad)
c TPSHGG      * TPS HIPER-GGD        N |  same as TOPHGG / TOPCON HIPER-GGD 
c TOPGB1      * TPS GB-1000          N |  Dual-frequency GPS/GLONAS 
c TPSNG3      * TPS NETG3            N |  GPS/GLONASS/Galileo 72 channel receiver   
c TPSGR3        TPS GR3              N | 

      else if (lowerc(rcvrsw).eq.lowerc('TOP')) then   

c           The newest Topcon/Javad receivers use firmware 2.1-2.5, which is 
c           a near conflict with the arbitrary codes we've assigned to the
c           early Topcon/Ashtech L12 and P12 receivers.  These 'TOP' receiver types,
c           however, are translated to 'ASH' in MAKEXP, so there should be no conflict. 
c           The NETG3 handles GLONASS and Galileo and hence may have 72 channels.
           
        if (      nint(100*swver) .ge. 210 .and. 
     .            nint(100*swver) .le. 260 ) then
cc             nchan = 20
             slop = 0.300d0    
        elseif(   nint(100*swver) .eq. 310 ) then
cc             nchan = 72
             slop = 0.300d0
        else
            write (uinfor,350) rcvrsw,swver
 350        format(1x,'SETTIM WARNING: unknown TOPCON version ',a3,f10.4
     .        , /,1x,'Current known versions: 2.1-2.5 (Javad), 3.1'
     .        ,' Assuming that data are sampled like Topcon 2.1')
            write (message,355) rcvrsw,swver
 355        format(1x,'Unknown TOPCON version ',a3,f10.4,1x
     .           ,'Assuming that data are sampled like Topcon 2.1')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
            slop = 0.300d0
cc            nchan = 20
        endif

c     Set parameters for the WM 102
c     *****************************

      else if (lowerc(rcvrsw).eq.lowerc('WM2')) then
cc         nchan = 8

c     Set parameters for the UNAVCO L1 volcano-monitoring system
c     ***********************************************************

      else if (lowerc(rcvrsw).eq.lowerc('CMC')) then
c        we'll tentatively reserve v1.0 for the early firmware with bogus time tags
         if( swver.ge.2.0 .and. swver.lt.3.0 ) then
cc            nchan = 12   
            slop = 0.1d0 
         else
            write (uinfor,356) rcvrsw,swver
 356        format(1x,
     . 'SETTIM WARNING: unknown version for CMC receiver ',
     . a3,f10.4,/,1x,'Current known versions: 2.00 ',/
     . ,1x,'Assuming that data are sampled like CMCA12 v 2.0')
            write (message,357) rcvrsw,swver
 357        format('Unknown version for CMC receiver ',a3,f10.4,
     .      ' Current known versions: 2.00',
     .      ' Assuming that data are sampled like CMCA12 v 2.0')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
cc            nchan = 12
            slop = 0.1d0
         endif


c     Set parameters for the Geotracer 
c     ********************************
           
      else if (lowerc(rcvrsw).eq.lowerc('GEO')) then
          slop = 0.1                       
         if( swver.ge.8.0 .and. swver.lt.9.0 ) then
c          arbitrarily assign 8.0 to single-frequency Geotracer 2000, vers 806-0200
cc           nchan = 9                                                  
         else if ( swver.ge.4.0 .and. swver.lt.4.1 .or.  
c          the Geotracer 3220 seems to use firmware C 04.05
     .             swver.ge.9.0 .and.swver.lt.9.1) then    
c          arbitrarily assign 9.0 to dual-frequency Geodtracer 2200, vers CU-MB1 0
cc           nchan = 12 
         else
            write (uinfor,358) rcvrsw,swver
 358        format(1x,
     . 'SETTIM WARNING: unknown version for Geotracer receiver ',
     . a3,f10.4,/,1x,'Current known versions: 4.05, 8.0, 9.0 ',/
     . ,1x,'Assuming that data are sampled like C 04.05 ')
            write (message,359) rcvrsw,swver
 359        format('Unknown version for Geotracer receiver ',a3,f10.4,
     .      ' Current known versions: 4.05, 8.0, 9.0',
     .      ' Assuming that data are sampled like C 04.05')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
cc            nchan = 12
            slop = 0.1d0
         endif
                                            

c     Set Parameters for the Rogue (except for Turbo-Rogue)
c     ****************************

      else if (lowerc(rcvrsw).eq.lowerc('ROG')) then
cc         nchan = 8

C Information from Ulf Lindqwister (10/29/92)

C Rogue
C -----
C
C Receiver Type     Antenna Type             Description
C -------------     -----------------------  -----------------------------
C SNR-8             Dorne-Margolin C146-6-1  2 unit rack-mounted (big Rogue)
C SNR-800           Dorne-Margolin C146-6-1  1 unit rack-mounted (big Rogue)
C SNR-8A            Dorne-Margolin C146-6-1  MiniRogue - not CONAN compatible
C SNR-8C            Dorne-Margolin C146-6-1  MiniRogue - CONAN compatible
C SNR-8000          Dorne-Margolin C146-6-1  TurboRogue

c Rogue v1.51 indicates UTC time tags for GOTEX  
c Rogue v2.20 seems to be UTC time tags (GEOMEX 89)
c Rogue v2.3 collected UTC time tags for a long time and then
c       GPS time tags (the version number was not changed)    
C       We will arbitrarily assign v2.30 to UTC time tags and v2.31
C       for GPS time tags   
c  rwk 96/9/23 : Made < 2.30 UTC
c  rwk 99/8/30 : Changed back to check specific firmware versions and
c                also made 0.0 default post-2.31 (GPST) versions
         if  ( nint(100*swver) .eq. 151 
     .   .or.  nint(100*swver) .eq. 220 
     .   .or.  nint(100*swver) .eq. 230 ) then  
c        Here for UTC time tags   
         write (message,360) swver
 360     format('ROGUE firmware version ',f10.4,' assume UTC time tag')
         call report_stat('WARNING','MAKEX','settim',' ',message,0)
c        add the GPST-UTC difference
c        calculate GPST - UTC
         igpsdow = int(sow0/86400.d0)
         jd  = julday(1,5,1980) + 7*iwkn0 + igpsdow
c        get the UTC offset according to JD in GPS time
         utcoff = taiutc(jd) - 19.d0
         call secsum( iwkn0,sow0,utcoff, iwkn0,sow0 )

         elseif (nint(100*swver) .eq. 231
     .       .or. nint(100*swver).eq. 240
     .       .or. nint(100*swver).eq. 150
     .       .or. nint(100*swver).eq. 111
     .       .or. nint(100*swver).eq. 550
     .       .or. nint(100*swver).eq. 560
     .       .or. nint(100*swver).eq. 561
     .       .or. nint(100*swver).eq. 611
     .       .or. nint(100*swver).eq. 700
     .       .or. nint(100*swver).eq. 710
     .       .or. nint(100*swver).eq. 720
     .       .or. nint(100*swver).eq. 730
     .       .or. nint(100*swver).eq. 740
     .       .or. nint(100*swver).eq. 750
     .       .or. nint(100*swver).eq. 760
     .       .or. nint(100*swver).eq. 770
     .       .or. nint(100*swver).eq. 780) then
          continue
         else
            write (uinfor,361) rcvrsw,swver
 361        format(1x,'SETTIM WARNING: unknown ROGUE version ',a3,f10.4
     .,/,1x,
     .'Current known versions: v.1.11,1.5,2.3,2.4,5.5,5.6,6.1,7.0-7.5'
     .,/,1x,'Assuming that data are sampled like ROGUE 7.5')
            write (message,362) rcvrsw,swver
 362        format('Unknown ROGUE version ',a3,f10.4,
     .' Current known versions: v.1.11,1.5,2.3,2.4,5.5,5.6,6.1,7.0-7.5',
     .' Assuming that data are sampled like ROGUE 7.5')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
         endif

c     Set Parameters for the Turbo-Rogue
c     **********************************

      else if (lowerc(rcvrsw).eq.lowerc('TRB')) then
cc         nchan = 8

         if (nint(100*swver) .ge. 100 .and.
     .       nint(100*swver) .lt. 320 ) then
cc             nchan = 8
         else if (nint(100*swver).ge.320 ) then
cc             nchan = 12
c            this not strictly true; some v3.20 still have only 8 channels
         else
            write (uinfor,370) rcvrsw,swver
 370        format(1x,
     .      'SETTIM WARNING: unknown Turbo-ROGUE version ',a3,f10.4
     .,/,1x,'Current known versions: v.1.00 - 3.30'
     .,/,1x,'Assuming that data are sampled like Turbo-ROGUE 3.20')
            write (message,375) rcvrsw,swver
 375        format('Unknown Turbo-ROGUE version ',a3,f10.4,
     .      ' Assuming that data are sampled like Turbo-ROGUE 3.20')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
         endif
                             
c     Set parameters for the Javad
c     ****************************

      else if (lowerc(rcvrsw).eq.lowerc('JPS').or.
     .         lowerc(rcvrsw).eq.lowerc('JNS') ) then
c        firmware encountered so far is 2.0, Aug 31, 1999, and 2.6 Apr 18, 2007
c        but assume all firmware is 20 channels and even-minute sampling
cc         nchan = 20
         slop = 0.3d0
         call secsum(iwkn0,sow0,-0.00d0,iwkn0,sow0)   
 

c     Set parameters for Sokkia receivers
c     ***********************************

      else if (lowerc(rcvrsw).eq.lowerc('SOK')) then 
c        The Sokkia receivers are products of Novatel, and contain the OEM3 or OEM4
c        GPS Cards.  The firmware listed on the one example RINEX file we have,
c        from a 'Sokkia Radian   ', is 4.502/2.03.  This is not listed specifically
c        in the Novatel manual, which has for the OEM3 4.503, and for the OEM4 2.140.
c        The earlier OEM cards seem to track 12 satellites, but the NOV OEMV3 receiver 
c        with firmware 3.20 seems to have the capability of tracking 14.  
c        rwk 110930: Allowed up to 24 channels to accommodate the GSR2700 RSX with 
c        firmware 3.701. This is just a guess, though. 
c         
      
cc         nchan = 24
         slop = 0.3d0
         call secsum(iwkn0,sow0,-0.00d0,iwkn0,sow0)     
                              

c     Set parameters for the NavCom 
c     ******************************

      else if (lowerc(rcvrsw).eq.lowerc('NAV')) then
c        The Navcom translator currently has firmware hardwired as '1_80',
c        so we'll code just this example for now.  R. King from e-mail of Kevin Dixon
c        rcved 17 Dec 2003.  All receivers are dual-frequency and track 10 GPS satellites
cc         nchan = 10
         slop = 0.3d0
         call secsum(iwkn0,sow0,-0.00d0,iwkn0,sow0)       
                      
c     Set parameters for CHC  
c     **********************

      else if (lowerc(rcvrsw).eq.lowerc('CHC')) then
c       Email from Alexander Bragin, Russian sales rep for CHC indicates a potential
c       firmware of 4.93/1.0.6. Put this into guess_rcvant.dat and allow any firmware here.
c       rwk 150422
cc         nchan = 10
         slop = 0.3d0
         call secsum(iwkn0,sow0,-0.00d0,iwkn0,sow0)       


c     Set parameters for the Septentrio 
c     **********************************

      else if (lowerc(rcvrsw).eq.lowerc('SEP')) then
c        So far we have the SEPT POLARX2 with firmware 2.5.0 (48 chns), 3.2.0-patch1,
c        the SEPT POLARX4 with firmware 2.3.4, and     
c        and the AsteRx2 with firmware 0.0-tst080930r20788 
cc         nchan = maxchn
         slop = 0.3d0
         call secsum(iwkn0,sow0,-0.00d0,iwkn0,sow0)     

c     Set parameters for the Novatel cards
c     ************************************
                                                   
c        The Novatel receivers contain the Euro OEM3 or OEM4  GPS Cards.  The firmware listed on 
c        the one example RINEX file we have, from a 'Sokkia Radian   ', is 4.502/2.03.  This is 
c        not listed specifically in the Novatel manual, which has for the OEM3 4.503, and for the 
c        OEM4 2.140.  The earlier OEM cards seem to track 12 satellites, but the NOV OEMV3 receiver 
c        with firmware 3.20 seems to have the capability of tracking 14.  
c        rwk 090628: An email from Juan Espinoza (see letters/nocquet) indicates that
c        the Novatel GSV 4004B (EURO3M card) has firmware 9.140S43; # of channel unknown. 
c   
      else if (lowerc(rcvrsw).eq.lowerc('NOV')) then  
cc       nchan = 12
       slop = 0.3d0      
c      known firmware versions for the OEM-3 'Millenium cards are 4.45, 4.50, 4.501, 4.503, 4.52 
       if( nint(100*swver).eq.320 ) then
cc         nchan =14
         slop = 0.3d0 
       elseif (nint(100*swver).ge.445.and.nint(100*swver).le.452) then
cc         nchan =12
         slop = 0.3d0  
       elseif( nint(100*swver).eq.914 ) then
cc         nchan =12
         slop = 0.3d0  
       else
          write (uinfor,381) rcvrsw,swver
 381      format(1x,
     . 'SETTIM WARNING: unknown version for Novatel  receiver '
     . , a3,f10.4,/,1x,'Current known versions: 3.20,4.45, 4.50,'
     , ,' 4.501, 4.503, 4.52, 9.14',/,1x
     . ,'Assuming that data are sampled like 4.45 ')
            write (message,382) rcvrsw,swver
 382        format('Unknown version for Novatel receiver ',a3,f10.4,
     .      ' Current known versions: 4.45, 4.50, 4.501, 4.503, 4.52',
     .      ' Assuming that data are sampled like 4.45')
            call report_stat('WARNING','MAKEX','settim',' ',message,0)
cc            nchan = 12
            slop = 0.1d0
        endif
                                                         

c     Set parameters for the Hemisphere 
c     *********************************

      else if (lowerc(rcvrsw).eq.lowerc('HEM')) then
c        The only example I have is a P320 Eclipse II with firmware MFA_1.2Qe
c        which I'll call 1.20 for now. 
cc         nchan = 24
         slop = 0.3d0
         call secsum(iwkn0,sow0,-0.00d0,iwkn0,sow0)      

c     Receiver software not recognized, print a warning
c     **************************************************

      else
         write (uinfor,400) rcvrsw,swver
 400     format(1x,'SETTIM WARNING: unknown software version ',a3,f10.4
     .    ,/,
     .    4x,'Current known versions: ',/,
     .    6x,'GESAR     (GES) v. 25.2, 1.0-1.5' ,/,
     .    6x,'CORE      (COR) v. 4.1,4.11,4.12,4.8',/,
     .    6x,'ROM (a.k.a. TI4100 NAVIGATOR, NAVDAPT',/,
     .    6x,'MACROMETER II (MAC)        ',/,
     .    6x,'Minimac   (MIN) v.1.49,1.50,1.61 - 1.64',/, 
     .    6x,'LEICA     (LEI)  v.2.0,2.5,3.5,4.0,4.5,5.0,5.5,6.0,6.5',/,
     .    6x,'Trimble   (TRM) v.1.xx,2.xx,3.12,3.22,3.25,4.10 - 7.29',/,   
     .    6x,'Ashtech   (ASH) v.1.00,2.00,3.00,6.00,
     .    7.10,7.20,7.30,7.80,7.90,7.92,8.00,8.01-8.05,8.10,8.20-8.29,
     .    8.30-8.36,8.40,8.50,8.70,8.72-8.74,8.80,8.85,
     .    8.90-8.92,8.95,9.00,9.10,9.20,9.40,9.45,9.50,9.70-9.73',/,
     .    6x,'ROGUE     (ROG) v.1.11,1.5,2.3,2.4,5.5,5.6,6.1,7.0-7.8',/,
     .    6x,'TurboRogue(TRB) v.1.00,2.00-3.30 ',/,
     .    6x,'Geotracer (GEO) v.8.0, 9.0 ',/,
     .    6x,'CMC Allstar 12 (CMC) v.1.0 ',/, 
     .    6x,'Javad (JPS or JNS ) all versions ',/,
     .    4x,'Assuming that data are sampled at even second')
         write (message,405) rcvrsw,swver
 405     format('Unknown firmware version in makex.batch ',a3,f10.4,
     .   ' Current known versions: GES COR ROM MAC MIN TRM ASH ROG LEI')
         call report_stat('WARNING','MAKEX','settim',' ',message,0)
         call secsum(iwkn0,sow0,0.0d0,iwkn0,sow0)
         slop = 0.25d0
cc         nchan = maxchn
      endif

cc      if( nchan.gt.maxchn ) then
cc         write(message,'(a,i3,a)') 
cc     .      'Number of channels requested for rcvr > maxchn (',maxchn
cc     .     ,') Change maxchn in makex.h and remake GAMIT'
cc         call report_stat('FATAL','MAKEX','settim',' ',message,0)
cc    endif
       
      if( debug ) then
cc        print *,'SETTIM: (2) nchan = ',nchan
        print *,'SETTIM: (2) iwkn0 = ', iwkn0
        print *,'SETTIM: (2) sow0  = ', sow0
        print *,'SETTIM: (2) slop   = ', slop
      endif

      return
      end
