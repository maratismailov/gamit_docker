      subroutine mapamp ( iaflag,rcvrsw,swver,l1l2, disnr,isnr,issi )

c     Map receiver signal amplitudes into RINEX signal strength indicators and/or
c     GAMIT X-file error flags


c     R. King -  16 January 1990

c        IAFLAG = 1  Input is receiver amplitude (FICA or raw data),
c                    stored in DISNR; output is GAMIT signal amplitude (ISNR)
c                    (temporary) and RINEX signal strength (ISSI)
c
c        IAFLAG = 2  Input is GAMIT amplitude (ISNR), output is RINEX
c                    signal strength indicator (ISSI)
c
c        IAFLAG = 3  Input is RINEX signal strength indicator (ISSI),
c                    output is GAMIT signal amplitude (ISNR)

c        L1L2  = 1   L1
c              = 2   L2
c
      implicit none

      include '../includes/makex.h'

      character*3 rcvrsw
      character*80 prog_name
      character*256 message

      integer iaflag,l1l2,isnr,issi,len,rcpar

      real*8 disnr

      real*4 swver

c    RCVRSW and SWVER indicate the receiver hardware and software
c      Currently supported:

c      GES x.x       GESAR versions 25.2, 1.0-1.5, 2.0 for the TI 4100
c                    or CORE versions 4.1,4.11,4.8
c
c      MAC x.x       MACROMETER II (MAC)
c
c      MIN x.x       Minimac version 1.49
c
c      TRM x.x       Trimble version 3.?

      logical known
      save known
      data known /.false./

      call uppers(rcvrsw)

c     Set parameters for the TI 4100
      if (rcvrsw.eq.'COR'. or.
     .    rcvrsw.eq.'GES'. or.
     .    rcvrsw.eq.'ROM' ) then

c       ** Provisional values from Werner Gurtner
        if( iaflag.eq.1 .or. iaflag.eq.2 ) then
           if( iaflag.eq.1 ) isnr = int(disnr)
c          if( isnr.eq.0 ) issi = 0
c          No SNR should be worst possible data rather than unknown
           if( isnr.eq.0 ) issi = 1
           if( isnr.gt.0  .and. isnr.lt.25 )  issi = 1
           if( isnr.ge.26 .and. isnr.le.28 )  issi = 2
           if( isnr.ge.29 .and. isnr.le.31 )  issi = 3
           if( isnr.ge.32 .and. isnr.le.34 )  issi = 4
           if( isnr.ge.35 .and. isnr.le.37 )  issi = 5
           if( isnr.ge.38 .and. isnr.le.40 )  issi = 6
           if( isnr.ge.41 .and. isnr.le.43 )  issi = 7
           if( isnr.ge.44 .and. isnr.le.46 )  issi = 8
           if( isnr.gt.46 )                   issi = 9
        endif

        if( iaflag.eq.3 ) then
           if( issi.eq.0 ) isnr =  0
           if( issi.eq.1 ) isnr = 20
           if( issi.eq.2 ) isnr = 26
           if( issi.eq.3 ) isnr = 30
           if( issi.eq.4 ) isnr = 33
           if( issi.eq.5 ) isnr = 36
           if( issi.eq.6 ) isnr = 39
           if( issi.eq.7 ) isnr = 42
           if( issi.eq.8 ) isnr = 45
           if( issi.eq.9 ) isnr = 49
        endif

      else if (rcvrsw.eq.'MAC') then
        if( iaflag.eq.1 .or. iaflag.eq.2 ) then
           if( iaflag.eq.1 ) isnr = int(disnr)
           if( l1l2.eq.1) then
              if( isnr.eq.0 ) issi = 0
              if( isnr.gt.0  .and. isnr.lt.7  )  issi = 1
              if( isnr.ge.7  .and. isnr.le.12 )  issi = 2
              if( isnr.ge.13 .and. isnr.le.18 )  issi = 3
              if( isnr.ge.19 .and. isnr.le.24 )  issi = 4
              if( isnr.ge.25 .and. isnr.le.30 )  issi = 5
              if( isnr.ge.31 .and. isnr.le.36 )  issi = 6
              if( isnr.ge.37 .and. isnr.le.42 )  issi = 7
              if( isnr.ge.43 .and. isnr.le.49 )  issi = 8
              if( isnr.gt.50 )                   issi = 9
           endif

           if( l1l2.eq.2) then
              if( isnr.eq.0 ) issi = 0
              if( isnr.gt.0  .and. isnr.lt.8  )  issi = 1
              if( isnr.ge.8  .and. isnr.le.21 )  issi = 2
              if( isnr.ge.22 .and. isnr.le.35 )  issi = 3
              if( isnr.ge.36 .and. isnr.le.49 )  issi = 4
              if( isnr.ge.50 .and. isnr.le.63 )  issi = 5
              if( isnr.ge.64 .and. isnr.le.77 )  issi = 6
              if( isnr.ge.78 .and. isnr.le.91 )  issi = 7
              if( isnr.ge.92 .and. isnr.le.105)  issi = 8
              if( isnr.gt.105)                   issi = 9
           endif
        endif

        if( iaflag.eq.3 ) then
           if( l1l2.eq.1 ) then
              if( issi.eq.0 ) isnr =  0
              if( issi.eq.1 ) isnr =  5
              if( issi.eq.2 ) isnr = 10
              if( issi.eq.3 ) isnr = 16
              if( issi.eq.4 ) isnr = 22
              if( issi.eq.5 ) isnr = 28
              if( issi.eq.6 ) isnr = 34
              if( issi.eq.7 ) isnr = 40
              if( issi.eq.8 ) isnr = 46
              if( issi.eq.9 ) isnr = 50
           endif
           if( l1l2.eq.2) then
              if( issi.eq.0 ) isnr =  0
              if( issi.eq.1 ) isnr =  5
              if( issi.eq.2 ) isnr = 15
              if( issi.eq.3 ) isnr = 30
              if( issi.eq.4 ) isnr = 42
              if( issi.eq.5 ) isnr = 58
              if( issi.eq.6 ) isnr = 70
              if( issi.eq.7 ) isnr = 85
              if( issi.eq.8 ) isnr = 98
              if( issi.eq.9 ) isnr =110
           endif
        endif

c     Minimac
      else if (rcvrsw.eq.'MIN') then
        if( iaflag.eq.1 .or. iaflag.eq.2 ) then
           if( iaflag.eq.1 ) isnr = int(disnr)
           if( l1l2.eq.1) then
              if( isnr.le.0 ) issi = 0
              if( isnr.gt.0  .and. isnr.lt.2  )  issi = 1
              if( isnr.ge.2  .and. isnr.le.3  )  issi = 2
              if( isnr.ge.4  .and. isnr.le.5  )  issi = 3
              if( isnr.ge.6  .and. isnr.le.8  )  issi = 4
              if( isnr.ge.9  .and. isnr.le.12 )  issi = 5
              if( isnr.ge.13 .and. isnr.le.16 )  issi = 6
              if( isnr.ge.17 .and. isnr.le.21 )  issi = 7
              if( isnr.ge.22 .and. isnr.le.30 )  issi = 8
              if( isnr.gt.30 )                   issi = 9
           endif

           if( l1l2.le.2) then
              if( isnr.eq.0 ) issi = 0
              if( isnr.gt.0  .and. isnr.lt.2  )  issi = 1
              if( isnr.ge.2  .and. isnr.le.3  )  issi = 2
              if( isnr.ge.4  .and. isnr.le.5  )  issi = 3
              if( isnr.ge.6  .and. isnr.le.8  )  issi = 4
              if( isnr.ge.9  .and. isnr.le.15 )  issi = 5
              if( isnr.ge.16 .and. isnr.le.26 )  issi = 6
              if( isnr.ge.27 .and. isnr.le.39 )  issi = 7
              if( isnr.ge.40 .and. isnr.le.55 )  issi = 8
              if( isnr.gt.55 )                   issi = 9
           endif
        endif

        if( iaflag.eq.3 ) then
           if( l1l2.eq.1 ) then
              if( issi.eq.0 ) isnr =  0
              if( issi.eq.1 ) isnr =  1
              if( issi.eq.2 ) isnr =  3
              if( issi.eq.3 ) isnr =  5
              if( issi.eq.4 ) isnr =  7
              if( issi.eq.5 ) isnr = 11
              if( issi.eq.6 ) isnr = 15
              if( issi.eq.7 ) isnr = 19
              if( issi.eq.8 ) isnr = 26
              if( issi.eq.9 ) isnr = 31
           endif
           if( l1l2.eq.2) then
              if( issi.eq.0 ) isnr =  0
              if( issi.eq.1 ) isnr =  1
              if( issi.eq.2 ) isnr =  3
              if( issi.eq.3 ) isnr =  5
              if( issi.eq.4 ) isnr =  7
              if( issi.eq.5 ) isnr = 11
              if( issi.eq.6 ) isnr = 20
              if( issi.eq.7 ) isnr = 32
              if( issi.eq.8 ) isnr = 48
              if( issi.eq.9 ) isnr = 56
           endif
        endif


c     Set parameters for the Trimble
      else if (rcvrsw.eq.'TRM') then
        if( iaflag.eq.1 .or. iaflag.eq.2 ) then
           if( iaflag.eq.1 ) isnr = int(disnr)

***********
c           if( l1l2.eq.1) then
           if( l1l2.eq.1 .or. l1l2.eq.2) then
              if( isnr.eq.0 ) issi = 0
              if( isnr.gt.0  .and. isnr.lt.3  )  issi = 1
              if( isnr.ge.3  .and. isnr.le.4  )  issi = 2
              if( isnr.ge.5  .and. isnr.le.6  )  issi = 3
              if( isnr.ge.7  .and. isnr.le.8  )  issi = 4
              if( isnr.ge.9  .and. isnr.le.11 )  issi = 5
              if( isnr.ge.12 .and. isnr.le.17 )  issi = 6
              if( isnr.ge.18 .and. isnr.le.27 )  issi = 7
              if( isnr.ge.28 .and. isnr.le.40 )  issi = 8
              if( isnr.gt.40 )                   issi = 9
           endif
*               ***************
cancel this for the time being
           if( l1l2.eq.2) then
              if( isnr.eq.0 ) issi = 0
              if( isnr.gt.0  .and. isnr.lt.10 )  issi = 1
              if( isnr.ge.10 .and. isnr.le.17 )  issi = 2
              if( isnr.ge.18 .and. isnr.le.24 )  issi = 3
              if( isnr.ge.25 .and. isnr.le.31 )  issi = 4
              if( isnr.ge.32 .and. isnr.le.42 )  issi = 5
              if( isnr.ge.43 .and. isnr.le.89 )  issi = 6
              if( isnr.ge.90 .and. isnr.le.149)  issi = 7
              if( isnr.ge.150.and. isnr.le.220)  issi = 8
              if( isnr.gt.220)                   issi = 9
           endif
        endif

        if( iaflag.eq.3 ) then
           if( l1l2.eq.1 ) then
              if( issi.eq.0 ) isnr =  0
              if( issi.eq.1 ) isnr =  1
              if( issi.eq.2 ) isnr =  4
              if( issi.eq.3 ) isnr =  6
              if( issi.eq.4 ) isnr =  8
              if( issi.eq.5 ) isnr = 10
              if( issi.eq.6 ) isnr = 15
              if( issi.eq.7 ) isnr = 23
              if( issi.eq.8 ) isnr = 34
              if( issi.eq.9 ) isnr = 41
           endif
           if( l1l2.eq.2) then
              if( issi.eq.0 ) isnr =  0
              if( issi.eq.1 ) isnr =  1
              if( issi.eq.2 ) isnr = 13
              if( issi.eq.3 ) isnr = 21
              if( issi.eq.4 ) isnr = 28
              if( issi.eq.5 ) isnr = 37
              if( issi.eq.6 ) isnr = 66
              if( issi.eq.7 ) isnr =120
              if( issi.eq.8 ) isnr =185
              if( issi.eq.9 ) isnr =225
           endif
        endif

******** modified mbo 5/13/1992
*** added rogue mapping : PROVISIONAL

c     Set parameters for the ROGUE
      else if (rcvrsw.eq.'ROG') then
        if( iaflag.eq.1 .or. iaflag.eq.2 ) then
           if( iaflag.eq.1 ) isnr = int(disnr)

***   The following values are obtained from july 1991 unavco RINEX translators
***  option files. I followed the logic as implemented for others and trimbles
***  since TRM numbers are very close to ROG's
***  ALso note that, issi = 5 is set to the threshold value

***MAXIMUM S/N FOR L1            --> :   255
***MAXIMUM S/N FOR L2            --> :   255
***THRESHOLD S/N VALUE FOR L1    --> :    40
***THRESHOLD S/N VALUE FOR L2    --> :    40
***MINIMUM S/N FOR L1            --> :     0
***MINIMUM S/N FOR L2            --> :     0

           if( l1l2.eq.1) then
              if( isnr.eq.0 ) issi = 0
              if( isnr.gt.0  .and. isnr.lt.10 )  issi = 1
              if( isnr.ge.10 .and. isnr.le.17 )  issi = 2
              if( isnr.ge.18 .and. isnr.le.24 )  issi = 3
              if( isnr.ge.25 .and. isnr.le.31 )  issi = 4
              if( isnr.ge.32 .and. isnr.le.40 )  issi = 5
              if( isnr.ge.41 .and. isnr.le.90 )  issi = 6
              if( isnr.ge.91 .and. isnr.le.169)  issi = 7
              if( isnr.ge.170.and. isnr.le.249)  issi = 8
              if( isnr.gt.250)                   issi = 9
           endif

           if( l1l2.eq.2) then
              if( isnr.eq.0 ) issi = 0
              if( isnr.gt.0  .and. isnr.lt.10 )  issi = 1
              if( isnr.ge.10 .and. isnr.le.17 )  issi = 2
              if( isnr.ge.18 .and. isnr.le.24 )  issi = 3
              if( isnr.ge.25 .and. isnr.le.31 )  issi = 4
              if( isnr.ge.32 .and. isnr.le.40 )  issi = 5
              if( isnr.ge.41 .and. isnr.le.90 )  issi = 6
              if( isnr.ge.91 .and. isnr.le.169)  issi = 7
              if( isnr.ge.170.and. isnr.le.249)  issi = 8
              if( isnr.gt.250)                   issi = 9
           endif
        endif

        if( iaflag.eq.3 ) then
           if( l1l2.eq.1 ) then
              if( issi.eq.0 ) isnr =  0
              if( issi.eq.1 ) isnr =  1
              if( issi.eq.2 ) isnr = 13
              if( issi.eq.3 ) isnr = 21
              if( issi.eq.4 ) isnr = 28
              if( issi.eq.5 ) isnr = 36
              if( issi.eq.6 ) isnr = 66
              if( issi.eq.7 ) isnr =131
              if( issi.eq.8 ) isnr =210
              if( issi.eq.9 ) isnr =250

           endif
           if( l1l2.eq.2) then
              if( issi.eq.0 ) isnr =  0
              if( issi.eq.1 ) isnr =  1
              if( issi.eq.2 ) isnr = 13
              if( issi.eq.3 ) isnr = 21
              if( issi.eq.4 ) isnr = 28
              if( issi.eq.5 ) isnr = 36
              if( issi.eq.6 ) isnr = 66
              if( issi.eq.7 ) isnr =131
              if( issi.eq.8 ) isnr =210
              if( issi.eq.9 ) isnr =250
           endif
        endif


c     Set parameters for the ASHTECH : PROVISIONAL  > MBO
      else if (rcvrsw.eq.'ASH' .or. rcvrsw.eq.'TOP') then
        if( iaflag.eq.1 .or. iaflag.eq.2 ) then
           if( iaflag.eq.1 ) isnr = int(disnr)

***   The following values are obtained from july 1991 unavco RINEX translators
***  option files. I followed the logic as implemented for others
***  ALso note that, issi = 5 is set to the threshold value

***MAXIMUM S/N FOR L1            --> :   130
***MAXIMUM S/N FOR L2            --> :    40
***THRESHOLD S/N VALUE FOR L1    --> :    40
***THRESHOLD S/N VALUE FOR L2    --> :    10
***MINIMUM S/N FOR L1            --> :     1
***MINIMUM S/N FOR L2            --> :     1


*** What do we do with MIN = 1
           if( l1l2.eq.1) then
              if( isnr.eq.0 ) issi = 0
              if( isnr.gt.0  .and. isnr.lt.10 )  issi = 1
              if( isnr.ge.10 .and. isnr.le.17 )  issi = 2
              if( isnr.ge.18 .and. isnr.le.24 )  issi = 3
              if( isnr.ge.25 .and. isnr.le.31 )  issi = 4
              if( isnr.ge.32 .and. isnr.le.40 )  issi = 5
              if( isnr.ge.41 .and. isnr.le.60 )  issi = 6
              if( isnr.ge.61 .and. isnr.le.109)  issi = 7
              if( isnr.ge.110.and. isnr.le.124)  issi = 8
              if( isnr.gt.125)                   issi = 9
           endif

           if( l1l2.eq.2) then
              if( isnr.eq.0 ) issi = 0
              if( isnr.gt.0  .and. isnr.lt.2  )  issi = 1
              if( isnr.ge.2  .and. isnr.le.3  )  issi = 2
              if( isnr.ge.4  .and. isnr.le.5  )  issi = 3
              if( isnr.ge.6  .and. isnr.le.7  )  issi = 4
              if( isnr.ge.8  .and. isnr.le.10 )  issi = 5
              if( isnr.ge.11 .and. isnr.le.16 )  issi = 6
              if( isnr.ge.17 .and. isnr.le.29 )  issi = 7
              if( isnr.ge.30 .and. isnr.le.37 )  issi = 8
              if( isnr.gt.38 )                   issi = 9
           endif
        endif

        if( iaflag.eq.3 ) then
           if( l1l2.eq.1 ) then
              if( issi.eq.0 ) isnr =  0
              if( issi.eq.1 ) isnr =  1
              if( issi.eq.2 ) isnr = 13
              if( issi.eq.3 ) isnr = 20
              if( issi.eq.4 ) isnr = 28
              if( issi.eq.5 ) isnr = 40
              if( issi.eq.6 ) isnr = 50
              if( issi.eq.7 ) isnr = 86
              if( issi.eq.8 ) isnr =117
              if( issi.eq.9 ) isnr =130
           endif
           if( l1l2.eq.2) then
              if( issi.eq.0 ) isnr =  0
              if( issi.eq.1 ) isnr =  1
              if( issi.eq.2 ) isnr =  2
              if( issi.eq.3 ) isnr =  4
              if( issi.eq.4 ) isnr =  6
              if( issi.eq.5 ) isnr =  9
              if( issi.eq.6 ) isnr = 14
              if( issi.eq.7 ) isnr = 23
              if( issi.eq.8 ) isnr = 34
              if( issi.eq.9 ) isnr = 40
           endif
        endif
************* end of modification


c       Receiver software not recognized, print a warning
c       only print once, though
      else
         if (.not. known) then
            write(message,'(a,a3,f10.4)') 'Rcvr/software not coded: '
     .                                  ,  rcvrsw,swver
c       Get calling program name for report_stat
            len = rcpar(0,prog_name)
            call report_stat('WARNING',prog_name,'lib/mapamp',' '
     .      ,message,0)
c            write (uscren,400) rcvrsw,swver
             write (uinfor,400) rcvrsw,swver
 400        format(1x,'WARNING: unknown software version',a3,f10.4,/,
     .            4x,'Current known formats: ',/,
     .            6x,'GESAR         (GES) v. 25.2, 1.0-1.5' ,/,
     .            6x,'CORE          (COR) v. 4.1,4.11,4.8',/,
     .            6x,'ROM (a.k.a. TI4100 NAVIGATOR, NAVDAPT',/,
     .            6x,'MACROMETER II (MAC)        ',/,
     .            6x,'Minimac       (MIN) v. 1.49, 1.61',/,
     .            6x,'Trimble       (TRM) v. 3.12' ,/,
**** added the next 2 lines
     .            6x,'ROGUE         (ROG) v. ?.??' ,/,
     .            6x,'TOPCON        (TOP) v. ?.??' ,/,
     .            6x,'Ashtech       (ASH) v. ?.??' )
            known = .true.
         endif
         isnr = issi
       endif

      return
      end
