c Copyright (c) Massachusetts Institute of Technology, 1988. All rights reserved.
      program ARGO2FIC
c*    formerly program ngs2fc
c
c
c     converts NGS ARGO format data to FICA format data.
c
c     Delaine Thompson, January 1988
c
c     roots:
c        makex                 Kurt Feigl et al.
c        sledge_hammer         peter morgan january 1987
c        sledge_hammer.test    Jerry Svarc on 4-30-87
c
c     Modifications:
c        version 1.5  New CIGNET ephemeris format: add 8 parameters,
c                     test date of data for output--mhm 880728
c        version 1.6  Correct v1.5 bugs--mhm 880729
c        version 1.7  Restructure the main program logic a bit to overcome the
c                      problem of misnamed files when days are missing.
c                     Add operator input of CORE version no.  and write BLK 101
c                         --rwk 880809
c        version 1.8  88/09/05 king
c                     Change logic to avoid losing phase recs
c                     after last BLK 9
c        version 1.9  89/05/02 kurt
c                     Adapt to change CIGNET: removal of frequency
c                     plan bias starting week 455
c                     Also begin testing Minimac data and
c                     general code-cleanup.
c                     Use standard gps library routines.
c        version 1.10 89/05/02 kurt
c                     Minimac FICA pseudoranges in kilometers.
c        version 1.11 89/07/01 kurt
c                     Incorporate Minimac data.
c        version 1.12 89/07/01 kurt
c                     Handle up to 100 read errors in .ORB file.
c        version 1.13 89/07/01 kurt
c                     Handle up to 100 read errors in .DAT file.
c        version 1.14 89/07/01 kurt
c                     Ooops. Last change didn't compile.
c        version 1.15 89/07/01 kurt
c                     This comes pretty close to working.
c        version 1.16 89/08/07 kurt
c                     BLOCK 1180s follow DOPPLER (not PSEUDORANGE)
c                     CONVENTION.
c                     Positive sign for increasing Doppler shift.
c                     Positive sign for decreasing Pseudorange
c        version 1.17 89/08/17 mhm
c                     Add 1 msec to Minimac time-tags in wb1180
c        version 1.18 90/07/05 mhm
c                     Adopt special CORE 4.12 for improper frequency
c                     plan implementation at Yellowknife in March 1990
c                     Pass variable swvrsn to WB55, which corrects
c                     the phases.
c        version 1.19 90/07/11 kurt
c                     Remove inline comments for Sun compatability.
c                     Send to Andrea for sun debugging.
c        version 1.20 90/11/01 kurt
c                     Repair bug found by Andrea on sun:
c                     set luday = 21
c        version 1.21 90/12/19 king
c                     Add Block 1201 and 1280 for CIGNET Trimble
c                     Add Block 1301 and 1380 for CIGNET Rogue
c                     91/06/18 king
c                     Correct changeover date for COR freq. plan in
c                        NGS data--from 455 to 456.
c        version 1.3  Change name to ARGO2FIC and put in /util directory
c                        King 92/06/22
c
c
c     This program takes ephemeris and phase data written in the NGS
c     format and rewrites it in the FICA format.  Since the NGS format
c     files contain a week's worth of data, the output from the program
c     is up to 7 FICA format day files.  The conversion requires some dummying
c     of variables which NGS does not include.  These are listed in
c     detail in the documentation.  Execution of this program is slow
c     because of the cumbersome rearrangement of data from two files
c     into one.
c
c     FILES:
c       input files include --
c           XXXX.orb     NGS file containing ephemeris and clock information
c           XXXX.dat     NGS file containing phase information
c         the user is asked for the names of both of these files
c
c       output files include --
c           [id]yr.day.fic  up to 7 FICA files, depending on the number
c                           of days of data
c
c  VARIABLE EXPLANATION
c
c     forb, fdat       character variables for NGS file names
c     lueph, luph       logical units associated with NGS files
c     luday             logical unit attached to a FICA day file
c     files             character array of FICA day file names
c     flag              logical variable; false after first record of
c                       phase has been read. (No longer used?)
c     ittag             ephemeris time tag info (yr,mo,dy,hr,min)
c     sec               ephemeris time tag -- seconds
c     ephr              ephemeris information (see subroutine)
c     gpswk             GPS week; returned by tconv
c     time              important sorting time tag in GPS sec of wk
c     soweph            same as above
c     sowold            second of week, buffers sowclk
c     locate            location array; holds positions of sorted time tags
c                       in time array
c     idoyp,idoyf       day of year; compared to open new day files
c     tphasr            phase time tag data
c     sowclk            GPS sec of week from phase data file
c     nsvtot            integer total of satellites
c     svid              integer array of the satellite ids
c     phasr            real array containing phase information
c  miscellaneous :
c     icnt,ncnt,ihold,ibot   integer counters
c     iplace                 location of correct ephemeris information
c
c  VARIABLE DECLARATION

      include '../includes/argo2fic.h'

c
c     L1,L2,P1,P2,SNR1,SNR2 for each SV
      real*8 phasr(maxsvs,6)
c     PRN numbers
      integer*4 svid(maxsvs)

      integer*4 icnt,ncnt,ihold,iplace,iflag
     .,         ittag(maxeph,7),idoyp,locate(maxeph)
     .,         nsvtot,idoyf,i,gpswk,ihnsec,ibot,j
     .,         iwknf,idoyd,idoyh
     .,         iyr,idoy,id,im,nerr
     .,         rt(6)
      integer*4 leng
      integer lueph,luph,luday,ios
      real*8 ephr(maxeph,27),soweph,sowclk,sowold,tphasr,time(maxeph)
     .,      sec,swvrsn,utcoff,argosv
      character   fyear*1
      character*3 fhold
      character*8 rcvrhw,vernum,filnam
      character*14 forb
      character*15 fdat
      character*16 dfname
      character*20 sitnam
      character*36 version
      character*80 buff80
      logical ludopn
c*      logical flag

      integer*2 tallyb(4,60)
      common /tally/ tallyb
c
c   DATA DECLARATIONS
c
      data version/'ARGO2FIC v. 1.3 of 92/06/22 16:00:00'/
c*      data flag /.true./
      data nerr /0/
c
c   RUN INITIALIZATION
c
c  get the run time and version number
      vernum = '1.17'
      call getdat(rt(1),rt(2),rt(3))
      call gettim(rt(4),rt(5),rt(6),ihnsec)
      write(6,20) version,(rt(i),i=1,6)
20    format (/,1x,a36,1x,' ON ',i4,'-',i2.2,'-',i2.2,
     $  1x,i2.2,':',i2.2,':',i2.2,//)
c
c  read the name of the NGS .orb and .dat files
      print 30
30    format(1x,'Enter the name of the CIGNET ephemeris file to be',
     $' translated.',/,1x,'nnnWWW.ORB e.g. MOJ477.ORB')
      read(5,'(a14)') forb
      print 35
35    format(/,1x,'Enter the name of the CIGNET phase data file to',
     $' translated.',/,1x,'nnnWWW.DAT e.g. MOJ477.DAT',/)
      read(5,'(a15)')fdat
c
c  open the NGS files
      lueph = 20
      luday = 31
      call openf(lueph,forb,'old','formatted','sequential')
c
      luph = 21
      call openf(luph,fdat,'old','formatted','sequential')
c
c  open an information file to count blocks
      i = index(forb,'.')
      filnam = forb(1:i-1)
      call openf(44,forb(1:i)//'inf','new',
     $          'formatted', 'sequential')
      call wb0(44,vernum,rt,filnam)
      write(44,*)'Summary of Blocks Written from NGS to FICA Format'
      call uptaly(0,0)
c
c   no FICA files are open at the beginning of processing
      ludopn = .false.
c
c  READ INITIAL EPHEMERIS FILE AND TITLE RECORD OF PHASE FILE
c  TO BEGIN PROCESSING
c
c      Pointers for data in ephemeris array:
c          ibot :  position in array of last record read
c          ncnt :  number of records in array
c      Day numbers used to match data with output files:
c          idoyp :  day-of-year of current phase record
c          idoyf :  day-of-year of open output FICA file
c          idoyh :  day-of-year from CIGNET header
c
c
c     Enter the CORE version number
      print 40
 40   format(1x,'Version number of reciever software?',/
     .      ,4x,'Reasonable choices include:',/
     .      ,4x,'CORE v. 4.1, 4.11, 4.12, 4.7, 4.8',/
     .      ,4x,'MINIMAC v. 1.49',/
     .      ,4x,'TRIMBLE v. 3.25, 4.1'
     .      )
      read(5,*) swvrsn
c
c     Read the header from the CIGNET phase file
      read(luph,'(a80)',end=56,err=56,iostat=ios) buff80
      read(buff80,50,end=56,err=56,iostat=ios)
     .        sitnam,iwknf,idoyh,iyr,im,id
 50   format(a20,5i5)
c     Receiver hardware, assumed TI 4100 unless otherwise stated.
      rcvrhw = buff80(51:58)
      if (rcvrhw(1:1) .eq. ' ') rcvrhw = 'TI4100  '
      read (buff80(71:76),'(f5.2)') argosv

 56   if(ios.eq.0) then
         write(44,55) sitnam,iwknf,idoyh,iyr,im,id,rcvrhw,swvrsn,argosv
         write(6,55)  sitnam,iwknf,idoyh,iyr,im,id,rcvrhw,swvrsn,argosv
 55      format(/,1x,'Title Record of Phase File - Station:  ',A20,/
     .         ,' Week:',I5,'  Day-of-yr:',I4,'   Date:',3I5,/
     .         ,' Hardware: ',a8,' Software version: ',f5.2,/
     .         ,' ARGOS version:  ',f5.2,//)
      else
         write(6,57)
 57      format(1x,'End or error reading title record of phase file')
         call ferror (ios,6)
         stop
      endif
c
c  Construct the root of the output FICA file names
      leng = index(forb,'.')
      write(fyear,'(i1)') mod(iyr,10)
c
c   Read in all the ephemeris and clock data and store it for later
c     interleaving with phase data
c
      icnt = 1
 100  continue
      call rdeph(lueph,ittag,sec,ephr,icnt,ios)
      if (ios .eq. -1) then
         ncnt = icnt - 1
         close(lueph)
         go to 200
      else
         icnt = icnt + 1
         go to 100
      end if

c     Convert CIGNET time tag for ephemeris (GPST yy mm dd hh mm ss.ssss) to
c     FICA time tag (GPST week number + second of week)
c     This time tag will be used for sorting
 200  do 210 i = 1,ncnt
c        get day of year
         idoyd = idoy(ittag(i,3),ittag(i,4),ittag(i,5))
c        convert to GPST week number and second of week
         iflag = -4
         call timcon(iflag,
     .              gpswk,soweph,
     .              ittag(i,3),idoyd,ittag(i,6),ittag(i,7),sec,
     .              utcoff)
c        if the data don't belong to this week flag time as negative
         if (gpswk .eq. iwknf) then
            time(i) = soweph
         else
            write (6,*) 'Found ephemeris block from week',gpswk
            time(i) = -soweph
         endif
 210  continue

c     build an array containing locations of sorted soweph time tags
c     for ephemeris blocks
      call tsort(time,ncnt,locate)
c
c  Initialize file day and ephemeris array pointers
c
      idoyf = idoyh - 1
      ibot = 1
c
c
c  LOOP OVER PHASE DATA, OPENING NEW FICA DAY FILES AS NECESSARY
c
c
 300  ihold = ibot
      call rdphas(luph,idoyp,tphasr,sowclk,nsvtot,svid,phasr,ios,nerr)
cd    write(6,8000) luph,idoyp,svid,ios
cd 8000  format(' at 300, luph,idoyp,svid,ios=',4i5)
      if (ios .eq. -1) goto 550
c
c  If the data day matches the currently open file day,
c  write FICA phase block (6 or 1080)
c  Then check for a (valid) current block 9 and write it.
c  Finally read another record
c
 400  if (idoyp.eq.idoyf) then
         if (rcvrhw .eq. 'TI4100  ') then
            call wb55(luday,iwknf,sowclk,nsvtot,svid,phasr,swvrsn)
         else if (rcvrhw .eq. 'MINIMAC ') then
            call wb1180(luday,iwknf,sowclk,nsvtot,svid,phasr)
         else if (rcvrhw .eq. 'TRIMBLE ') then
            call wb1280(luday,iwknf,sowclk,nsvtot,svid,phasr)
         else if (rcvrhw .eq. 'ROGUE   ') then
            call wb1380(luday,iwknf,sowclk,nsvtot,svid,phasr)

         else
            write (6,403) rcvrhw
 403        format (1x,'ERROR: unrecognized hardare = ',a8)
         endif
         do 420 i = ihold,ncnt
            iplace = locate(i)
c***            if (time(iplace) .ge. 0 .and. time(iplace) .le. sowclk) then
            if (  time(iplace) .ge. 0
     .      .and. dabs(time(iplace)-sowclk) .lt. 1.d0 ) then
               call wb9(luday,iplace,time(iplace),ephr,ittag)
               ibot = ibot + 1
            end if
 420     continue
         ihold = ibot
         sowold = sowclk
         go to 300
      else
c        If the days don't match increment the file day, close the old file,
c        open a new one, and repeat the test
 500     idoyf = idoyf + 1
CD        write(6,8001) idoyh,idoyp,idoyf
CD 8001     format(' at 500, idoyh,idoyp,idoyf=',3I5)
         if( idoyf.gt. (idoyh+6) ) goto 600
         if( idoyf.ne.idoyp ) goto 500
         if( ludopn ) then
            call wtime ( 6,iwknf,sowold,'gps',dfname//' ends at   ')
            call wtime (44,iwknf,sowold,'gps',dfname//' ends at   ')
            close( luday )
            do 510 i = 1,60
               if (tallyb(2,i).ne.0 .or. tallyb(3,i).ne.0 ) then
                  write ( 6,507) dfname,(tallyb(j,i),j=2,1,-1)
                  write (44,507) dfname,(tallyb(j,i),j=2,1,-1)
 507              format (1x,a15,'has ',i4,' FICA block ',i4,'s')
               endif
 510        continue
            call uptaly(0,0)
         endif 
          write(fhold,'(i3.3)') idoyf  
c        The following code used until 971024 will produce a name of the 
c        form sssswwwwy.ddd.  This was always incorrect, I think.
c         dfname = forb(1:leng-1)//fyear//'.'//fhold 
c        In any case, the apparent standard is ssssyddd.fic, which is
c        now coded below 
         dfname = forb(1:4)//fyear//fhold//'.fic'
         call openf(luday,dfname,'new','formatted','sequential')
         write (44,*) ' File name: ',dfname
         ludopn= .true.

c        tell the user about the start time of the file
         call wtime ( 6,iwknf,sowclk,'gps',dfname//' starts at ')
         call wtime (44,iwknf,sowclk,'gps',dfname//' starts at ')
c
c        WRITE HEADERS
         call wb0a(luday,argosv,rcvrhw,sitnam)
         call wb0(luday,vernum,rt,filnam)
         if (rcvrhw .eq. 'TI4100  ') then
            call wb101 (luday,swvrsn,sitnam)
         else if (rcvrhw .eq. 'MINIMAC ') then
            call wb1101 (luday,swvrsn,sitnam )
         else if (rcvrhw .eq. 'TRIMBLE ') then
            call wb1201 (luday,swvrsn,sitnam )
         else if (rcvrhw .eq. 'ROGUE   ') then
            call wb1301 (luday,swvrsn,sitnam )

         else
            write (6,513) rcvrhw
 513        format (1x,'ERROR: unrecognized hardare = ',a8)
         endif

c        write a BLK 9 first if one exists at the time of the first phase obs
         do 520 i = ihold,ncnt
           iplace = locate(i)
           if (time(iplace).ge.0.0d0 .and. time(iplace).le.sowclk) then
              call wb9(luday,iplace,time(iplace),ephr,ittag)
              ibot = ibot + 1
           end if
 520     continue
         ihold = ibot

c        then go write the ephemeris and phase records
         goto 400
      endif

c     all phase records are read and written -- see if there are more
c     orbit records to write

 550  do 560 i = ihold,ncnt
         iplace = locate(i)
         call wb9(luday,iplace,time(iplace),ephr,ittag)
 560  continue

C     write the final information on the last day file to XXXX.inf
c     close the files if all the records have been written

 600  continue
      call wtime ( 6,iwknf,sowold,'gps',dfname//' ends at   ')
      call wtime (44,iwknf,sowold,'gps',dfname//' ends at   ')
      do 610 i = 1,60
         if (tallyb(2,i).ne.0) then
            write ( 6,607) dfname,(tallyb(j,i),j=2,1,-1)
            write (44,607) dfname,(tallyb(j,i),j=2,1,-1)
 607        format (1x,a15,'has ',i4,' FICA block ',i4,'s')
         endif
 610  continue

      write ( 6,900)
      write (44,900)
 900  format (/,/,1x,'Normal stop on ARGO2FIC')


      close(44)
      close(luday)
      close(luph)
c
      stop
      end

c
      subroutine openf(lunit,nam,how,frm,acs)
c     open a file

      character*(*)   nam,how,frm,acs
      integer ios,lunit

         open (unit   =  lunit,
     &         file   =  nam,
     &         status =  how,
     &         form   =  frm,
     &         access =  acs,
     &         iostat =  ios)

      if (ios .ne. 0) then
         print 100,nam
 100     format (1x,'ERROR opening file: ',/,1x,a)
         call ferror(ios,6)
         stop
      else
 200     format (1x,'Opened: ',a)
         print 200,nam
      endif
      return
      end
c
      subroutine tsort(time,ncnt,loc)


      include '../includes/argo2fic.h'

      integer*4 loc(maxeph),ncnt,temp,i,j,jmin
      real*8 time(maxeph)

c     Sort array of time tags via a pointer array.
c     The time array is not actually changed.
c     The loc array contains the pointers to the sorted time tag data.

c     Load pointers
      do 10 i = 1,ncnt
         loc(i) = i
   10 continue
c
      do 30 i = 1,ncnt-1
         jmin = i
         do 20 j = i+1,ncnt
            if (time(loc(j)) .lt. time(loc(jmin)))jmin=j
 20      continue
         temp = loc(i)
         loc(i) = loc(jmin)
         loc(jmin) = temp
30    continue
      return
      end
c
      subroutine rdeph(lu,ittag,sec,ephr,icnt,ios)
c     Read CIGNET ephemeris data file
c     Note that the CIGNET format changed on week 424
c

      include '../includes/argo2fic.h'

      integer*4 ittag(maxeph,7),icnt,igpswk,i,iwke
      integer lu,ios,nerr,idoy,idoyd,iflag
      real*8 sec
      real*8 ephr(maxeph,27),soweph,utcoff

      save nerr
      data nerr/0/
c
      ios = 0
c
c     Read the first line containing the GPS week number
c     Do this only once.
      if (icnt .eq. 1) then
         read(lu,100,end=600,err=600,iostat=ios) ittag(1,1)
 100     format(i4)
c        New format blank line
         if (ittag(1,1).ge.424) read(lu,*)
      end if
c
c     Read the date information
c     Different format begins during gps week 424 (M. Chin)
c
 150  continue
      if (ittag(1,1).lt.424) then
         read(lu,200,end=600,err=500,iostat=ios)(ittag(icnt,i),i=2,7),
     $       sec,(ephr(icnt,i),i=1,19)
 200     format(/,i2,5i3,f5.1,3d19.12,4(/,3x,4d19.12))
      else
         read(lu,201,end=600,err=500,iostat=ios)(ittag(icnt,i),i=2,7),
     $       sec,(ephr(icnt,i),i=1,27)
 201     format(i2,5i3,f5.1,3d19.12,6(/,3x,4d19.12))
      endif

c     protect from bad week numbers
c     if bad, read another record, overwriting the bad one

      if (ittag(icnt,3) .gt. 0) then
         idoyd = idoy(ittag(icnt,3),ittag(icnt,4),ittag(icnt,5))
c        convert to GPST week number and second of week
         iflag = -4
         call timcon(iflag,
     .              igpswk,soweph,
     .              ittag(icnt,3),idoyd,ittag(icnt,6),ittag(icnt,7),sec,
     .              utcoff)

         if (ittag(1,1).lt.424) then
            iwke = igpswk
         else
            iwke = ephr(icnt,22)
         endif
      else
         iwke = 0
      endif

      if (iwke .ne. ittag(1,1)) then
         if (mod(nerr,20) .eq. 0) then
            print *,'Week number mismatch: ',igpswk,iwke
         endif
         nerr = nerr + 1
         goto 150
      endif

500   if (ios .ne. 0) then
        print *,'RDEPH: file error.'
        nerr = nerr + 1
        call ferror (ios,6)
        if (nerr .lt. 1000) then
           goto 150
        else
           call suicid('RDEPH: too many errors.')
        endif
      end if
600   return
      end
c

      subroutine rdphas(lu,idoyp,ttagr,sowclk,nsvtot,svid,phasr
     .  ,ios,nerr)

c     Read phase data one block at a time

c     modified by rwk 880809:  remove read of title record to main program
c
c     since CIGNET tapes come from NGS blocked as 3200 bytes, the
c     maximum number of skipped 80-byte lines after a bad block is 40.
c     We will allow 5 bad blocks, so 200 bad lines.
c
c
      include '../includes/argo2fic.h'

c     L1,L2,P1,P2,SNR1,SNR2 for each SV
      real*8 phasr(maxsvs,6)
c     PRN numbers
      integer*4 svid(maxsvs)

      integer*4 iyr,imo,iday,ihr,imn,nsvtot
      integer*4 i,j,idoy,idoyp
      integer lu,ios,nerr
      real*8 ttagr,sowclk
      character*80 buff80
c
c     read time tag line
  10  read(lu,'(a80)',end=1000,err=500,iostat=ios) buff80
      read(buff80,100,err=110,iostat=ios)
     .   iyr,imo,iday,ihr,imn,ttagr,sowclk,nsvtot,(svid(i),i=1,nsvtot)
      call check_y2k(iyr)
 100  format(i4,4i3,f11.7,f16.7,i4,11i3)
 110  if (ios.ne.0 .and. nerr .lt. 200) then
         write (6,112) nerr
         write (44,112) nerr
 112     format (1x,'ERROR # ',i4,' Bad line in phase file:')
         write (6,'(1x,a80)') buff80
         write (44,'(1x,a80)') buff80
         call ferror (ios,6)
         call ferror (ios,44)
         write (6,114)
         write (44,114)
 114     format (1x,'Continuing to next line....',/)
         nerr = nerr+1
         goto 10
      endif

c     calculate the day-of-year for the phase data
      idoyp = idoy(iyr,imo,iday)

c     read 1 line per SV
      do 300 j = 1,nsvtot
         read(lu,'(a80)',end=1000,err=500,iostat=ios) buff80
         read(buff80,400,err=410,iostat=ios) (phasr(j,i),i=1,6)
 400     format(4f17.5,2f5.1)
 410     if (ios.ne.0 .and. nerr .lt. 200) then
            write (6,412) nerr
            write (44,412) nerr
 412        format (1x,'ERROR # ',i4,' Bad line in phase file:')
            write (6,'(1x,a80)') buff80
            write (44,'(1x,a80)') buff80
            call ferror (ios,6)
            call ferror (ios,44)
            write (6,414)
            write (44,414)
 414        format (1x,'Setting data to 0.0 and continuing...',/)
            do 450 i=1,6
               phasr(j,i) = 0.0d0
 450        continue
            nerr = nerr+1
         endif
  300 continue

c     Handle nasty file errors.
 500  if (ios .ne. 0) then
        if (nerr .lt. 200) then
           write ( 6,512) nerr
           write (44,512) nerr
 512       format (1x,'RDPHAS: # ',i4,' Cannot read phase data file.')
           call ferror (ios,6)
           call ferror (ios,44)
           nerr=nerr+1
         else
           write ( 6,*) 'RDPHAS: Too many errors in phase file.'
           write (44,*) 'RDPHAS: Too many errors in phase file.'
           call suicid ('RDPHAS: Too many errors.')
         endif
      end if
1000  return
      end
c
      subroutine wb0(lu,vernum,rt,ficfil)

c     write a fica block zero with header info

      integer iblkid
      integer nf,ni,nc,lu,ioerr
C
c     run time yr,mo,da,hr,min,sec
      integer*4 rt(6),i
c
c     the FICA arrays
      real*8      ff(2)
      integer*4   fi(2)
      character*8 fc(8)

c     username
      character*16 uname
      character*8 ficfil,vernum

c     number of elements in each FIC array
      ioerr = 0
      iblkid = 0
      nf = 0
      ni = 0
      nc = 8

c     program name
      write(fc(1),1)
   1  format('ARGO2FIC')

c     program version
      fc(2) = vernum

c     institution
      write(fc(3),2)
   2  format('MIT GPS')

c     date of conversion
      write(fc(4),'(i4.4,2i2.2)') (rt(i),i=1,3)

c     time of conversion
      write(fc(5),10) (rt(i),i=4,5)
  10  format (1x,i2.2,':',i2.2,1x)

c     name of this file
      write(fc(6),'(a8)') ficfil

c     name of the user
      call getusr(uname)
      write(fc(7),'(a8)') uname

c     source format
      write (fc(8),'(a8)') 'CIGNET'

c     Write the header info to information file or
c     the current output day file
      if (lu .eq. 44) then
         write (lu,20) (fc(i),i=1,8)
  20     format (/,1x,'FICA header information:',/,8(1x,a8),/)
      else
         call wfica (lu,iblkid,ff,fi,fc,nf,ni,nc)
      end if

      return
      end
c
      subroutine wb0a(lu,argosv,rcvrhw,sitnam)

c     write a fica block zero with header info
c     for the NGS people

      integer iblkid
      integer nf,ni,nc,lu,ioerr
C
c     run time yr,mo,da,hr,min,sec
      integer*4 i
c
c     the FICA arrays
      real*8      ff(2)
      integer*4   fi(2)
      character*8 fc(8)

c     username
      character*16 uname
      character*8 rcvrhw,sitnam
      real*8      argosv

c     number of elements in each FIC array
      ioerr = 0
      iblkid = 0
      nf = 0
      ni = 0
      nc = 8

c     program name
      write(fc(1),1)
   1  format('ARGOS   ')

c     program version
      write(fc(2),'(f5.2)') argosv

c     institution
      write(fc(3),2)
   2  format(' NGS   ')

c     date of conversion
      write(fc(4),'(8X)')

c     time of conversion
      write(fc(5),'(8X)')

c     name of this file
      write(fc(6),'(a8)') sitnam

c     name of the user
      call getusr(uname)
      write(fc(7),'(8X)')

c     source format
      write (fc(8),'(a8)') rcvrhw

c     Write the header info to information file or
c     the current output day file
      if (lu .eq. 44) then
         write (lu,20) (fc(i),i=1,8)
  20     format (/,1x,'FICA header information:',/,8(1x,a8),/)
      else
         call wfica (lu,iblkid,ff,fi,fc,nf,ni,nc)
      end if

      return
      end
c


      subroutine wb9(lu,it,sowclk,ephr,ittag)
c
c  Write a FICA block 9 with ephemeris information.
c
c  EXPLANATION
c
c  Because the data in the NGS format has been decimated and 'filtered'
c  certain pieces of information that the FICA format utilizes are not
c  present.  These are the following:
c
c        TLM word (preamble)
c        TLM word (message
c        HOW word (time)
c        DATA id
c        C/A and/or P flag, L2 flag
c        SV accuracy
c        SV health
c        Age of Data clock
c        L2 P data flag
c        Group delay differential
c        Tracker
c        Inclination time derivative
c
c  DUMMIED VARIABLES
c
c  FICA BLOCK 9
c  map of the floating point items:
c    ff(i)      is                                             DUMMIED TO :
c    (sub frame 1)
c      1     TLM word (preamble)            dimless                 0.0d0
c      2     TLM word (message)             dimless                 0.0d0
c      3     HOW word (time)                GPS sec of week         sowclk (clock epoch time tag)
c      4     Data ID                        dimless                 0.0d0
c      5     Sub-frame ID (1)               dimless                 same
c      6     full week number (10 bit)      GPS weeks               same
c      7     C/A and/or P flag, L2 flag     dimless                 0.0d0
c      8     SV accuracy                    dimless                 0.0d0
c      9     SV health                      dimless                 0.0d0
c      10    Age of Data clock              sec                     0.0d0
c      11    L2 P data flag                 dimless                 0.0d0
c      12    group delay differential       sec                     0.0d0
c      13    clock epoch                    GPS sec of week         same
c      14    clock drift rate               sec/(sec**2)            same
c      15    clock drift                    sec/sec                 same
c      16    clock bias                     sec                     same
c      17    (not used)                                             same
c      18    (not used)                                             same
c      19    Tracker                        dimless                 0.0d0
c      20    SV PRN                         dimless                 same
c     (subframe 2)
c      21    TLM word preamble              dimless                 0.0d0
c      22    TLM word message               dimless                 0.0d0
c      23    HOW word (time)                GPS sec of week         sowclk (clock epoch time tag)
c      24    Data id                        dimless                 0.0d0
c      25    Sub-frame id should be 2       dimless                 same
c      26    Age of data (ephemeris)        GPS sec of week         same
c      27    Radial sine correction (CRS)   meters                  same
c      28    Correction to mean motion      radians/sec             same
c      29    Mean anomaly at epoch          radians                 same
c      30    in-track cosine amp. (CUC)     radians                 same
c      31    eccentricity                  dimless                  same
c      32    in-track sine amp (CUS)        radians                 same
c      33    square root of sem-major axis  (meters)**1/2           same
c      34    time of epoch                  GPS sec of week         same
c      35    fit interval flag            mless                     same
c      36    unused (set to 0.0d0)                                   "
c      37    unused (set to 0.0d0)                                   "
c      38    unused (set to 0.0d0)                                   "
c      39    unused (set to 0.0d0)                                   "
c      40    unused (set to 0.0d0)                                   "
c     (subframe 3)
c      41    TLM word preamble              dimless                 0.0d0
c      42    TLM word message               dimless                 0.0d0
c      43    HOW word (time)                dimless                 0.0d0
c      44    Data id                        dimless                 0.0d0
c      45    Sub-frame id (should be 3)     dimless                 same
c      46    inclination cosine correction CIC rads                 same
c      47    right ascension of ascending node  rads                same
c      48    incl. sine correction (CIS)    rads                    same
c      49    inclination                    rads                    same
c      50    radial cosine adj              rads                    same
c      51    argument of perigee            rads                    same
c      52    right ascension of ascending node  rads/sec            same
c      53    age of data (ephemeris)       GPS sec of week          same
c      54    inclination time derivative    rads/sec                0.0d0
c      55    unused (set to 0.0d0)                                   "
c      56    unused (set to 0.0d0)                                   "
c      57    unused (set to 0.0d0)                                   "
c      58    unused (set to 0.0d0)                                   "
c      59    unused (set to 0.0d0)                                   "
c      60    unused (set to 0.0d0)                                   "
c

      include '../includes/argo2fic.h'
      integer*4 ittag(maxeph,7),it,i,fi(2)
      integer*4 nf,ni,nc,iblkid
      integer lu
      real*8 ff(60),ephr(maxeph,27),sowclk
      character*8 fc(2)
c
c  begin rearrangement and cramming of data
c  only deal with non-zeroed data
      do 10 i = 1,60
        ff(i) = 0.0d0
 10   continue

c  NOTE: the important time tag throughout ARGO2FIC's code is going to be
c        clock epoch time tag ff(13); this will cause confusion
c        unless comparisons of data are done on this time tag, NOT the HOW
c        word time (ff(3),ff(23),ff(43))
      ff(3) = sowclk
      ff(5) = 1.0d0
      ff(6) = dble(ittag(1,1))
      ff(13) = sowclk
      ff(14) = ephr(it,3)
      ff(15) = ephr(it,2)
      ff(16) = ephr(it,1)
      ff(20) = dble(ittag(it,2))
      ff(25) = 2.0d0
      do 20 i=26,34
         ff(i) = ephr(it,i-22)
 20   continue
      ff(45) = 3.0d0
      ff(46) = ephr(it,13)
      ff(47) = ephr(it,14)
      ff(48) = ephr(it,15)
      ff(49) = ephr(it,16)
      ff(50) = ephr(it,17)
      ff(51) = ephr(it,18)
      ff(52) = ephr(it,19)
      ff(53) = ephr(it,4)
c
c     New CIGNET format (including more ephemeris parameters)
c     was implemented begininging GPS week 426
c     Handle time-dependent file format here.
c
      if (ittag(1,1).ge.424) then
         ff(54) = ephr(it,20)
         ff(7 ) = ephr(it,21)
         ff(6 ) = ephr(it,22)
         ff(11) = ephr(it,23)
         ff(8 ) = ephr(it,24)
         ff(9 ) = ephr(it,25)
         ff(12) = ephr(it,26)
         ff(10) = ephr(it,27)
      endif
c
c     write the fica block data to the current output day file
      iblkid = 9
      nf = 60
      ni = 0
      nc = 0
      call wfica (lu,iblkid,ff,fi,fc,nf,ni,nc)

      return
      end
c
      subroutine wb55(lu,iwknf,sowclk,nsvtot,svid,phasr,swvrsn)
c
c  Write a FICA block 55 with TI4100 (CORE) phase data
c
c
c  EXPLANATION
c
c  NGS format data does not include some of the information that
c  FICA format uses.  These variables must be dummied.  They
c  are the following:
c
c        Receiver internal temperature
c        Tracker mode
c        L1, L2 quality vectors (tracker, frequency)
c
c  DUMMIED VARIABLES
c
c  FICA BLOCK 55
c  map of the floating point items:                           DUMMIED TO :
c  ff(i)    is
c    1      pseudorange of FTF validity   sec                  0.0D0
c    2      FTF bias offset               sec                  0.0D0
c    3      user epoch time of pseudorng  GPS sec of week      same
c    4-11   L1,L2 carrier signal to noise db-Hz                same
c           tracker, frequency
c    12-15  L1 pseudorange                kilometers           same
c    16-19  L2 pseudorange                kilometers           same
c    20-23  L1 carrier Doppler phase at   cycles               same
c           code FTF
c    24-27  L2 carrier Doppler phase at   cycles               same
c           code FTF
c    28-31  L1 carrier velocity at code   Hz                   0.0d0
c           FTF
c    32-35  L2 carrier velocity at code   Hz                   0.0d0
c           FTF
c    36-39  averaged line of site accel   m/s**2               0.0d0
c    40     recvr internal temperature    deg C                0.0d0
c  map of the integer items:
c  fi(i)  is:                        units
c    1-4    SV PRN of each tracker        dimless              same
c    5-8    tracker mode                  dimless               -1
c    9-16   L1,L2 quality vector          dimless           16#FFFF0000 (hex)
c           (tracker,frequency)
c

      include '../includes/argo2fic.h'

c     L1,L2,P1,P2,SNR1,SNR2 for each SV
      real*8 phasr(maxsvs,6)
c     PRN numbers
      integer*4 svid(maxsvs)

      integer*4 nsvtot,fi(16),nf,ni,nc,iblkid,iwknf
      integer*4 i
      integer lu
      real*8 ff(40),sowclk,swvrsn
      character*8 fc(2)

c     Begin rearrangement of the data
c     only deal with the non-zeroed (non-dummied) data
c     first the integer data
      do 10 i=1,16
         fi(i) = 0
 10   continue
      do 20 i = 1,nsvtot
         fi(i) = svid(i)
 20   continue
c
c     dummy the tracker mode
c     NOTE: -1 will ensure that the tracker mode is
c     interpreted to be safe
      do 30 i = 5,8
         fi(i) = -1
 30   continue

c     dummy the L1 and L2 quality vectors (used to filter the data)
c     NOTE: these flags are bit-mapped, thus
c        -1 is not a good dummy value since it will
c        turn on all error bit flags.
c        To denote health, this is dummied 16#FFFF0000 = -65536
c        so that the significant byte is mapped to zero
c        and the remaining two bytes are used as the flag for health
      do 33 i = 9,16
         fi(i) = -65536
 33   continue
c     zero out the floating point array
      do 35 i = 1,40
         ff(i) = 0.0d0
 35   continue

c     Floating point data.
c     In late 1988, NGS changed the CIGNET format.
c     Data in the new format begins with week 456.
c     (**Note: This changed from 455 by rwk based on Burc Oral's
c              data for week 455;  it may be later still, since our
c              next data set is GOTEX--week 459.)
c     So, we asssume:
c     GPS week .lt. 456: CIGNET format has biased phases.
c           L1 + 6000 * dt
c           L2 - 7600 * dt
c     GPS week .ge. 456: CIGNET format has negative unbiased phases.
c           -L1
c           -L2
c     Since, however, FICA block 55's are defined as
c     being biased phases, we will write (in all cases):
c           L1 + 6000 * dt
c           L2 - 7600 * dt
c     Note, that the origin time for calculating dt
c     is not really defined.  So, for pre-week 456, we
c     keep whatever TI uses for dt, but for post-456,
c     we assume that dt is from the beginning of the GPS week.

      ff(3) = sowclk

c     Signal to noise ratios
      do 43 i = 1,nsvtot
         ff(i+3) = phasr(i,5)
         ff(i+7) = phasr(i,6)
 43   continue

c     phases, old format
      if (iwknf .lt. 456) then
         do 45 i = 1,nsvtot
            ff(i+19) = phasr(i,1)
            ff(i+23) = phasr(i,3)
 45      continue
c     CORE 4.12 phases, an incorrect implementation of the
c     frequency plan which applies to Yellowknife data circa
c     March 1990 (week 533)
      else if (swvrsn .eq. 4.12) then
         do 46 i = 1,nsvtot
            ff(i+19) = phasr(i,1) - 6000.0D0 * sowclk
            ff(i+23) = phasr(i,3) + 7600.0D0 * sowclk
 46      continue
c     phases, new format
      else
         do 47 i = 1,nsvtot
            ff(i+19) = -1.0d0 * phasr(i,1) - 6000.0D0 * sowclk
            ff(i+23) = -1.0d0 * phasr(i,3) + 7600.0D0 * sowclk
 47      continue
      endif

c     Pseudoranges must be output in kilometers ff(12-19)
      do 50 i = 1,nsvtot
         ff(i+11) = phasr(i,2)/1000.
         ff(i+15) = phasr(i,4)/1000.
 50   continue
c
c     write the block 55 data to the current data file
      iblkid = 55
      nf = 40
      ni = 16
      nc = 0
      call wfica (lu,iblkid,ff,fi,fc,nf,ni,nc)

      return
      end

      subroutine wb101( lu,swvrsn,sitnam )
c     R.W. King from Kurt Feigl code in X2FIC -- 880809

c     FICA arrays
      real*8 ff(5)
      integer*4 fi(230)
      character*8 fc(14)
      character*20 sitnam
      character*4 blank4
      integer*4 nf,ni,nc,itype,i,lu
      real*8 swvrsn             
      logical debug 

      data blank4/'    '/,debug/.false./

c     write some things to a FICA block 101, which usually
c     contains GESAR input data
      nf = 5
      ni = 230
      nc = 14

c     zero everything, overwrite some entries
      do i=1,ni
         fi(i) = 0
      enddo
      do i=1,nf
         ff(i) = 0.0d0
      enddo
      do i=1,nc
         fc(i) = blank4//blank4
      enddo

c     collection rate now is X-file interval
c     ff(4) = dble(inter)

c     version number *10
      fi(1) = nint(swvrsn*10)

c     'current date and time' is initial epoch
c     fi(4) = iy
c     fi(5) = im
c     fi(6) = id
c     fi(7) = ihr
c     fi(8) = min

c     site name      
c     fc(4) = sitnam(1:8)
c     fc(5) = sitnam(9:12)//blank4
      if(debug) then
        write(*,'(a20)') sitnam
      endif

c     user code
      fc(1) = 'ARGO2FIC'

c     write the block 101 to the FICA file
      itype = 101
      call wfica(lu,itype,ff,fi,fc,nf,ni,nc)

      return
      end

c
      subroutine wb1101(lu,swvrsn,sitnam )

c
c     Write block 1101
c
c     Translates Minimac station information into FICA block 1101 format
c
c--------------------------------------------------------------------------------
c
c    Provisional MIT Definition of reduced FICA MINIMAC Data Block  Feigl 89-05-01
c
c   FIC Block Name         : Station information
c   FIC Block Number       : 1180
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


c     FICA arrays
      real*8 ff(7)
      integer*4 fi(2)
      character*8 fc(8)
      character*(*) sitnam
      character*32 buff32
      integer*4 nf,ni,nc,itype,lu
      real*8 swvrsn


c     write some things to a FICA block 101, which usually
c     contains GESAR input data
      nf = 7
      ni = 0
      nc = 8

c     allow 32 characters of site name
      write (buff32,'(a32)') sitnam
      fc(1) = buff32(1:8)
      fc(2) = buff32(9:16)
      fc(3) = buff32(17:24)
      fc(4) = buff32(25:32)

c     Operator name?
      write (buff32,'(a32)') 'WARNING: TEST BLOCK.....'
      fc(5) = buff32(1:8)
      fc(6) = buff32(9:16)
      fc(7) = buff32(17:24)
      fc(8) = buff32(25:32)

c     Software version number
      ff(1) = swvrsn

c     L1 phase center above mark?
      ff(2) = 0.0d0
c     L2 phase center above mark?
      ff(3) = 0.0d0
c     antenna base above mark?
      ff(4) = 0.0d0
c     Antenna east offset
      ff(5) = 0.0d0
c     Antenna north offset
      ff(6) = 0.0d0

c     Data collection interval
      ff(7) = 30.0d0

c     write the block 101 to the FICA file
      itype = 1101
      call wfica(lu,itype,ff,fi,fc,nf,ni,nc)

      return
      end
      subroutine wb1180(lu,iwknf,sowclk,nsvtot,svid,phasr)
c
c     Write block 1180
c
c     Translates Minimac phase data into FICA block 1180 format
c
c--------------------------------------------------------------------------------
c
c    Provisional MIT Definition of reduced FICA MINIMAC Data Block  Feigl 89-05-01
c
c   FIC Block Name         : Tracking data
c   FIC Block Number       : 1180
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
c                          (positive sign for decreasing Pseudorange)
c   6 (10, 14, 18, ...)    L2 carrier Doppler phase                half cycles
c                          (positive sign for increasing Doppler shift)
c                          (positive sign for decreasing Pseudorange)
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


      include '../includes/argo2fic.h'

c     L1,L2,P1,P2,SNR1,SNR2 for each SV
      real*8 phasr(maxsvs,6)
c     PRN numbers
      integer*4 svid(maxsvs)
c     The FICA arrays
      real*8 ff(4+4*maxsvs)
      integer*4 fi(6+3*maxsvs)
      character*8 fc(2)

      integer*4 nsvtot,nf,ni,nc,iblkid,iwknf
      integer*4 i,iflag
c     Logical unit and I/O error (may be system dependent)
      integer*4 lu
      real*8      sowclk,utcoff

c     number of integer items
      ni = 6 + 3 * nsvtot

c     number of floating point items
      nf = 4 + 4 * nsvtot

c     number of character items
      nc = 0

c     zero the arrays
      do 10 i=1,ni
         fi(i) = 0
 10   continue
      do 20 i = 1,nf
         ff(i) = 0.0d0
 20   continue

c     number of channels
      fi(1) = nsvtot

c     TIME TAG CORRECTION
c     The GPS bulletin (July-August 1989 p. 14) says:
c     Starting GPS week 496 (July 9, 1989) at 00:00 Minimac 2816 AT
c     recievers have installed the latest processor firmware (version
c     RP 1908 and NP 1906) and MMAT software (version 1.61)
c     Measurement time tag (still in UTC system) has moved by 1 millisecond.
c     Previous data required one millisecond to be added to align to the
c     GPS second

c     BUT: The NGS translation program (ARGOS) converts from UTC to GPST.
c
c     SO: The old data (before week 496) were really collected at
c     times like xx.001 seconds
c     Here, we add 1 millisecond for data before week 496
      if (iwknf .lt. 496) then
         call secsum (iwknf,sowclk,1.0d-03,iwknf,sowclk)
      endif

c     second of GPS week and GPS week number (GPST)
      ff(1) = sowclk
      fi(2) = iwknf

c     convert the time tag from week number and second of week (GPST)
c     into year, day-of-year, hours, minutes, seconds (GPST)

      iflag = 4
      call timcon (iflag,
     .             iwknf,sowclk,
     .             fi(3),fi(4),fi(5),fi(6),ff(2),
     .             utcoff)


      do 30 i=1,nsvtot
c        PRN number
         fi(7+3*(i-1)) = svid(i)
c        amplitudes
         fi(8+3*(i-1)) = nint(phasr(i,5))
         fi(9+3*(i-1)) = nint(phasr(i,6))

c        Phases
c        BLOCK 1180s follow DOPPLER (not PSEUDORANGE) CONVENTION.
c        Positive sign for increasing Doppler shift
c        Positive sign for decreasing Pseudorange
c        L1 phase is in full cycles in both CIGNET and FICA
         ff(5+4*(i-1)) = -1.0d0 * phasr(i,1)
c        L2 phase is in full cycles in CIGNET,
c        should be in half-cycles in FICA
         ff(6+4*(i-1)) = -1.0d0 * phasr(i,3)*2.0d0

c        L1 pseudorange in kilometers
         ff(7+4*(i-1)) = phasr(i,2)/1000.0d0
c        L2 pseudorange is not observed
         ff(8+4*(i-1)) = phasr(i,4)
 30   continue
c
c     write the FICA block to the current file
      iblkid = 1180
      call wfica (lu,iblkid,ff,fi,fc,nf,ni,nc)

      return
      end
c

      subroutine wb1201(lu,swvrsn,sitnam )

c
c     Write block 1201
c
c     Translates Trimble station information into FICA block 1201 format
c
c--------------------------------------------------------------------------------
c
c    Provisional MIT Definition of reduced FICA Trimble Data Block  King 90-12-19
c    (exactly the same as the Minimac 1101 block of Feigl 89-05-01)
c
c   FIC Block Name         : Station information
c   FIC Block Number       : 1201
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


c     FICA arrays
      real*8 ff(7)
      integer*4 fi(2)
      character*8 fc(8)
      character*(*) sitnam
      character*32 buff32
      integer*4 nf,ni,nc,itype,lu
      real*8 swvrsn


c     write some things to a FICA block 101, which usually
c     contains GESAR input data
      nf = 7
      ni = 0
      nc = 8

c     allow 32 characters of site name
      write (buff32,'(a32)') sitnam
      fc(1) = buff32(1:8)
      fc(2) = buff32(9:16)
      fc(3) = buff32(17:24)
      fc(4) = buff32(25:32)

c     Operator name?
      write (buff32,'(a32)') 'WARNING: TEST BLOCK.....'
      fc(5) = buff32(1:8)
      fc(6) = buff32(9:16)
      fc(7) = buff32(17:24)
      fc(8) = buff32(25:32)

c     Software version number
      ff(1) = swvrsn

c     L1 phase center above mark?
      ff(2) = 0.0d0
c     L2 phase center above mark?
      ff(3) = 0.0d0
c     antenna base above mark?
      ff(4) = 0.0d0
c     Antenna east offset
      ff(5) = 0.0d0
c     Antenna north offset
      ff(6) = 0.0d0

c     Data collection interval
      ff(7) = 30.0d0

c     write the block 1201 to the FICA file
      itype = 1201
      call wfica(lu,itype,ff,fi,fc,nf,ni,nc)

      return
      end
      subroutine wb1280(lu,iwknf,sowclk,nsvtot,svid,phasr)
c
c     Write block 1280
c
c     Translates Trimble phase data into FICA block 1280 format
c
c--------------------------------------------------------------------------------
c
c    Provisional MIT Definition of reduced FICA Trimble Data Block  King 90-12-19
c      (exactly the same as Minimac 1101 of Feigl 89-05-01)
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
c                          (positive sign for decreasing Pseudorange)
c   6 (10, 14, 18, ...)    L2 carrier Doppler phase                half cycles
c                          (positive sign for increasing Doppler shift)
c                          (positive sign for decreasing Pseudorange)
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


      include '../includes/argo2fic.h'

c     L1,L2,P1,P2,SNR1,SNR2 for each SV
      real*8 phasr(maxsvs,6)
c     PRN numbers
      integer*4 svid(maxsvs)
c     The FICA arrays
      real*8 ff(4+4*maxsvs)
      integer*4 fi(6+3*maxsvs)
      character*8 fc(2)

      integer*4 nsvtot,nf,ni,nc,iblkid,iwknf
      integer*4 i,iflag
c     Logical unit and I/O error (may be system dependent)
      integer*4 lu
      real*8      sowclk,utcoff

c     number of integer items
      ni = 6 + 3 * nsvtot

c     number of floating point items
      nf = 4 + 4 * nsvtot

c     number of character items
      nc = 0

c     zero the arrays
      do 10 i=1,ni
         fi(i) = 0
 10   continue
      do 20 i = 1,nf
         ff(i) = 0.0d0
 20   continue

c     number of channels
      fi(1) = nsvtot

c     second of GPS week and GPS week number (GPST)
      ff(1) = sowclk
      fi(2) = iwknf

c     convert the time tag from week number and second of week (GPST)
c     into year, day-of-year, hours, minutes, seconds (GPST)

      iflag = 4
      call timcon (iflag,
     .             iwknf,sowclk,
     .             fi(3),fi(4),fi(5),fi(6),ff(2),
     .             utcoff)


      do 30 i=1,nsvtot
c        PRN number
         fi(7+3*(i-1)) = svid(i)
c        amplitudes
         fi(8+3*(i-1)) = nint(phasr(i,5))
         fi(9+3*(i-1)) = nint(phasr(i,6))

c        Phases
c        BLOCK 1280s follow DOPPLER (not PSEUDORANGE) CONVENTION.
c        Positive sign for increasing Doppler shift
c        Positive sign for decreasing Pseudorange
c        L1 phase is in full cycles in both CIGNET and FICA
         ff(5+4*(i-1)) = -1.0d0 * phasr(i,1)
c        L2 phase is in full cycles in CIGNET,
c        should be in half-cycles in FICA
         ff(6+4*(i-1)) = -1.0d0 * phasr(i,3)*2.0d0

c        L1 pseudorange in kilometers
         ff(7+4*(i-1)) = phasr(i,2)/1000.0d0
c        L2 pseudorange is not observed
         ff(8+4*(i-1)) = phasr(i,4)
 30   continue
c
c     write the FICA block to the current file
      iblkid = 1280
      call wfica (lu,iblkid,ff,fi,fc,nf,ni,nc)

      return
      end

      subroutine wb1301(lu,swvrsn,sitnam )

c
c     Write block 1301
c
c     Translates Rogue station information into FICA block 1301 format
c
c--------------------------------------------------------------------------------
c
c    Provisional MIT Definition of reduced FICA Rogue Data Block  King 90-12-28
c    (exactly the same as the Minimac 1101 block of Feigl 89-05-01)
c
c   FIC Block Name         : Station information
c   FIC Block Number       : 1301
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
c   1                      Rogue   software version                none
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


c     FICA arrays
      real*8 ff(7)
      integer*4 fi(2)
      character*8 fc(8)
      character*(*) sitnam
      character*32 buff32
      integer*4 nf,ni,nc,itype,lu
      real*8 swvrsn


c     write some things to a FICA block 1301, which usually
c     contains GESAR input data
      nf = 7
      ni = 0
      nc = 8

c     allow 32 characters of site name
      write (buff32,'(a32)') sitnam
      fc(1) = buff32(1:8)
      fc(2) = buff32(9:16)
      fc(3) = buff32(17:24)
      fc(4) = buff32(25:32)

c     Operator name?
      write (buff32,'(a32)') 'WARNING: TEST BLOCK.....'
      fc(5) = buff32(1:8)
      fc(6) = buff32(9:16)
      fc(7) = buff32(17:24)
      fc(8) = buff32(25:32)

c     Software version number
      ff(1) = swvrsn

c     L1 phase center above mark?
      ff(2) = 0.0d0
c     L2 phase center above mark?
      ff(3) = 0.0d0
c     antenna base above mark?
      ff(4) = 0.0d0
c     Antenna east offset
      ff(5) = 0.0d0
c     Antenna north offset
      ff(6) = 0.0d0

c     Data collection interval
      ff(7) = 30.0d0

c     write the block 1301 to the FICA file
      itype = 1301
      call wfica(lu,itype,ff,fi,fc,nf,ni,nc)

      return
      end


      subroutine wb1380(lu,iwknf,sowclk,nsvtot,svid,phasr)
c
c     Write block 1380
c
c     Translates Rogue phase data into FICA block 1380 format
c
c--------------------------------------------------------------------------------
c
c    Provisional MIT Definition of reduced FICA Rogue Data Block  King 90-12-28
c      (same as Minimac 1180 of Feigl 89-05-01 except that P2 added)
c
c   FIC Block Name         : Tracking data
c   FIC Block Number       : 1380
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
c                          (positive sign for decreasing Pseudorange)
c   6 (10, 14, 18, ...)    L2 carrier Doppler phase                half cycles
c                          (positive sign for increasing Doppler shift)
c                          (positive sign for decreasing Pseudorange)
c   7 (11, 15, 19, ...)    L1 pseudorange                          kilometers
c   8 (12, 16, 20, ...)    L1 pseudorange                          kilometers
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


      include '../includes/argo2fic.h'

c     L1,L2,P1,P2,SNR1,SNR2 for each SV
      real*8 phasr(maxsvs,6)
c     PRN numbers
      integer*4 svid(maxsvs)
c     The FICA arrays
      real*8 ff(4+4*maxsvs)
      integer*4 fi(6+3*maxsvs)
      character*8 fc(2)

      integer*4 nsvtot,nf,ni,nc,iblkid,iwknf
      integer*4 i,iflag
c     Logical unit and I/O error (may be system dependent)
      integer*4 lu
      real*8      sowclk,utcoff

c     number of integer items
      ni = 6 + 3 * nsvtot

c     number of floating point items
      nf = 4 + 4 * nsvtot

c     number of character items
      nc = 0

c     zero the arrays
      do 10 i=1,ni
         fi(i) = 0
 10   continue
      do 20 i = 1,nf
         ff(i) = 0.0d0
 20   continue

c     number of channels
      fi(1) = nsvtot

c     second of GPS week and GPS week number (GPST)
      ff(1) = sowclk
      fi(2) = iwknf

c     convert the time tag from week number and second of week (GPST)
c     into year, day-of-year, hours, minutes, seconds (GPST)

      iflag = 4
      call timcon (iflag,
     .             iwknf,sowclk,
     .             fi(3),fi(4),fi(5),fi(6),ff(2),
     .             utcoff)


      do 30 i=1,nsvtot
c        PRN number
         fi(7+3*(i-1)) = svid(i)
c        amplitudes
         fi(8+3*(i-1)) = nint(phasr(i,5))
         fi(9+3*(i-1)) = nint(phasr(i,6))

c        Phases
c        BLOCK 1380s follow DOPPLER (not PSEUDORANGE) CONVENTION.
c        Positive sign for increasing Doppler shift
c        Positive sign for decreasing Pseudorange
c        L1 phase is in full cycles in both CIGNET and FICA
         ff(5+4*(i-1)) = -1.0d0 * phasr(i,1)
c        L2 phase is in full cycles in CIGNET,
c        should be in half-cycles in FICA
         ff(6+4*(i-1)) = -1.0d0 * phasr(i,3)*2.0d0

c        L1 pseudorange in kilometers
         ff(7+4*(i-1)) = phasr(i,2)/1000.0d0
c        L2 pseudorange is not observed
         ff(8+4*(i-1)) = phasr(i,4)/1000.0d0
 30   continue
c
c     write the FICA block to the current file
      iblkid = 1380
      call wfica (lu,iblkid,ff,fi,fc,nf,ni,nc)

      return
      end
c
      subroutine wfica(lun,itype,fa,ia,ca,nf,ni,nc)
c
c     write a file in fica (float, int, char (ascii)) format
c     modified from judah levine's wtfica by kf mar 5, 87
c     modified by k feigl, jan 17, 1988 to deal with
c     variable logical unit number (lun)
c
c     this subroutine writes a fica record.  itype is the
c     type, nf,ni and nc are the number of floating variables,
c     integer variables and character variables, respectively.
c
c     fa, ia and ca are the arrays of floating variables, integer
c     variables and character variables, respectively.
c
      integer*4 lun,i
      integer*4 itype,nf,ni,nc
      real*8 fa(nf)
      integer*4 ia(ni)
      character*8 ca(nc)

      integer*2 tallyb(4,60)
      common /tally/ tallyb

c
c
c
      write(lun,1) itype,nf,ni,nc
    1 format('BLK  ',4i5)
      if(nf .gt. 0) write(lun,2)(fa(i),i=1,nf)
    2 format(4(1pd20.13))
      if(ni .gt. 0) write(lun,3)(ia(i),i=1,ni)
    3 format(6i12)
      if(nc .gt. 0) write(lun,4)(ca(i),i=1,nc)
    4 format(10a8)

      call uptaly(itype,1)

      return
      end


      subroutine uptaly(blocknum, itype)
c     update the read, write tally
      integer itype
C      0 to initialize the counts
c
C      1 to increment seen count
c
C      2 to increment read count
c
C      3 to increment written count
c
      integer blocknum,i,nmax,j
      integer*2 tallyb(4,60)
      common /tally/ tallyb

      nmax = 55

      if (itype.eq.-1) then
         do i = 1,nmax
            if (tallyb(1,i) .eq. blocknum) then
               itype = i
               return
            endif
         enddo
C     initialize everything
      else if (itype.eq.0) then
         tallyb(1,1) = 1
         tallyb(1,2) = 2
         tallyb(1,3) = 3
         tallyb(1,4) = 6
         tallyb(1,5) = 7
         tallyb(1,6) = 8
         tallyb(1,7) = 9
         tallyb(1,8) = 10
         tallyb(1,9) = 11
         tallyb(1,10) = 12
         tallyb(1,11) = 13
         tallyb(1,12) = 50
         tallyb(1,13) = 51
         tallyb(1,14) = 52
         tallyb(1,15) = 53
         tallyb(1,16) = 54
         tallyb(1,17) = 55
         tallyb(1,18) = 56
         tallyb(1,19) = 57
         tallyb(1,20) = 58
         tallyb(1,21) = 59
         tallyb(1,22) = 62
         tallyb(1,23) = 63
         tallyb(1,24) = 1001
         tallyb(1,25) = 1002
         tallyb(1,26) = 1003
         tallyb(1,27) = 1004
         tallyb(1,28) = 1005
         tallyb(1,29) = 1006
         tallyb(1,30) = 1007
         tallyb(1,31) = 1008
         tallyb(1,32) = 1009
         tallyb(1,33) = 101
         tallyb(1,34) = 102
         tallyb(1,35) = 109
         tallyb(1,36) = 162
         tallyb(1,37) = 400
         tallyb(1,38) = 401
         tallyb(1,39) = 402
         tallyb(1,40) = 403
         tallyb(1,41) = 404
         tallyb(1,42) = 405
         tallyb(1,43) = 410
         tallyb(1,44) = 411
         tallyb(1,45) = 420
         tallyb(1,46) = 421
         tallyb(1,47) = 423
         tallyb(1,48) = 424
         tallyb(1,49) = 425
         tallyb(1,50) = 426
         tallyb(1,51) = 70
         tallyb(1,52) = 80
         tallyb(1,53) = 1180
         tallyb(1,54) = 1101
         tallyb(1,55) = 0
         do i = 2,4
             do j = 1,nmax
                tallyb(i,j) = 0
             enddo
         enddo
      else
         do i = 1,nmax
            if (tallyb(1,i) .eq. blocknum) then
               tallyb(itype+1,i) = tallyb(itype+1,i) + 1
               return
            endif
         enddo
      endif
      return
      end



