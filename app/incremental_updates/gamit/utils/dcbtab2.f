      Program DCBTAB2

c     Create or update a Version 2 GAMIT dcb.dat file

c     Command-line arguments control the mode:
c       First argument is the input file, either merged CODE dcb file
c         (special format) to be translated to GAMIT Version 2 dcb.dat,
c         or an existing Version 2 dcb.dat to be updated
c       Second argument, if given, is the COD P1C1 file to be added to get
c         dcb.dat.gnss.new.  if omitted, the program will translate the special
c         -format CODE file, with the output file named by appending 'new' to
c         the first argument.
c

c     R. King 19 May 2015
c     Modified to make input file names more general. 8 August 2016

      implicit none

* nl    - number of total values in the arrays/dcb.dat
* axnl - dimension of the arrays -- change also in get_times, insert_arrays, sortsvn

      integer*4 maxnl,nl,ng,nr,nc,ne
      parameter(maxnl= 30000)

* Array values associated with each entry

      character*1 gnss(maxnl)
      integer*4 svn(maxnl),prn(maxnl),start(5,maxnl),stop(5,maxnl)
     .        , svng(maxnl),prng(maxnl),startg(5,maxnl),stopg(5,maxnl)
     .        , svnr(maxnl),prnr(maxnl),startr(5,maxnl),stopr(5,maxnl)
     .        , svne(maxnl),prne(maxnl),starte(5,maxnl),stope(5,maxnl)
     .        , svnc(maxnl),prnc(maxnl),startc(5,maxnl),stopc(5,maxnl)
*     Although we write only yr doy hr min, make the time arrays 5 for
*     for consistency with subroutine itimdif (and other GAMIT routines)

      real*4 dcb(maxnl),rms(maxnl)
     .     , dcbg(maxnl),rmsg(maxnl)
     .     , dcbr(maxnl),rmsr(maxnl)
     .     , dcbe(maxnl),rmse(maxnl)
     .     , dcbc(maxnl),rmsc(maxnl)


* Pointers to SVN within the arrays

      integer*4 indx(maxnl)

* Input and output files

* infile - current dcb.dat file, either version 1 of version 2
* monfile - monthly COD P1C1 file to be added to the input
* outfile - new version 2 dcb.dat file, named by appending 'new' to enfile
      character*20 infile,monfile,outfile
      integer*4 luin,lumon,luout
      parameter(luin=1,lumon=2,luout=3)

* Other variables

      character*1 ag
      character*80 line
      character*256 message
      integer*4 iprn,isvn,year,month,day,doy
     .        , startyr,startdoy,dcbstart(5)
     .        , stopyr,stopdoy,dcbstop(5),svnlast
     .        , irunt(3),iarg,iclarg,ioerr,i,j,nblen
      logical eof,endsvs,debug/.false./

* Function
      integer*4 idoy
      integer*8 itimdif
      logical leapyr


c  Read the command-line to get the input file name

      iarg = iclarg(1,infile)
      if( iarg.le.0 ) then
         write(*,'(a)') 'Missing arguments for dcbtab2 '
         write(*,'(a)') '  dcbtab2 dcb.dat [COD-file] '
         stop
      endif
      monfile = ' '
      iarg = iclarg(2,monfile)
      if( iarg.le.0 ) then
         write(*,'(a )') 'Translating a v1 file to v2'
      endif

*  Open the input and output files

      open(unit=luin,file=infile,status='old',iostat=ioerr)
      if( ioerr.ne.0 )  then
         call report_stat('FATAL','DCBTAB2','Main',infile
     .        ,'Cannot open input file',ioerr)
      else
          write(*,'(2a)') 'Opened input file ',infile
      endif
      if( monfile(1:1).ne.' ') then
        open(unit=lumon,file=monfile,status='old',iostat=ioerr)
        if( ioerr.ne.0 ) then
            call report_stat('FATAL','DCBTAB2','Main',monfile
     .              ,'Cannot open COD monthly dcb file',ioerr)
        else
           write(*,'(2a)') 'Opened monthly file ',monfile
        endif
      endif
      outfile = infile(1:nblen(infile))//'.new'
      print *,'outfile ',outfile
      open(unit=luout,file=outfile,status='unknown'
     .      ,iostat=ioerr)
      if( ioerr.ne.0 )  call report_stat('FATAL','DCBTAB','Main',outfile
     .       ,'Failure in opening the output dcb.dat file',ioerr)

*  Read all the values from the input dcb file into storage

      if( monfile(1:1).eq.' ' ) then
c-----translating a special-format CODE file (1991-2015)
        eof = .false.
        nl = 0
        do while( .not.eof)
          read(luin,'(a80)',iostat=ioerr) line
          if( ioerr.eq.-1.or.line(1:3).eq.'   '  ) then
            eof = .true.
          elseif( ioerr.ne.0 ) then
            call report_stat('FATAL','DCBTAB','Main',infile,
     .            'Error reading Version 1  DCB file ',ioerr)
          elseif( line(1:1).eq.'*' ) then
            continue
          elseif( line(14:17).eq.'CODE') then
            if(debug) print *,'Epoch line ',line
            read(line,'(2i4,37x,2i3)') stopyr,stopdoy,month,day
            if(debug) print *,'new span stopyr stopdoy ',stopyr,stopdoy
            dcbstop(1) = stopyr
            dcbstop(2) = stopdoy
            dcbstop(3) = 23
            dcbstop(4) = 59
            startyr = stopyr
            startdoy = stopdoy - day + 1
            dcbstart(1) = startyr
            dcbstart(2) = startdoy
            do i=3,5
              dcbstart(i) = 0
            enddo
            dcbstop(5) = 0
            endsvs = .false.
            do while( .not.endsvs )
              read(luin,'(a80)',iostat=ioerr) line
              if(debug)  print *,' line ',line
              if( ioerr.eq.-1 ) then
                eof =.true.
                endsvs =.true.
              elseif( line(14:17).eq.'CODE') then
                endsvs =. true.
                backspace(luin)
*** Temporarily skip the Glonass entries since svnav.dat doesn't yet include the entries
c**           elseif( line(1:1).eq.'G' ) then
              else
                nl = nl + 1
                call checkmax(nl,maxnl)
                read(line,'(a1,i2,20x,2f12.0)',iostat=ioerr)
     .            gnss(nl),prn(nl),dcb(nl),rms(nl)
                  if(debug)  print *,'READ nl prn dcb rms '
     .                , nl,prn(nl),dcb(nl),rms(nl)
                if(ioerr.ne.0) then
                  write(message,'(a,i6)')
     .                  'Error reading input file dcb line ',nl
                  call report_stat('FATAL','DCBTAB2','Main',infile
     .                            , message,ioerr)
                endif
c               get the SVN and check for a mid-month PRN change
                call get_times( gnss(nl),prn(nl),dcbstart,dcbstop
     .                        , svn(nl),start(1,nl),stop(1,nl) )
                if(debug) write(*,'(a,3i4,f8.3,4i5)')
     .            'nl iprn isvn dcb start stop'
     .            , nl,prn(nl),svn(nl),dcb(nl)
     .            ,(start(i,nl),i=1,2),(stop(i,nl),i=1,2)
              endif
            enddo
          endif
        enddo

      else
c-------incrementing an existing V2 file
c       read the input dcb.dat values into storage
        eof = .false.
        do while( .not.eof )
          read(luin,'(a80)',iostat=ioerr) line
          if( ioerr.eq.-1 ) then
             eof = .true.
          elseif( ioerr.ne.0 ) then
            call report_stat('FATAL','DCBTAB2','Main',infile
     .             , 'Error reading original DCB file ',ioerr)
          elseif( line(1:1).ne.' ' ) then
c           comment
            continue
          else
            nl = nl + 1
            call checkmax(nl,maxnl)
            read(line,'(1x,a1,i6,i4,2(i6,i4,2i3),2f10.3)')
     .         gnss(nl),svn(nl),prn(nl)
     .        ,(start(i,nl),i=1,4),(stop(i,nl),i=1,4),dcb(nl),rms(nl)
              start(5,nl) = 0
              stop(5,nl) = 0
          endif
        enddo
        print *,'after reading input file, nl = ',nl
c       read the new monthly values and add them to the arrays
c       read the header
        read(lumon,'(a80)',iostat=ioerr) line
cd        print *,'CODE header line ',line
        if( ioerr.ne.0 ) then
           call report_stat('FATAL','DCBTAB2','Main',monfile
     .           ,'Error reading CODE monthly DCB file ',ioerr)
        else
c         There are at least two versions of the header for GNSS files:
cCODE'S MONTHLY GNSS P1-C1 DCB SOLUTION, YEAR 2015, MONTH 07      03-AUG-15 08:47
cCODE'S 30-DAY GNSS P1-C1 DCB SOLUTION, ENDING DAY 019, 2016      20-JAN-16 07:33
cd          print *,'line 8:10 ',line(8:10)
          if(line(8:10).eq.'MON') then
            read(line,'(45x,i4,8x,i2)',iostat=ioerr) year,month
            if(debug) print *,'read year month ',year,month
c           set the start/stop times as the first and last days of the month
            dcbstart(1) = year
            dcbstart(2) = idoy(year,month,1)
            dcbstop(1) = year
            if(debug) print *,'dcbstart ',dcbstart
            if( month.eq.12 ) then
              dcbstop(2) = idoy(year,month,31)
            else
              month = month + 1
              dcbstop(2) = idoy( year,month,1) -1
            endif
            write(*,'(a,i5,i3)') 'CODE P1-C1 for ',year,month
          elseif(line(8:10).eq.'30-') then
            read(line,'(50x,i3,2x,i4)',iostat=ioerr) day,year
c           set the start time 30 days prior to the end day
            if(debug) print *,'day year ',day,year
            dcbstop(1) = year
                 dcbstop(2) = day
            dcbstart(1) = year
            dcbstart(2) = dcbstop(2) - 29
                 if( dcbstart(2).le.0 ) then
                   dcbstart(1) = year - 1
                   if(leapyr(year) ) then
                     dcbstart(2) = dcbstart(2) + 366
                   else
                     dcbstart(2) = dcbstart(2) + 365
                   endif
                 endif
            if(debug)  print *,'dcbstop ',dcbstop
            if(debug)  print *,'dcbstart',dcbstart
            write(*,'(a,i5,i4)') 'CODE P1-C1 for ',year,day
          else
            call report_stat('FATAL','DCBTAB2','Main',monfile,
     .         'Unrecognized header on CODE monthly RINEX.DCB file ',0)
          endif
          if( ioerr.ne.0 )
     .      call report_stat('FATAL','DCBTAB2','Main',monfile,
     .      'Error decoding date on CODE monthly RINEX.DCB file ',ioerr)
          dcbstop(3) = 23
          dcbstop(4) = 59
          do i=3,5
            dcbstart(i) = 0
          enddo
          dcbstop(5) = 0
c         skip the rest of the header
          do i=1,6
            read(lumon,'(a)')
          enddo
        endif
        if(debug) print *,'finished header '
c       now read the values
        eof = .false.
        do while(.not.eof )
          read(lumon,'(a80)',iostat=ioerr) line
          if(debug)  print *,'LINE ',line
          if( ioerr.eq.-1.or.line(1:3).eq.'   '.or.
     .        line(6:9).ne.'    ' ) then
            eof = .true.
            if(debug) print *,'EOF ',eof
          elseif( ioerr.ne.0 ) then
            call report_stat('FATAL','DCBTAB2','Main',infile,
     .            'Error reading Version 2  DCB values ',ioerr)
c** temporary to skip Glonass
* MOD MAF 20210224: Commented out to include GLONASS
c         elseif( line(1:1).eq.'R') then
c           continue
          else
            nl = nl + 1
            if(debug) print *,'CODE nl ',nl
            call checkmax(nl,maxnl)
            read(line,'(a1,i2,20x,2f12.0)',iostat=ioerr)
     .           gnss(nl),prn(nl),dcb(nl),rms(nl)
            if(debug) print *,'nl gnss prn dcb ',nl,gnss(nl),dcb(nl)
            if( ioerr.ne.0 ) then
              write(message,'(a,i6)') 'Error decoding dcb line ',nl
              call report_stat('FATAL','DCBTAB2','Main',monfile
     .                          , message,ioerr)
            else
c             get the SVN and check for a mid-month PRN change
              if(debug) print *
     .              ,'calling get_times nl gnss prn dcbtart dcbstop '
     .      ,nl,gnss(nl),prn(nl),(dcbstart(i),i=1,5),(dcbstop(i),i=1,5)
              call get_times( gnss(nl),prn(nl),dcbstart,dcbstop
     .                      , svn(nl),start(1,nl),stop(1,nl) )
            endif
          endif
        enddo
c     endif on whether translation or incrementing
      endif

* Populate sub-arrays for each system

      ng = 0
      nr = 0
      ne = 0
      nc = 0
      do i=1,nl
        if( gnss(i).eq.'G') then
           ng = ng+1
           call insert_arrays(
     .          gnss(i),prn(i),svn(i),dcb(i),rms(i),start(1,i),stop(1,i)
     .        , ng,prng,svng,dcbg,rmsg,startg,stopg )
        elseif( gnss(i).eq.'R') then
          nr = nr + 1
          call insert_arrays(
     .          gnss(i),prn(i),svn(i),dcb(i),rms(i),start(1,i),stop(1,i)
     .        , nr,prnr,svnr,dcbr,rmsr,startr,stopr )
        elseif( gnss(i).eq.'E') then
          ne = ne + 1
          call insert_arrays(
     .          gnss(i),prn(i),svn(i),dcb(i),rms(i),start(1,i),stop(1,i)
     .        , ne,prne,svne,dcbe,rmse,starte,stope )
        elseif( gnss(i).eq.'C') then
          nc = nc + 1
          call insert_arrays(
     .          gnss(i),prn(i),svn(i),dcb(i),rms(i),start(1,i),stop(1,i)
     .        ,  nc,prnc,svnc,dcbc,rmsc,startc,stopc )
        else
          call report_stat('FATAL','DCBTAB2','Main',infile
     .                  , 'Only gnss G R E C coded',0)
        endif
      enddo
      if(debug) then
        print *,'G values'
        do i=1,ng
          print *,svng(i),prng(i),dcbg(i)
     .      ,(startg(j,i),j=1,2),(stopg(j,i),j=1,2)
        enddo
      endif

*  Within each system, sort by SVN

      call sort_svn( ng,svng,prng,dcbg,rmsg,startg,stopg )
* MOD MAF 20210224: Uncommented to include non-GPS systems
      call sort_svn( nr,svnr,prnr,dcbr,rmsr,startr,stopr )
      call sort_svn( ne,svne,prne,dcbe,rmse,starte,stope )
      call sort_svn( nc,svnc,prnc,dcbc,rmsc,startc,stopc )

      if(debug) then
        print *,'Sorted G values'
        do i=1,ng
          print *,svng(i),prng(i),dcbg(i)
     .     ,(startg(j,i),j=1,2),(stopg(j,i),j=1,2)
        enddo
      endif

* Change the last two entries so that the stop time of the 2nd-to-last one,
* previously  2100 001 is set to the start of the last one, and the last
* is set to  2100 001.  Do not change if the last entry was not 2100 since
* this indicates a no-long-active SVN.

c** GPS
      svnlast = svng(1)
      do i=2,ng
        if( svng(i).eq.svnlast ) then
          continue
        else
          if(debug) then
            print *,'i ',i,'svng ',svng(i),' ne svnlast ',svnlast
            print *,'start i-1 ',startg(1,i-1),' stop i-1 ',stopg(1,i-1)
            print *,'start i-2 ',startg(1,i-2),' stop i-2 ',stopg(1,i-2)
          endif
          if( stopg(1,i-2).eq.2100 ) then
            stopg(1,i-2) = startg(1,i-1)
            stopg(2,i-2) = startg(2,i-1) - 1
            if(debug) then
              print *,'stop i-2 = 2100 '
              print *,'set stop i-2 = start i-1 ',stopg(1,i-2)
              print *,' start stop i ',startg(1,i),stopg(1,i)
            endif
            stopg(3,i-2) = 23
            stopg(4,i-2) = 59
            stopg(1,i-1) = 2100
            stopg(2,i-1) = 1
            stopg(3,i-1) = 0
            stopg(4,i-1) = 0
c**           svnlast = svng(i)
            if(debug) print *,'i-1 i-2 startgi1 stopgi1 i2  '
     .         ,i-1,i-2,startg(1,i-1),stopg(1,i-1),stopg(1,i-2)
          endif
c** moved rwk 190304
          svnlast = svng(i)
        endif
      enddo
      if( stopg(1,ng-1).eq.2100 ) then
        stopg(1,ng-1) = startg(1,ng)
        stopg(2,ng-1) = startg(2,ng) - 1
        stopg(3,ng-1) = 23
        stopg(4,ng-1) = 59
        stopg(1,ng) = 2100
        stopg(2,ng) = 1
        stopg(3,ng) = 0
        stopg(4,ng) = 0
        if(debug) print *,'2nd i-1 i-2 startgi1 stopgi1 i2  '
     .         ,i-1,i-2,startg(1,i-1),stopg(1,i-1),stopg(1,i-2)
      endif

* MOD MAF 20210224: Added GLONASS (copy of GPS block, above,
*                   with variable names changed to "r" from "g")
      svnlast = svnr(1)
      do i=2,nr
        if( svnr(i).eq.svnlast ) then
          continue
        else
          if(debug) then
            print *,'i ',i,'svnr ',svnr(i),' ne svnlast ',svnlast
            print *,'start i-1 ',startr(1,i-1),' stop i-1 ',stopr(1,i-1)
            print *,'start i-2 ',startr(1,i-2),' stop i-2 ',stopr(1,i-2)
          endif
          if( stopr(1,i-2).eq.2100 ) then
            stopr(1,i-2) = startr(1,i-1)
            stopr(2,i-2) = startr(2,i-1) - 1
            if(debug) then
              print *,'stop i-2 = 2100 '
              print *,'set stop i-2 = start i-1 ',stopr(1,i-2)
              print *,' start stop i ',startr(1,i),stopr(1,i)
            endif
            stopr(3,i-2) = 23
            stopr(4,i-2) = 59
            stopr(1,i-1) = 2100
            stopr(2,i-1) = 1
            stopr(3,i-1) = 0
            stopr(4,i-1) = 0
c**           svnlast = svnr(i)
            if(debug) print *,'i-1 i-2 startri1 stopri1 i2  '
     .         ,i-1,i-2,startr(1,i-1),stopr(1,i-1),stopr(1,i-2)
          endif
          svnlast = svnr(i)
        endif
      enddo
      if( stopr(1,nr-1).eq.2100 ) then
        stopr(1,nr-1) = startr(1,nr)
        stopr(2,nr-1) = startr(2,nr) - 1
        stopr(3,nr-1) = 23
        stopr(4,nr-1) = 59
        stopr(1,nr) = 2100
        stopr(2,nr) = 1
        stopr(3,nr) = 0
        stopr(4,nr) = 0
        if(debug) print *,'2nd i-1 i-2 startri1 stopri1 i2  '
     .         ,i-1,i-2,startr(1,i-1),stopr(1,i-1),stopr(1,i-2)
      endif

* MOD MAF 20210224: Added Galileo (copy of GPS block, above,
*                   with variable names changed to "r" from "e")
      svnlast = svne(1)
      do i=2,ne
        if( svne(i).eq.svnlast ) then
          continue
        else
          if(debug) then
            print *,'i ',i,'svne ',svne(i),' ne svnlast ',svnlast
            print *,'start i-1 ',starte(1,i-1),' stop i-1 ',stope(1,i-1)
            print *,'start i-2 ',starte(1,i-2),' stop i-2 ',stope(1,i-2)
          endif
          if( stope(1,i-2).eq.2100 ) then
            stope(1,i-2) = starte(1,i-1)
            stope(2,i-2) = starte(2,i-1) - 1
            if(debug) then
              print *,'stop i-2 = 2100 '
              print *,'set stop i-2 = start i-1 ',stope(1,i-2)
              print *,' start stop i ',starte(1,i),stope(1,i)
            endif
            stope(3,i-2) = 23
            stope(4,i-2) = 59
            stope(1,i-1) = 2100
            stope(2,i-1) = 1
            stope(3,i-1) = 0
            stope(4,i-1) = 0
c**           svnlast = svne(i)
            if(debug) print *,'i-1 i-2 startei1 stopei1 i2  '
     .         ,i-1,i-2,starte(1,i-1),stope(1,i-1),stope(1,i-2)
          endif
          svnlast = svne(i)
        endif
      enddo
      if( stope(1,ne-1).eq.2100 ) then
        stope(1,ne-1) = starte(1,ne)
        stope(2,ne-1) = starte(2,ne) - 1
        stope(3,ne-1) = 23
        stope(4,ne-1) = 59
        stope(1,ne) = 2100
        stope(2,ne) = 1
        stope(3,ne) = 0
        stope(4,ne) = 0
        if(debug) print *,'2nd i-1 i-2 startei1 stopei1 i2  '
     .         ,i-1,i-2,starte(1,i-1),stope(1,i-1),stope(1,i-2)
      endif

* MOD MAF 20210224: Added Beidou (copy of GPS block, above,
*                   with variable names changed to "r" from "c")
      svnlast = svnc(1)
      do i=2,nc
        if( svnc(i).eq.svnlast ) then
          continue
        else
          if(debug) then
            print *,'i ',i,'svnc ',svnc(i),' ne svnlast ',svnlast
            print *,'start i-1 ',startc(1,i-1),' stop i-1 ',stopc(1,i-1)
            print *,'start i-2 ',startc(1,i-2),' stop i-2 ',stopc(1,i-2)
          endif
          if( stopc(1,i-2).eq.2100 ) then
            stopc(1,i-2) = startc(1,i-1)
            stopc(2,i-2) = startc(2,i-1) - 1
            if(debug) then
              print *,'stop i-2 = 2100 '
              print *,'set stop i-2 = start i-1 ',stopc(1,i-2)
              print *,' start stop i ',startc(1,i),stopc(1,i)
            endif
            stopc(3,i-2) = 23
            stopc(4,i-2) = 59
            stopc(1,i-1) = 2100
            stopc(2,i-1) = 1
            stopc(3,i-1) = 0
            stopc(4,i-1) = 0
c**           svnlast = svnc(i)
            if(debug) print *,'i-1 i-2 startci1 stopci1 i2  '
     .         ,i-1,i-2,startc(1,i-1),stopc(1,i-1),stopc(1,i-2)
          endif
          svnlast = svnc(i)
        endif
      enddo
      if( stopc(1,nc-1).eq.2100 ) then
        stopc(1,nc-1) = startc(1,nc)
        stopc(2,nc-1) = startc(2,nc) - 1
        stopc(3,nc-1) = 23
        stopc(4,nc-1) = 59
        stopc(1,nc) = 2100
        stopc(2,nc) = 1
        stopc(3,nc) = 0
        stopc(4,nc) = 0
        if(debug) print *,'2nd i-1 i-2 startci1 stopci1 i2  '
     .         ,i-1,i-2,startc(1,i-1),stopc(1,i-1),stopc(1,i-2)
      endif


* Write out the merged file

c     get the date of the update and write the headers
      call getdat( irunt(1),irunt(2),irunt(3) )
      write(luout,'(a)')  '* dcb.dat Version 2.0 - units are ns '
      write(luout,'(a,i4,a,i2,a,i2 )')
     .   '* Last updated ',irunt(1),'-',irunt(2),'-',irunt(3)
      write(luout,'(a)')
     .  '*SYS SVN PRN  Start           Stop               P1-C1    rms'
c     write the values
      if( ng.ne.0 ) then
        do i=1,ng
          write(luout,'(1x,a1,i6,i4,2(i6,i4,2i3),2f10.3)')
     .           'G',svng(i),prng(i)
     .           ,(startg(j,i),j=1,4),(stopg(j,i),j=1,4),dcbg(i),rmsg(i)
        enddo
      endif
** RWK temporary: don't write out the GLonass values since svnav.dat does not yet have the SVN
* MOD MAF: Commented out to include GLONASS
c     nr = 0
**
      if( nr.ne.0 ) then
        do i=1,nr
          write(luout,'(1x,a1,i6,i4,2(i6,i4,2i3),2f10.3)')
     .           'R',svnr(i),prnr(i)
     .           ,(startr(j,i),j=1,4),(stopr(j,i),j=1,4),dcbr(i),rmsr(i)
        enddo
      endif
      if( ne.ne.0 ) then
        do i=1,nc
          write(luout,'(1x,a1,i6,i4,2(i6,i4,2i3),2f10.3)')
     .           'E',svne(i),prne(i)
     .           ,(starte(j,i),j=1,4),(stope(j,i),j=1,4),dcbe(i),rmse(i)
        enddo
      endif
      if( nc.ne.0 ) then
        do i=1,nc
          write(luout,'(1x,a1,i6,i4,2(i6,i4,2i3),2f10.3)')
     .           'C',svnc(i),prnc(i)
     .           ,(startc(j,i),j=1,4),(stopc(j,i),j=1,4),dcbc(i),rmsc(i)
        enddo
      endif

      stop
      end

c------------------------------------------------------------------------------

      Subroutine get_times( ag,iprn,dcbstart,dcbstop
     .                    , isvn,newstart,newstop )

*     Check the iprn for a change midmonth and reset the start and
*     stop times to match the correct SVN

      integer*4 maxnl
      parameter(maxnl= 30000)

      integer*4 iprn,dcbstart(5),dcbstop(5),isvn,svnstart(5),svnstop(5)
     .        , newstart(5),newstop(5),startdoy,idum
      real*8 rdum
      character*1 ag
      character*256 message

cd    print *,'ag iprn dcbstart dcbstop ',ag,iprn,dcbstart,dcbstop
      if( iprn.le.32 ) then
* MOD TAH 190702: Added place holder for antpwr to snav_read call
        call svnav_read( -1,dcbstart(1),dcbstart(2),0,0,ag,iprn,isvn
     .                 , idum,rdum,rdum,rdum,rdum, rdum
     .                 , svnstart,svnstop )
c       if svn is returned as zero, PRN not active at start, search through
c       the month until it's valid:
        if( isvn.eq.0 ) then
          startdoy = dcbstart(2)
          do while( isvn.eq.0 )
            startdoy = startdoy + 1
* MOD TAH 190702: Added place holder for antpwr to snav_read call
            call svnav_read( -1,dcbstart(1),startdoy,0,0,ag,iprn,isvn
     .                 , idum,rdum,rdum,rdum,rdum, rdum
     .                 , svnstart,svnstop )
cd            print *,'DEBUG iprn isvn,yr doy '
cd     .            ,iprn,isvn,dcbstart(1),startdoy
            if( (startdoy-dcbstart(2)).gt.32 ) then
cd              print *,'DEBUG iprn isvn dcbstart startdoy '
cd     .           , iprn,isvn,dcbstart(2),startdoy
              write(message,'(a,i3,a,4i5)') 'No valid SVN entry for PRN'
     .             , iprn,' in span '
     .             , dcbstart(1),dcbstart(2),dcbstop(1),dcbstop(2)
              call report_stat('FATAL','DCBTAB2','get_times',' '
     .                  ,message,0)
            endif
          enddo
          do i=1,5
            newstart(i) = svnstart(i)
            newstop(i) =  dcbstop(i)
          enddo
        else
          do i=1,5
            newstart(i) = dcbstart(i)
            if( svnstop(1).ne.0 ) then
              if( itimdif(svnstop,dcbstop).lt.0 ) then
                 newstop(i) = svnstop(i)
              else
                newstop(i) = dcbstop(i)
              endif
            else
              newstop(i) = dcbstop(i)
            endif
          enddo
        endif
cd  DEBUG
cd       write(*,'(a,2i4,4i5)') 'iprn isvn svnstart svnstop '
cd     .    ,iprn,isvn,(svnstart(i),i=1,2),(svnstop(i),i=1,2)
      else
        print *,'** NEW iprn ',iprn
        iprn = iprn - 50
* MOD TAH 190702: Added place holder for antpwr to snav_read call
        call svnav_read( -1,dcbstop(1),dcbstop(2),0,0,ag,iprn,isvn
     .                 , idum,rdum,rdum,rdum,rdum, rdum
     .                 , svnstart,svnstop )
        do i=1,5
          newstart(i) = svnstart(i)
          newstop(i) = dcbstop(i)
        enddo
cd        write(*,'(a,2i4,4i5)') 'iprn isvn svnstart svnstop '
cd     .    ,iprn,isvn,(svnstart(i),i=1,2),(svnstop(i),i=1,2)
      endif

      return
      end


c---------------------------------------------------------------------

      Subroutine checkmax(nl,maxnl)

      integer*4 nl,maxnl
      character*80 message

      if( nl.gt.maxnl ) then
         write(message,'(a,i7)') 'Number of dcb entries > maxnl ',maxnl
         call report_stat('FATAL','DCBTAB2','checkmax',' ',message,0)
      endif
      end

c-----------------------------------------------------------------------

      Subroutine insert_arrays( cgnss,iprn,isvn,rdcb,rrms,istart,istop
     .                        , ng,prng,svng,dcbg,rmsg,startg,stopg )

*     Insert one set of values from the master arrays into the GNSS-
*     specific arrays.

      integer*4 maxnl
      parameter(maxnl= 30000)

*     Counter for the GNSS-specific array
      integer*4 ng

*     Values to be inserted
      integer*4 iprn,isvn,istart(5),istop(5)
      real*4 rdcb,rrms
      character*1 cgnss

*     GNSS-specific arrays
      integer*4 prng(maxnl),svng(maxnl),startg(5,maxnl),stopg(5,maxnl)
      real*4 dcbg(maxnl),rmsg(maxnl)

*     Local
      integer*4 i

cd      print *,'INSERT ng iprn isvn rdcb rrms istart istop '
cd     .        ,ng,iprn,isvn,rdcb,rrms,istart,istop
      prng(ng) = iprn
      svng(ng) = isvn
      dcbg(ng) = rdcb
      rmsg(ng) = rrms
      do i=1,5
        startg(i,ng) = istart(i)
        stopg(i,ng) = istop(i)
      enddo

      return
      end

c-----------------------------------------------------------------------
      Subroutine sort_svn( n,svn,prn,dcb,rms,start,stop )

*     Sort the DCB arrays by SVN

      implicit none

      integer*4 maxnl
      parameter(maxnl= 30000)


      integer*4 n,i,j,k
     .        , svn(maxnl),prn(maxnl),start(5,maxnl),stop(5,maxnl)
     .        , svnbuf1,svnbuf2,prnbuf1,prnbuf2
     .        , startbuf1(5),stopbuf1(5),startbuf2(5),stopbuf2(5)

      real*4 dcb(maxnl),rms(maxnl),dcbbuf1,dcbbuf2,rmsbuf1,rmsbuf2

c     Sort by SVN

      do i = 1,n-1
        do j = 1,n-i
           svnbuf1 = svn(j)
           svnbuf2 = svn(j+1)
           prnbuf1 = prn(j)
           prnbuf2 = prn(j+1)
           dcbbuf1 = dcb(j)
           dcbbuf2 = dcb(j+1)
           rmsbuf1 = rms(j)
           rmsbuf2=  rms(j+1)
           do k=1,5
             startbuf1(k) = start(k,j)
             startbuf2(k) = start(k,j+1)
             stopbuf1(k) = stop(k,j)
             stopbuf2(k) = stop(k,j+1)
           enddo
           if( svnbuf1.le.svnbuf2 ) then
             svn(j) = svnbuf1
             svn(j+1) = svnbuf2
             prn(j) = prnbuf1
             prn(j+1) = prnbuf2
             dcb(j) = dcbbuf1
             dcb(j+1) = dcbbuf2
             rms(j) = rmsbuf1
             rms(j+1) = rmsbuf2
             do k=1,5
               start(k,j)   = startbuf1(k)
               start(k,j+1) = startbuf2(k)
               stop(k,j)   = stopbuf1(k)
               stop(k,j+1) = stopbuf2(k)
             enddo
           else
             svn(j) = svnbuf2
             svn(j+1) = svnbuf1
             prn(j) = prnbuf2
             prn(j+1) = prnbuf1
             dcb(j) = dcbbuf2
             dcb(j+1) = dcbbuf1
             rms(j) = rmsbuf2
             rms(j+1) = rmsbuf1
             do k=1,5
               start(k,j)   = startbuf2(k)
               start(k,j+1) = startbuf1(k)
               stop(k,j)   = stopbuf2(k)
               stop(k,j+1) = stopbuf1(k)
             enddo
           endif
         enddo
      enddo

      return
      end

c------------------------------------------------------------------

      LOGICAL FUNCTION LEAPYR ( IYEAR )

c     Determine if the input year is a leap year
c     R. King   2 January 1992
c     based on S.Shimada's integer function LEAP

      integer*4 iyear,iyear1

      leapyr = .false.
      iyear1 = iyear
      if( iyear1.lt.1900 ) iyear1 = iyear1 + 1900
      IF( MOD(IYEAR1,  4) .NE. 0 )  GO TO 9
      IF( MOD(IYEAR1,100) .NE. 0 )  GO TO 1
      IF( MOD(IYEAR1,400) .NE. 0 )  GO TO 9
    1 leapyr = .true.
    9 RETURN
      END



