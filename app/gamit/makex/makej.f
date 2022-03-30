Copyright (c) Massachusetts Institute of Technology,1989. All rights reserved.
      Program MAKEJ
c
c     Make J files from an SP3 file, a RINEX navigation, or C-files

c     R.W. King December 1989 from code in MAKEK
c     Last modified by R. King July 2015

c     Input files:
c
c       Navigation file and optionally an SP3 
c
c     Output file:
c
c       J-file   SV CLOCK file, assumed created.
         
      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'

c                        DECLARATIONS
c                        ************
                       
c     In makex.h: 
c         File names  : character*80 fnav,fsvclk,fsp3
c         Unit numbers: integer unav,usvclk,usp3
c         Dimensions  : maxepc,maxseg 

      logical          lask,batch,fcheck

      character*1      gnss 
      character*3      type
      character*16     uname          
      character*80     clkfile, wildcard,pickfn
      character*109    version
      character*154    header
      character*60     afmt
      character*256    message

      integer*4        irunt(6),i,ihnsec
     .               , iprns(maxsat),nprns
     .               , icalc,ioerr
     .               , orbflen,iarg,iclarg,len


c     return non-blank length of string
      integer nblen


c     here is the format statement for the J-file 
c        (year increased to I4 by rwk 990803)
c        (gnss added by rwk 151120)
      data afmt
     ./'(i4,1x,i4,2i3,1x,f10.7,2x,i4,1x,f14.7,1x,a1,i2.2,2x,3d19.11)'/


c     satellite to debug
      integer*4 iprndb/1/
      logical  debug/.false./
       

c     RUN INITIALIZATION

c     the format for user questions:
 1000 format (/,1x,a,1x,$)

C     unit numbers
      uscren = 6 
      unav = 10
      usvclk = 11
      usp3   = 12 

c     tell the user who we are
      version = ' '
      call mversn(version)
      call getdat(irunt(1),irunt(2),irunt(3))
      call gettim(irunt(4),irunt(5),irunt(6),ihnsec )
      call getusr(uname)

c     write makej status line
      WRITE(message,5)version
    5 FORMAT('Started MAKEJ ',a109)
      call report_stat('STATUS','MAKEJ','makej',' ',message,0)

c Read the command-line input

c   Keep the original numerical flags even through the c-file and nav-file options are rarely used
c      icalc = 1   use the broadcast ephemieris
c      icalc = 2   use c-files (SA on)  (works only with interactive input)
c      icalc = 3   use the sp3 file with broadcast as fall-back for missing SVs (invoke3 if 3 arguments)

c     If there are command-line arguments, use them and skip the interactive questions   
      iarg = iclarg(1,fnav)
      if( iarg.gt.0 ) then
        batch = .true.
c       Four arguments: nav-file j-file  sp3-file [optional] gnss [optional] 
        iarg = iclarg(2,fsvclk)  
        if( iarg.le.0 ) call report_stat('FATAL','MAKEJ','makej',' ',
     .     'Missing command-line argument for J-file',0) 
        iarg = iclarg(3,fsp3)
        if( .not.fcheck(fsp3)) then
cc        if( iarg.le.0.or.fsp3(1:1).eq.' ' ) then
           call report_stat('WARNING','MAKEJ','makej',' ',
cc    .      'No command-line argument for SP3 file, use the nav-file',0)
     .    'SP3 file not specified or not available, use the nav-file',0)
           icalc = 1                  
        else
          icalc = 3 
        endif   
        iarg = iclarg(4,gnss)
cd        print *,'read gnss ',gnss
        if( iarg.le.0 ) gnss = 'G'
      else   
c       no command-line arguments, enter info interactively             
        write(6,*) 
        write(6,*) 
     .     'To run with command-line arguments, abort and enter: '
        write(6,*) 
     .    ' makej [nav-file] [j-file] [sp3-file] [gnss] '
     .      ,' (last two optional)'
        write(6,*) ' '
        write(6,*) 
     .    ' Command-line input works only for use of broadcast message'
     .   ,' You must run interactively to use C-files '
     .   ,' An existing J-file will be overwritten'
        write(6,*) ' '
        batch = .false.
        write (*,*)
     .' Choose source of SV oscillator frequency corrections:'
        write (*,*) '  1  ',
     .'E-file broadcast message. [OK for simultaneous sampling]'
        write (*,*) '  2  ',
     .'Second order fit to C-file from site with H-maser [see Chap 4.6]'
        call imenu (icalc,2)
      endif
        
C Open the input files

       if(icalc.ne.2 ) then
c        open the navigation file 
         if( .not.batch) then
           fnav = ' ' 
           wildcard = '[ae]*.???   ' 
cd           print *,'wildcard ',wildcard
           len = 12
           write (6,1000) 'Choose an ephemeris file from'
           fnav(1:len) = pickfn (wildcard,len) 
c          print *,'len,fnav ',len,fnav(1:len)
c          also check for RINEX N-files   - this doesn't work; pickfn hangs up until it is satisfied
           if (fnav(1:1) .eq. ' ') then    
             wildcard = '????????.??n'  
c            print *,'wildcard ',wildcard
             len = 12
             fnav(1:len) = pickfn ( wildcard,len )
c            print *,'len,fnav ',fnav(1:len)
           endif  
         endif
         open(unit=unav,file=fnav,status='old',iostat=ioerr)
         if (ioerr .eq. 0) then
           call report_stat('STATUS','MAKEJ','makej',fnav,
     .      'Opened navigation file: ',ioerr)
         else
           call report_stat('FATAL','MAKEJ','makej',fnav,
     .       'Error opening navigation file: ',ioerr)
         endif 
       endif
       if(icalc.eq.3 ) then
c        open the SP3 file
         open(unit=usp3,file=fsp3,status='old',iostat=ioerr)
         if (ioerr .eq. 0) then
           call report_stat('STATUS','MAKEJ','makej',fsp3,
     .      'Opened SP3 file: ',ioerr)
         else
           call report_stat('FATAL','MAKEJ','makej',fsp3,
     .       'Error opening SP3 file: ',ioerr)
         endif 
       endif

c Open the output J-file (overwriting any existing file)

      if( .not.batch ) then 
        write(*,1000) 'Enter output J-file name >: '
        read(5,'(a)') fsvclk
      endif
      open(unit=usvclk,file=fsvclk,status='unknown',
     .     form='formatted',iostat=ioerr)
      if (ioerr .eq. 0) then
c        file is open, write header info
         call report_stat('STATUS','MAKEJ','makej',fsvclk
     .                   , 'Opened J-file: ',ioerr)
c        Write the header records of the J-file
         header = ' '
         if( icalc.eq.1 ) clkfile = fnav
         if( icalc.eq.2 ) clkfile = 'C-file'
         if( icalc.eq.3)  clkfile = fsp3 
         header ='SV clock terms from '//clkfile(1:nblen(clkfile))//
     .            ' '//uname(1:7)//' MAKEJ '//version
         write(usvclk,70) header,afmt(1:nblen(afmt))
   70    format(a154,/,'YEAR  DOY HR MN SEC(UTC)    WKNO  SOW(GPST)   '
     .          ,'  PRN        XEAF0            XEAF1             XEAF2'
     .          ,'              rms(ns)  #val source'
     .          ,/, a)
      else
         call report_stat('FATAL','MAKEJ','makej',clkfile
     .                   , 'Error opening J-file: ',ioerr)
      endif
                                
c Get the clock terms
         
      if (icalc .eq. 1) then
c        use the navigation file
         call j_from_nav(batch,afmt,gnss,nprns,iprns,debug,iprndb)
      else if (icalc .eq. 2) then
c        use second-order fit to C-file
c        this option relevant only to GPS
         if(gnss.ne.'G') call report_stat('FATAL','MAKEJ','makej'
     .       ,' ','Cannot use C-file clocks for GNSS other than GPS',0)
         call j_from_c (afmt,nprns,iprns,debug,iprndb)
      elseif (icalc.eq.3) then
        call j_from_sp3(afmt,gnss,nprns,iprns,debug,iprndb )
      else
         call report_stat('FATAL','MAKEJ','makej',' '
     .                   , 'Unknown icalc value ',0) 
      endif

      call report_stat('STATUS','MAKEJ','makej',' ',
     .'Normal end in MAKEJ ',0)

      stop
      end



