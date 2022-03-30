Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1995.   All rights reserved.

      PROGRAM  FIXDRV
C     Make GAMIT batch files and I-file (station clock file)
C
C     S.Shimada              original at NRCDP
C     S.Shimada   90/02/07   modified at IGPP
C     Y. Bock     May 1991   modified for multi-session mode
C     Y. Bock & R.W. King    1992 for kinematic version of GAMIT
C
C     Files
c       screen    : Terminal                  LUN= LSCRN (6)  output
C       scratch file                          LUN=       (8)  output   
c       fixdrv.out: Clock log file            LUN=LPRNT  (9)  ouput
C       DFILE     : D-file                    LUN=       (10  input (possible output)
c       DFILE1    : Backup D-file             LUN=       (11) input/output
C       TFILE     : T-file                    LUN=LTFIL  (12) input
C       SITTBL    : sittbl.                   LUN=LSITE  (13) input
C       SESTBL    : sestbl.                   LUN=LSESS  (14) input  
c       SESSFO    : session.info              LUN=LSESFO (15) input
C       IFILE     : I-file                    LUN=LIFIL  (16) output
C       BFILE     : main batch file           LUN=       (17) output  
c       SFILE     : simulation control file   LUN=SFILE  (18) input 
c       STNFO     : station.info              LUN=LSTNFO (19) input
C       KFILE     : K-file                    LUN=       (20) input (assigned in CLKPRM/SESCHK)
C       BFIL2     : secondary batch file      LUN=     (20,21)output
C       XFILE     : beginning X- (or C-) file LUN=LXFIL  (22) input
C       GFILE     : G-file                    LUN=LGFIL  (23) input
C       CFILE     : C-file                    LUN=LCFIL  (24) input
C       LFILE     : L-file or apr file        not read        input
c       JFILE     : J-file                    not read        input
c       IANT      : antmod.dat                LUN=IANT   (32) input 
c       autcln.cmd: autcln commmand file      LUN=IAUT   (33) input (assinged in SOMAKE)   
c       hi.dat    : antenna mechanical offsts LUN=IHI    (34) input 
                                              

C     subroutines and functions
C        FIXDRV
c          + RDSEST
C          +-BMAKE  +-ARMAKE +-NBLEN .
C          |        |        +-RDSEST --UPPERS .
C          |        +-ACMAKE +-NBLEN .
C          |        |        +-UPNAM1 --NEWCHR .
C          |        +-CHRNUM .
C          |        +-CMMAKE .
C          |        +-EILOAD +-NBLEN .
C          |        |        +-UPNAM1 --NEWCHR .
C          |        +-GETDAT
C          |        +-GETTIM
C          |        +-GETUSR
C          |        +-LOWERS .
C          |        +-MDMAKE +-LFNTOP .
C          |        |        +-LOWERS .
C          |        |        +-NBLEN .
C          |        |        +-RDSITT --UPPERS .
C          |        |        +-UPNAM1 --NEWCHR .
C          |        +-NBLEN
C          |        +-RDSEST --UPPERS .
C          |        +-SATNMR .
C          |        +-XCSATS +-CHDRED +-READC1 +-FERROR .
C          |        |        |        +-READC2 +-FERROR .
C          |        |        |        +-READC3 +-FERROR .
C          |        |        +-COPENS --FERROR .
C          |        |        +-GETMD  .
C          |        |        +-LFNTOP .
C          |        |        +-LOWERS .
C          |        |        +-NBLEN .
C          |        |        +-NEXTT1 --LEAPYR .
C          |        +-THDRED +-JULDAY
C          |        |        +-LOWERS .
C          |        |        +-RJUST  .
C          |        |        +-TOPENS --FERROR .
C          |        |        +-UPPERC .
C          |        +-SOMAKE +-LFNTOP .
C          |        |        +-NBLEN .
C          |        |        +-NXEPOC +-CHDRED +-READC1 +-FERROR .
C          |        |        |        |        +-READC2 +-FERROR .
C          |        |        |        |        +-READC3 +-FERROR .
C          |        |        |        +-COPENS --FERROR .
C          |        |        |        +-LFNTOP .
C          |        |        |        +-LOWERS .
C          |        |        |        +-NBLEN .
C          |        |        +-NXSAT  +-CHDRED +-READC1 +-FERROR .
C          |        |        |        |        +-READC2 +-FERROR .
C          |        |        |        |        +-READC3 +-FERROR .
C          |        |        |        +-COPENS --FERROR .
C          |        |        |        +-LFNTOP .
C          |        |        |        +-LOWERS .
C          |        |        |        +-NBLEN .
C          |        |        +-RDSEST --UPPERS.
C          |        |        +-RDSITT --UPPERS.
C          |        |        +-SOLVRN +-GETDAT
C          |        |                 +-GETTIM
C          |        |                 +-GETUSR
C          |        |                 +-LJUST  .
C          |        |                 +-SOLATM .
C          |        |                 +-SOLCLK .
C          |        |                 +-SOLCRD .
C          |        |                 +-SOLORB .
C          |        |                 +-UPNAM1 --NEWCHR .
C          |        +-UPNAM1 --NEWCHR .
C          |        +-UPNAM2 --UPCHR  .
C          +-GETDAT
C          +-GETTIM
C          +-GETUSR
C          +-IMAKEF +-CHDRED +-READC1 +-FERROR .
C          |        |        +-READC2 +-FERROR .
C          |        |        +-READC3 +-FERROR .
C          |        +-CLKPRM +-DAYNUM .
C          |        |        +-JULDAY
C          |        |        +-MONDAY .
C          |        +-COPENS --FERROR .
C          |        +-LFNTOP .
C          |        +-LOWERS .
C          |        +-NBLEN .
C          |        +-UPPERS .
C          |        +-XHDRED +-MONDAY .
C          |                 +-UPPERC .
C          +-LFNTOP .
C          +-LOWERS .
C          +-NBLEN  .
C          +-SORT_STRING.
C          +-RDSITT --UPPERS.
C          +-SESCHK +-RSESFO

      implicit none

      include '../includes/dimpar.h'     
      include '../includes/global.h'

c** FOR DEBUG
      logical ex,iop
      integer num

      LOGICAL  IFLSW, fexist, fcheck, reqd, simulation, old_stinf
     .            , warnings

      integer  *4   ldfil, ltfil, lsite, lsess, lstnfo, lsesfo, lifil
     .            , lxfil, lgfil, lcfil, lsfil, lscrn, iarg, iant,lprint
     .            , nblen, year, doy,  istat1, jstat 
     .            , klock(maxsit),klocks
     .            , iy,imonth, ihr, imn, ihnsec, isec, isessn
     .            , nepoch, inter
     .            , ndummy, iday, ispan
     .            , ians,  nout, nk, len
     .            , icall, istarts(5), istops(5), nsod, ill, lcmd, ic
     .            , nstat, iclarg, ioerr, ierr_sestbl, ierr_sittbl
     .            , count_arg,nchr,trimlen,ihi,i
                                         
      real     * 4  swver 
      real     * 8  anth,antn,ante,offarp(3)
      real*8 antdaz  ! Antenna aligment from True N (deg).
 
      CHARACTER* 1  XORC, KPICK, CLKFT, use_ifile
      CHARACTER* 3  proc
      CHARACTER* 4  SNAMES(maxsit), SNAME,  sitei,  sitcod
      character* 5  htcod, fixdrv_vers, radome
      character* 6  rcvcod,antcod
      CHARACTER*16  DFILE, TFILE, GFILE, IFILE, IFILEB, MFILE
      CHARACTER*16  LFILE, KFILE, SITTBL, SESTBL, sfile
      CHARACTER*16  XFILE(maxsit)
     .           ,  xjunk(maxsit), JFILE,kfiles(100)
      CHARACTER*16  UNAME, stanam, wildcard, pickfn  
      character*20  rcvers,rcvrsn,antsn
      character*40  vers
      character*48  formti
      character*80  LINE,BUFF80
      character*120 wcmd     
      character*80 scratch_dir
      character*256 message
   
      parameter ( lprint=9, ldfil=10, ltfil=12, lsite=13, lsess=14
     .          , lsesfo=15, lifil=16, lsfil=18, lstnfo=19, lxfil=22
     .          , lgfil=23, lcfil=24, lscrn=6, iant=32, ihi=33)

      data  formti /'(a4,1x,i4,1x,i3,1x,i1,1x,2i3,1x,f7.4,3x,4d16.8) '/

c Remove old versions of the status, warning, and error files
c*  rwk/scm/tah 060628: No longer clear since now written to GAMIT.status/warning/fatal
c      call report_stat('CLEAR','FIXDRV',' ',' ', ' ',0)
c      call report_stat('CLEAR','LIB',' ',' ', ' ',0)
                        
c Start the status file and write the Header and version number to the screen

      call fversn( vers )
      fixdrv_vers = vers(1:5)
              
c Get the user, date, and time

      call getusr ( uname )
      call getdat ( iy, imonth,iday )
      call gettim ( ihr,imn,isec,ihnsec )
                
c Check for a prvious failure (in MAKEXP or MAKEX)

      if( fcheck('GAMIT.fatal') )
     .  call report_stat('FATAL','FIXDRV','fixdrv',' '
     .                  ,'GAMIT.fatal exists: FIXDRV not executed',0)

c=========================================================================
 

c  Open the session and site tables

         sestbl = 'sestbl.'
         if (fcheck(sestbl)) then
            open( lsess, file=sestbl, status='old' )
         else
            call report_stat('FATAL','FIXDRV','fixdrv',' '
     .                      ,'sestbl. not available',0)
         endif
c        initialize the counter for sestbl errors
         ierr_sestbl = 0

         sittbl = 'sittbl.'
         if (fcheck(sittbl)) then
            open( lsite, file=sittbl, status='old' )
         else
            call report_stat('FATAL','FIXDRV','fixdrv',' '
     .                      ,'sittbl. not available',0)
         endif                                   
c        initialize the counter for sittbl errors
         ierr_sittbl = 0

         if (fcheck('session.info')) then
            open( lsesfo,file='session.info', status='old',iostat=ioerr)
            if( ioerr.ne.0 ) 
     .      call report_stat('FATAL','FIXDRV','fixdrv',' '
     .                            , 'Error opening session.info',ioerr)
         endif

         if (fcheck('station.info')) then
            OPEN( lstnfo,file='station.info', STATUS='OLD' )    
 
         else
           call report_stat('WARNING','FIXDRV','fixdrv',' '
     .                        , 'No station.info ',0)
           call report_stat('WARNING','FIXDRV','fixdrv',' '
     .         , '--MODEL will use X-file values for antenna offsets',0)
         endif    

c          Open the height-of-instrument table for converting field measurements
         if( fcheck('hi.dat')) then
           ioerr = 0
           open(unit=ihi,file='hi.dat',status='old',iostat=ioerr)
           if( ioerr.ne.0 ) then
              call report_stat('FATAL','FIXDRV','fixdrv','hi.dat '
     .                  ,'Error opening antenna dimension table',ioerr)
           endif 
         else 
             call report_stat('WARNING','FIXDRV','fixdrv',' '
     .                      ,'No hi.dat file--MODEL will use hisub',0)
         endif

c          Open the antenna phase centre model file
         if (fcheck('antmod.dat') .or. fcheck('station.info') ) then
           ioerr = 0
           open(unit=iant,file='antmod.dat',status='old',iostat=ioerr)
           if (ioerr .ne. 0) then
              call report_stat('FATAL','FIXDRV','fixdrv','antmod.dat '
     .                      ,'Error opening antenna model table',ioerr)
           endif
         else
              call report_stat('WARNING','FIXDRV','fixdrv',' '
     .                      ,'No antmod.dat file--required by MODEL',0)
         endif
c
c     check for Processing Agency (this is mandatory)
      proc = '   '
      reqd = .true.
      call rdsest( 17, 'Processing Agency', 3, proc, lsess
     .            , reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill   

c     check for override of processing (day) directory for scratch (orbit, c-files, normal eqs)
      scratch_dir = ' '
      reqd = .false.
      call rdsest( 17, 'Scratch directory', 80, scratch_dir, lsess
     .            , reqd, ill ) 
      scratch_dir = scratch_dir(1:trimlen(scratch_dir))    
      call lowers(scratch_dir) 
      if( ill.ne.0 ) ierr_sestbl = ill   


c======================================================================================

c  Get D-file name from the command-line argument if present
                      
      iarg = iclarg(1,dfile)
      if  (iarg .le. 0) then
c     If no D-file on the command line, use pickfn to prompt the user
c        machin = getmac(1)
c        if (index(machin,'apollo') .gt. 0) then 
           wildcard = 'd*.*'
           len = 4
           dfile = pickfn(wildcard,len) 
           dfile = dfile(1:len)
      else
c       stop if the d-file given on the command line does not exist
        if( .not.fcheck(dfile) ) 
     .   call report_stat('FATAL','FIXDRV','fixdrv',dfile
     .                   ,'D-file not found:',0)
      endif  

c   Make sure the D-file name has the proper format--else problems later

      if( dfile(7:7).ne.'.' )
     .   call report_stat('FATAL','FIXDRV','fixdrv',' '
     .                   ,'D-file name must have . as 7th character',0)

c   Name of M-file and Q-file is same as D-file except for 1st and 6th characters

c     6th character is '1' for quick solution, 'a' for regular
      MFILE = DFILE
      MFILE(1:1) = 'm'
      MFILE(6:6) = '1'

c   See if simulation

      reqd = .false.
      simulation = .false.                      
      sfile = ' '
      call rdsest(14,'Simulation con',16,sfile,lsess,reqd,ill)
      if( ill.ne.0 ) ierr_sestbl = ill
      if( sfile(1:1).ne.' ') then  
        call lowers(sfile)
        simulation = .true. 
        xorc = 's'   
      endif 
  
c   Open the D-file and check for the presence of all the X-files 
c   --if one or more missing, remove them from the D-file, issuing
c     a warning and saving the old file with extent .saved_original
      
      if (simulation) then 
        open( ldfil, file=dfile, iostat=ioerr,status='old')
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','FIXDRV','fixdrv',dfile
     .                     ,'Error opening D-file',ioerr)

      else   
        call dcheck(dfile,ldfil)
      endif

c==========================================================================

c   --- Read the D-file-----------------------------------------------

c   Prior to release 10.62 the first two entries are the now-defunct 
c   solution and session number.   Replace the first with the GNSS system
c   and temporarily keep the second one = 1 
c   
c   Read the now-dummy solution and session number    
      read( ldfil,'(a)',iostat=ioerr) gnss
      if( gnss.eq.'1') then 
         call report_stat('WARNING','FIXDRV','fixdrv','  '
     .         ,'Old-style d-file, set gnss = G',0)
         gnss = 'G'
      elseif ( gnss.ne.'G'.and.gnss.ne.'R'.and.gnss.ne.'C'.and.
     .    gnss.ne.'E'.and.gnss.ne.'I'.and.gnss.ne.'J') then
         call report_stat('FATAL','FIXDRV','fixdrv','  '
     .         ,'First entry in the d-file must be a valid GNSS code',0)
      endif
      read( ldfil, * )  ndummy
      if(ndummy.ne.1)
     .  call report_stat('FATAL','FIXDRV','fixdrv','  '
     .             ,'FIXDRV no longer supports multiple sessions',0)

C   Read the L- or apr file name

      read( ldfil, '(a)'  )  lfile
      call  lowers( lfile )       
c     make sure the file naming conventions have been followed: lxxxxy.ddd or *.apr (<=16 characters)
      nchr = nblen(lfile)
      if( nchr.le.16 .and.
     .   (lfile(1:1).eq.'l' .or. lfile(nchr-2:nchr).eq.'apr') ) then
         if( .not.fcheck(lfile))
     .      call report_stat('FATAL','FIXDRV','fixdrv',lfile
     .                 ,'Coordinate file not available',0)
      else
        call report_stat('FATAL','FIXDRV','fixdrv',lfile,
     .   'First entry in D-file should be L- or apr file; check name',0) 
      endif

C   Read the T- (G-) file name

      READ( ldfil, '(A)'  )  TFILE
      CALL  LOWERS( TFILE )
      IF( TFILE(1:1) .EQ. 't' )  THEN
         GFILE = TFILE
         GFILE(1:1) = 'g'
      ENDIF
      IF( TFILE(1:1) .EQ. 'g' )  THEN
         GFILE = TFILE
         TFILE(1:1) = 't'
      ENDIF
      if( .not.fcheck(tfile) .and. .not.fcheck(gfile))
     .    call report_stat('FATAL','FIXDRV','fixdrv',tfile
     .                    ,'Neither T- nor G-file available',0)

C   Read the I-file name and determine if it needs to be created
                
      iflsw = .false.
      read( ldfil, '(a)'  )  ifile
      if( ifile(1:4).eq.'none' ) ifile(1:4)='NONE' 
      reqd = .false.
      call rdsest(10,'Use I-file',1,use_ifile,lsess,reqd,ill )
      if( ill.ne.0 ) ierr_sestbl = ill
      if( use_ifile.eq.'N' ) then
c       set filename to signal skipping in MDMAKE
        ifile = 'NONE'  
      endif
      if( ifile(1:1).ne.' ' .and.ifile(1:4).ne.'NONE' ) then
c       set default IFLSW to indicate make new I-file
        IFLSW = .TRUE.
        CALL  LOWERS(IFILE)
        inquire( file=ifile(1:nblen(ifile)),exist=fexist)
        if( fexist ) then
c         I-file exists
          if (iarg.gt.0) then
c         Command-line argument was used, so use the old I-file
          iflsw = .false.
        else
           WRITE( 6, '(/,3A,3(/,10X,A))' )
     .          '   File ', IFILE(1:NBLEN(IFILE)),
     .          ' already exists.',
     .              '1. use the old one',
     .              '2. save old and create new',
     .              '3. overwrite old I-file'
           call imenu (ians,3)
           IF( IANS .EQ. 1 )  IFLSW = .FALSE.
           IF( IANS .EQ. 2 ) THEN
             CALL UPNAM1(IFILE,IFILEB)
             IFILE=IFILEB
           ENDIF
        endif
      else
c       No I-file, so make a new one
        iflsw = .true.
      endif

C   Open the I-file
                   
      if( .not.iflsw .and. ifile(1:4).ne.'NONE' ) then
         open( lifil,file=ifile,status = 'OLD' )
      else if  ( iflsw )  then
        OPEN( lifil, FILE=IFILE, STATUS = 'UNKNOWN' )
          call report_stat('STATUS','FIXDRV','fixdrv',' ',
     . 'New Clock-polynomial (I-) file being written--see fixdrv.out',0)
        WRITE( lifil, '(A,A8,A,2(I2,A),I4,2X,2(I2,A),I2,a,a5)' )
     .       ' Site clock polynomials by ',UNAME,
     .       ' on  ', IMONTH, '/', IDAY, '/', IY,
     .       IHR, ':', IMN, ':', ISEC,'  FIXDRV Vers '
     .       , fixdrv_vers
c        second line is titles
        write(lifil,'(a)') 
     .     'SITE YEAR DOY S HR MN SEC(GPST)   EPOCH (SEC)        RATE    
     .        ACCEL (1/SEC)   CUBIC (1/SEC**2)'
c        third line is format
          write(lifil,'(a)') formti
        endif
      endif
                 
c   Open the print file only if calculating clock values (so don't overwrite an old one)
             
      if( iflsw ) then
        open (unit=lprint,file='fixdrv.out',form='formatted'
     .       , status='unknown',iostat=ioerr)
        if(ioerr .ne. 0 ) then
          call report_stat('FATAL','FIXDRV','fixdrv','fixdrv.out',
     .    'Error opening fixdrv print file: ',ioerr)
        endif
        write(lprint,*) ' FIXDRV v.'//vers
        write(lprint,*) ' '
      else
        call report_stat('STATUS','FIXDRV','fixdrv',' '
     .  ,'Old I-file used, print file fixdrv.out not written',0)
      endif


C   Read J-file name

      READ( ldfil, '(A)'  )  JFILE
      CALL  LOWERS(JFILE)
      if( jfile(2:5) .ne. 'none' ) THEN
        IF( .NOT. FCHECK(JFILE))
     .      call report_stat('FATAL','FIXDRV','fixdrv',jfile
     .                      ,'J-file not available',0)
      endif                          

c   Read the number of stations (X-files) 

      READ( ldfil, * )  nstat        


C   Read and accumulate the X-file names for this session

        DO  JSTAT = 1, nstat
          READ( ldfil,'(A)',iostat=ioerr ) XFILE(JSTAT)
          if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV','fixdrv'
     .     ,' ','Error or EOF reading D-file--check # of X-files',ioerr)
          CALL LOWERS( XFILE(JSTAT) )
        enddo        


c   Reorder the X-files for this session alphabetically, removing duplicates
 
        call sort_string( maxsit,nstat,xfile,0,xjunk
     .                  , istat1,xjunk) 
        nstat = istat1
        do i=1,nstat
           xfile(i) = xjunk(i)
        enddo
                       

c   Read session.info to get the year, doy, and span, and check these 
c   against the x-file header
       
       if( simulation ) then 
c        for simulation, we must get the year and day from S-file
         open( unit=lsfil,file=sfile,form='formatted',iostat=ioerr
     .        ,status='old')   
         if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV','fixdrv'
     .             ,sfile,'Error opening simulation-control file',ioerr)
         call getcmd(lsfil,'spanid',wcmd,lcmd,1)
         if (lcmd.le.0)  call report_stat('FATAL','FIXDRV','fixdrv',' '
     .       ,'Missing spanid in simulation-control file',0)
         ic =count_arg(wcmd)
         if( ic.lt.3) call report_stat('FATAL','MODEL','simred',' '
     .     ,'Spanid error in simulation-control file',0 )
          read(wcmd,*) year,doy,isessn     
       endif  
       xorc = xfile(1)(1:1)
       call lowers( xorc ) 
       call seschk( xorc,xfile(1),lxfil,lsesfo
     .           , year,doy,isessn,nsod,nepoch,inter,ispan )
     

c     Loop thru the stations for this session, checking for the validity of
c     station.info values and writing I-file records if requested

          DO  200  JSTAT = 1, nstat

            sname = xfile(jstat)(2:5) 
            call uppers(sname)   
            if( iflsw) write(lprint,'(//,80("_"),/,a,a4,a,i3,/)')
     .          'Site code ', sname, '   Day ', doy
c           save the station list for CFMRG and SOLVE batch files
            snames(jstat) = sname

c              Check only: read station.info and call hisub to avoid a later stop in MODEL
c              Must have valid start time from x- or c-file header (seschk returns 99999
c              if not available)

            if( fcheck('station.info') .and. nsod.ne.99999 ) then
c             isessn from seschk (=1 if session.info not available)
* MOD TAH 200203: Added AntDAZ to list of values from station.info
              call rstnfo( lstnfo,sname,year,doy,nsod
     .                    , ispan,stanam,anth,antn,ante,antdaz
     .                    , rcvcod,antcod
     .                    , htcod, radome,swver,rcvers,rcvrsn,antsn
     ,                    , istarts, istops )  
c             check hi.dat only if station.info available
              if( fcheck('hi.dat') ) then 
                 call hisub( ihi,anth,antn,ante,antcod,htcod,sname
     .                   , year,doy,isessn,offarp,warnings) 
              endif   
            endif


c             Get the clock-model type from the receiver type or sittbl. 
c        
c            KLOCK value   Clock epoch     Phase-clock model     Receiver example
c            -----------   -----------     -----------------     ----------------
c               1           0.0              none                MiniMac
c               2           polynomial       polynomial          Atomic clock or good crystal (TI4100)
c               3           epoch-by-epoch   polynomial          All modern receivers

    
            if( rcvcod.eq.'MIN6AT' ) then
               klock(jstat) = 1
            else
               klock(jstat) = 3
            endif 
            reqd = .false.
            call rdsitt (sname,5,'KLOCK',nout,line,lsite,reqd,ill)
            if( ill .lt. 0 ) then
                ierr_sittbl = ill
            else
              if( line(1:5).ne.'     ' ) then
                 read (line,'(i5)',iostat=ioerr)  klocks  
                 if( ioerr.ne.0 ) then
                    write(message,'(2a)') 
     .                 'Error reading klock from sittbl.; line=',line
                    call report_stat('FATAL','FIXDRV','fixdrv',' '
     .                              , message,ioerr)  
                 endif
                 if( klocks.ne.klock(jstat) ) then                         
                   write(message,'(a4,a,i1,a,i1,a,a6)') sname 
     .                           ,' klock from sittbl. ( ',klocks
     .                           ,' ) overrides station.info value ( '
     .                           ,klock(jstat),' ) for ',rcvcod 
                   if( iflsw) write(lprint,'(a)') message
                   call report_stat('WARNING','FIXDRV','fixdrv',' '
     .                           , message,0) 
                   klock(jstat) = klocks
                 endif
               endif
            endif

C             Check for stations on old I-file or make a new one 

            if( .not.iflsw .and. ifile(1:4).ne.'NONE' ) then
               rewind( lifil )
               read(lifil,'(//)',iostat=ioerr) 
               if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV'
     .            ,'fixdrv',' ','Error reading start of i-file',ioerr)
 
               ioerr = 0
               do while (ioerr.eq.0 )
                 read(lifil,'(a4)',err=195,end=191,iostat=ioerr) sitei
                 if( sname.eq.sitei ) goto 200
               enddo   
  191          write(message,'(a,a4,a)') 'D-file site ',sname
     .            ,' missing from I-file; clock polynomials set = 0.'
               call report_stat('WARNING','FIXDRV','fixdrv',' '
     .          , message,0)
               call report_stat('WARNING','FIXDRV','fixdrv',' '
     .     ,'(cont) possible numerical problem in SOLVE if large gap',0)
               goto 200
  195          call report_stat('FATAL','FIXDRV','fixdrv',' '
     .                      ,'Error reading site on I-file',ioerr)
            elseif ( iflsw) then
c             need K-file for all clock models except klock = 1
              IF(KLOCK(jstat) .GT. 1) THEN
                 KFILE = XFILE(JSTAT)
                 KFILE(1:1) = 'K'
                 CALL  LOWERS( KFILE )
                 KPICK = 'Y'
c                check for existence of K-file
c                if we cannot find the right one, use the most recent one available
                 IF( .NOT. fcheck (kfile) ) then
                   write(message,'(a,a16)') 'K-file does not exist:  '
     .                                     , kfile
                   call report_stat('WARNING','FIXDRV','fixdrv',' '
     .                              ,message,0)
                   buff80 = ' ' 
                   buff80(1:16) = kfile
                   buff80(6:6) = '?'
                   call getdir (buff80,100,kfiles,nk)
                   if (nk .gt. 0) then
                     kfile = kfiles(nk)
                     write(message,'(a,a16)')'FIXDRV using K-file '
     .                                      , kfile
                     call report_stat('WARNING','FIXDRV','fixdrv',' '
     .                               ,message,0)
                   else
                     call report_stat('FATAL','FIXDRV','fixdrv',' '
     .                              ,'Cannot find a K-file',0)
                   endif
                 ENDIF
               ELSE
                 KPICK = 'N'
               ENDIF
C
C             cubic or linear fitting
              IF( KPICK .EQ. 'Y' )  THEN
                reqd = .false.
                CALL RDSITT(SNAME,5,'CLKFT',NOUT,LINE,LSITE,reqd,ILL)  
                if( ill.eq.0 ) then
                  IF ( INDEX(LINE(1:NOUT),'C') .NE. 0 )  THEN
                    CLKFT = 'C'
                  ELSEIF( INDEX(LINE(1:NOUT),'L') .NE. 0 )  THEN
                    CLKFT = 'L'
                  ELSE
                    CLKFT = ' '
                  ENDIF
                else
                  ierr_sittbl=ill
                endif
              endif  
                      
C             Write the session records of the I-file
              CALL  imakef( XFILE(JSTAT), XORC, KFILE, KPICK
     .                   , CLKFT, LXFIL, LCFIL, formti, lprint )
            ENDIF
                
  200     CONTINUE
C        ---End of loop over stations--


C     Close the I-file

      close( lifil )
                 
C     Now call BMAKE to create the batch files for this solution
       
      if (iflsw ) then
        write(lprint,'(/,80("_"),/)')
        write( lprint, '(/,A,A16,A,2(I2.2,A),I4,2X,2(I2.2,A),I2.2)' )
     .     '   Batch file created by ', UNAME,
     .     ' on  ', IMONTH, '/', IDAY, '/', IY,
     .     IHR, ':', IMN, ':', ISEC
      endif    
      Call  BMAKE( nstat, snames, year, doy, nepoch, inter
     .            , dfile, mfile, tfile, gfile, ifile, jfile, lfile
     .            , simulation, sfile, xorc, xfile, klock 
     .            , ltfil, lsite, lsess,lstnfo, lsesfo, lxfil
     .            , lgfil, lcfil, proc
     .            , ierr_sestbl, ierr_sittbl, fixdrv_vers, scratch_dir )
      close( lsess )
      close( lsite )

c     End of FIXDRV

      if (ierr_sestbl.ne.0 .or. ierr_sittbl.ne.0 ) then
        if ( ierr_sestbl.ne.0 )  call report_stat('WARNING'
     .            ,'FIXDRV','fixdrv',' ','Errors in sestbl entries',0)
        if ( ierr_sittbl.ne.0 )  call report_stat('WARNING'
     .            ,'FIXDRV','fixdrv',' ','Errors in sittbl entries',0)
         call report_stat('FATAL','FIXDRV','fixdrv',' '
     .        ,'Sestbl or sittbl errors--see GAMIT.warning',0)
      else
         call report_stat('STATUS','FIXDRV','fixdrv',' ','Normal end',0)
      endif
      stop
      end

