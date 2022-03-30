       SUBROUTINE BMAKE( nstat, snames, year, doy, nepoch, inter
     .              , dfile, mfile, tfile, gfile, ifile, jfile, lfile
     .              , simulation, sfile, xorc, xfile, klock
     .              , ltfil, lsite, lsess, lstnfo, lsesfo, lxfil
     .              , lgfil, lcfil, proc
     .              , ierr_sestbl, ierr_sittbl, fixdrv_vers,scratch_dir)

C     program to make batch file for GAMIT runs
C
C     S.Shimada              original at NRCDP
C     S.Shimada   90/02/08   modified at IGPP
C     S.Shimada   90/09/21   modified at IGPP
C     Y. Bock, R. King, K. Feigl : 91/03/15 ff. - see FVERSN   
C
C     arguments               
c        nstat  : number of stations
C        SNAMES : 4-chracter site codes      
C        year   : 4-digit integer year of session 
c        doy    : integer day-of-year of session
c        nepoch : number of epochs
c        inter  : sampling interval (sec) 
C        File names written to B-file :
C          DFILE  D-file
c          MFILE  M-file
C          TFILE  T-file
C          GFILE  G-file
C          IFILE  I-file
C          JFILE  J-file
C          LFILE  L-file or apr file
C        XFILE  : Starting X- (or C-) file name list
C        XORC   : First character of X-,C-files
c        SESTBL : session control file name (sestbl.)
c        SITTBL : site control file name (sittbl.)
c        Logical unit numbers for files read :
c           LTFIL    T-file         11
c           LSITE    sittbl.        13
c           LSESS    sessbl.        14
c           LSESFO   session.info   15
c           LSTNFO   station.info   19
c           LXFIL    X-file         22
c           LGFIL    G-file         23
c           LCFIL    C-file         24
C        IERR_SESTBL, IERR_SITTBL - return codes for reading sestbl and sittbl (0 if ok)
c
C     Local files
C        scratch file                    LUN= 8
c        BFILE  : primary batch file     LUN=LBFIL=17 (also assigned in ARMAKE
C        BFIL2  : secondary batch file   LUN=20,21    (MDMAKE, SOMAKE, etc)

C     Dimension parameters
c        MAXSIT : maximum number of sites 
C        MAXSAT : maximum number of satellites
C        MAXCFL : maximum X-(C-)files in a solution

C     variables
C        NTSAT  : number of satellites in G- and T-file
C        ITB, ITSTP : T-file begin and stop time
C                     (last two digit of year, month, day)
C        TBB, TSTP  : T-file begin and stop time (hour,minute,second)

C     subroutines and functions
C        BMAKE
C          +-ARMAKE +-NBLEN .
C          |        +-RDSEST --UPPERS .
C          +-DBMAKE .
C          +-ACMAKE +-NBLEN .
C          |        +-UPNAM1 --NEWCHR .
C          +-CHRNUM .
C          +-C          +-GETDAT --PM_$GET_SID_TXT .
C          +-GETTIM --CAL_$DECODE_LOCAL_TIME .
C          +-GETUSR --PM_$GET_SID_TXT .
C          +-LFNTOP .
C          +-LOWERS .
c          +-GRDMAKE .
C          +-MDMAKE +-LFNTOP .
C          |        +-LOWERS .
C          |        +-NBLEN .
C          |        +-RDSITT --UPPERS.
C          |        +-UPNAM1 --NEWCHR .
C          +-NBLEN
C          +-RDSEST --UPPERS.
C          +-SATNMR .
C          +-XCSATS +XHDRED
C          |        +CHDRED +-READC1 +-FERROR .
C          |        |        |
C          |        |        +-READC2 +-FERROR .
C          |        |        +-READC3 +-FERROR .
C          |        +-COPENS --FERROR .
C          |        +-LFNTOP .
C          |        +-LOWERS .
C          |        +-NBLEN .
C          |        +-NEXTT1 --LEAPYR .
C          +-TGSATS +-SATNMR
C          |        +-THDRED
C          +-SOMAKE +-LFNTOP .
C          |        +-NBLEN .
C          |        +-NXEPOC +-CHDRED +-READC1 +-FERROR .
C          |        |        |        +-READC2 +-FERROR .
C          |        |        |        +-READC3 +-FERROR .
C          |        |        +-COPENS --FERROR .
C          |        |        +-LFNTOP .
C          |        |        +-LOWERS .
C          |        |        +-NBLEN .
C          |        +-NXSAT  +-CHDRED +-READC1 +-FERROR .
C          |        |        |        +-READC2 +-FERROR .
C          |        |        |        +-READC3 +-FERROR .
C          |        |        +-COPENS --FERROR .
C          |        |        +-LFNTOP .
C          |        |        +-LOWERS .
C          |        |        +-NBLEN .
C          |        +-RDSEST --UPPERS.
C          |        +-RDSITT --UPPERS.
C          |        +-SOLRUN +-GETDAT --PM_$GET_SID_TXT .
C          |                 +-GETTIM --CAL_$DECODE_LOCAL_TIME .
C          |                 +-GETUSR --PM_$GET_SID_TXT .
C          |                 +-LJUST  .
C          |                 +-SOLATM .
C          |                 +-SOLCLK .
C          |                 +-SOLCRD .
C          |                 +-SOLORB .
C          |                 +-UPNAM1 --NEWCHR .
C          +-ACMAKE
C          +-NGMAKE .
C          +-UPNAM1 +-NEWCHR .
C          +-UPNAM2 +-UPCHR  .
c          +-XCHECK +-XHDRED
c                   +-RSTNFO

c      IMPLICIT  REAL*8 (A-H,O-Z)
      implicit none

      include '../includes/dimpar.h'
      include '../includes/global.h'

      LOGICAL     found, fcheck, debug, reqd, missing_sv, simulation
     .          , kbit,zenest,gradest,eof
     .          , otll,otlg,atmll,atmlg,atll,atlg,metl,metg,mapl,mapg

      CHARACTER* 1  XORC,  EXPO
     .             , first_arc,delall, delmod,  delaut
     .             , finarc, xcompress, lowerc, run_ctox
     .             , autcln_post,autinchr,autoutchr,outchr,yawmod,buf1
     .             , rewgt,autcln_redo,dquote,skd,atmlflag 
     .             , ans,autwl,znthdl,gradient,autcln_clk
      CHARACTER* 2  CYR,atide,buf2
      CHARACTER* 3  CHRNUM, cdoy,avedit_opt,rawcln,proc,buf3
     .             , zenmod,gradmod
      character* 4  snames(maxsit),delecl,pre_post
     .              ,dmap,wmap,stnerr,saterr,ionsrc,ansnum,anshr,sitnam
c**   probably temporary
      character* 6  magfield
      CHARACTER* 5  typana, EXPRMT,scan_control,fixdrv_vers,saterropt
     .             ,buf5
c**   character* 5  geoddt, datum - uncomment if 'GEODETIC DATUM' restored
      character* 6  iterate,etidemod
      character* 7  tlupdt,stnerropt                 
      character*10  grdfiles(10)
      character*12  xtemplate
      character*15  observ,observp
      CHARACTER*16  BFILE, BFIL2, TFILE, GFILE, IFILE, DFILE, LFILE
     .           ,  xfile(maxsit), CFILE(maxsit)
     .           ,  JFILE, CFILE1, MFILE, qfile, qfilep
     .           ,  cfile2,yfile,nfile,sfile,ffile,satnam(maxsat)
c     .          ,  GDATUM       
      character*20  qfileopt,metsource
      character*25  autclnopt
      character*30   eclopt
      character*80  autcln_cmd,message,line,scratch_dir
      character*256 blnk256

      integer*4  iscrn,lbfil,lxfil,lsite,lstnfo,lsesfo,nstat
     .       , jstat,iday,kcol,nsat,ntsat,check_flg
     .       , lsess,idrv,lexp,year,doy,imn,idatum
     .       , iystp,month,idoystp,itb(3),itstp(3),ixb(3),ixstp(3)
     .       , ill,ierr_sestbl,ierr_sittbl,itflag
     .       , ihr,lcfil,lgfil,ltfil
     .       , klock(maxsit)
     .       , nintrs,iwkn,iwknstp,icheck(maxsat),isessn
     .       , nepoch,inter,idec,idecp
     .       , itsat(maxsat),ixsat(maxsat),iarray,nblen,nchr
     .       , norb,glen,iper,nzen,ngrad,ietide
     .       , ioerr,i,nyr,ndoy,nsod,nyr2,min,ngfiles,i0,nout

      real* 8  sec,tbb(3),tstp(3),txb(3),txstp(3),sow,sowstp,tspan
     .       , spanhr,spanday,utcoff,sod,zenint,gradint

c       External functions
      integer*4 mchkey

c       Local files
      PARAMETER ( ISCRN=6, LBFIL=17 )
c      DATA  GDATUM / 'gdetic.dat' /

c=============================================================================


c    ** Initialization

      blnk256 = ' '
c     Set double quote for batch-file shell scripts (only way to get printed correctly on all systems)
      dquote=char(34)
      idatum = 0

C     Initialize export orbits for first calls to ARMAKE
      EXPO = 'N'

C     Create the batch file name from the D-file name
              
      bfile = dfile
      BFILE(1:1) = 'b'
      KCOL = INDEX( DFILE, '.' )
      BFILE(KCOL+1:KCOL+3) = 'bat'
      BFIL2 = BFILE
      OPEN( lbfil, FILE=BFILE, FORM='FORMATTED',status='unknown'
     .     ,iostat=ioerr)
      if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV','bmake',bfile
     .  ,'Error opening primary batch file',ioerr)  


c      Make the batch file a UNIX shell script and remove old status, error, and warning file
               
      write(17,'(a)',iostat=ioerr) '#!/bin/csh'  
      if( ioerr.ne.0 ) call report_stat('FATAL','FIXDRV','bmake',bfile
     .  ,'Cannot write to primary batch file--check permissions',ioerr)
c** rwk/tah/scm: No longer remove existing status/warning/fatal files since these
c**              written into now by makexp, makej, makex, fixdrv.  Remove in sh_gamit
c      write(17,'(2a)') '# Remove any existing '
c     .                ,'GAMIT.status, .warning, .fatal files'
c      write(17,'(a)') "\\rm *.status"
c      write(17,'(a)') "\\rm *.warning"
c      write(17,'(a)') "\\rm *.fatal"
c      write(17,'(a)') "/bin/rm GAMIT.status >& /dev/null"
c      write(17,'(a)') "/bin/rm GAMIT.warning >& /dev/null"
c      write(17,'(a)') "/bin/rm GAMIT.fatal >& /dev/null"
      write(17,'(a)') '#'
               

c     Create the file-name ids from the year and day-of-year
      
      cyr = '  '
      cdoy = '   '   
      call even_minute(year,doy,sod,nyr,ndoy,nsod)
      nyr2 = mod(nyr,100)
      write(buf2,'(i2)') nyr2
      read(buf2,'(a2)') cyr
      write(buf3,'(i3)') ndoy
      read(buf3,'(a3)') cdoy
      if(cdoy.ne.'   ') then
        if(cdoy(1:1).eq.' ') cdoy(1:1) = '0'
        if(cdoy(2:2).eq.' ') cdoy(2:2) = '0'
      endif                 


C      Initialize the counter for module batch files

      idrv = 1

   
c========================================================================================

C     ** Read the sestbl to get program controls
                                          

C     Choice of experiment (parameters estimated)

c       BASELINE (default) : no orbits or EOP estimated
c       RELAX.             : estimate orbits and EOP along with coordinates
c       ORBIT              : estimate orbits and EOP only, no coordinates

      reqd = .false.   
      CALL  RDSEST( 20,'Choice of Experiment',5,EXPRMT,LSESS,reqd,ILL )
      if( exprmt.eq.'     ' ) exprmt = 'BASEL'
      if ( ill.ne.0 ) ierr_sestbl = ill
      LEXP = 0
      IF( EXPRMT .EQ. 'BASEL' )  LEXP = 1
      IF( index(EXPRMT,'RELAX') .gt. 0 )  LEXP = 2
      IF( EXPRMT .EQ. 'ORBIT' )  LEXP = 3
      IF( LEXP .EQ. 0 ) then
        write(message,'(a,a5)') 'Improper Choice of Experiment',exprmt
        call report_stat('WARNING','FIXDRV','bmake',' ',message,0)
        ierr_sestbl = -1
      endif


c     Choice of observable and ambiguity resolution 

c      Final solution ('a' q-file)
c        LC_AUTCLN   : dual-freq with good P-code; WLs resolved in autcln  (default)
c        LC_HELP     : dual-freq with codeless L2; WLs resolve with ion constraint    
c        LC_ONLY     : dual-freq with no ambiguity resolution
c        L1_ONLY     : dual-freq cleaning but L1-only solution
c        L2_ONLY     : dual-freq cleaning but L2-only 
c        L1L2_INDEP  : dual-freq cleaning but L1 + L2 solution (no ion constraint)
c        L1&L2       : dual-freq cleaning but L1 + L2 solution with ion constraint
c        L1_RECEIVER : single-freq cleaning(must comment out 'noL1only in autcln.cmd)

c      Preliminary solution ('p' q-file)  
c        LC_ONLY     : dual-frequency with no ambiguity resolution (default for 
c                      all final-solution observables except L1_RECEIVER)
c        L1_RECEIVER : single-freq cleaning 

      observ = 'LC_AUTCLN'                                              
      call rdsest( 20,'Choice of Observable',15,observ,lsess,reqd,ill )
      if( ill .ne. 0 ) ierr_sestbl = ill                                   
      if( observ(1:7).eq.'L1_SING' ) observ = 'L1_RECEIVER    '
      if( observ(1:7).eq.'L1_ONLY' ) observ = 'L1_ONLY        '
      if( observ(1:7).eq.'L2_ONLY' ) observ = 'L2_ONLY        '
      if( observ(1:6).eq.'L1&L2 ' .or. observ(1:6).eq.'L1,L2 ' )
     .                               observ = 'L1&L2          '
      if( observ(1:7).eq.'L1,L2_I'.or.observ(1:7).eq.'L1L2_I' ) 
     .    observ = 'L1,L2_INDEPEND.'
      if( observ(1:3).eq.'LC ' )     observ = 'LC_ONLY        ' 
      if( observ(1:4).eq.'LC_A')     observ = 'LC_AUTCLN      '
      if( observ(1:4).eq.'LC_H')     observ = 'LC_HELP        '
      if( observ(1:5).eq.'LC_HE')    observ = 'LC_HELP        '
c     if LC_AUTCLN, check that this set in the autcln command file
      if( observ(1:4).eq.'LC_A') then
        if( fcheck('autcln.cmd') ) then  
          open(unit=33,file='autcln.cmd',status='old')
          found = .false.              
          eof = .false.
          do while ( .not.found .and. .not.eof ) 
            read(33,'(a)',iostat=ioerr) line  
            if( ioerr.eq.-1 ) then
              eof =.true.
            elseif( ioerr.lt.0 ) then
              call report_stat('FATAL','FIXDRV','bmake',' '
     .        ,'Error reading autcln.cmd file',ioerr)
            endif                            
            if( line(1:1).eq.' '.or. line(1:4).eq.'POST') then  
               call uppers(line)
               if(mchkey(line,'LC_AUTCLN',256,9).gt.0) found=.true.
            endif
          enddo
          if( .not.found ) then
            call report_stat('FATAL','FIXDRV','bmake',' '
     .       ,'LC_AUTCLN in sestbl. but missing from autcln.cmd file',0)
          endif   
          rewind(33)
        else
          call report_stat('WARNING','FIXDRV','bmake',' '
     .       ,'LC_AUTCLN in sestbl. autcln.cmd file not present',0)
        endif
      endif    
      observp = 'LC_ONLY' 
      if( observ(1:8).eq.'L1_RECEI') observp = observ
      call rdsest( 20,'Quick-pre Observable',7,observp,lsess,reqd,ill)  
      if( ill.ne.0 ) ierr_sestbl = ill  
      if( observp(1:7).eq.'L1_SING' ) observp = 'L1_RECEIVER    '
      if( observp(1:3).eq.'LC ' )     observp = 'LC_ONLY        '  
       

c     Type of analysis (number of iterations and raw/clean data)
             
c**    rwk 060818:  The old 0-ITER, 1-ITER, 2-ITER, QUICK, 1-CLEAN, 2-CLEAN, 3-CLEAN
c**    have been removed, leaving only full solution w/optional iterations for autcln.
c      New scheme:
c         0-ITER : single full solution 
c         1-ITER : invokes autcln postfit option 'R' unless 'AUTCLN  Postfit = Y' (no redo).
c                  1-ITER with 'AUTCLN Postfit = N' not allowed
c         For backward compatibility, 0-ITER  paired with AUTCLN Postfit = Y or R becomes 1-ITER

      reqd = .true.                  
      call rdsest( 16,'Type of Analysis',5,typana,LSESS,reqd,ill )  
      reqd = .false.
      call  rdsest( 11, 'Data Status', 3, rawcln, LSESS, reqd, ILL)
      if( simulation ) then
        rawcln = 'CLN'  
        typana = '0-ITE' 
        autcln_post = 'N'  
        autcln_redo = 'N'
      else  
c       this just to cover case of typing 'CLEAN' instead of 'CLN'
        if( rawcln(1:2).eq.'CL' ) rawcln = 'CLN'
        if( rawcln.eq.'   ' ) then
           if( typana(3:5).eq.'CLE' ) then
              rawcln = 'CLN'
           else
              rawcln = 'RAW'
           endif
        endif    
        if( typana.eq.'     ') typana = '1-ITE' 
        if( typana.eq.'PREFI'.or.rawcln.eq.'CLN' ) then    
          autcln_post = 'N'
          autcln_redo = 'N'
        elseif( typana.eq.'0-ITE' ) then  
          autcln_post = 'N'
          autcln_redo = 'N'
        elseif( typana.eq.'1-ITE' ) then              
          autcln_post = 'Y'
          call rdsest( 11,'AUTCLN Redo',1,autcln_redo,lsess,reqd,ill) 
          if( autcln_redo.eq.' ' ) autcln_redo = 'Y'
        endif  
c**rwk 090528  temporary trap old combination and set to new    
        if( typana.eq.'0-ITE' ) then
           call rdsest( 14,'AUTCLN Postfit',1,buf1,lsess,reqd,ill)  
           if( buf1.eq.'Y'.or.buf1.eq.'R' ) typana = '1-ITE'
           if( buf1.eq.'Y' ) then
              autcln_post = 'Y'
              autcln_redo = 'N'
           elseif( buf1.eq.'R' ) then
              autcln_post = 'Y'
              autcln_redo = 'R'
           endif             
        endif
c----- end temporary trap
      endif
   
c     See if an extra autcln run is needed for clocks or emipircal antenna model

      reqd = .false.                  
      call rdsest( 12,'AUTCLN clock',1,autcln_clk,LSESS,reqd,ill ) 
 

c     Trap incompatibilities
   
c     1-ITER now implies use of autcln 'POST' options though possibly w/o 'use_postfit' 
      if( typana.eq.'1-ITE' .and. autcln_post.eq.'N' )
     .   call report_stat('FATAL','FIXDRV','bmake',' '
     .       ,'1-ITER with Autcln Postfit = N not allowed',0)

c    Require clean data if C-files input, to avoid file-naming problems
      if( xorc.eq.'c' .and. rawcln. ne. 'CLN' ) then
            call report_stat('FATAL','FIXDRV','bmake',' '
     .         ,'C-files in D-file allowed only with clean data',0)
      endif
                                                       
c     if LC_AUTCLN, check that this set in the autcln command file
      if( observ(1:4).eq.'LC_A') then
        if( fcheck('autcln.cmd') ) then  
          open(unit=33,file='autcln.cmd',status='old')
          found = .false.              
          eof = .false.
          do while ( .not.found .and. .not.eof ) 
            read(33,'(a)',iostat=ioerr) line  
            if( ioerr.eq.-1 ) then
              eof =.true.
            elseif( ioerr.lt.0 ) then
              call report_stat('FATAL','FIXDRV','bmake',' '
     .        ,'Error reading autcln.cmd file',ioerr)
            endif                            
            if( line(1:1).eq.' '.or. line(1:4).eq.'POST') then  
               call uppers(line)
               if(mchkey(line,'LC_AUTCLN',256,9).gt.0) found=.true.
            endif
          enddo
          if( .not.found ) then
            call report_stat('FATAL','FIXDRV','bmake',' '
     .       ,'LC_AUTCLN in sestbl. but missing from autcln.cmd file',0)
          endif   
          rewind(33)
        else
          call report_stat('WARNING','FIXDRV','bmake',' '
     .       ,'LC_AUTCLN in sestbl. autcln.cmd file not present',0)
        endif
      endif


c    Read the yaw control and get the y-file name
       
c     set default
      yawmod = 'Y'
      reqd = .false.
      call rdsest( 9,'Yaw Model',1,buf1,lsess,reqd,ill )
      if( ill.ne.0 ) ierr_sestbl = ill
      if( buf1.ne.' ' ) yawmod = buf1
c     if yaw model, construct y-file names from t-file name
c        yfile is a binary file of attitude written by yawtab and used  MODEL.  
c        It has 't' as the sixth character rathen than the year. 
      yfile = '                '
      if( yawmod.eq.'Y') then
        yfile = tfile 
        yfile(1:1) = 'y'   
        yfile(6:6) = 't'
      endif    
                

c    Set name of autcln base command file (could move to later)

      autcln_cmd = 'autcln.cmd'
                   

c    Edit autcln command file for bad clocks and eclipse data?

      reqd = .false.
      avedit_opt = 'NO '
      call rdsest
     .   ( 24, 'Edit AUTCLN Command File',3,avedit_opt,lsess,reqd,ill)
      if ( ill.ne.0 ) ierr_sestbl = ill  
c     check if eclipse or post-eclipse data edited (scan arcout)
      glen = nblen( gfile )
      reqd = .false.
      call rdsest( 19,'Delete eclipse data',4,delecl,lsess
     .           , reqd,ill )
      if ( ill.ne.0 ) ierr_sestbl = ill  
      eclopt = ' ' 
      if( delecl(1:1).ne.'n' .and. delecl(1:1).ne.'N' ) then 
        read ( cdoy,'(i3)') doy
        iper = index( tfile,'.')
        if( delecl(1:3).eq.'all' .or. delecl(1:3).eq.'ALL' ) then
            write(eclopt,'(2a,1x,i3.3,1x,a4)')
     .           '-autecl ',yfile(1:nblen(yfile)),doy,'Y 30'
        elseif( delecl(1:4).eq.'post'.or.delecl(1:4).eq.'POST') then
            write(eclopt,'(2a,1x,i3.3,1x,a4)')
     .           '-autecl ',yfile(1:nblen(yfile)),doy,'N 30'
        endif
      endif


c    Read the C, L, and T-file update controls

c     Defaults set, so no entry is required: RDSEST will return blank
      reqd = .false.

      call rdsest(11, 'Initial ARC', 1, first_arc, lsess, reqd, ill )
      if ( ill.ne.0 ) ierr_sestbl = ill
      if ( first_arc.eq.' ' ) then
        first_arc = 'Y'  
      endif  
c      print *,'lexp first_arc ',lexp,first_arc
      call rdsest( 16, 'Update T/L files', 7, tlupdt, lsess, reqd, ill )
      if ( ill.ne.0 ) ierr_sestbl = ill
      if ( tlupdt.eq.'       ')   tlupdt = 'L_ONLY '
c     Answer to delxxx is 'Y' (delete) or 'N' (do not delete)
      call rdsest( 24, 'Delete all input C-files'
     .           ,  1, delall, lsess, reqd,ill )
      if ( ill.ne.0 ) ierr_sestbl = ill
      if( delall.eq.' ' ) delall = 'N'
      call rdsest( 26, 'Delete MODEL input C-files'
     .            , 1, delmod, lsess, reqd, ill )
      if ( ill.ne.0 )  ierr_sestbl = ill
      if( delmod.eq.' ' ) then
        if( delall.eq.'Y' ) then
          delmod = 'Y'
        else
          delmod = 'N'
        endif
      endif   


c     C-file removal control
                       
      call rdsest( 27, 'Delete AUTCLN input C-files'
     .             ,1, delaut, lsess, reqd, ill )
      if ( ill.ne.0 ) ierr_sestbl = ill    
      if( delaut.eq.' ') then
        delaut = 'Y'
      endif


c     Final ARC ?

      call rdsest( 9, 'Final ARC', 1, finarc, lsess, reqd, ill )
      if ( ill.ne.0 ) ierr_sestbl = ill
      if( finarc.eq.' ') finarc = 'N'


c     Iteration control command line  (does c-file iteration still work?)

c     iterate=XFILES/CFILES
      reqd= .false.
      call  rdsest( 9,'Iteration',6,iterate,lsess,reqd,ill )
      if ( ill.ne.0 ) ierr_sestbl = ill
      if( iterate.eq.' ') iterate = 'XFILES'

c    Read the control for executing SCANDD

      reqd = .false.
      call rdsest( 14,'SCANDD control',5,scan_control,lsess,reqd,ill )
      if ( ill.ne.0 ) ierr_sestbl = ill
      if( scan_control(1:1).eq.' ') scan_control = 'IFBAD'


c    Read the first X- or  C-file to get the satellite numbers and integration span

      if( .not.simulation .and. fcheck(xfile(1)) ) then
         call xcsats( xfile(1), xorc, lxfil, nsat, ixsat
     .              , ixb, txb, ixstp, txstp )
      else
         call report_stat('WARNING','FIXDRV','bmake',' '
     . ,'No X- or C-file available; SVs determined from G- or T-file',0)
      endif
c     The first X- or C-file array determines the total satellites at
c     present since SOLVE cannot handle mismatched satellite numbers or order
c     X- or C-files must all agree in multi-session runs
                       

     
c    Get the total span in hours for setting number of zenith and gradient-delay parameters
         
      spanhr = (inter*nepoch)/3600.d0


c    Read the G- or T-file to check for satellite consistency
               
      call tgsats( tfile, gfile, ltfil, lgfil, iscrn 
     .           , ntsat, itsat, satnam, nintrs, norb
     .           , itb, tbb, itstp, tstp )  
cd      print *,'BMAKEE after tgsats ntsat itsat satnam ',ntsat
cd      do i=1,ntsat
cd        print *,itsat(i),satnam(i)
cd      enddo
      if( first_arc.eq.'Y' .and. nintrs.gt.3 )
     .   call report_stat('WARNING','FIXDRV','bmake',' '
     .          ,'A T-file with partials will be overwritten ',0)  
      if( fcheck(xfile(1)) ) then
        missing_sv = .false.
        do i=1,nsat
          icheck(i) = iarray( ixsat(i),itsat,ntsat )
          if( icheck(i).eq.0 ) missing_sv = .true.
        enddo
        if( missing_sv ) then
          write(iscrn,'(a)')
          do i=1,nsat
            if( icheck(i).eq.0 ) then
              write(message,'(a,i2,a)')
     .             'PRN ',ixsat(i),' in X-file not in G- or T-file'
              call report_stat('WARNING','FIXDRV','bmake',' ',message,0)
            endif
          enddo            
          call report_stat('WARNING','FIXDRV','bmake',' '
     .   ,'SVs missing from G- or T-file; see list in GAMIT.warning',0)
        endif
      else
        nsat = ntsat
        do i=1,ntsat
          ixsat(i) = itsat(i)
        enddo
      endif

c    If all else fails, try getting the satellite list from session.info
          
      if( .not.fcheck(xfile(1)) .and. .not.fcheck(tfile) .and.
     .    .not.fcheck(gfile) ) then
          if ( fcheck('session.info') ) then
            debug = .false.
            check_flg = 2.
            call rsesfo( lsesfo,debug,check_flg,year,doy,isessn
     .                 , ihr,min,inter,nepoch,nsat,ixsat,found )
            if (found )  then
              call report_stat('WARNING','FIXDRV','bmake',' '
     .      ,'No X, C, G, or T-file: SV list taken from session.info',0)
              call report_stat('WARNING','FIXDRV','bmake',' '
     .      ,'No X, C, G, or T-file: Integration span not known',0)
            else
              call report_stat('FATAL','FIXDRV','bmake',' '
     .    ,'Cannot get SV list from X, C, G, T, or session.info file',0)
            endif
          else
            call report_stat('FATAL','FIXDRV','bmake',' '
     .    ,'Cannot get SV list from X, C, G, T, or session.info file',0)
          endif
      endif                               
             
                       
c    Read the zenith delay and gradient information

      reqd = .true.
      zenest = .false.
      call rdsest( 23, 'Zenith Delay Estimation', 1, znthdl,
     .             lsess, reqd, ill )
      if( ill.ne. 0 ) ierr_sestbl = ill
      if( znthdl.eq.'Y' ) then          
        zenest = .true.   
        reqd = .false.
        call rdsest( 10,'Zenith Mod',3,zenmod,lsess,reqd,ill )
        if( ill.ne.0 ) ierr_sestbl = ill
        if( zenmod(1:1).eq.' ') zenmod = 'PWL' 
c       set default and see if nzen explicit in the sestbl
        call rdsest( 10,'Number Zen',4,ansnum,lsess,reqd,ill )
        if( ill.ne. 0 ) ierr_sestbl = ill
        if( ansnum.ne.'    ' ) read(ansnum,'(i4)') nzen     
        if( nzen.gt.maxatm ) then
          write(message,'(a,i3,a,i3,a)') 'Number Zen from sestbl. ('
     .          ,nzen,') > maxatm in dimpar.h (',maxatm,')'
          call report_stat('FATAL','FIXDRV','bmake',' ',message,0)
        endif
        if( nzen.eq.1 ) then
          zenmod = 'con'
        else
          reqd = .false.
          call rdsest( 12,'Interval Zen',4,anshr,lsess,reqd,ill ) 
          if( ill.ne. 0 ) ierr_sestbl = ill
          if( anshr.ne.'   ' ) then 
            read(anshr,'(f4.0)') zenint  
            if (spanhr.le. 0.d0 ) call report_stat('FATAL','FIXDRV'
     .         ,'bmake',' ','Span length zero for numzen calculation',0)
            if( zenmod.eq.'PWL' ) then
              nzen = spanhr/zenint + 1                              
            else
              nzen = spanhr/zenint
            endif
            write(message,'(a,i3,a,f4.1,a)') 'Setting numzen = ',nzen
     .           ,' from zenint = ',zenint,' hr'
           call report_stat('STATUS','FIXDRV','bmake',' ',message,0)
         endif
        endif
c        reading the number from the sittbl. no longer supported
         sitnam = cfile(i)(2:5)
         call rdsitt( sitnam,4,'NZEN',nout,line,lsite,reqd,ill )
         if( line(1:4).ne.'    ' ) then 
             call report_stat('WARNING','FIXDRV','bmake',' '
     .        , 'NZEN from sittbl. now ignored',0)
         endif
      else
        zenest = .false.
c       zenith delay partials always on c-file even if not estimated
        nzen = 1
        zenmod = 'con'
      endif
      reqd = .false.
      gradest = .false.
      call rdsest( 20, 'Atmospheric gradient', 1, gradient,
     .             lsess, reqd, ill )
      if( ill.ne. 0 ) ierr_sestbl = ill
      if( gradient.eq.'Y' ) then     
        gradest = .true.
        call rdsest( 12,'Gradient Mod',3,gradmod,lsess,reqd,ill )
        if( ill.ne.0 ) ierr_sestbl = ill
        if( gradmod(1:1).eq.' ') gradmod = 'PWL'

c       set default and see if ngrad explicit in the sestbl
        ngrad = 1
        call rdsest( 11,'Number grad',4,ansnum,lsess,reqd,ill )
        if( ill.ne. 0 ) ierr_sestbl = ill
        if( ansnum.ne.'    ' ) read(ansnum,'(i4)') ngrad
        if( ngrad.gt.maxgrad ) then
          write(message,'(a,i3,a,i3,a)') 'Number grad from sestbl. ('
     .          ,ngrad,') > maxgrad in dimpar.h (',maxgrad,')'
          call report_stat('FATAL','FIXDRV','bmake',' ',message,0)
        endif
        if( ngrad.eq.1 ) then
          gradmod = 'con'
        else
          reqd = .false.
          call rdsest( 13,'Interval grad',4,anshr,lsess,reqd,ill ) 
          if( ill.ne. 0 ) ierr_sestbl = ill
          if( anshr.ne.'   ' ) then 
            read(anshr,'(f4.0)') gradint
            if (spanhr.le. 0.d0 ) call report_stat('FATAL','FIXDRV'
     .       ,'bmake',' ','Span length zero for numgrad calculation',0)
            if( gradmod.eq.'PWL' ) then
              ngrad = spanhr/gradint + 1                              
            else
              ngrad = spanhr/gradint
            endif
            write(message,'(a,i3,a,f4.1,a)') 'Setting numgrad = ',ngrad
     .           ,' from gradint = ',gradint,' hr'
            call report_stat('STATUS','FIXDRV','bmake',' ',message,0)
          endif
         endif
      else
        ngrad = 0 
        gradest = .false.
      endif


c    Read the control for source of (2nd & 3rd order) ionospheric corrections
      
      reqd = .false.
      call rdsest( 9,'Ion model',4,ionsrc,lsess,reqd,ill )
c       only blank (no model) and GMAP allowed for now   
c       if GMAP, construct the f-file name 
      if( ionsrc.eq.'GMAP' ) then
        ffile = lfile
        ffile(1:1) = 'f'  
        ffile(6:6) = cyr(2:2)      
        call rdsest(9, 'Mag field',6,magfield,lsess,reqd,ill ) 
* MOD TAH 200125: Make IGRF123 default 
        if( magfield(1:2).eq.'  ' ) magfield = 'IGRF13'
        if( magfield.ne.'DIPOLE'.and.magfield.ne.'IGRF10'.and.
     .      magfield.ne.'IGRF11'.and.magfield.ne.'IGRF12'.and.
     .      magfield.ne.'IGRF13' ) then 
           call report_stat('FATAL','FIXDRV','bmake',' '
     .        ,'Mag field entry not recognized',0)
        endif
      else
        ionsrc = 'NONE'
        ffile = ' ' 
      endif    
      

c    Read the controls for Earth tides, tidal and non-tidal loading, and
c    meterological inputs, needed for both GRDTAB and MODEL
                                    
c     ----solid-earth and pole tides
c *   default changed to include ocean tides by rwk/tah 021029 
c *   default changed to include removal of mean pole of 2000 for pole tide by rwk/tah 051101
c *     ietide = 7
c *     ietide = 15
      ietide = 31
      reqd = .false.
c     'Tide model' changed to 'Tides applied' but check both for backward 
c      compatibility.  RWK 050217  
      call rdsest( 13,'Tides applied',2,atide,lsess,reqd,ill )
      if( atide.eq.'  ' ) then  
        call rdsest( 10,'Tide Model',2,atide,lsess,reqd,ill ) 
      endif            
      call rdsest( 11,'Etide model',6,etidemod,lsess,reqd,ill )   
      if( etidemod(1:1).eq.' ' ) then
         etidemod = 'IERS03'
      endif            
      if( ill.ne.0 ) ierr_sestbl = ill  
      if( atide.ne.'  ' ) read(atide,'(i2)') ietide
      otll = .false.
      otlg = .false. 
c     ----ocean loading
      if( kbit(ietide,4) ) then  
        call rdsest(12,'Use otl.list',1,ans,lsess,reqd,ill )
        if( ans.eq.'Y' ) then
          if( .not.fcheck('otl.list') ) 
     .      call report_stat('FATAL','FIXDRV','bmake','olt.list'
     .        ,'Ocean loading list file requested but not available',0)
          otll = .true.
        endif
        if( .not.fcheck('otl.grid') ) then     
          if( .not.otll )  then
            call report_stat('FATAL','FIXDRV','bmake',' '
     .          ,'Ocean loading requested no list or grid file',0)
          else
            call report_stat('WARNING','FIXDRV','bmake',' '
     .      ,'No ocean loading grid file--list file must be complete',0)
          endif
        else
          otlg = .true. 
        endif     
      endif   
c     -----atmospheric tidal loading
      atll = .false.
      atlg = .false.
      if( kbit(ietide,6) ) then    
        call rdsest(12,'Use atl.list',1,ans,lsess,reqd,ill )
        if( ans.eq.'Y' ) then
          if( .not.fcheck('atl.list') ) 
     .      call report_stat('FATAL','FIXDRV','bmake','atl.list'
     .     ,'Atm tidal loading list file requested but not available',0)
c         atll = .true.
        endif
        if( .not.fcheck('atl.grid')) then   
          if( .not.atll ) then
            call report_stat('FATAL','FIXDRV','bmake',' '
     .  ,'Atm tidal loading requested but no list or grid file',0)
          else
            call report_stat('WARNING','FIXDRV','bmake',' '
     .  ,'No atm tidal loading grid file--list file must be complete',0)
          endif
        else
          atlg = .true.
        endif     
      endif            
c     ------non-tidal atmospheric loading
      atmll = .false.
      atmlg = .false.   
      atmlflag = 'N'   
      call rdsest(9,'Apply atm',1,ans,lsess,reqd,ill)
      if( ans.eq.'Y') then
        atmlflag = 'Y'
        call rdsest(13,'Use atml.list',1,ans,lsess,reqd,ill )
        if( ans.eq.'Y' ) then
          if( .not.fcheck('atml.list') ) 
     .      call report_stat('FATAL','FIXDRV','bmake','atml.list'
     .       ,'Atm loading list file requested but not available',0)
          atmll = .true.
        endif
        if( .not.fcheck('atml.grid') ) then   
          if( .not.atmll ) then  
          call report_stat('FATAL','FIXDRV','bmake',' '
     .           ,'Atm loading requested but no list or grid file',0)
          else
            call report_stat('WARNING','FIXDRV','bmake',' '
     .      ,'No atm loading grid file--list file must be complete',0)
          endif
        else
          atmlg = .true.
        endif   
      endif
c     ------met data and mapping functions
      metl = .false.
      metg = .false.    
      mapl = .false.
      mapg = .false.
      call rdsest(14,'met obs source',20,metsource,lsess,reqd,ill )    
c       Tokens are hierrarchical, with GPT or STP followed by a humidity 
c       value; e.g. 'RNX ufile GPT 50'.  Searched here to for files 
c       needed for grdtab call;  passed to mdmake to write lines for model.
      if( metsource(1:3).eq.'   ') then
        metsource = 'GPT 50.             '
      endif          
      if( metsource(1:3).eq.'UFL' ) then
        if( fcheck('map.grid') ) mapg = .true.
        if( fcheck('map.list') ) mapl = .true.
      endif
      call rdsest(12,'Use met.list',1,ans,lsess,reqd,ill )
      if( ans.eq.'Y' ) then
        if( .not.fcheck('met.list') ) 
     .    call report_stat('FATAL','FIXDRV','bmake','met.list'
     .        ,'Met list file requested but not available',0)
c       metl = .true.
      endif        
      call rdsest(12,'Use met.grid',1,ans,lsess,reqd,ill )
      if( ans.eq.'Y' ) then
        if( .not.fcheck('met.grid') ) 
     .    call report_stat('FATAL','FIXDRV','bmake','met.grid'
     .        ,'Met grid file requested but not available',0)
c       metg = .true.
      endif        
      call rdsest(4,'dmap',4,dmap,lsess,reqd,ill ) 
      call rdsest(4,'wmap',4,wmap,lsess,reqd,ill )   
      if( dmap.eq.'IMFH'.or.dmap.eq.'VMFH'.or.dmap.eq.'VMF1' .or.
     .    wmap.eq.'IMFW'.or.wmap.eq.'VMFW'.or.wmap.eq.'VMF1' ) then
        call rdsest(12,'Use map.list',1,ans,lsess,reqd,ill )
        if( ans.eq.'Y' ) then
          if( .not.fcheck('map.list') ) 
     .      call report_stat('FATAL','FIXDRV','bmake','map.list'
     .      ,'Mapping function list file requested but not available',0)
c         mapl = .true.
        endif
        if( .not.fcheck('map.grid')) then   
          if( .not.mapl ) then 
            write(message,'(2a)') 'Site-dependent mapping function '
     .          ,'requested but no list or grid file'
            call report_stat('FATAL','FIXDRV','bmake',' ',message,0)
          else
            call report_stat('WARNING','FIXDRV','bmake',' '
     .   ,'No mapping function grid file--list file must be complete',0)
          endif
        else
          mapg = .true.
        endif
      endif       
c     if u-file requested for met data, must have a met or map grid or list file
      if( mchkey(metsource,'ufile',20,5).gt.0 .and. 
     .    .not.metl .and. .not.metg .and. .not.mapl .and. .not.mapg  )
     .   call report_stat('WARNING','FIXDRV','bmake',' '
     . ,'U-file requested for met data but no met or map file available'
     . , 0)
                          

c    Decimation factor in SOLVE 

      idec = 1    
      call rdsest( 17, 'Decimation Factor', 30,line, lsess, reqd, ill)
      if( ill.ne.0 ) ierr_sestbl = ill
      i0 = index( line, '=' ) + 1
      if( line.ne.blnk256 ) read(line(i0:nblen(line)),*) idec  
c     allow different decimation factor for quick iterations or autcln pre-fit
      idecp = idec 
      call rdsest( 27, 'Quick-pre Decimation Factor', 30,line, lsess
     .              , reqd, ill )  
      if( ill.ne.0 ) ierr_sestbl = ill
      i0 = index( line, '=' ) + 1
      if( line.ne.blnk256 ) read(line(i0:nblen(line)),*) idecp    


c     Set SOLVE file names and controls for initial and final solutions

       qfile = mfile
       qfile(1:1) = 'q'  
       qfile(6:6) = 'a'
       qfilep = qfile
       qfilep(6:6) = 'p'   
       if( typana.eq.'0-ITE' ) then
         qfilep = qfile
         idecp = idec
         observp = observ
       endif       
c** comment this temporarily to allow 0-iter to work for Tom
c       if( observp(1:6).eq.'LC_AUT' ) 
c     .    call report_stat('FATAL','FIXDRV','bmake',' '
c     .        ,'LC_AUTCLN not allowed for pre-fit solution',0)


c===============================================================================


c    **Generate the first ARC batch file**
                         
      if( first_arc.eq.'Y' ) then
             
         write(17,'(a,/,a)') '#','# Initial orbital integration '
         BFIL2(KCOL+1:KCOL+3) = CHRNUM(IDRV)

c        give the X-file limits priority in setting the start and stop times
         if( simulation ) then  
           tspan = dfloat(nepoch*inter)     
           call monday(doy,month,iday,year)
           itb(3) = year
           itb(1) = month
           itb(2) = iday 
           itb(3) = year
           tbb(1) = dfloat(ihr)
           tbb(2) = dfloat(imn)
           tbb(3) = sec   
           itflag = -4                       
           call timcon( itflag,iwkn,sow,year,doy,ihr,imn,sec,utcoff)   
           call secsum(iwkn,sow,tspan,iwknstp,sowstp )                
           itflag = 4
           call timcon( itflag,iwknstp,sowstp,iystp,idoystp,ihr,imn,sec
     .                , utcoff)
           call monday(idoystp,month,iday,iystp)  
           itstp(1) = month
           itstp(2) = iday   
           itstp(3) = iystp
           tstp(1) = dfloat(ihr)
           tstp(2) = dfloat(imn)
           tstp(3) = sec           
         else
           if( fcheck(xfile(1)) )  then
             do i=1,3
              itb(i) = ixb(i)
              tbb(i) = txb(i)
              itstp(i) = ixstp(i)
              tstp(i) =  txstp(i)
             enddo
           endif
         endif
c        if no T-file exists, check for the existence of the G-file
         if( .not.fcheck(tfile) ) then
            if( .not.fcheck(gfile) ) then
              write(iscrn,'(2a,/,2a)')
     .              ' Could not find G-file: ', gfile
     .             ,' Enter name of the G-file associated with: '
     .             , tfile
               read ( 5, '(a16)' )  gfile
               if( .not.fcheck(gfile) ) then
                  write(iscrn,'(2a)') gfile, 'does not exist : Stop'
                  stop
               endif
            endif
         endif

c        if a T-file exists, check for orbit partials if needed
         if( fcheck(tfile) ) then
           if( (nintrs.eq.3 .and. lexp.eq.2)  .or.
     .         (nintrs.eq.3 .and. lexp.eq.3) ) then
               call report_stat('WARNING','FIXDRV','bmake',' '
     .          , 'No orbit partials in T-file',0)
           endif
         endif    

         call armake( bfil2, tfile, gfile, ntsat, itsat, satnam
     .              , itb, tbb, itstp, tstp, lsess, ierr_sestbl
     .              , lexp, expo , norb )
         idrv = idrv + 1

      endif
          


c=================================================================================

c    ** Conversion of ARC-generated y-file to an attitude file read by MODEL
                                                        
                                    
      if( yawmod.eq.'Y' ) then
        if( inter.eq.0 ) call report_stat('FATAL','FIXDRV','bmake',' '
     .    ,'Yaw-table interval is zero',0)
        write(17,'(a,/,a)') '#','# Generation of yaw file ' 
        write(17,'(a,2(1x,a16),i5)')'yawtab',tfile,yfile,inter
      endif

c=================================================================================

c    ** Repeat integration using the yaw file to determine the spacecraft attitude
c       if needed for the radiation-pressure models (possibly make this automatic
c       in the future)

c** here read the sestbl for the radiation models to see if the repeat arc is needed
       
      reqd = .false.
      call  rdsest( 23, 'Radiation Model for ARC', 5, buf5, lsess,
     .              reqd, ill )
      if( ill.ne.0 ) ierr_sestbl = ill
c     buffer will be blank if nothing found
      if( buf5(1:3).eq.'UCL') then
         write(17,'(a,/,a)') 
     .       '#','# Repeat orbital integration using the actual yaw ' 
           call armake( bfil2, tfile, gfile, ntsat, itsat, satnam
     .              , itb, tbb, itstp, tstp, lsess, ierr_sestbl
     .              , lexp, expo , norb )
         idrv = idrv + 1
      endif

                                      
c==================================================================================

c     Invoke GRDTAB to create the u-file for ocean loading, atm loading, met values, 
c     mapping function from external files (station list or grid)
              
      ngfiles = 0
      if( otll ) then
        ngfiles = ngfiles + 1
        grdfiles(ngfiles) = 'otl.list'
      endif
      if( otlg ) then
        ngfiles = ngfiles + 1
        grdfiles(ngfiles) = 'otl.grid'
      endif
      if( atmll ) then
        ngfiles = ngfiles + 1
        grdfiles(ngfiles) = 'atml.list'
      endif
      if( atmlg ) then
        ngfiles = ngfiles + 1
        grdfiles(ngfiles) = 'atml.grid'
      endif
      if( atll ) then
        ngfiles = ngfiles + 1
        grdfiles(ngfiles) = 'atl.list'
      endif
      if( atlg ) then
        ngfiles = ngfiles + 1
        grdfiles(ngfiles) = 'atl.grid'
      endif
      if( metl ) then
        ngfiles = ngfiles + 1
        grdfiles(ngfiles) = 'met.list'
      endif
      if( metg ) then
        ngfiles = ngfiles + 1
        grdfiles(ngfiles) = 'met.grid'
      endif
      if( mapl ) then
        ngfiles = ngfiles + 1
        grdfiles(ngfiles) = 'map.list'
      endif
      if( mapg ) then
        ngfiles = ngfiles + 1
        grdfiles(ngfiles) = 'map.grid'
      endif 
c      print *,'ngfiles ',ngfiles
      if( ngfiles.gt.0) then 
c       need to cover from start of day through one 6-hr pt beyond the end of the span
        spanday = tbb(1)/24.d0 + tbb(2)/1440.d0 
     .            + tbb(3)/86400.d0 
     .            + (nepoch*inter)/86400.d0 + 0.2510  
c       now round to the nearest 0.25 day
cd        print *,'nepoch inter tbb spanday ',nepoch,inter,tbb,spanday  
        spanday = dint(4.d0*spanday)/4.d0
        write(17,'(a,/,a)') '#','# Generation of u-file ' 
        write(17,'(1x,a,1x,a,1x,2i4,f5.2,1x,10a10)') 
     .     'grdtab',dfile,year,doy,spanday,(grdfiles(i),i=1,ngfiles)
      endif
       

c===================================================================================


c     ** Write out batch files for the initial (possibly final) solution **
              
       write(17,'(a,/,a)') '#','# Initial solution '

c       Uncompress X-files if requested and not SOLVE-only
        reqd = .false.
        call rdsest(10,'X-compress',1,xcompress,lsess,reqd,ill )
        if ( ill.ne.0 ) ierr_sestbl = ill
        if( lowerc(xcompress).eq.'y' ) then
           xtemplate = 'x*.'//cdoy//'.Z'   
           write( 17, '(a,a12)' )  'uncompress ',xtemplate
        endif  
        
C       Station loop for MODEL

        do  jstat = 1, nstat

c         Check each X-file for format type and agreement with SV list
c         (if X-file not available, guess format from station.info)  
          if( xorc.eq.'x' )  then
            call xcheck( xfile(jstat), lxfil, lstnfo, exprmt
     .          , nsat, ixsat,  year, doy, skd )
          endif

c         Set the output C-file name

c          The 6th character for the C-file is always the year
c          unless the input is a C-file (allowed only for 0-ITER
c          with clean data), in which case the 6th character is
c          incremented

          cfile(jstat) = xfile(jstat)
          cfile(jstat)(1:1)='c'
          cfile(jstat)(6:6) = cyr(2:2) 
          if( xfile(jstat)(1:1).eq.'c' .and.
     .        typana.eq.'0-ITE' .and. rawcln.eq.'CLN' ) 
     .       call upnam1( cfile(jstat),cfile(jstat) )   
                     
C         Write the MODEL batch file
                  
            bfil2(kcol+1:kcol+3) = chrnum(idrv)
             call MDMAKE( 0, bfil2, tfile, ifile, lfile, jfile
     .                  , yfile, ffile, xfile(jstat), cfile(jstat)
     .                  , klock(jstat), lsite, idatum, ietide, year,lexp
     .                  , delmod, etidemod, ionsrc, magfield, gnss
     .                  , lsess, sfile, ierr_sestbl, ierr_sittbl
     .                  , atmlflag, scratch_dir, metsource )

            idrv = idrv + 1

c       End of station loop
        enddo


C          AUTCLN command line

        if ( rawcln.eq.'RAW' ) then

c         save C-file version for post-fit
          autinchr= cfile(1)(6:6)
          pre_post = 'PRE '
	  call uppers(pre_post)
          call ACMAKE ( nstat,dfile,cfile,autcln_cmd,delaut
     .                , typana,avedit_opt,eclopt,pre_post
     .                , autcln_clk,mfile,ierr_sestbl ) 
          autoutchr= cfile(1)(6:6)
          idrv = idrv + 1

        endif

c      --- Exit here if no solution
      if(typana .eq. 'PREFI') goto 900


c          CFMRG batch file

C        Change M-file 6th character from '1' to 'a' for regular solution
         CALL  UPNAM2( MFILE, MFILE )
         BFIL2(KCOL+1:KCOL+3) = CHRNUM(IDRV)
C        WRITE( lbfil, '(A)' ) 'if existf abort  then dlf abort'
c        Use D-file inputs if SOLVE-only solution (can now be x or c files)
         if( rawcln.eq.'CLN' ) then
               do i = 1, nstat
                  cfile(i)= xfile(i)
                  cfile(i)(1:1)= 'c'
               enddo
         endif        
         call CMMAKE( nstat, cfile, mfile, snames, bfil2
     .               , exprmt, nsat, ixsat, lsite, lsess 
     .               , nzen, gradest, ngrad,ierr_sestbl, ierr_sittbl )  
         idrv = idrv + 1


C           SOLVE batch file

         bfil2(kcol+1:kcol+3) = chrnum(idrv)
C        write( lbfil, '(A)' ) 'if exist abort  return'
c        no noise file until final solution
         nfile = '          ' 
              
         call SOMAKE(mfile,qfilep,lfile,gfile,nfile,xfile,bfil2   
     .              ,lsess,lsite,lsesfo,lxfil,lcfil
     .              ,xorc,nstat,snames,nsat,norb,observp,idecp
     .              ,zenest,nzen,zenmod,gradest,ngrad,gradmod
     .              ,year,doy,idatum,proc,ierr_sestbl,ierr_sittbl
     .              ,fixdrv_vers, scratch_dir )     
         idrv = idrv + 1                       
         write(17,'(a)') 'sh_chksolve '
         write(17,'(a)') 'if( -e GAMIT.fatal ) exit '

    
c------------------------------------------------------------------------------     

c        Repeat solution with updated coordinates for AUTCLN using post-fit residuals
         
         if( autcln_post.eq.'Y' ) then
                       
           write(17,'(a,/,a)') '#','# Post-fit editing and solution'
c          if input AUTCLN C-files are not kept, erase the output and recreate the input
           if( delaut.eq.'Y' ) then 
              write(17,'(3a)') "/bin/rm c????",autoutchr,".???"
c             update L-file name
              if( tlupdt.eq.'T_AND_L' .or. tlupdt.eq.'L_ONLY ' ) then
                 nchr = nblen(lfile)
                 if( lfile(nchr-2:nchr).eq.'apr' ) then
                   call upnam3( lfile, lfile )
                 else
                   call upnam1( lfile,lfile )
                 endif
              endif     

c             Station loop for MODEL 

              do jstat = 1, nstat
                bfil2(kcol+1:kcol+3) = chrnum(idrv) 
                call upnam1( cfile(jstat),cfile2 )
                if( iterate.eq.'XFILES' ) then
                   cfile1 = xfile(jstat)
                else
                   call report_stat('FATAL','FIXDRV','bmake',' '
     . ,'Iterate = C incompatible with Delete AUTCLN input C-files = Y'
     .                ,0)
                endif  
                call  MDMAKE( 1, bfil2, tfile, ifile, lfile, jfile
     .                      , yfile, ffile, cfile1, cfile2 
     .                      , klock(jstat),lsite,idatum,ietide,year,lexp
     .                      , delmod, etidemod, ionsrc, magfield, gnss
     .                      , lsess, sfile,ierr_sestbl, ierr_sittbl  
     .                      , atmlflag, scratch_dir, metsource )
                cfile(jstat) = cfile2  
                idrv = idrv + 1
              enddo  

c          otherwise, move the last set of AUTCLN input and erase the last set of output 

           else 
              do jstat = 1, nstat
                call upnam1(cfile(jstat),cfile2)
                cfile(jstat) = cfile2
              enddo
              outchr = cfile2(6:6)  
              write(17,'(4a)') 'mvcf ',autinchr,' ',outchr 
              if( delaut.ne.'N') then
                write(17,'(3a)') "/bin/rm c????",autoutchr,".???"
              endif
           endif 
           pre_post = 'POST' 
	   call uppers(pre_post)
           autinchr = cfile(1)(6:6)
           call ACMAKE( nstat,dfile,cfile,autcln_cmd
     .                , delaut,typana,avedit_opt,eclopt,pre_post
     .                , autcln_clk,mfile,ierr_sestbl ) 
           autoutchr= cfile(1)(6:6)
           bfil2(kcol+1:kcol+3) = chrnum(idrv)  
                      autoutchr= cfile(1)(6:6)    
           call CMMAKE( nstat, cfile, mfile, snames, bfil2  
     .                , exprmt, nsat, ixsat, lsite, lsess 
     .                , nzen, gradest, ngrad,ierr_sestbl, ierr_sittbl )
           idrv = idrv + 1   

c  RWK 070117: Now always use an N-file for the second solution since we will 
c              always use it to loosen bad contraints; optionally use for 
c              AUTCLN WL biases and/or reweighting by station and possibly satellites.     
c          --check if WL biases from autcln
           if( observ(1:6).eq.'LC_AUT') then
             autwl = 'y'   
           else
             autwl = 'n'
           endif
c          --check if noise reweighting    
           reqd = .false.
           call rdsest( 15,'AUTCLN reweight',1,rewgt,lsess,reqd,ill )   
c          for backward compatibility, old command is proxy for the new      
           if( rewgt.eq.' ' ) then
             call rdsest( 10, 'Use N-file', 1 ,ans, lsess, reqd,ill )
             write(message,'(2a)') 'Use N-file is obsolete sestbl entry'
     .            ,'; ok now but use AUTCLN reweight' 
             call report_stat('WARNING','FIXDRV','bmake',' ',message,0)
             rewgt = ans
           endif
           if( rewgt.eq.' ' ) rewgt = 'N'
           call lowers(rewgt)  
           call rdsest( 13, 'Station Error',4,stnerr,lsess,reqd,ill)    
           call rdsest( 15, 'Satellite Error',4,saterr,lsess,reqd,ill) 
           call lowers(stnerr)
           call lowers(saterr)                                         
           if( rewgt.eq.'y' .and. stnerr.ne.'elev' )    
     .         call report_stat('FATAL','FIXDRV','bmake',' '
     .      ,'Station error model must be ELEVATION for reweight option'
     .         ,0)
           if( ill.ne.0 ) ierr_sestbl = ill 
c          --construct the sigelv command line
           nfile = dfile  
           nfile(1:1) = 'n'
           nfile(6:6) = cfile(1)(6:6)    
           qfileopt = ' -qfile '//qfilep  
           autclnopt = ' '
           if( autwl.eq.'y' .or. rewgt.ne.'n' ) 
     .        autclnopt = ' -acmd autcln.cmd.postfit'
           stnerropt = ' '
           saterropt = ' '    
           if (rewgt.eq.'n' ) then  
              stnerropt = ' -noelv '  
            else
              if( rewgt.eq.'s' ) saterropt = ' -sv ' 
           endif
           write(17,'(a,a10,4a)')  'sh_sigelv -nfile ',nfile
     .                  ,qfileopt,saterropt,stnerropt,autclnopt

           bfil2(kcol+1:kcol+3) = chrnum(idrv)
           call SOMAKE(mfile,qfile,lfile,gfile,nfile,xfile,bfil2 
     .              ,lsess,lsite,lsesfo,lxfil,lcfil
     .              ,xorc,nstat,snames,nsat,norb,observ,idec
     .              ,zenest,nzen,zenmod,gradest,ngrad,gradmod
     .              ,year,doy,idatum,proc,ierr_sestbl,ierr_sittbl
     .              ,fixdrv_vers,scratch_dir )
           idrv = idrv + 1     
           write(17,'(a)') 'sh_chksolve'
           write(17,'(a)') 'if( -e GAMIT.fatal ) exit '

         endif 

c----------------------------------------------------------------------------

c         Optional second AUTCLN run to get clock estimates or residuals for 
c         empirical antenna modeling

         if( autcln_clk.eq.'Y') then
	   write(17,'(a,/,a)') '#','# Extra autcln for clocks'
c
	   pre_post = 'CLK ' 
	   call uppers(pre_post)
c          use the model-created c-files 
	   do i=1,nstat
	     cfile(i)(6:6) = autinchr
           enddo
           call ACMAKE( nstat,dfile,cfile,autcln_cmd
     .                , delaut,typana,avedit_opt,eclopt,pre_post
     .                , autcln_clk,mfile,ierr_sestbl ) 
           autoutchr= cfile(1)(6:6)
	 endif
  
c----------------------------------------------------------------------------

c         Optional second post-fit solution using post-fit M-file
c         if the adjustments high enough to have compromised cleaning
         if( rawcln.ne.'CLN' .and. autcln_redo.eq.'Y' ) then
               
           write(17,'(a,/,2a)') '#','# Re-do AUTCLN and SOLVE with '
     .                                              ,'updated M-file'
           write(17,'(2a)') '# if the ratio of the pre-fit to '
     .                         ,'the post-fit nrms exceeds 1.5' 
           qfileopt = ' -qfile '//qfile  
           write(17,310) dquote,dquote,qfilep
  310      format
     .      (17hset prms = `grep ,a1,12hPostfit nrms,a1,1x,a10,
     .       30h | head -1 | awk '{print $6}'` )    
           write(17,311) dquote,dquote,qfile
  311      format
     .      (17hset arms = `grep ,a1,12hPostfit nrms,a1,1x,a10,
     .       30h | head -1 | awk '{print $6}'` )
           write(17,312) dquote,dquote
  312      format
     .      (55hset redo = `echo $prms $arms | awk '$1/$2 > 1.5 {print 
     .      ,a1,3hyes,a1,8h;exit}'` ) 
           write(17,313) dquote,dquote
  313      format(13hif( $redo == ,a1,3hyes,a1,7h ) then )
c          if input AUTCLN C-files are not kept, erase the output and recreate the input
c          (necessary here with both 'Delete' and 'Intermediate' options for AUTCLN input C-files)
           if( delaut.eq.'Y' .or. delaut.eq.'I' ) then 
              write(17,'(3a)') "/bin/rm c????",autoutchr,".???"
              do jstat = 1, nstat
                bfil2(kcol+1:kcol+3) = chrnum(idrv) 
                call upnam1( cfile(jstat),cfile2 )
                if( iterate.eq.'XFILES' ) then
                   cfile1 = xfile(jstat)
                else
                   call report_stat('FATAL','FIXDRV','bmake',' '
     . ,'Iterate = C incompatible with Delete AUTCLN input C-files = Y'
     .                ,0)
                endif  
                call  MDMAKE( 1, bfil2, tfile, ifile, lfile, jfile
     .                      , yfile, ffile, cfile1, cfile2 
     .                      , klock(jstat),lsite,idatum,ietide,year,lexp
     .                      , delmod, etidemod, ionsrc, magfield, gnss
     .                      , lsess, sfile, ierr_sestbl, ierr_sittbl  
     .                      , atmlflag, scratch_dir, metsource )
                cfile(jstat) = cfile2  
                idrv = idrv + 1
              enddo
c          otherwise, move the last set of AUTCLN input and erase the last set of output 
           else   
               do jstat = 1, nstat
                call upnam1(cfile(jstat),cfile2)
                cfile(jstat) = cfile2
              enddo
              outchr = cfile2(6:6)  
              write(17,'(a)') 'if( -e GAMIT.fatal ) exit '
              write(17,'(4a)') 'mvcf ',autinchr,' ',outchr 
              if( delaut.ne.'N') then
                write(17,'(3a)') "/bin/rm c????",autoutchr,".???"
              endif
           endif
           pre_post = 'POST'
	   call uppers(pre_post) 
	   autinchr = cfile(1)(6:6)  
           call ACMAKE( nstat,dfile,cfile,autcln_cmd
     .                , delaut,typana,avedit_opt,eclopt,pre_post
     .                , autcln_clk,mfile,ierr_sestbl )
           bfil2(kcol+1:kcol+3) = chrnum(idrv) 
           call CMMAKE( nstat, cfile, mfile, snames, bfil2   
     .                 , exprmt, nsat, ixsat, lsite, lsess 
     .                 , nzen,gradest,ngrad,ierr_sestbl, ierr_sittbl )
           idrv = idrv + 1     
c          create an N-file for constraint-checking and optionally lc_autcln and/or reweighting  
           nfile(6:6) = cfile(1)(6:6)    
           qfileopt = ' -qfile '//qfile 
           write(17,'(a,a10,4a)')  'sh_sigelv -nfile ',nfile 
     .                   ,qfileopt,saterropt,stnerropt,autclnopt
           bfil2(kcol+1:kcol+3) = chrnum(idrv) 
           call SOMAKE(mfile,qfile,lfile,gfile,nfile,xfile,bfil2 
     .              ,lsess,lsite,lsesfo,lxfil,lcfil
     .              ,xorc,nstat,snames,nsat,norb,observ,idec
     .              ,zenest,nzen,zenmod,gradest,ngrad,gradmod
     .              ,year,doy,idatum,proc,ierr_sestbl,ierr_sittbl
     .              ,fixdrv_vers, scratch_dir )
           idrv = idrv + 1                     
           write(17,'(a)') 'sh_chksolve '
           write(17,'(a)') 'if( -e GAMIT.fatal ) exit '
           write(17,'(a)') 'endif'
	 
c          Optional second AUTCLN run to get clock estimates or residuals for 
c          empirical antenna modeling

           if( autcln_clk.eq.'Y') then
             write(17,'(a,/,a)') '#','# Extra autcln for clocks'
	     pre_post = 'CLK ' 
	     call uppers(pre_post)
c            use the model-created c-files 
	     do i=1,nstat
	       cfile(i)(6:6) = autinchr
             enddo
             autinchr = cfile(1)(6:6)
             call ACMAKE( nstat,dfile,cfile,autcln_cmd
     .                  , delaut,typana,avedit_opt,eclopt,pre_post
     .                  , autcln_clk,mfile,ierr_sestbl ) 
             autoutchr= cfile(1)(6:6)
	   endif
	 endif
                   
c------------------------------------------------------------------------------------------

c        Final check on quality of full solution to determine if SCANDD should be run
                       
           if( .not.simulation ) then
             write(17,'(a,/,a)') 
     .          '#','# Check the quality of the final solution'
             write(17,320) dquote,dquote,qfile
  320        format
     .      (16hset rms = `grep ,a1,12hPostfit nrms,a1,1x,a10,
     .       35h | awk '$6 > 1 {print "no";exit}'` )
             write(17,'(a)')
     .         'if ( $rms == "no" ) echo Full solution rms is too high'
             if( rawcln.eq.'RAW' ) then
               if(scan_control.eq.'NONE '.or.scan_control.eq.'FIRST') 
     .             then
                 write(17,'(a,a16)') '#scandd ',mfile
              elseif(scan_control.eq.'BOTH '.or.scan_control.eq.'FULL ')
     .           then   
                 write(17,'(a,a16)') 'scandd ',mfile
               else
                 write(17,'(a,a16)') 'if ( $rms == "no" ) scandd ',mfile  
               endif
             endif     
             write(17,'(a)') 'if ( $rms == "no" ) exit'
           endif

                                              
c-----------------------------------------------------------------------

c           CTOX run to save clean X-files

         reqd = .false.
         call rdsest(8, 'Run CTOX', 1, run_ctox, lsess, reqd, ill )
         if ( ill.ne.0 ) ierr_sestbl = ill
         if ( run_ctox.eq.' ' ) run_ctox = 'N'
         if ( run_ctox.eq.'Y') then
            write(17,'(a,a16,1x,a1,1x,a1)')
     .                'ctox ',dfile,cfile(1)(6:6),'a'
         endif


c----------------------------------------------------------------------------

c           One more ARC run to get best T-file

c         Default set to 'N' at beginning of BMAKE, no entry required
          reqd = .false.
          if( exprmt.ne.'BASEL' .and. exprmt.ne.'KINEM' .and.
     .        finarc.eq.'Y' )  then
               call  rdsest( 13, 'Export Orbits', 1, expo, lsess,
     .              reqd, ill )
               if ( ill.ne.0 ) ierr_sestbl = ill
               if( expo.eq.' ') expo = 'N'
               call  upnam1( gfile, gfile )
               call  upnam1( tfile, tfile )
               bfil2(kcol+1:kcol+3) = chrnum(idrv)
               call  ARMAKE( BFIL2, TFILE, GFILE, ntsat, itsat, satnam
     .                     , ITB, TBB, ITSTP, TSTP, LSESS, ierr_sestbl
     .                     , lexp, EXPO, norb )
               idrv = idrv + 1
C              Convert T-file to NGS format
               if(expo.eq.'Y') then
                 bfil2(kcol+1:kcol+3) = chrnum(idrv)
                 call  NGMAKE( bfil2,tfile,cyr,cdoy,lsess
     .                       , ierr_sestbl )
                 idrv = idrv + 1
               endif
          endif


c============================================================================


C     Close the batch files and print the run message

      if( lowerc(xcompress).eq.'y' ) then
        xtemplate = 'x*.'//cdoy//'  '
        write( 17, '(a,a12)' )  'compress ',xtemplate
      endif
  900 close( lbfil ) 
      write(message,'(a,a16)') 'Created GAMIT batch file ',bfile
      call report_stat('STATUS','FIXDRV','bmake',' ',message,0)
      write(6,901) bfile,bfile,bfile
  901 format(/,
     .     ' To run this job in the foreground enter: csh ', a16,/,
     .     ' To run his job and minimize output enter: csh ',
     .          a16,' > /dev/null', /,
     .     ' To run this job in the background enter: gbat ',a16,/)

      RETURN
      END

