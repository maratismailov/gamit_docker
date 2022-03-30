       Program MAKEXP
C     ---------------
C
C  Purpose:    Create input files for programs
C              HI, MAKEJ, MAKEX, BCTOT, AND FIXDRV
c
c
C  Written by Werner Gurtner and Yehuda Bock  5 November 1990

C  Calls  :    RSTNFO, GETDIR, and NBLEN from GAMIT Library


C  version 1.0
C  Changes:    Modify for GAMIT release 8.2
C              Yehuda Bock 910710
C  version 2.0 Modified for GAMIT release 9.1 (kinematic option)
C              Implement master STATION.INFO
C              Yehuda Bock 10/17/91
C  version 2.1 Add option to read in X-files instead of RINEX files
C              (in this case bypass MAKEX)
C              Yehuda Bock 11/06/91
C         2.11 Modify HI.RAW format YB 11/08/91
C         2.12 Allow station.info to have antenna history YB 12/22/91
C              Check year and day properly YB & MB (Mike Bevis) 01/25/92
c         2.13 Use RSTNFO to read station.info file.
c              Create a scenario file if it doesn't exist.
c              Add session number and antenna type to input.
c              Add dimpar.fti
c                 King 92/2/10
C         2.14 Change "s" to "i" file YB 02/27/92
c              Modify initial output to screen YB 03/02/92
c         9.1  Copied from MIT, modified Makefile YB 05/05/92
c         9.11 Add L-file to D-file   King 6/4/92
c              Changed maxfil to 1000 Bock 6/6/92
c              Fixed infinite loop problem Bock 6/10/92
c         9.12 Add calling argument in CHECKE for library READE.
c              Translate TI400 code into GES or COR for MAKEX batch file.
c                King 6/24/92
c         9.13 Minor change in E-file message.   King 92/9/5
c              Change calling arguments for READE to match lib change. King 92/9/21
c         9.14 Add Kurt's user friendly suggestion. Bock 93/03/15
c         9.15 Add WM-102 ==> WM2 translation for MAKEX batch file.  King 93/4/5
c         9.16 Increase formats for 32 satellites.   King 93/10/18
c         9.17 Output to screen, increased format to 32 satellites.  Bock/Fang 94/01/03
c              Do not update session.info if more than one entry.  Bock 94/01/04
c         9.18 Check for either 0 or 1 for RINEX session number if either input.
c                King 94/03/26
c              Change HP compiler flags from +e +E1 to +U77.  Fang/King 94/05/06
c              Add screen message when searching for RINEX files.  King 94/05/06
c         9.19 Add calling argument for library RSTNFO.   King 94/05/11
c              Makefile:  Shift ranlib to execute only once.  Herring/King 94/06/24
c         9.20 Read sestbl. to get the inertial frame for BCTOT.  King/McClusky 95/10/03
c         9.21 Makefile:  variable ranlib for compatibility with Solaris 2.  Fang/King 951208
c         9.22 Add calls to report_stat to MAKEXP and CHECKE.   McClusky 960223
c         9.23 Change structure to compute site lists from rinex/xfiles and then call
c               station.info using icall=3   Tregoning/King 960508
c              Fix bug in writing out *.makex.batch file  Tregoning 960528
c         9.24 Use 'MAKEXP' rather than 'GAMIT' for the module name for fatal, warning,
c                and status files (matching change in /lib/report_stat).  King 960610
c              Allow running without navigation file if session.info exists; add
c                comment of sp3_to_gt as alternative to BCTOT.  King 960611
c              Set session times for rstnfo call to bypass check (not needed here
c                and awkward because session.info not yet read.  King 960614
c         9.25 Remove extraneous warning message Bock/King 960610/960710
c              Add 'checke' message Bock/King 960623/960710; change to 'checkj' King 960805
c              Remove extraneous comma.  King 960808
c              Change name of SP3 script and add copy g-file option in message displayed
c                  for user. King 960825
c         9.26 Fix typo in declaring 'rinex1' leading to abort with session=0. King 961016
c              Put unit number in write at 999 (caught by gcc compiler.  Tregoning/King 961017/18
c         9.27 CHECKE:  Change calling arguments for READE for orbits/BCCHECK additions. King 961107
c              Change 'rcvers' => 'swvers' to match lib/rstnfo.  King 970111
c         9.28 G77 compiler changes:
c              MAKEXP: Initialised previously uninitialised variables, removed unused variables,
c              and explicitly typed all unexplicitly typed variables. McClusky 970113
c         9.29 Update instruction messages at end to use sh_check_sess.  King 970129
c              Set TR8000 and TR8100 to TRB as receiver software type.  King 970219
c         9.30 Allow search on multiple RINEX files per day, removing duplicates.  King 970306
c              Allow search on multiple days.  King 970324
c              Reorder printed instructions; remove sort_string (now in /lib).  King 970326  
c              Define print unit for lversn and proper calls.  Calais/King 970414
c         9.31 Set SR299 or SR399 to LEI as receiver software type.  King 970717 
c         9.32 Add yr and day-of-year to first-line input for BCTOT, to allow use of session.info
c                 rather than an X-file.  King 971028
c              Use session.info rather than X-file to get session span for sh_bctot.  King 971208
c              Add SR9500 to translation for Leica firmare ids.  King 980105   
c         9.33 Fix bug in comparing nrsess with 99 (was .99!).  Herring 980205
c         9.34 Translate all SR and CR rcvr codes to rcver LEI.  King 980902  
c         9.35 Fix format statement in printout of satellites available.  King 981112  
c         9.36 Change instructions for running FIXDRV.  King 990201
c         9.37 Change instructions for running makej and makex.  King 990324
c              No longer create input files for makej, makex, and fixdrv.  King 990504   
c         9.38 Add 'LC' to translation for Leica firmware id.  King 990529
c         9.49 Changes for 4-digit years: new formats for MAKEX batch file and session.info,
c               allow 2- or 4-digit inputs.  King 990802
c              Change default inertial frame to J2000 (bug--changed in FIXDRV in 1996). King 990816
c              Change session.info to free-format with comments.  McClusky/King 990816 
c              Add 'AO' to translation for AOA firmware id's.  Bock 990820
c              Change name maxfil->maxrxfil to avoid conflict with maxfil in makex.h,
c                now used to get maxepc for check in rsesfo.  Herring/King 990824 
c              Add translation 'ROG_12'-> firmware 'TRB'.   King 990826
c              Add translation 'ROG800'-> firmware 'TRB'.   King 990830
c         9.50 Fix bug in year digits for bctot input.  Sugimoto/King 991129
c              Fix comment line in MAKEX batch file.  King 991201
c         9.51 Line-length for session.info buffer too short for 28 SVs.  Fang 000104
c         9.52 Add translation 'TR4' and 'TR7' --> firmware 'TRM'.  King 010131  
c         9.60 Allow use of command-line arguments, give 'nav' file priority over
c               'e' file, and read new station.info. King 010516   
c              Add x-file version to command-line arguments.  King010710 
c              Fix bug in getting nav-file name in interactive mode.  King 010719  
c         9.61 Rearrange order so that j-file always gets named by nav-file name, not
c               project name; fix bug in getting x-file 6th character when run interactively.  King 010809 
c              Fix problem with naming j-file when nav file is e????y.ddd.  King 010912   
c         9.62 Allow new-style station.info; change 'swver' from R*8 to R*4 to match rest of GAMIT.  King 020311
c              Print session project,orbit,day-search info to makexp.out.  King020318   
c         9.63 Add 'TR5' (for 5700) to alias list for 'TRM' firmware.  King 020531
c         9.64 Fix improper re-initialization of expt and orbt after input.  Herring/McClusky/King 020710 
c         9.65 Read the l- or apr file name from the command line.  King 020806 
c              Remove the feature of updating the day in session.info.  King 020806 
c              Fix bug in 020806 changes.  King 020923   
c         9.66 Change *,* to *,[fmt] for DEC alpha.  Laurain/King 030305    
c              Translate Topcon receiver codes into Ashtech, Trimble, or AOA as appropriate.  King 030313   
c         9.67 Fix bug in calling rstnfo2 (add span).  McClusky/King 030502   
c         9.68 Fix bug in failing to accept interactive session input. King 030527
c         9.69 Fix syntax error in report_stat call for error reading session.info. King 031027
c         9.70 Add 'NAV' for all NavCom receivers.  King 031218
c              Fix problem in interactive mode for default L-file name.  King 031230     
c              Translate TOPSSI to TRMSSI.  King 040131  
c              Initialize 'isod' for new-format station.info.  Elosegui/King 040205        
c         9.71 Reorder creating session.info so that the start/stop times are known by rstnfo. King 051229
c         9.72 Change lib/rstnfo2 so that if no entry is found for a linked RINEX file and the calling program
c              is makexp, only a warning is issued, leaving to makex the detection of a real problem  
c              (no change to code in makexp.f)    
c         9.73 No longer erase the status/warning/fatal files before running; this
c               necessary to accommodate change to write to GAMIT.status/warning/fatal. King 060626
c         9.74 Change arguments for rsesfo (see lib/lversn.f).  King 080815   
c         9.75 Change arguments for rsesfo (see lib/lversn.f).  King 061207
c              CHECKE: Trap # SVs > maxsat on nav file.  King 061226
c              Add trap for nsite > 99, maximum allowed by cfmrg and solve.  King 070104
c         9.76 MAKEXP: Remove print unit number from 'lversn' call. King 070416
c         9.77 MAKEXP: Fix problem when an empty session.info exists.  King 080403
c              MAKEXP: Add translation of rcvcod TRNTRS to swver TRM.  King 080429  
c              MAKEXP: Remove old-style station.info option.  King 080509
c         9.78 MAKEXP: Remove debug output for translating rcvcod.  King 080811   
c         9.79 MAKEXP: Remove call to obsolete check_oldstnfo; remove 'trkcod' from stnfo call.  King 100209
c         9.80 MAKEXP: Add rcvers, rcvrsn, antsn to rstnfo calling arguments. King 100908
c         9.81 MAKEXP: change swver format from 4.2 to 5.2 to accomodate versions > 10.00. King 141007
c         9.82 MAKEXP: Fix long-standing format error trapped by gfortran 4.6 and later. King 141111
c         9.83 MAKEXP: Add gnss option. King 141208 
c         9.84 MAKEXP: Add j-file name to the command-line options.  King 150806 
c              MAKEXP: Fix bug in 150806 change.
c              MAKEXP: Trap a missing j-file argument. King 150821 
c              MAKEXP: Change new argument from j-file name to 'jclock' (brdc or sp3) so that the 
c                      program is self-contained.  King 150824
c         9.85 CHECKE: Match calling arguments for modified lib/reade.  King 151119
c              CHECKE: Match GNSS code as well as PRN. King 151130
c              MAKEXP: Set i-file to NONE for Glonass.  King 151209
c              MAKEXP: Write an sp3-file into the makex batch file if jclock is not brdc;
c                      for Glonass, set i-file to NONE if no sp3 file used.  King 151222
c              CHECKE: Skip record nav-file record if bad.  King 160218
c         9.86 MAKEXP: Remove obsolete 'skd' from call to rstnfo. King 180322 
c              MAKEXP: Replace the obsolete solution number with the GNSS code in the d-file. King 180427
c              MAKEXP: Add command-line argument for the sp3 file since constructing it is dangerous. King 180428 
c        [9.87 MAKEXP: Assign sp3 j-file name when running manually. Zou Xuan 180904 (revoked; Floyd 190121/20200928)]
c         9.88 MAKEXP: Allow passing of lower GNSS frequemcy to over-ride default. Herring 200511.
c         9.89 MAKEXP: Assign sp3 name when running manually. Floyd 20200928 (originally 190121)
c         9.90 MEXEXP: Added reading satellites from sp3 file is jclock=sp3 and no satellites are found
c                      in nav-file. TAH 201215.
c 
      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
               
      integer*4 maxrxfil
      parameter(maxrxfil=1000)

      character*1  xorx,xver,ans,gnss
* MOD TAH 200511: Added gnsslf
      character*2  gnsslf  ! GNSS designation with optional lower frequency 
      character*8  buff    ! Small buffer for readinf runstring (e.g., gnsslf)   
      character*2  kyear  
      character*3  aday
      character*4  savcod(maxrxfil),scodes(maxrxfil),new_codes(maxrxfil)
     .          ,  scodes1(maxrxfil),sitcod,blank4
     .          ,  orbt,project,navf,buf4,orbt1,project1,jclock 
      character*5  buf5,frame,hgtcod(maxrxfil),radome(maxrxfil)    
      character*6  rcvcod(maxrxfil),antcod(maxrxfil)
      character*10 getmac,machin,sestbl,jfile
      character*12 sp3file
      character*16 stanam(maxrxfil)
      character*14 wild
      character*20 rcvers,rcvrsn,antsn
      character*80 navfile,aprfile,makexbatch,dfile
     .           , filcod(maxrxfil),ifile,tfile
      character*45 libver     
      character*50 vers
      character*80 rinex(maxrxfil),rinex1(maxrxfil)*80,xfile(maxrxfil)  
      character*256 message
      character*1024 message1
                                
      real*4   swvers(maxrxfil)
      real*8   anthgt(maxrxfil)
     .       , offstn(maxrxfil),offste(maxrxfil)
      real*8 antdaz(maxrxfil)  ! Antenna aligment from True N (deg).

      integer*4 
c                 unit numbers  
     .          lusest,lustnfo,lusinfo,lunav,lumakexb,luprint  
c                 data file numbers
     .        , nobs_files,nobs_files1  
c                 days and sessions to search for RINEX files
     .        , iday_search,nrsess_search
c                 primary session descriptors
     .        , iyear,iyr2,iday,nrsess,inter,ihr,min,nepoch
     .        , nprn,iprn(maxsat)
c                 for/from station.info
     .        , jstart(5),jstop(5),isod,jsessn,span  
c                 d-file numbers
     .        , nsolu,nsess,nsites
c                 GPS week, day-of-week
     .        , gpsw,gpsd
c                  other integers
     .        , nblen,is,ioerr,iclarg,iarg,ixver
     .        , ill,i,j,k

      logical  batch,fcheck,reqd,nav_found,cmdline_sess
     .       , old_stinf

* MOD TAH 201215: Added to read the SP3 file header.
      integer*4 lusp3     ! Unit number for SP3 file.
     .,         iyrsp3, imosp3, idysp3, ihrsp3, imnsp3 ! YMDHM for sp3 file
     .,         nepsp3    ! Number of epochs
     .,         numsp3sv  ! Total number svs in SP#
     .,         numsp3    ! Number of satellites of GNSS type
     .,         itsat(maxsat)   ! List of PRSN
      real*8    delt      ! Not described in rsp3hd (spacing) 
     .,         secsp3    ! Seconds tag. 
     .,         mjd, fmjd ! MJD and fraction of day                                            
     .,         accsat(maxsat)  ! Accuracy values 
     .,         pvsigb     ! Base for position and velocity
     .,         clksigb    ! Base for clocks

      character*8 otlmod    ! Ocean tide model


               
c Setup defaults and read the input parameters
c ---------------------------------------------

c   Remove old versions of the status, warning, and error files
c*    rwk/scm/tah 060628: No longer clear since now written to GAMIT.status/warning/fatal
c      call report_stat('CLEAR','MAKEXP',' ',' ', ' ',0)
c      call report_stat('CLEAR','LIB',' ',' ', ' ',0)
           
c   Set the unit numbers and initialize literals

      lunav  =  11
      lusest =  12 
      lustnfo = 13
      lusinfo = 14
      luprint=  19
      lumakexb= 20
      lusp3   = 21   ! New unit for readinf SP3 file
      message  = ' '
      message1 = ' '
      wild = ' '   
      navfile = ' ' 

c   Open the MAKEXP print file
       
      open (unit=luprint,file='makexp.out',form='formatted'
     .     ,status='unknown',iostat=ioerr)
      if( ioerr.ne.0 ) call report_stat('FATAL','MAKEXP','makexp'
     .  ,'makexp.out','Error opening MAKEXP print file: ',ioerr)

c   Get and write the program header and version number

      MACHIN = GETMAC(1)
      WRITE (VERS,1) MACHIN(1:nblen(machin))
    1 format (' MAXEXP Ver. 9.90 2020/12/15 22:57 UTC (',a,')')
c     get library version
      CALL LVERSN(libver)
      WRITE(6,'(a)')' '
      WRITE(message,2) vers,libver
    2 FORMAT('Started ',a44,' Library Ver. ',a45)
      call report_stat('STATUS','MAKEXP','makexp',' ',message,0)
      CALL PROPER(luprint) 
      write(luprint,*) message
   
c   Initialize the program and session parameters
                       
      batch = .true.
      cmdline_sess = .false. 
      iday_search = 999              
      iday = 0 
      nrsess_search = 99 
      nrsess = 0
      navf = 'auto'   
      aprfile = ' '

c   Get the session parameters from the command-line if present
   
c     makexp [project] [orbt] [gnsslf] [nav-file] [jclock] [sp3file]  [year] [doy] [xver] [apr-file] [inter] [start hr] [start min] ] [epochs]      
c                1       2      3      4          5         6         7      8    9        10        11       12          13          14
c
c       project   = 4-char project code
c       orbt      = 4-char orbit code              
c       gnsslf    = GNSS system ( G R C E J I ) with option lower frequency e.g., C7 G5 E7 to use
c                   lower frequencies different from default in GAMIT.
c       nav-file  = full navigation file name or e-file name, or possibly 4-character code;
c                   'none' if x-, j-, and k-files available and nav-file not needed
c       jclock    = 4-letter code indicating whether SV clocks and orbits are to be extracted
c                   from an sp3 file or navigation file: 'sp3 ' or 'brdc' [default]
c       sp3file   = sp3-file name (cccwwwd.sp3) 
c       year      = 4-digit year
c       doy       = day-of-year       
c       xver      = 6th character of X-file (default is single-digit year)   
c       apr-file  = name of l- or apr file
c       inter     = sampling interval (seconds)
c       start     = hr min of start   
c       nepochs   = number of epochs
c         These last three can be omitted if session.info exists or the default (30 2880 0 0) is to be used

      iarg = iclarg(1,project)
      if ( iarg.le.0 ) then
        batch = .false.  
      endif
      if( batch ) then
         iarg = iclarg(2,orbt)
         if( iarg.le.0 ) call report_stat('FATAL','MAKEXP','makexp',' '
     .     ,'Missing orbit code in command line ',0) 
* MOD TAH 200511: Read into longer buffer
         iarg = iclarg(3,buff)   
         if( iarg.le.0 ) call report_stat('FATAL','MAKEXP','makexp',' '
     .     ,'Missing GNSS system requested in command line ',0)
* MOD TAH 200511: Read gnsslf
         read(buff,'(a2)')  gnsslf   
         gnss = gnsslf(1:1)  
         if( gnss.eq.' ' ) gnss = 'G'    
         iarg = iclarg(4,navfile)   
         if( iarg.le.0 ) call report_stat('FATAL','MAKEXP','makexp',' '
     .     ,'Missing nav-file name or code in command line ',0)  
         iarg = iclarg(5,jclock)
c        trap a user mistake of omitting the jclock variable (value will be numerical year)
         read(jclock,'(i4)',iostat=ioerr) iyear
         if(ioerr.eq.0 )  call report_stat('FATAL','MAKEXP','makexp',' '
     .            ,'Missing jclock argument, read year instead',ioerr)
         iarg = iclarg(6,sp3file)
         if( iarg.le.0 ) call report_stat('FATAL','MAKEXP','makexp',' '
     .     ,'Missing sp3-file name in command line ',0) 
         iarg = iclarg(7,buf4)
         if( iarg.le.0 ) call report_stat('FATAL','MAKEXP','makexp',' '
     .     ,'Missing year in command line ',0)   
         read(buf4,'(i4)') iyear 
         iyr2=mod(iyear,100) 
         iarg = iclarg(8,buf4)
         if( iarg.le.0 ) call report_stat('FATAL','MAKEXP','makexp',' '
     .     ,'Missing day-of-year in command line ',0)   
         read(buf4,'(i4)')  iday 
         iarg = iclarg(9,buf4)
         if( iarg.le.0 ) call report_stat('FATAL','MAKEXP','makexp',' '
     .     ,'Missing x-file 6th character in command line ',0)   
         read(buf4,'(a1)')  xver 
         if( xver.eq.' ' ) then
           ixver = mod(iyear,10)
           write(xver,'(i1)') ixver 
         endif  
         iarg = iclarg(10,aprfile)     
         iarg = iclarg(11,buf4)  
         if( iarg.gt.0 )  then
            read(buf4,'(i4)') inter   
            cmdline_sess = .true.
         else
            inter = 30
         endif   
         iarg = iclarg(12,buf4)     
         if( iarg.gt.0 )  then
            read(buf4,'(i4)') ihr 
         else
            ihr = 0
         endif
         iarg = iclarg(13,buf4)     
         if( iarg.gt.0 )  then
            read(buf4,'(i4)') min 
         else
            min = 0
         endif   
         iarg = iclarg(14,buf4)     
         if( iarg.gt.0 )  then
            read(buf4,'(i4)') nepoch 
         else
            nepoch = 2880
         endif            
cd         print *,'DEBUG makexp1 iyr2 iday xver aprfile sp3file '
cd     .       ,iyr2,iday,xver,aprfile,sp3file
cd         print *,'inter ihr min nepoch ',inter,ihr,min,nepoch
              
c If not on the command-line, read the parameters interactively

      else 
        write(*,3)
    3   format(/,'    Create Input File for Programs',
     .         /,'   MAKEJ, MAKEX, BCTOT, and FIXDRV',
     .         /,'   --------------------------------',//)  
                                       
c **  rwk 020710: These next four lines were commented out with the last major revision
c **              but I no longer remember why.  We will need this input with the
c **              new station.info format because the new file does not have these parameters
        write(*,*) 'Enter 4-character project code'
        read(*,'(a4)') project
        write(*,*) 'Enter 4-character orbit code'
        read(*,'(a4)') orbt
* MOD TAH 200511: Allow lower frequency
        write(*,*) 'Enter GNSS code (G R C E J I, C7 G5 E6/7/8)'   
        read(*,'(a2)') gnsslf
        gnss = gnsslf(1:1)
        write(*,*) 'Enter year'  
        read(*,'(i4)') iyear 
        call check_y2k(iyear)
        iyr2=mod(iyear,100) 
        write(*,*) 'Enter day of the year (999 to search all)'
        read(*,'(i3)') iday_search
        write(*,*) 'Enter session number (99 to search all)'
        read(*,'(i2)') nrsess  
c       set the x-file version number  (change to allow interactive input?)
        ixver = mod(iyear,10)
        write(xver,'(i1)') ixver  
        write(*,*) 'Enter the l-file or apr file name (CR to use l-file'
     .              ,' default)'    
        read(*,'(a80)',iostat=ioerr) aprfile  
       if( aprfile(1:1).ne.' '.and. .not.fcheck(aprfile) ) then
         call report_stat('WARNING','MAKEXP','makexp',aprfile
     .        ,'L-file not found',0)
         write(*,*) 'Continue (y/n) ?'
         read(*,'(a)') ans
         if( ans.ne.'y' ) stop
       endif
c       MAF (2020-09-28; originally 2019-02-15), in replacement of incorrect and revoked fix by Zou Xuan (2018-09-04), below
        write(*,*) 'Enter SP3 file name or "none" to skip'
        read(*,'(a12)') sp3file
      endif      


c  Check for the existence of the navigation file and get the available SVs
c  ------------------------------------------------------------------------

      nav_found = .false.  
      if( navfile(1:4).ne.'none' ) then
   55   open(unit=lunav,file=navfile,status='old',err=60,iostat=ioerr)
        write(luprint,'(1x,a,a20,/)' ) 'Opened navigation file ',navfile
   60   if (ioerr .eq. 0) then    
           call checke ( lunav,nprn,gnss,iprn )
cd           print *,'MAKEXP after checke ',gnss,nprn
cd           write(*,'(200i4)') (iprn(i),i=1,nprn)
           nav_found = .true.
        else   
          if ( batch ) then  
            call report_stat('FATAL','MAKEXP','makexp',navfile
     .                     ,'Error opening navigation file ',ioerr)
          else
            write (*,*) ' MAKEXP: could not open navigation file: '
     .               ,   navfile
            write (*,*) " Enter name of RINEX nav file or e-file "
     .               ,"or 'none' to skip: "
            read  (*,'(a20)') navfile  
            if( navfile(1:4).eq.'none' ) goto 70
            goto 55 
          endif
        endif  

* MOD TAH 201215: Check to see if zero satellites, if so and jclock=sp3
*       get satellite llist form SP2 file
        if( nprn.eq.0 ) then
           call report_stat('WARNING','MAKEXP','makexp',navfile
     .                     ,'No satellites found ',0)
           if( jclock(1:1).eq.'s' .or.jclock(1:1).eq.'S' ) then
*             read SP3 file for satellites
              open(lusp3,file=sp3file,status='old',iostat=ioerr)
              if( ioerr.ne.0 ) call report_stat('FATAL','MAKEXP',
     .                     'makexp',sp3file,
     .                     'Error opening sp3file ',ioerr)
              call rsp3hd( lusp3,gnss, iyrsp3,imosp3, idysp3 
     .,                 ihrsp3, imnsp3, secsp3, delt,mjd,fmjd, nepsp3 
     .,                 numsp3sv,numsp3,itsat, accsat 
     .,                 pvsigb,clksigb,otlmod )
              close(lusp3)
*             Now copy to arrays from navfile
              nprn = numsp3
              iprn = itsat(1:nprn)
           endif
         endif 
* MOD TAH 201215: Finish new code; continue as if list from navfile.
c          
 
c       arrange in ascending order
cd        print *,'nprn ',nprn
        do j=2,nprn
          is = iprn(j)
          do i= j-1,1,-1
            if( iprn(i).le.is ) goto 61
            iprn(i+1)=iprn(i)
          enddo
          i=0
   61     iprn(i+1) = is
        enddo 
      endif
          

               
C  Create the session.info file if necessary
c  ------------------------------------------
c
c     No longer allow changing only the date of an existing file (too dangerous if SV list changes)
c     but retain checks to avoid overwriting a potentially valuable multi-day session.info
            
      if( cmdline_sess .or. .not.batch ) then

        if( fcheck('session.info') ) then
          open( unit=lusinfo,file='session.info',status='old'
     .        , iostat=ioerr)           
          if(ioerr .eq. 0) call sinfo_check( lusinfo )
        endif
        write(luprint,'(/,2a)') ' Command-line or interactive entries'
     .                     ,' override one-line session.info'   
        close(unit=lusinfo)
        open( unit=lusinfo,file='session.info',status='unknown'
     .         , iostat=ioerr)         
        if(ioerr .ne. 0)  call report_stat('FATAL','MAKEXP','makexp'
     .            ,' ',  'Error opening new session.info ',ioerr)    
        if( .not.nav_found )   
c          must have nav file to create session.info anew
     .     call report_stat('FATAL','MAKEXP','makexp',' '
     .        ,'Cannot create session.info without navigation file',0)
        if ( .not.batch ) then
           write(*,'(a)') 'Create the scenario file'
           write(*,'(/,2a)')
     .      'Enter sampling interval (s)   Start time (hh mm)'
     .     ,'   and Number of epochs'
           read(*,*) inter,ihr,min,nepoch
        endif
c       write the output session.info file
        if( iday.eq.0 ) iday = iday_search
        if( nrsess_search.eq.99 ) nrsess = 1
c         write two comment lines  (optional)
       write(lusinfo,'(a,/,a)') 
     . '# Session.info : free format, non-blank first column is comment'
     . ,'#Year Day  Sess#  Interval  #Epochs  Start hr/min  Satellites'
c       file is free-format, but put in blanks here for aesthetics (to match comment)
* MOD TAH 200610: Increased the number of PRN allowed in format from 32 to 50.
       write(lusinfo,'(1x,i4,1x,i3,4x,i1,6x,i4,4x,i5,4x,2i3,4x,50i3)')
     .      iyear,iday,nrsess,inter,nepoch,ihr,min,(iprn(i),i=1,nprn)
       close (unit=lusinfo)    

      else   
        if( fcheck('session.info') ) then
          write(luprint,'(/,a,/)') 
     .       ' No command-line session input, use existing session.info'
          if ( .not.batch ) write(*,'(/,a,/)') 
     .       ' No command-line session input, use existing session.info' 
          open( unit=lusinfo,file='session.info',status='old'
     .         , iostat=ioerr)  
          if(ioerr .ne. 0)  call report_stat('FATAL','MAKEXP','makexp'
     .            ,' ',  'Error opening old session.info ',ioerr)
          call rsesfo (lusinfo,.false.,0,iyear,iday,jsessn
     .                ,ihr,min,inter,nepoch,nprn,iprn )
          rewind(unit=lusinfo)   
        else
          call report_stat('FATAL','MAKEXP','makexp',' '
     .        ,' No session.info and no command-line session entries',0)
        endif
  
      endif 


c  Print the experiment summary 
c -----------------------------
* MOD TAH 200511: Report gnsslf
         write(luprint,'(/,a,a4,a,a4,a,a2,/,a,i4,1x,i3,1x,i2)')
     .      ' Project: ',project,'   Orbit: ',orbt,'  GNSS ',gnsslf
     .     ,' Search year,day,session: ',iyear,iday,nrsess_search
                                             

c  Get the list of RINEX or X files from the directory
c  ---------------------------------------------------

c     set RINEX/X flag to RINEX
      xorx = 'r'                            
cd      print *,'DEBUG makexp iday_search ',iday_search
      if( iday_search.eq.999 ) then 
        if( .not.batch ) then
          write(*,*) 'Global search on RINEX files selected'
          write(*,*) '  --enter day of the year for processing'
          read(*,*) iday  
        endif
        write(wild,'(a2,i2.2,a1)') '*.',iyr2,'o'
        call getdir(wild,maxrxfil,rinex,nobs_files) 
      else                
        iday = iday_search
        if( nrsess_search.eq.99) then
          write(wild,'(a4,i3.3,a2,i2.2,a1)') '????',iday,'?.',iyr2,'o'
          call getdir(wild,maxrxfil,rinex,nobs_files)
        else if (nrsess_search.eq.0 .or. nrsess_search.eq.1 ) then
c         this case is redundant, but GAMITeers will be used to entering
c         0 or 1 and expecting to get both, so keep the code until they
c         get used to entering 99  
          write(wild,'(a4,i3.3,a2,i2.2,a1)') '????',iday,'0.',iyr2,'o'
          call getdir(wild,maxrxfil,rinex,nobs_files) 
          write(wild,'(a4,i3.3,a2,i2.2,a1)') '????',iday,'1.',iyr2,'o'
          call getdir(wild,maxrxfil,rinex1,nobs_files1)  
          if( nobs_files1.gt.0 ) then                        
            do i=1,nobs_files1
               rinex(nobs_files + i) = rinex1(i)
            enddo
          endif
          nobs_files = nobs_files + nobs_files1  
        else     
          write(wild,'(a4,i3.3,i1,a1,i2.2,a1)') 
     .                '????',iday,nrsess_search,'.',iyr2,'o'
          call getdir(wild,maxrxfil,rinex,nobs_files)     
        endif
      endif 
c     create station-code names
      if( nobs_files.gt.0 ) then
        do i = 1,nobs_files
          scodes(i) = rinex(i)(1:4)
        enddo
      endif            


c   If no RINEX files, look for X-files

      if( nobs_files.eq.0 ) then 
         if( nrsess_search.eq.99 ) then
           write(message,'(a,i4,a,i4)')
     .       'No RINEX files found for year ',iyear,' day ',iday
         elseif ( nrsess_search.eq.0 .or. nrsess_search.eq.1 ) then
           write(message,'(a,i4,a,i4,a)')
     .        'No RINEX files found for year ',iyear,' day ',iday
     .      ,' session 0 or 1' 
         else
           write(message,'(a,i4,a,i4,a,i4)')
     .       'No RINEX files found for year ',iyear,' day ',iday
     .      ,' session',nrsess_search
         endif
         call report_stat('WARNING','MAKEXP','makexp',' ',message,0)
         call report_stat('WARNING','MAKEXP','makexp',' ',
     .      '---using X-files (and ignoring MAKEX.BATCH',0)
         write(kyear,'(i2)') iyr2
         write(wild,'(a6,a1,i3.3)') 'x?????','.',iday
         call getdir(wild,maxrxfil,xfile,nobs_files)
         if( nobs_files.le.0 ) then
           call report_stat('FATAL','MAKEXP','makexp',' '
     .                     ,'No RINEX or X-files found',0)
         else
           xorx = 'x'
           do i=1,nobs_files
              scodes(i) = xfile(i)(2:5)
           enddo
         endif
      endif
             
c     Sort the station codes alphabetically, removing duplicates

      call sort_string( maxrxfil,nobs_files,scodes,0,new_codes
     .                , nsites,scodes1 )
c     check now to make sure < 99 sites, the most allowed by cfmrg and solve
      if( nsites.gt.99 )  call report_stat('FATAL','MAKEXP','makexp',' '
     .                        ,'Number of x-files > 99 ',0)
      do i=1,nsites
         scodes(i) = scodes1(i)
      enddo     
                                                                

c     Create the X-file names and print the session summary 
                    
      do i=1,nsites
          write(xfile(i),'(a1,a4,a1,a1,i3.3)')
     .       'x',scodes(i),xver,'.',iday
      enddo     
      write(luprint,'(/,1x,i3,a,/)') nsites
     .      ,' X-files to be used or created:'
      write(luprint,'(1x,a20)') (xfile(i),i=1,nsites) 
      write(luprint,'(1x)')
      if (.not.batch ) then
       write(*,'(/,1x,i3,a,/)') nsites,' X-files to be used or created:'
       write(*,'(1x,a20)') (xfile(i),i=1,nsites) 
       write(*,'(1x)')   
      endif


c  Read the station.info file to check for an entry and get the required information
c-----------------------------------------------------------------------------------

      if( fcheck('station.info') ) then
        open(unit=lustnfo,file='station.info',status='OLD',iostat=ioerr)
        if (ioerr .ne. 0) then
           call report_stat('FATAL','MAKEXP','makexp','station.info'
     .                     ,'Error opening  station.info file',ioerr)
        endif
      else
         call report_stat('WARNING','MAKEXP','makexp',' '
     .        ,'File station.info not found',0)
      endif   

      write(luprint,'(/,a,/)') ' Reading station.info '
      do i = 1,nsites             
 
c       session 0 always convert to 1 internally in GAMIT
        if (nrsess.eq.0 .or. nrsess.eq.99 ) nrsess=1  
        sitcod = scodes(i)   
        isod = ihr*3600 + min*60
        span = inter*nepoch                               
* MOD TAH 200203: Added AntDAZ to list of values from station.info
       call rstnfo(lustnfo, sitcod, iyear, iday, isod
     .               , span, stanam(i), anthgt(i), offstn(i), offste(i)
     .               , antdaz(i), rcvcod(i), antcod(i)
     .               , hgtcod(i), radome(i), swvers(i)
     .               , rcvers, rcvrsn, antsn, jstart, jstop)   
        savcod(i) = sitcod     
        write(luprint,'(1x,a4,1x,a6,1x,a6,2(i5,4i3))') 
     .      sitcod,rcvcod(i),antcod(i),jstart,jstop
      enddo

C  Create other file names 
c  -----------------------
                  
      call lowers(navfile)
      if( navfile(1:4).ne.'none' ) then 
        if( navfile(7:7).eq.'.' ) then
c         filename is e????y.ddd'
          navf = navfile(2:5)
        else   
c         navfile is presumed to be ????ddd?.yyn
          navf = navfile(1:4)       
        endif
        if( navfile(5:6).eq.' ') then
          navf = project   
          write(aday,'(i3)') iday
          if( aday(1:1).eq.' ' ) aday(1:1) = '0'
          if( aday(2:2).eq.' ' ) aday(2:2) = '0'
          write(navfile,'(a4,a3,a2,i2.2,a1)') navf,aday,'0.',iyr2,'n'
         endif 
      endif
      if( jclock.eq.'brdc' ) then
        write(jfile,'(a1,a4,i1,a,i3.3)') 'j',navf,mod(iyr2,10),'.',iday 
cd        print *,'MAKEXP DEBUG jfile ',jfile
        sp3file = 'none        '
      else
        write(jfile,'(a1,a4,i1,a,i3.3)') 'j',orbt,mod(iyr2,10),'.',iday 
cd        print *,'MAKEXP DEBUG jfile ',jfile
        call doygwk(iday,iyear,gpsw,gpsd)
C       Added by Zou Xuan 2018/9/4 (revoked; Floyd 190121/20200928)
C       write(sp3file,'(a3,i4,i1,a4)')
C    .                orbt(1:3),gpsw,gpsd,'.sp3'
      endif
      write(tfile,'(a1,a4,i1,a,i3.3)') 't',orbt,mod(iyr2,10),'.',iday
      write(dfile,'(a1,a4,i1,a,i3.3)') 'd',project,mod(iyr2,10),'.',iday
      if( gnss.eq.'R'.and.jclock.eq.'brdc' ) then
c       Glonass clocks cannot be computed from the nav-file since the orbits 
c       are not accurate enough without an integration
        write(ifile,'(a4)') 'NONE'
      else
        write(ifile,'(a1,a4,i1,a,i3.3)') 
     .     'i',project,mod(iyr2,10),'.',iday
      endif
      if( aprfile(3:3).eq.' ')  write(aprfile,'(a1,a4,i1,a,i3.3)')
     .                                 'l',project,mod(iyr2,10),'.',iday
      write(luprint,'(2a)') ' Coordinate file is ',aprfile
  

C  If RINEX input, create the batch and user input files for program MAKEX
c  -----------------------------------------------------------------------

   70 if( xorx.eq.'r' ) then
                 
        makexbatch = project//'.'//'makex.batch'
        open(unit=lumakexb,file=makexbatch,status='unknown')
         
        write(lumakexb,'(a)')  'infor 1'
        write(lumakexb,'(a)')  'sceno 1 session.info'
        write(lumakexb,'(a)')  'rinex 1 ./' 
        write(lumakexb,'(a)')  'fica  0 '
        write(lumakexb,'(2a)') 'coord 1 ',aprfile
        write(lumakexb,'(a)')  'stnfo 1 station.info'
        write(lumakexb,'(a)')  'xfile 1 ./x'
        write(lumakexb,'(2a)') 'svclk 1 ',jfile
        write(lumakexb,'(a)')  'clock 1 ./k' 
cd        print *,'DEBUG makexp sp3file ',sp3file 
        if( sp3file(1:4).eq.'none' ) then                
          write(lumakexb,'(a)')  'sp3   0 '
        else
          write(lumakexb,'(2a)') 'sp3   1 ',sp3file
        endif
        write(lumakexb,'(2a)')  'rdorb 1 ',navfile
* MOD TAH 200511: Write the gnsslf version
        write(lumakexb,'(2a)')  'gnss  1 ',gnsslf
        write(lumakexb,'(a)') 'site year doy sn  sw  ver'
        write(lumakexb,'(a)') '(a4,1x,a4,1x,a3,1x,a1,2x,a3,1x,f4.2)'
              
cd        print *,'before makex codes nsites ',nsites
        do k = 1,nsites

c         Create the MAKEX firmware codes
c--old format
c          write(filcod(k),'(a4,i2.2,i1,a1,i3.3)')
c     .               scodes(k),iyr2,nrsess,'.',iday
c--new format
           write(filcod(k),'(a4,1x,i4,1x,i3.3,1x,i1)')
     .                scodes(k),iyear,iday,nrsess
c         Translate TI4 to COR or GES for TI400
          if( rcvcod(k)(1:3).eq.'TI4') then
             rcvcod(k)(1:3)='COR'
             if (swvers(k).lt.4.1) rcvcod(k)(1:3)='GES'
          endif
c         Translate WM-102 to WM2
          if( rcvcod(k).eq.'WM-102') then
             rcvcod(k)(1:3)='WM2'
          endif   
c         Translate 'AT' receiver codes to ASH firmware
          if( rcvcod(k)(1:2).eq.'AT' ) then
             rcvcod(k)(1:3)='ASH'
          endif      
c         Translate Topcon receiver codes to the original manufacturers name
          if( rcvcod(k).eq.'TOPRID' ) rcvcod(k) = 'ASHL12'
          if( rcvcod(k).eq.'TOPRIP' ) rcvcod(k) = 'ASHP12'
          if( rcvcod(k).eq.'TOPRIY' ) rcvcod(k) = 'ASHZ12' 
          if( rcvcod(k).eq.'TOPTRB' ) rcvcod(k) = 'AO_8RS'  
          if( rcvcod(k).eq.'TOPDX1' ) rcvcod(k) = 'TR4700' 
          if( rcvcod(k).eq.'TOPLEG' ) rcvcod(k) = 'JPSLEG'    
          if( rcvcod(k).eq.'TOPODY' ) rcvcod(k) = 'JPSODY'
          if( rcvcod(k).eq.'TOPHIP' ) rcvcod(k) = 'JPSHIP'  
          if( rcvcod(k).eq.'TOPSSI' ) rcvcod(k) = 'TRMSSI'  

c         Translate Trimble receiver codes to TRM firmware
cd         print *,'k orig rcvcod ',k,rcvcod(k)
          if( rcvcod(k)(1:3).eq.'TR4' .or. rcvcod(k)(1:3).eq.'TR7' .or. 
     .        rcvcod(k)(1:3).eq.'TR5'.or. rcvcod(k)(1:3).eq.'TRN' ) then
             rcvcod(k)(1:3)='TRM'
          endif                               
cd          print *,'translated rcvcod ', rcvcod(k)
c         Translate TR8000, TR8100, ROG800, ROG_12, and 'AO' rcvr codes to TRB firmware
          if( rcvcod(k).eq.'TR8000'.or. rcvcod(k).eq.'TR8100' .or.
     .        rcvcod(k).eq.'ROG800'.or. rcvcod(k).eq.'ROG_12' .or. 
     .        rcvcod(k)(1:2).eq.'AO') then
             rcvcod(k)(1:3)='TRB'
          endif 
c         Translate SRxxxx, LCxxxx, and CRxxxx to LEI 
          if( rcvcod(k)(1:2).eq.'SR' .or. rcvcod(k)(1:2).eq.'LC' .or.
     .        rcvcod(k)(1:2).eq.'CR' ) then
             rcvcod(k)(1:3)='LEI'
          endif         
c         Translate NTCxxx, SF2xxx, or RT3xxx to NAV for NavCom receivers
          if( rcvcod(k)(1:3).eq.'NCT' .or. rcvcod(k)(1:3).eq.'SF2' .or.
     .        rcvcod(k)(1:3).eq.'RT3' ) then
             rcvcod(k)(1:3)='NAV'
          endif       


c         write the line for each station  
c--       old format
c--          write(lumakexb,'(a11,2X,A3,1X,F4.2)')
c--     .       filcod(k),rcvcod(k)(1:3),swvers(k)
c--       new format
             write(lumakexb,'(a15,2x,a3,1x,f5.2)')
     .          filcod(k),rcvcod(k)(1:3),swvers(k)
        enddo
c       close the MAKEX batch file
        close(unit=lumakexb)
c       No longer create makexp.inp -- use command-line args instead:  makexp [expt].makex.batch
c       open the makex.inp file, write a single line ('1'), and close it
c        open(unit=lumakexb,file='makex.inp',status='unknown')
c        write(lumakexb,'(I1)') 1
c        close(unit=lumakexb)

      endif

C  Create user input file for program BCTOT
C  ----------------------------------------

c     read the sestbl to get the inertial frame
      sestbl = 'sestbl.'
      if (fcheck(sestbl)) then
         open(lusest,file=sestbl,status='old',iostat=ioerr)
         if(ioerr .ne. 0) then
           call report_stat('FATAL','MAKEXP','makexp',' ',
     .     'Error opening file sestbl. ',0)
         endif
         reqd = .false.
         call rdsest( 14, 'Inertial frame', 5, buf5, lusest, reqd, ill )
         if( ill .ne. 0 ) then
            call report_stat('FATAL','MAKEXP','makexp',' ',
     .     'Error while reading sestbl. ',0)
         endif
c        buffer will be blank if nothing found
         if( buf5.ne.'     ' )  frame = buf5
         if( frame.ne.'B1950' .and. frame.ne.'J2000' ) then
            write(message,80)
   80       format('Valid inertial frame not found in sestbl. '
     .      ,'Assuming inertial frame is J2000')
            call report_stat('WARNING','MAKEXP','makexp',' ',message,0)
            frame = 'J2000'
         endif
      else
         call report_stat('WARNING','MAKEXP','makexp',' ',
     1   'Could not open sestbl. Assuming inertial frame is J2000',0)
cd         print *,' '
         frame = 'J2000'
      endif
      open(unit=lumakexb,file='bctot.inp',status='unknown')
      write(lumakexb,'(A,1x,i4,1x,i3)') 'b',iyear,iday
      write(lumakexb,'(A)') navfile
C     Use blank instead of x-file name to force use of session.info
      write(lumakexb,'(A)') '                 '
      write(lumakexb,'(A)') tfile
      write(lumakexb,'(A)') 'y'
      write(lumakexb,'(a)') frame
      close(unit=lumakexb)


C  Create the D-file
C  -----------------
                             
cd      print *,'before d-file nsites ',nsites
      open(unit=lumakexb,file=dfile,status='unknown')
      nsess=1  
* NO MOD TAH 200511: D-file only needs gnss not the lower frequency
      write(lumakexb,'(a1)') gnss
      write(lumakexb,'(i2)') nsess
      write(lumakexb,'(a)')  aprfile
      write(lumakexb,'(a)')  tfile
      write(lumakexb,'(a)')  ifile
      write(lumakexb,'(a)')  jfile
      write(lumakexb,'(i2)') nsites
      do i=1,nsites
        write(lumakexb,'(a)') xfile(i)
      enddo
      close(unit=lumakexb)


c  Create the user input file for program FIXDRV
c  ---------------------------------------------
c                                        
c     No longer create fixdrv.inp -- use command-line args instead:  fixdrv [d-file]
c     open(unit=lumakexb,file='fixdrv.inp',status='unknown')
c     write(lumakexb,'(a20)') dfile
c     close(unit=lumakexb)
   

c  Write the summary to makexp.out
c  -------------------------------

* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
      write(luprint,'(a,i4,/,a,i5,/,a,2i3,/,a,50i3)')
     .        ' Sampling interval =',inter
     .      , ' Number epochs     =',nepoch
     .      , ' Start time (HHMM) =',ihr,min
     .      , ' Satellites        =',(iprn(i),i=1,nprn)
           if( .not.batch ) write(*,'(a,i4,/,a,i5,/,a,2i3,/,a,50i3)')
     .        ' Sampling interval =',inter
     .      , ' Number epochs     =',nepoch
     .      , 'Start time (HHMM) =',ihr,min
     .      , 'Satellites        =',(iprn(i),i=1,nprn)

      call report_stat('STATUS','MAKEXP','makexp',' ',
     .'Normal end in Program MAKEXP ',0)   

      if( .not. batch ) then  
        write(*,'(1x)')
        write(*,'(/,a)') '---------------------------------------------'
        write(*,'(/,a,/)') 'Now run, in order:'  
        write(*,'(2a)') ' sh_sp3fit -f <sp3 file> OR sh_bcfit bctot.inp'
     .                  ,'     OR copy a g-file from SOPAC'
        write(*,'(/,a,i3,a)') ' sh_check_sess -sess ',iday
     .                        ,' -type gfile -file <g-file>'
        if( jclock.eq.'brdc' ) then
          write(*,'(/,a,a15,a10)') ' makej ',navfile(1:15),jfile
        else
          write(*,'(/,a,a15,1x,a10,1x,a12)') 
     .              ' makej ',navfile(1:15),jfile,sp3file
        endif
        write(*,'(/,a,i3,a,a20)') ' sh_check_sess -sess ',iday
     .                            ,' -type jfile -file ',jfile
        write(*,'(/,a,a20)') ' makex ',makexbatch
        write(*,'(/,a,a20,a)') ' fixdrv ',dfile,' OR  run interactively'
        write(*,'(/,a)') '--------------------------------------------'
      endif
      stop
      end


      Subroutine sinfo_check ( lusinfo )   

      implicit none

      integer*4 lusinfo,icnt,iys,ids,ioerr

      character*256 line

      icnt=0
   90 read(lusinfo,'(a)',iostat=ioerr,end=91) line  
        if( ioerr.gt.0 ) call report_stat('FATAL','MAKEXP','sinfo_check'
     .     ,' ','Error reading session.info',ioerr)
c       a non-valid data line is any that doesn't have a reasonable year and
c       day-of-year (no need to check explicitly on non-blank first column
        read(line,*,iostat=ioerr) iys,ids
c       no iostat check since non-integer tokens will give an error
cd        print *,'iys ids ioerr ',iys,ids,ioerr
        if( iys.lt.80 .or. iys.gt.2100   .or.
     .      ids.le.0 .and. ids.gt.366 ) then
          goto 90
        else
          icnt = icnt + 1
        endif   
        goto 90
   91   if( icnt.eq.0 ) then
          call report_stat('WARNING','MAKEXP','sinfo_check',' '
     .                  , 'No valid records on old session.info ',ioerr)
        elseif ( icnt.gt.1 ) then
          call report_stat('FATAL','MAKEXP','sinfo_check',' '
     .     ,'Stop to avoid overwriting a multiline session.info ',ioerr)
        endif
           
        rewind ( lusinfo )
        return
        end
