Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego, 1996.  All rights reserved.      

      Program MAKEX
c
c     Make x files from RINEX files
c
c    (FICA no longer supported: use fic2rx to translate.) 
c
c     Runs from command line arguments:
c
c       makex <batch file> <rx_doy_minus> <rx_doy_plus>
c
c        where <batch file> is the MAKEX batch file required; if the name is
c                           debug.makex.batch, debug print will be invoked
c              <rx_doy_minus> is the number of days before the present to search for
c                             data in RINEX files (optional, default = 7)
c              <rx_doy_plus>  is the number of days after the present to search for 
c                             data in RINEX files (optional, default = 1)
c
c     Kurt Feigl June, 1987
c     Extensive mods by Kurt Feigl and Bob King  Oct 88 - Jan 90
c     Mods for kinematics by Jeff Genrich and Bob King, Nov 91 - Feb 92 
c     Mods for pre-searching RINEX files for data.  Bob King Feb 97
c     Roots:  sledgehammer         peter morgan january 1987
c         /rhea    sledgehammer.test    Jerry Svarc on 4-30-87
c     Last modified by R. King Aug 2002

c     This program opens a lot of files.  Each file is designated in this
c     program by 3 variables: a file name, unit number, and a logical flag.
c     The file name is a string containing the path name of the file.
c     The unit number gives the unit where the file is connected
c     The logical flag is .true. if we are using this file.
c     These 3 variables have a naming convention fXXXXX, uXXXXX, qXXXXX,
c     respectively, where XXXXX is gipwdven in the table below.  For example,
c     The batch input file is described by fbatch and ubatch. All these
c     variables are declared and included in a common block in makex.h.
c     Please try to maintain this convention!
c
c     The files are:
c
c     * batch lu  4    Command file, assumed to exist.
c     *       lu  5    Terminal keyboard.
c     * scren lu  6    Terminal screen.
c     * infor lu  8    Log file, created.
c     * sceno lu 10    Scenario (session.info) file, assumed to exist.
c     # rinex lu 11    RINEX data file.
c     * coord lu 12    Coordinates (L-) file, assumed to exist.
c     * sited lu 13    Station information (station.info)file, assumed to exist.
c     #       lu 14    (formerly FICA file - NO LONGER USED)
c       xfile lu 20    X-file, created.
c     * svclk lu 21    SV clock (J-) file, assumed to exist
c       clock lu 22    Station clock (K-) file, created.
c       sp3   lu 23    Orbit sp3 file 
c     * nav   lu 24    Orbit RINEX navigation file    
c
c
c     The starred files must exist, plus a RINEX observation file;
c          the others are optional.
c     The subroutine OPENF opens all files, and stops if it can't.
c
      implicit none
c
      include '../includes/dimpar.h'  
      include '../includes/units.h'
      include '../includes/makex.h' 
      include '../includes/global.h'
      include '../includes/errflg.h'
      include '../includes/model.h'
c
c                        DECLARATIONS
c                        ************
c
      logical          timeset,debug,opnerr,fcheck
      logical          epochok,late_epoch,satok,pastend_epoch
      logical          fend,ferr,newfile,ant_event,ircint_msg
      logical          clkwarnings,antwarnings
      logical          old_stinf 
      logical          dcb_override                             
c     temporary while old-style station.info supported
      character*256 line

c     session variables
      character*4 upperc  
      integer*4  interval,ideci,nsats
              
c     site info
      character*4 sites(maxsit),site
      character*3 rcvrs(maxsit) 
      real*4      vers(maxsit)

c     counter and list of RINEX files to be searched
      integer*4 infilcount,rx_doy_minus,rx_doy_plus
      character*80 infiles(maxfil)

C     run identifiers
      character        runnam*10
      character*1      dots(60)  
      character*3      rcvrsw 
      character*16     uname
C     the version number of this prog  (short: makex-only, long: w/ lib version)
      character*40      version40
      character*109     version
C     header lines for X-file
      integer*4        nxhdlin
      character*80     xhdlin(maxlin)
c     prefixes for file names
      character*80     pclock,pxfile
c     K-file headers and format statment (new 160831)
      character*105 kheader1
      character*126 kheader2 
      character*12  ksource
      character*51 afmt
      character*4  kversion

C     standard GAMIT error code
c     integer*4        ier(maxchn)  (model.h)
                       
c     new code to indicate need for C1/P1 corrections
      integer*4        idcb(maxchn)
                                           
C     run time yr,mo,da,hr,min,sec,fraction
      integer*4        irunt(6),ihnsec

c     iflag=0 for scanning only (get_rxfiles; 1 to return observables)
      integer*4        iflag
C     for broadcast ephem
      integer*4        nprn
C     sat id character and number
      character*1      asvid(maxchn)
      integer*4        isvid(maxchn)
C     tracking mode
      integer*4        tmode(maxchn)
C     gps week number for epoch
      integer*4        igpswk(maxchn)
C     quality vector
      integer*4        iqvec(maxchn,2)
c     site number and loop variable
      integer*4        numsit,is
c     dumb programming variables
      integer*4        i,j,k,idots
C     how many data records?
      integer*4        iread

C     true if the week number is OK
      logical          weekset

      integer*4        iyr,id,im,iy1,
     1                 ibadtag,idflag,ireject

c     station coordinates
      real*8           coords(3)
      real*8           seclat,seclon
      integer*4        dlat,dlon,mlat,mlon
      character*1      latflag,lonflag

C     Quantities for receiver clock correction
      real*8           tau
c     rclock, svclock declared in model.h 
c
      real*8           dsecs,deltak,ltvel,deltatag

c     values for x-file
c     receiver software version number (GAMIT)
c      real*4 swver  - in includes/model.h
c     station offset variables
      real*8           offset_l1(3),offset_l2(3)

c     integer*4 dattyp(maxdat),lambda(maxsat,maxdat) --in includes/model.h
C     L1, L2 doppler phase in cycles
      real*8      dofl1(maxchn), dofl2(maxchn)
C     L1, L2  pseudorange
      real*8      prgl1(maxchn),prgl2(maxchn)
      
C     quality for L1,L2
      real*8      disnr(maxchn,2) 

c     SV coordinates from SP3 or nav file
      real*8 xsat(3)

c     Primary time quantities 
C     gps sec of week for epoch
      real*8      gpssec(maxchn)
c       selection of data will always have: start < open < tag < close < stop
c     tag on current epoch
      integer*4   iwkntag,iyrtag,idoytag,isodtag,iwknchck
C     first epoch to write
      integer*4   iwknstart,iyrstart,idoystart,isodstart
C     last epoch to write
      integer*4   iwknstop
C     beginning of window for current epoch
      integer*4   iwknopen
C     end of window for current epoch
      integer*4   iwknclose
C     week number in the data file
      integer*4   iwknfile
c     week number of ephemeris
      integer*4   iwkne
c     week number for K-file
      integer*4   iwknk
C     tag on current epoch
      real*8      sowtag,sodtag    
c     tag for ephemeris
      real*8      sowe 
c     tag corrected by PR for checking epoch alignment
      real*8      sowchck
C     first epoch to write
      real*8      sowstart,sodstart
C     last epoch to write
      real*8      sowstop
C     beginning of window for current epoch
      real*8      sowopen
C     end of window for current epoch
      real*8      sowclose
c     K-file epoch
      real*8      sowk,sodk 
      integer*4   nsowk,jdk
c     j-file sv clock quantities 
      integer*4 jdtoc
      real*8    toc,svepc,svcrat,svcacc,valid

c     currently processing this epoch  (epoch1 = epoch -1 to check for 0 in rrinex)
      integer          epoch,epoch1   

c     number of good channel this epoch
      integer          ngood_chan

C     keep track of how many full epochs
      integer          ifull(maxchn+1)  

C     window define +/- this many seconds
      real*8      slop
C     offset from GPS time to UTC time
      real*8      utcoff
C     number of seconds in scenario
      integer*4 ispan

c     functions
C     returns difference between 2 epochs
      real*8        secdif
C     removes trainling blanks
      character     unpad*80
C     return the channel number
      integer*4     inscen
C     true if the data are acceptable
      logical       isgood,iyuk
c     returns the Julian day
      integer*4 julday  
c     return difference in seconds between two integer epochs y/doy/h/m/s
      integer*8 itimdif
c     return TAI-UTC
      real*8 taiutc
                      
      character*3     buf3
      character*40    buff40

c     prefixes for file names
      save             pclock,pxfile

c     things for RINEX
      character*3  rxobtyp(maxobt)
      character*20 rxpgm     
c     indices for observable types to be used
      integer*4 iobtypx(6)
      real*4 rxver  
      character*3 rxtime 
      integer*4    nobtyp,issi(maxchn,4),illi(maxchn,4)
c     wavelength factors
      integer*4 nwave1,nwave2
      real*8 anth, ante, antn

c     variables for testing for gaps
      logical gap
      integer obsmat(maxsat,maxepc),itimel
     .      , ioktot(maxsat,maxsit),maxtot(maxsat,maxsit)
     .      , maxgap,itotal,igap,ibad

c     variables for station.info 
      character*4 sitcod
      character*5 hgtcod
      character*9 exptyp
      character*16 stanam,amodel
      integer*4 icall,itimes(5)  
c     kstarts(5),kstops(5) in ..includes/model.h
      real*4 swver_stnfo
      real*8    anthgt,offstn,offste,dhpab
      character*512 message   
                     
c     obsolete variables for kinematic 
      integer*4 kflag
c     temporary and general variables  
      integer*4  min,iserr,idum,check_flg
     .        ,  iclarg,len,iarg
     .        ,  number,ihr,itflag,jyr,jdoy
     .        ,  sats(maxsat),isnrx(maxchn,2)
     .        ,  month,iday
      real*8      rdum,rdum3(3),site_epoch,decyrs
      character*5 char5
      character*6 char6
      character*20 char20
      character*80 wildcard,pickfn
            
c     k-file format 
      data afmt/'(2i4,2i3,f11.7,i6,f15.7,2x,a1,i2.2,3f17.9,i6,f15.7)'/ 

c     light velocity (m/s)
      data ltvel/2.99792458d8/

c     Initialize file unit numbers
c     ****************************

      uinfor = 8
      uscren = 6
      usceno = 10
      urinex = 11
      ucoord = 12
      usited = 13
c*     uficaf = 14 no longer used
      uxfile = 20
      usvclk = 21
      uclock = 22
      usp3 = 23
      unav   = 24
      uanthi = 32
c    kin change
      iul = ucoord
c    kin end change

c     Get the run time,user name and version number
c     *********************************************

      call mversn(version)
      version40 = version(1:40)
      call getdat(irunt(1),irunt(2),irunt(3))
      call gettim(irunt(4),irunt(5),irunt(6),ihnsec )
      call getusr(uname)
                            
c     write makex status line
      WRITE(uscren,'(a)')' '
      WRITE(message,5)version
    5 FORMAT('Started Makex ',a109)
      call report_stat('STATUS','MAKEX','makex',' ',message,0)

c     write (uscren,10) version,uname,(irunt(i),i=1,6)
   10 format (/,1x,'MAKEX ',a,/
     .   ,1x,'Run by ',a,1x,'on',1x,
     .   i4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2)
c      call proper(uscren)

c     Print warnings for missing antenna info
c     ***************************************
      antwarnings = .true.
                           

c     Read the command line
c     *********************
                                 
c     allow command-line arguments
      wildcard = '*.makex.batch'
      len = 13
      if( iclarg(1,fbatch).eq.0 ) then 
        write(uscren,'(//,a)') 'To run with command-line arguments: '
        write(uscren,'(/,a)') 
     .  ' makex  <batch-file>  <rx_doy_minus>   <rx_doy_plus>'
        write(uscren,20)
   20   format(/
     .   ,4x,'where <batch_file> is the batch file name (required)'
     ,   ,/,10x,'<rx_doy_minus> is days backward to search (default 7)'
     ,   ,/,10x,'<rx_doy_plus>  is days forward to search (default 1)'
     .   , //)
        write (uscren,*) 'Choose a batch file.'
        fbatch = pickfn(wildcard, len )
        fbatch = fbatch(1:len)
      endif
      i = index(fbatch,'.')
      runnam = fbatch(1:i)
      if (runnam(1:5) .eq. 'debug') then
         debug = .true.
      else
         debug = .false.
      endif
      ubatch = 4          

c     Read in the optional limits to RINEX file-days to be searched
c        allow 7 days in the past to account for week-long RINEX files
c        named for the start day; allow 1 day in the future to allow for
c        mis-named files.  Will need to input arguments for case of 
c        multiday files named with the stop day
      iarg = iclarg(2,buf3)   
      if( iarg.ne.0 ) then 
        read(buf3,'(i3)') rx_doy_minus  
      else
        rx_doy_minus =7  
      endif
      iarg = iclarg(3,buf3)    
      if( iarg.ne.0 ) then
        read(buf3,'(i3)') rx_doy_plus
      else
        rx_doy_plus = 1  
      endif
c     if( debug)  print *,'rx_plus rx_minus: ',rx_doy_plus, rx_doy_minus
            

c     Open the MAKEX log file
c     ***********************
                   
      finfor=unpad(runnam,'makex.infor')
      call openf (uinfor,finfor,'unknown','formatted',
     .  'sequential',.true. )
      if( opnerr ) then
         call report_stat('FATAL','MAKEX','makex',' ',
     .   'Error, cannot open information file',0)
      endif    
cd      print *,'MAKEX finfor uinfor version ',finfor,uinfor,version
      write (uinfor,10) version,uname,(irunt(i),i=1,6)
      call proper(uinfor)
             


c     Read the Input Batch File
c     *******************************************

      call openf (ubatch,fbatch,'old','formatted','sequential',
     .   .true. )
      if( opnerr ) then
         call report_stat('FATAL','MAKEX','makex',' ',
     .   'Error, cannot open batch file',0)
      endif 
c     rbatch reads the file names into common and returns the session parameters
* MOD TAH 200511: Return gnsslf with optional lower frequency.
      call rbatch
     .    ( debug,iyrstart,idoystart,gnsslf,numsit,sites,rcvrs,vers )
      gnss = gnsslf(1:1)
                                        
c     Begin loop on sites
c     *******************

      do is=1,numsit    
c     (this loop terminated at the very end of the routine)

        site = sites(is)
        rcvrsw = rcvrs(is)
        swver = vers(is)

        if(debug) write(uscren,110) upperc(site),iyrstart,idoystart
        write(uinfor,110) upperc(site),iyrstart,idoystart
  110   format(//,1x,'*********************************',/,
     .            1x,'BEGIN PROCESSING: ',a4,1x,i4,1x,i3,/,
     .            1x,'*********************************',/)
        write(message,115) upperc(site),iyrstart,idoystart
  115   format('**Begin processing: ',a4,1x,i4,1x,i3)
        call report_stat('STATUS','MAKEX','makex',' ',message,0) 
c       write into warning file to provide session id for warnings
        call report_stat('WARNING','MAKEX','makex',' ',message,0)
c         X- and K-file names for GAMIT are xssssy.ddd
c         ssss = 4 char station id
c         yy = last digit of year.
c         ddd = day of year 
        iy1 = mod(iyrstart,10)
        if(qxfile) write(fxfile,'(a,a4,i1,a1,i3.3)') 
     .                         'x',site,iy1,'.',idoystart
        if(qclock) write(fclock,'(a,a4,i1,a1,i3.3)') 
     .                         'k',site,iy1,'.',idoystart
           
c     Record the GNSS in the print file
c     *********************************

      write(uinfor,'(a,a2)') 'Collecting data for GNSS ',gnsslf
                         
c     Open all files except for input RINEX and output X
c     **************************************************

        call openf (usceno,fsceno,
     .   'old','formatted','sequential',qsceno )
        call openf (ucoord,fcoord,
     .    'old','formatted','sequential',qcoord )
      call openf (usited,fsited,
     . 'old','formatted','sequential',qxfile )
c rwk 971215: don't open until we know there are data 
c      call openf (uxfile,fxfile,
c     . 'unknown','formatted','sequential',qxfile )
      call openf (usvclk,fsvclk,
     . 'old','formatted','sequential',qsvclk )
      read(usvclk,'(a)') line  
      backspace (usvclk) 
      call openf (uclock,fclock,
     .  'unknown','formatted','sequential',qclock )
      call openf(usp3,fsp3  ,
     .  'unknown','formatted','sequential',qsp3 )
      call openf (unav  ,fnav  ,
     .  'old','formatted','sequential',qnav   )
c     define the file name for the hi.dat file
      fanthi = 'hi.dat'
      qanthi = .true.
      call openf (uanthi,fanthi,
     . 'old','formatted','sequential',qanthi )  
        read(uanthi,'(a)') line
c     hold the prefixes for the file names
      pclock = fclock
      pxfile = fxfile

      clkwarnings = .false.

c     say all epochs are bad for gap testing
      do i = 1,maxepc
         do j=1,maxsat
           obsmat(j,i) = 1
         enddo
      enddo

c     initialize the counter for data rejected for unreasonableness
      ireject = 0   


c     Read the scenario (session.info) file to get the times and satellites
c     *********************************************************************
           
cd     print *,'calling RSESFO '
       check_flg = 1    
       call rsesfo( usceno,debug,check_flg,iyrstart,idoystart,isessn
     .            , ihr,min,interval,nepoch,nsats,sats ) 
       if( nepoch.gt.maxepc ) then
       write(message,'(a,i5,a,i5,a)')  '# epochs (',nepoch
     .    ,') gt maxepc {',maxepc
         call report_stat('FATAL','MAKEX','makex',' ',message,0)  
       endif
       dsecs = 0.d0
cd       print *,'after RSESFO iyrstart,idoystart,ihr,min,dsecs,interval '
cd     .        ,              iyrstart,idoystart,ihr,min,dsecs,interval  
cd       print *,'             isessn,nsats,sats '
cd     .        ,              isessn,nsats,sats

c     Convert the start time 
c     ***********************
 
c** RWK 060815
c*    get sec-of-day in GPST, avoiding day boundaries by adding 0.5 sec
c*      itflag = +4             
c*      call timcon (itflag,iwknstart,sowstart,iyrstart,jdoy,ihr,min,dsecs
c*     .            , utcoff ) 
      itflag = -4 
      call timcon ( itflag,iwknstart,sowstart,iyrstart,idoystart
     .            , ihr,min,dsecs,utcoff ) 
      isodstart = 3600*ihr + 60*min + idint(dsecs+0.1)   
      sodstart = dble(isodstart)   
cd      print *,'AFTER timcon' 
cd      print *,'             iwknstart,sowstart,interval,nepoch '
cd     .                    , iwknstart,sowstart,interval,nepoch


c     Read the SV clock entries into memory
c     **************************************
      icall = 0 
      if (qsvclk) then
        icall = 0      
        call readj( usvclk,sats,nsats,idum,idum,rdum,rdum
     .            , icall,idum,rdum,rdum,rdum,rdum,rdum )
      endif
                  
c     Read the SP3 file or nav-file file into memory
c     **********************************************
cd      print *,'ICALL 0 qsp3 qnav ',qsp3,qnav
      if( qsp3 ) then 
c       requires an exact time match since position-only
        call getsp3(debug,icall,idum,rdum,gnss,idum,xsat,rdum,satok)
      elseif ( qnav ) then     
c       takes the nearest previous point (Keplerian elements)
        call getnav ( debug,icall,idum,rdum,gnss,idum,idum,sats
     .              , rdum3,rdum,satok,idum,rdum )
      else
        call report_stat('FATAL','MAKEX','makex',' '
     .                  , 'Neither SP3 nor RINEX nav file available',0)
      endif                                                                      

c     Write the header lines for the K-file
c     *************************************
      kheader1 = ' ' 
      kheader2 = ' ' 
      if( qsp3 ) then
        ksource = fsp3(1:12) 
      elseif( qnav ) then
        ksource = fnav(1:12)
      endif                  
      kversion = ' 2.0'                                                
      kheader1= 'Ver '//kversion//' Station clock values from '//ksource
     .//fsvclk(1:12)//'  MAKEX '//version40
      if( qsp3 ) then
      kheader2 = 'YEAR DOY HR MN  SEC(UTC)   WKNO   SOW(GPST)     PRN
     .OBSERVED PR(sec)    SV CLOCK    SITE CLOCK ' 
      else
      kheader2 = 'YEAR DOY HR MN  SEC(UTC)   WKNO   SOW(GPST)     PRN
     .OBSERVED PR(sec)    SV CLOCK    SITE CLOCK   NAV-FILE WKNO SOW(GP
     .ST)' 
      endif
      write(uclock,'(a,/,a,/,a)') kheader1,kheader2,afmt 
              

c     Get the site coordinates
c     ************************

c     spherical coordinates stored in 'modkin.h' common/kinpar1/
      sitecd = site
      iul = ucoord   
c     determine whether l-file or apr file 
      call crd_file_type(fcoord,kfflg)
      site_epoch = decyrs( iyrstart,idoystart,sodstart ) 
      call lread ( site,site_epoch ) 
      do i=1,3
        coords(i)= kpos(i) + kvel(i)*(site_epoch-kepoch0)
      enddo


c     Set the critical time variables
c     *******************************

c     The variables iwknstart and sowstart contain in GPS week
c     number and seconds of week the required starting epoch.
c     Determine these values from the receiver and software version
c     While we are there, decide the number of channels in the receiver

      call settim( debug,rcvrsw,swver,iwknstart,sowstart,slop )
      if(debug) print *,'After SETTIM rcvrsw swver iwknstart sowstart '
     .                 ,              rcvrsw,swver,iwknstart,sowstart

c     The interval for the K-file depends on the source of ephemerides
c     and receiver  
      call uppers (rcvrsw)   
      if( qsp3 ) then
c       if an sp3 file, then must match the 15-minute intervals exactly
        deltak = 900.d0
      elseif (rcvrsw .eq. 'TRM' .or. rcvrsw. eq. 'ASH'
     .   .or. rcvrsw. eq. 'TOP') then
c        120 seconds for Trimbles and Ashtechs
         deltak = 120.d0
      else
c        900 seconds for everybody else.
         deltak = 900.d0
      endif
      if (debug) print *,'RCVRSW, DELTAK = ',rcvrsw,deltak

c     calculate the span and stop time for the scenario
      ispan = nepoch*interval 
      write(message,'(a,i5,a,i3,a,f4.1)')  'Epochs',nepoch
     .     ,'  X-file interval ',interval
     .     ,'  Length of session (hrs) ',dble(ispan)/3600. 
      call report_stat('STATUS','MAKEX','makex',' ',message,0)  
      call secsum(iwknstart,sowstart,dble(ispan),iwknstop, sowstop)

c     translate starting epoch for information file and screen  
      call wtime (uinfor,iwknstart,sowstart,'GPST','Starting epoch ')
      call wtime (uinfor,iwknstop ,sowstop ,'GPST','Stopping epoch ') 
                                                                                  

c     Get the list of RINEX files that have data within the requested span
c     ****************************************************************************
                               
c     search ranges in days set in get_rxfiles
              
      call get_rxfiles( debug,site,iwknstart,sowstart,iwknstop
     .                , sowstop,rx_doy_minus,rx_doy_plus,infiles ) 

c     skip this site if no data found on the RINEX  file
      if( infiles(1)(1:1) == ' ') goto 4000
      
        
c     Data found on at least one RINEX file:  Open the X-file
c     ***************************************************************

      call openf (uxfile,fxfile,
     . 'unknown','formatted','sequential',qxfile )


c     Open the first RINEX file and get the xfile header information
c     **************************************************************

      infilcount = 1      
      frinex = infiles(1)                                                                  
      call openf (urinex,frinex,
     .             'old','formatted','sequential',qrinex )
      if(debug) print *,'MAKEX calling RHEAD '
      call rhead( debug,gnss,rcvrsw,swver,rxver,rxpgm,rxtime
     .          , nobtyp,rxobtyp,nwave1,nwave2,ircint
     .          , nxhdlin,xhdlin )               

      if(debug)  print *
     .    ,'MAKEX aft 1st rhead rxver nwave1 nwave2 nobtyp rxobtyp'
     .       ,                 rxver,nwave1,nwave2,nobtyp,rxobtyp


c     Determine the data types available and which ones should be used on the x-file
c     ******************************************************************************
c MOD TAH 200511: Passed gnsslf into sel_obtyp                             
      call sel_obtyp(gnsslf,nobtyp,rxobtyp,iobtypx)
c Changed by MAF following RWK's advice in 2020-05-04 email
c     write(uinfor,'(a,4(1x,a3))') 'Observation types used: '
c    .        ,(rxobtyp(iobtypx(i)),i=1,4)
      write(uinfor,'(a)')  'Observation types used: '
      write(uinfor,'(a,2(1x,a3))')  
     .      'Higher frequency: ',rxobtyp(iobtypx(1)),rxobtyp(iobtypx(3))
* MOD TAH 200526: Fixed format (missing 2 repeat).
      if(iobtypx(2).ne.0) write(uinfor,'(a,2(1x,a3))') 
     .      'Lower frequency : ',rxobtyp(iobtypx(2)),rxobtyp(iobtypx(4))
      if( iobtypx(5).ne.0.or.iobtypx(6).ne.0) then
        write(uinfor,'(a,2(1x,a3))') 
     .       'Alternate observation types available: '
     .       ,rxobtyp(iobtypx(5)),rxobtyp(iobtypx(6))
      endif
      if(debug) print *
     .      ,'MAKEX aft sel_obtyp gnss nobtyp rxobtyp iobtypx '
     .      ,gnss,nobtyp,rxobtyp,iobtypx
c     set iflag for rrinex call
      iflag = 1 

      call settyp ( debug,rxver,nobtyp,rxobtyp,iobtypx
     .            , nwave2,nsats,rcvrsw,swver,ndat,dattyp,lambda )
      if(debug) print *,'Aft settyp nwave2 dattyp lambda ',nwave2,dattyp
                                     
      
c     Read the the station.info file to determine occupation data
c     *****************************************************************
      
c     Find the requested receiver-track 

      exptyp = 'STATIC   '
      sitcod = site
* MOD TAH 200205: Added antenna azimuth to call. (antdaz in model.h)
      call rstnfo( usited
     .        , sitcod, iyrstart, idoystart, isodstart, ispan
     .        , stanam,  anthgt, offstn, offste, antdaz, rcvcod, antcod
     .        , hgtcod, radome_in, swver_stnfo, rcvers, rcvrsn, antsn
     .        , kstarts, kstops )  

      if( nint(100*swver).ne.nint(100*swver_stnfo) ) then  
        write(uinfor,200) swver,swver_stnfo
  200   format(1x,'WARNING: MAKEX input software version (',f5.2
     .        ,' ) different from station.info version (',f5.2
     .        ,' )',/,3x,'MAKEX value written on X-file header')
        write(message,205) swver,swver_stnfo
  205   format('Input software version (',f5.2
     .        ,') different from station.info version (',f5.2
     .        ,'). MAKEX value written on X-file header')
        call report_stat('WARNING','MAKEX','makex',' ',message,0)  
      endif
cd     print *,usited,icall,sitcod
cd    .    ,   iyrstart,idoystart,isessn
cd    .    ,   stanam,anthgt,offstn,offste,rcvcod,antcod,hgtcod
cd    .    ,   swver_stnfo,kstarts,kstops
                         
                
      if( debug ) write(uscren,'(a,f8.3,15i8)' ) 
     .   'Initial anthgt,itimes,kstarts,kstops'
     .   ,        anthgt,itimes,kstarts,kstops

c     Convert the raw antenna offst position to ARP-above-mark
      call hisub( uanthi,anthgt,offstn,offste,antcod,hgtcod,site
     .           , iyrstart,idoystart,isessn,offarp,antwarnings )

c     Get the PCN code from rcvant.dat indicating what differential code biases are to be applied

      call read_rcvant( 1,2,char6,char20,char5,rcvcod,char20,pcncod )
  

c     Write the X-file header
c     ***********************
            
                                
      call wxhead (uxfile,stanam,coords,offarp,
c**     2           l1z,l1n,l1e,l2z,l2n,l2e,
     3           gnss,nsats,sats,ndat,dattyp,lambda,
     .           rxobtyp,iobtypx,
     4           iwknstart,sowstart,interval,
     5           nepoch,isessn,version40,irunt,xhdlin,nxhdlin,
     6           uname,rcvrsw,swver,rcvcod,ircint,exptyp,antcod,
     7           dcb_override ) 


c     Initialize flags and arrays
c     ***************************

      timeset =   .false.
      late_epoch = .false.
      epochok   = .false.
      pastend_epoch = .false.
      newfile = .false.
      ant_event = .false. 
      iread = 0
      epoch = 1  
      if( ircint.ne.0.and.ircint.gt.interval ) then
        ideci = ircint/interval
      else
        ideci = 1   
      endif  
      ircint_msg = .false.
      ibadtag = 0

c     zero a few things for safety
      do i=1,maxchn
         prgl1(i) = 0.0d0
         prgl2(i) = 0.0d0
         dofl1(i) = 0.0d0
         dofl2(i) = 0.0d0
         gpssec(i) = 0.0d0
         disnr(i,1) = 0.0d0
         disnr(i,2) = 0.0d0
      enddo

c     initialize counts
      do j=1,maxchn+1
         ifull(j) = 0 
      enddo 
      do j=1,maxsat
         do i=1,maxepc
            obsmat(j,i) = 1
         enddo
      enddo


c     Setup the time span
c     *******************

c     define a window +/- slop second long, depending on receiver type
      call secsum
     .   (iwknstart,sowstart,-slop,iwknopen,sowopen)
      call secsum
     .   (iwknstart,sowstart,+slop,iwknclose,sowclose)
      write(message,'(a3,f6.2,a,f5.3,a)') rcvrsw,swver
     .             ,': accept data within +-',slop,'s of nominal epochs'
      write(uinfor,'(a)') message
      call report_stat('STATUS','MAKEX','makex',' ',message,0)
      if( debug ) then
        call wtime (uscren,iwknopen,sowopen,'GPST','Opening at start')
        call wtime (uscren,iwknclose,sowclose,'GPST','Closing at start')
        call wtime (uinfor,iwknopen,sowopen,'GPST','Opening at start')
        call wtime (uinfor,iwknclose,sowclose,'GPST','Closing at start')
      endif

c     Set the initial K-file epoch at the first X-file epoch for nav-files,
c     but force to an even 15-minute spacing for SP3 files since these 
c     need to match exactly
      iwknk= iwknstart
      sowk = sowstart              
cd      print *,'orig sowk ',sowk                       
      if( qsp3 ) then
        nsowk = nint(sowk/900.d0)
        sowk = dfloat(nsowk)*900.d0
      endif                    
cd      print *,'mod sowk ',sowk

c     Give the file reading routine the week number

      iwknfile = iwknstart
      weekset = .true.

c          Begin the Loop that Reads to the End of the File
c----------************************************************--------------

  400 continue
c     ---come back here when from several points whenever a new record needed
c        The data file is open: read it.
c          iflag defines what flavor data it is
c           -1 is a big problem
c            0 is documentation
c            1 is phase data
c            2 is an ephemeris
c
      if (debug) then
        call wtime (uscren,iwknopen, sowopen,'GPST', 'Opening now     ')
        call wtime (uscren,iwknclose,sowclose,'GPST','Closing now     ')
      endif

c       Skip the read if the last epoch was late compared with the window

      if ( late_epoch ) goto 580


c     Read a RINEX
c     **********************

      if (qrinex) then
                          
c        this used to check for bogus times on RINEX files
         epoch1 = epoch -1   
         call rrinex ( debug,iflag,rxver,gnss,nobtyp,rxobtyp,iobtypx
     .               , nprn,isvid,rxtime,igpswk,gpssec(1),epoch1
     .               , dofl1,dofl2,prgl1,prgl2,illi,issi
     .               , anth,ante,antn
     .               , fend,ferr )             
         if( debug ) then
            do i=1,nprn
              print *,'i illi ',i,(illi(i,j),j=1,4) 
            enddo
         endif
c        Set a flag if an antenna event record encountered
         if( anth.ne.0.d0 .or. ante.ne.0.d0 .or. antn.ne.0d0 ) 
     .       ant_event = .true.
c        Fix for bad TI ROM translator - round to even second
         if( rxpgm.eq.'TISTRX 10/90       ') gpssec(1)=dnint(gpssec(1))
c        Fix for bad MiniMac translator - 1 ms offset not in time tag
c          swver 1.89 gets changed to swver 1.59 in XTORX since the time
c          tag will have been fixed
         if( rcvrsw.eq.'MIN' .and. nint(100.d0*swver).eq.189 )
     .        gpssec(1) = gpssec(1) + 0.001d0
      else
         call report_stat('WARNING','MAKEX','makex',' ',
     .    'No input file--must have rinex = 1 in batch file',0)
           write(uinfor,'(a)') 'MAKEX: No RINEX input file  '
      endif
                                                               

c     See if an error is encountered reading the file 
c     ************************************************
                                                                                 
      if (ferr) then 
         write (message,'(a,i5)') 'Error in RINEX file at epoch '
     .                          ,  epoch
         write(uinfor,'(a)') message
         call report_stat('WARNING','MAKEX','makex',' ',message,0)
      endif
     

c     See if end-of-file and more data needed
c     **************************************     

c     if the file ends before the epochs are filled, open the next file
c     --if no more in list, goto 4000 to finish writing the X-file with zeroes
c     don't bother if within 10 epochs of the end of the requested span  
      if (fend) then
         if(epoch.gt.(nepoch-10)) goto 4000
         infilcount = infilcount + 1
         if( infilcount.gt.maxfil ) then
            call report_stat('WARNING','MAKEX','makex',' '
     .        ,'More than MAXFIL RINEX files requested',0) 
            goto 4000
         endif 
         frinex = infiles(infilcount)
         if( frinex(1:1).eq.' ' ) goto 4000                                                                  
         call openf (urinex,frinex,
     .               'old','formatted','sequential',qrinex )
         fend = .false.  
c        skip over the header and then go read the data record
         call rhead( debug,gnss,rcvrsw,swver,rxver,rxpgm,rxtime
     .             , nobtyp,rxobtyp,nwave1,nwave2,ircint
     .             , nxhdlin,xhdlin )
         if( debug ) 
     .  print *,'MAKEX aft 1st rhead rxver nwave1 nwave2 nobtyp rxobtyp'
     .       ,                 rxver,nwave1,nwave2,nobtyp,rxobtyp
         newfile = .true. 
c        read another data record
         goto 400
      endif        
      

c     Report the observation record read
c     **********************************

      iread = iread+1

c     Only consider this time tag if it is in the scenario
c     and less than a day after the last good time tag.

      call tagchk ( timeset,weekset,ibadtag,epoch
     .             , iwkntag,sowtag,igpswk(1),gpssec(1)
     .             , iwknstart )
      if( .not.timeset .or. .not.weekset ) epochok = .false.
      if ( iread .lt. 5) then
        call wtime (uinfor,iwkntag,sowtag,'GPST','Read epoch ')
        if (debug ) call wtime (uscren,iwkntag,sowtag,'GPST'
     .                             ,'Read epoch ')
      endif        
      itflag = 4
      call timcon( itflag,iwkntag,sowtag,iyrtag,idoytag,ihr,min,dsecs
     .              , utcoff )    
c     get the time seconds of day and y/d/h/m/s for checking station.info entries
      itimes(1) = iyrtag
      itimes(2) = idoytag
      itimes(3) = ihr
      itimes(4) = min          
      itimes(5) = int(dsecs)   
      isodtag = ihr*3600 + min*60 + int(dsecs)  
      sodtag = dble(isodtag)

 580  continue
        

c     Test time tag on epoch
c     **********************

c         Correct the nominal time tag by the pseudorange to avoid
c         problem with nominal tags running several tenth of a second
c         off true time.                             
      deltatag = (prgl1(1)-26000.d3)/ltvel  
c     if GLONASS, remove the leap seconds for checking the window: no, not now needed
c*      if( gnss.eq.'R') then 
c*        call monday(idoytag,month,iday,iyrtag)
c*        jdobs = julday(month,iday,iyrtag)
c*        utcoff = taiutc(jdobs) - 19.d0       
c*        deltatag = (prgl1(1)-utcoff*ltvel-26000.d3)/ltvel
c*      endif
      call secsum(iwkntag,sowtag,-deltatag,iwknchck,sowchck)  
      if( debug ) write(uscren,*) 'Window correction: ',
     . 'sowtag prgl1 deltatag sowchck:',sowtag,prgl1(1),deltatag,sowchck
      
      if (secdif(iwknopen,sowopen,iwknchck,sowchck).ge.0.d0) then
c          too soon, go read another record
        if (debug) call wtime (6,iwknchck,sowchck,'GPST','Too early   ')
        late_epoch = .false.
        goto 400
      elseif
     . (secdif(iwknclose,sowclose,iwknchck,sowchck).le.0.d0 .and.
     .  secdif(iwknclose,sowclose,iwknchck,sowchck).gt.-7.0d5) then
c          too late --set a flag to skip reads until the window catches up
       if (debug) call wtime (6,iwknchck,sowchck,'GPST','Too late    ')
       epochok =.false.
       late_epoch = .true.
c         if the record read is past the end of the scenario, set flag to avoid reading station.info
      if(secdif(iwknchck,sowchck,iwknstop,sowstop).gt.0.d0)  
     .        pastend_epoch = .true.
      else if (secdif(iwknchck, sowchck,
     .                iwknopen,sowopen) .gt. 0.d0 .and.
     .         secdif(iwknclose,sowclose,
     .                iwknchck, sowchck ) .gt. 0.d0) then
c         within the window: good epoch and the clock is presumably OK 
c         but if input data sampling is less dense than x-file output, skip the epoch
c         (this added to avoid confusing autcln, solve, and cview if data are unexpectedly
c         in the expected gaps--should not occur except for short periods within a RINEX file)
        if( mod((epoch-1),ideci).ne.0 ) then
           epochok = .false. 
           if( .not.ircint_msg ) then 
             write(message,'(2a,i4,a)') ' RINEX has data '
     .       ,'more frequently than header sampling (ircint=',ircint,')'
              call report_stat('WARNING','MAKEX','makex',' ',message,0)   
              write(uinfor,'(a)')  message  
              write(message,'(a,i5,a)') 
     .         ' First occurs at epoch ',epoch,';  message not repeated'     
              call report_stat('WARNING','MAKEX','makex',' ',message,0)   
              write(uinfor,'(a)')  message    
              ircint_msg = .true. 
            endif
        else
          timeset = .true.
          epochok = .true.
          late_epoch = .false.
          call timcon( itflag,iwkntag,sowtag,iyrtag,idoytag,ihr,min
     .                , dsecs,utcoff ) 
        endif
      else
         write (message,*) ' Warning, bad time tag: ',iwkntag,sowtag
         write (uinfor,*) ' Warning, bad time tag: ',iwkntag,sowtag 
         call report_stat('WARNING','MAKEX','makex',' ',message,0)
         if (debug) call wtime(uscren,iwkntag,sowtag,'GPST'
     .                        ,'Epoch mismatch')
         epochok = .false.
      endif
                                                   


c     See if new antenna offsets are indicated
c     ****************************************         

c       Key by   1) New RINEX file opened
c                2) Event flag encountered in RINEX file
c                3) Current time more than 1 minuate past end-time from station.info entry
                    
c     this logical variable can go away when old-style station.info is removed 
      newant = .false.
      if( .not.pastend_epoch .and. ( newfile .or. ant_event .or.
     .  ( itimdif(itimes,kstops).gt.60 ) ) ) then  
        newant = .true.
        sitcod = sitecd  
        if( debug ) write(uscren,'(a)') 'Calling rstnfo for new file'
* MOD TAH 200203: Added AntDAZ to list of values from station.info
        call rstnfo( usited
     .         , sitcod, iyrtag, idoytag, isodtag, 0
     .         , stanam,  anthgt, offstn, offste, antdaz,rcvcod, antcod
     .         , hgtcod, radome_in, swver_stnfo, rcvers, rcvrsn, antsn
     .         , kstarts, kstops ) 
          if( debug ) write(uscren,'(a,f8.3,3(2x,5i5))' )  
     .      'New anthgt,itimes,kstarts,kstops'
     .        , anthgt,itimes,kstarts,kstops  
      endif
      if( newant ) then    
c       convert the raw antenna offsets to ARP-above-mark
        call hisub( uanthi,anthgt,offstn,offste,antcod,hgtcod,site
     .            , iyrstart,idoystart,isessn,offarp,antwarnings )
c       reset the newfile flag
        newfile = .false.
      endif
                                    


c     Update coordinates and antenna offsets for kinematic or pseudokinematic survey
c     ******************************************************************************

        if( exptyp .ne.'STATIC   ' ) then

          sitecd = sitcod   
          site_epoch = decyrs( iyrtag,idoytag,sodtag ) 
          call lread ( site,site_epoch )
        endif


c     Write a record on the X-file and K-file
c     ***************************************

      number = 0
      if(epochok) then
c       count how many good ones we've got:
c       ISGOOD checks for adequate SNR and reasonableness
c       If not good enough, then iyuk is returned true.
        ngood_chan = 0                  
cd       print *,'MAKEX nsats sats ',nsats,(sats(i),i=1,nsats)
cd        print *,'MAKEX write x nprn  isvid '
cd     .     ,nprn,(isvid(i),i=1,nprn)
        do j=1,nprn
c         is the satellite in the scenario ?                            
cd          print *,'j isvid inscen '
cd     .         ,j,isvid(j),inscen(isvid(j),nsats,sats)
          if (inscen(isvid(j),nsats,sats).gt.0) then  
            ngood_chan = ngood_chan+1   
            if(debug) write(*,*) 'j ndat isgood ',j,ndat
     .               ,isgood(ndat,dattyp,dofl1(j),issi(j,1),dofl2(j)
     .              , issi(j,2),prgl1(j),prgl2(j),iyuk)
            if (isgood(ndat,dattyp,dofl1(j),issi(j,1),dofl2(j)
     .           , issi(j,2),prgl1(j),prgl2(j),iyuk) ) then
              number = number + 1
            else
              if(iyuk) then
                ireject = ireject + 1
                if( ireject.lt.1000 ) then
                  if ( debug ) write(uscren,585) 
     .                        epoch,isvid(j),dofl1(j),issi(j,1)
     .                      , dofl2(j),issi(j,2),prgl1(j),prgl2(j)
                  write(uinfor,585) epoch,isvid(j),dofl1(j),issi(j,1)
     .                        , dofl2(j),issi(j,2),prgl1(j),prgl2(j)
  585             format(1x,'Data rejected  (epoch,PRN) : ',
     .                     i4,i3,' L1 ',1pe8.1,
     .                     i2,   ' L2 ',1pe8.1,
     .                     i2,   ' P1 ',1pe8.1,
     .                           ' P2 ',1pe8.1)
                endif
              endif
            endif
          else
            if (debug .and. isvid(j) .ne. 0) then
             write (uinfor,*) 'PRN ',isvid(j),' not in scenario.'
            endif
          endif        
c       enddo on loop over channels checking whether in scenario
        enddo  
c     endif on epochok 
      endif

c     Write the time line
c     -------------------            
      if( .not.epochok. or. ngood_chan.le.0 ) then
c       no valid data at this epoch - write an empty epoch
        if (debug) print *,'Empty epoch',epoch
        write(uxfile,592) epoch, 0
        ifull(1) = ifull(1)+1
      else 
        itflag = 4
        call timcon(itflag,iwkntag,sowtag,iyr,jdoy,ihr,min,dsecs,utcoff)
        call xyz2sph(coords,latr,lonr,radius)
        call raddms( latr,latflag,dlat,mlat,seclat )
        if( latflag.eq.'-') then
          latflag = 'S'
        else
         latflag = 'N'
        endif
        call raddms( lonr,lonflag,dlon,mlon,seclon )
        if( lonflag.eq.'-' ) then
          lonflag = 'W'
        else
          lonflag = 'E'
        endif             
c       this now dummy
        kflag = 0 
        write(uxfile,592) epoch,number,iyr,jdoy,ihr,min,dsecs
     .                  , kflag,upperc(sitecd),latflag,dlat,mlat,seclat
     .                  , lonflag,dlon,mlon,seclon,radius
     .                  , offarp
  592    format(/,2I4,I5,I4,2I3,F11.7,1x,i2,1x,a4,1x,
     .          a1,i2,1x,i2,1x,f8.5,1x,a1,i3,1x,i2,1x,f8.5,f13.4,
     .          1x,3F8.4,3X,3F8.4)
        if (debug) then 
          write(*,*) 'MAKEX: nprn ',nprn
          do j = 1,2
            write (uscren,*) 'MAKEX: ISSIs:',(issi(k,j),k=1,nprn)
          enddo               
          write(*,*) 'MAKEX: nprn ',nprn
        endif

c     write a data line for each channel
c     ----------------------------------   
        do j=1,nprn
          if( inscen(isvid(j),nsats,sats).gt.0 ) then
            if ( isgood( ndat,dattyp,dofl1(j),issi(j,1),dofl2(j)
     .       , issi(j,2),prgl1(j),prgl2(j),iyuk) ) then
c            What is the number of the X-file channel?
             i = inscen(isvid(j),nsats,sats)
cd              print *,'ngood loop j isvid i ',j,isvid(j) 
c             Set the X-file error flag according to the RINEX ISSI
c             0    unknown, assumed good for now
c             1-2  low amplitude
c             3-9  good
c             NOTE that that the order here is significant.
              ier(j) = iggood
c              Do not flag issi 1 or 2 as marginal for Trimble/Ashtech codeless data
              if( (rcvrsw.eq.'TRM'.and.swver.le.550).or.
     .           (rcvrsw.eq.'ASH' .and. swver.le.800).or.
     .           ( rcvrsw. eq. 'TOP') ) then
                 continue
              else
               if ( issi(j,1) .eq. 1 .or. issi(j,2) .eq. 1 .or.
     .              issi(j,1) .eq. 2 .or. issi(j,2) .eq. 2) then
                 ier(j) = iglamp                                
              endif
            endif
c           For all receivers and zero amplitude
            if(issi(j,1) .eq. 0 .or. issi(j,2) .eq. 0) then
              ier(j) = iglamp
              if( debug)  write(*,*) 
     .         'zero amplitude j issi ',issi(j,1),issi(j,2)
            endif
c            Patch for JPL RINEX screw up
            if(issi(j,1) .eq. 0 .or. issi(j,2) .eq. 0 .and.
     .         rcvrsw.eq.'ROG') then
               ier(j) = 0
            endif
c             If the RINEX file indicates loss-of-lock, put a bias flag in the X-file
            if (illi(j,1) .eq. 1 ) then
              ier(j) = igbias
            else
              ier(j) = 0
            endif 
c             Get the (new) flag to indicate whether corrections needed for cross-correlating rcvrs
            call set_dcb_flag( sitcod,rcvcod, pcncod,rxver,illi(j,1)
     .                     , dcb_override,idcb(j))
c            Record the flag for gap testing
            obsmat(i,epoch) = ier(j)
c           Write a record to the X-file
            write(uxfile,594) idcb(j),ier(j),i,dofl1(j),issi(j,1)
     .            , dofl2(j),issi(j,2),prgl1(j),prgl2(j)
 594        format(8x,i1,1x,2i2,1x,d22.15,1x,i3,1x,d22.15,1x,i3,
     .              2x,d22.15,2x,d22.15)
c           endif on valid data
            endif
c         endif on in-scenario
          endif
c       loop on writing channels
        enddo 
c     endif on non-empty epoch
      endif  


c     Tell us what you have done with this epoch
c     ------------------------------------------
      if (debug) then
        write (buff40,'(a,i2)') 'Wrote. number = ',number
        call wtime (uscren,iwkntag,sowtag,'GPST',buff40)
      endif
      ifull(number+1) = ifull(number+1)+1

      if( debug ) then
        if (epoch .eq. 1 .or. epoch. eq. nepoch .or.
     .   mod(epoch,200).eq.0) then
         write(message,'(a,i5)') ' Wrote epoch ',epoch
         write (uscren,'(a)') message
         write (uinfor,'(a)') message
        endif
      endif
      

c     Write a K-file Record
c     *********************
                                       
      if ( secdif(iwkntag,sowtag,iwknk,sowk). ge. 0.d0
     .    .and. nprn.gt.0 .and. epochok) then
c       The requested epoch has been found.  Get the ephemeris and
c       clock information for the latest epoch before the requested
c       one, for each satellite.  Then compute a clock correction.
        do j=1,nprn
          if ( inscen(isvid(j),nsats,sats).le.0 ) then
            continue 
          else
            if ( isgood( ndat,dattyp,dofl1(j),issi(j,1),dofl2(j)
     .         , issi(j,2),prgl1(j),prgl2(j),iyuk) ) then
             icall = 1 
             if(debug) print *,'MAKEX getting ephem and clock,'
     .             , 'iprn qsp3 qnav ',isvid(j),qsp3,qnav
             if( qsp3 ) then
               call getsp3( debug,icall,iwkntag,sowtag,gnss,isvid(j)
     .                    , xsat,svclock(1,j),satok )
c              set this to omit the nav-file epoch entries in wclock
               iwkne = 0 
             else
               call getnav ( debug,icall,iwkntag,sowtag,gnss,isvid(j)
     .                     , nsats,sats 
     .                     , xsat,svclock(1,j),satok,iwkne,sowe )  
               if(debug) print *,'SV clock from nav file: ',svclock(1,j)
             endif  
             call timcon(4,iwkntag,sowtag,iyr,jdoy,ihr,min,dsecs,utcoff)
             call monday(jdoy,month,iday,iyr)
             jdk = julday(month,iday,iyr)              
             sodk = 3600.d0*ihr + 60.d0*min + dsecs
             call readj(usvclk,sats,nsats,isvid(j),jdk,sodk,svclock(1,j)
     .                 ,icall,jdtoc,toc,svepc,svcrat,svcacc,valid )
             if(debug) print *,'SV clock from j-file: ',svclock(1,j)
             if (satok) then
               call stnclk (debug,gnss,iwkntag,sowtag,isvid(j),prgl1(j)
     .                     ,coords,xsat,svclock(1,j)
     .                     ,rclock,iserr,clkwarnings )          
               if( debug ) print *
     .            ,'MAKEX gnss isvid iserr rclock iserr'
     .           ,     gnss,isvid(j),iserr,rclock,iserr
               if (iserr .eq. 0)
     .           call wclock ( uclock,afmt,gnss,isvid(j),iwkntag,sowtag
     .                       , prgl1(j),svclock(1,j),rclock,iwkne,sowe )
             endif
            endif 
          endif
        enddo
c         Update the K-file epoch to the next "even" epoch to have the best
c         chance of matching the (usually) hourly epochs of the broadcast ephemeris
        call secsum(iwknk,sowk,deltak,iwknk,sowk)
        sowk= deltak*dint( sowk/deltak )
      endif
         

c     Update the Window Opening and Closing
c     *************************************

      call secsum
     .  (iwknopen,sowopen,dble(interval),iwknopen,sowopen)
      call secsum
     .  (iwknclose,sowclose,dble(interval),iwknclose,sowclose)

      epoch = epoch +1
cd DEBUG print stop
cd      if ( epoch.gt.54 ) stop 
      if(epoch.gt.nepoch) go to 4100
c     loop back for another record
      go to 400               
c-----end of of loop over data records
c**************************************************************

C     come here at end of input data file
 4000 continue
c     fill remaining epochs with zeros
      number = 0
      do  i=epoch, nepoch
         write(uxfile,592)i,number
         ifull(1)=ifull(1)+1
      enddo
              

c     Write screen and info file messages
c     ************************************

       call report_stat('STATUS','MAKEX','makex',' ',
     .'End of file encountered in reading RINEX input file.',0)
 4001 format(1x,'MAKEX: End of file encountered in reading RINEX'
     .        ,' input data file.')
      write(uinfor,4001)
      go to 4105
   
 4100 write(message,'(a)') ' Wrote all the epochs requested'   
      call report_stat('STATUS','MAKEX','makex',' ',message,0)
      write(uinfor,'(a)') message

 4105 continue
c     print out what you did
      if( debug ) write (uscren,4107)
      write (uinfor,4107)
 4107 format(1x,/,2x,'Sats  Epochs ',/,2x,'----  -----')
      do i=1,nprn+1
        idots = 60*ifull(i)/nepoch
        do j=1,60
          if (j .lt. idots) then
            dots(j) = '!'
          else
            dots(j) = ' '
          endif
        enddo
        write(message,'(1x,i5,1x,i5,1x,60a)')  
     .            i-1,ifull(i),(dots(j),j=1,60)
        write(uinfor,'(a)') message
        if(debug) write(uscren,'(a)') message
      enddo


c     Do gap analysis
c     ***************

      itotal = 0
      do i = 1, nsats
        gap = .false.
        ibad = 0
        igap = 0
        maxgap = 0
        do itimel = 1, nepoch-1
          if( obsmat(i,itimel).eq.0 .and.
     .        obsmat(i,itimel+1).ne.0) then
             gap = .true.
          endif
          if( obsmat(i,itimel).ne.0) then
             ibad = ibad+1
             if( gap ) igap = igap+1
          endif
          if( obsmat(i,itimel).ne.0 .and.
     .        obsmat(i,itimel+1).eq.0) then
            maxgap = max(maxgap,igap)
            gap = .false.
            igap = 0
          endif
        enddo
        if( obsmat(i,nepoch).ne. 0) ibad = ibad+1
        ioktot(i,1) = itimel-ibad
        maxtot(i,1) = maxgap
        itotal = itotal + ioktot(i,1)
      enddo
      write(uinfor,'(A)') ' '
      if( clkwarnings ) then
        call report_stat('WARNING','MAKEX','makex',' ',
     .  '**Warnings issued in .infor file for bad clocks',0)
      endif
      write(uinfor,'(A)')
     . ' Good observations per channel: Total Number and Maximum Gap'
      write(uinfor,905) 'PRNs ',(sats(I),I=1,nsats)
      write(uinfor,'(A)') ' '
      write(uinfor,905) 'OBS  ',(ioktot(i,1),i=1,nsats)
      write(uinfor,905) 'GAPS ',(maxtot(i,1),i=1,nsats)  
* MOD TAH 200618: Updated 32I to 50I to allow for 35 Beidou satellites
  905 format(1x,a6,1x,50i5)
      write(message,907) itotal,ireject 
      write(uinfor,907) itotal,ireject 
  907 format(i7,' observations written to xfile ',
     .       i5,' observations rejected as unreasonable')
      call report_stat('STATUS','MAKEX','makex',' ',message,0) 
      write(uinfor,'(A)') ' '
      if( itotal.le.0 ) then
        write(message,920)
  920   format('No observations written to X-file: '
     .     ,'Check requested and read times for timetag mismatch '
     .     ,'Use RXSCAN to check for observations '
     .     ,'within the requested span') 
        write(uinfor,'(a)') message
        call report_stat('WARNING','MAKEX','makex',' ',message,0)
      endif
      write (uinfor,1000) upperc(site),iyr,idoystart,isessn
 1000 format (1x,'END PROCESSING: ',a4,1x,i4,1x,i3,1x,i2/,
     .        1x,'..........................',//)
      write (message,1001) upperc(site),iyr,idoystart,isessn
 1001 format ('End processing: ',a4,1x,i4,1x,i3,1x,i2)
      call report_stat('STATUS','MAKEX','makex',' ',message,0)
      call closem     

      enddo
c-----End of loop on sites---------------------------------------------

      call report_stat('STATUS','MAKEX','makex',' ',
     .  'Normal End of MAKEX',0)
      stop
      end

