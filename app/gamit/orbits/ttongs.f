Copyright (c) Massachusetts Institute of Technology and the University of
California, 1994,2000.  All rights reserved.

      Program ttongs   

C Write an NGS Standard Product #1 or #3 emphemeris file (earth-fixed)
C from a GAMIT T-file (inertial or earth-fixed)

C Written by R. King November 1991 from Program ngstot  King & Bock  March 1988
c Options added for output SP3 format by P. Fang March 1993
c Option for command-line arguments and clock-file by R. King March 2000   
c Modified to write sp3-C by R King September 2009
c Modified to accept a sample_rate value TAH 210122: Passed in runstring.

      implicit none

      include '../includes/dimpar.h'
                                                          
      integer*4 iut,iungs,iuclkfile,iscrn,iprnt,idir
     .        , jdb,jdstp,jde,jdn,jdn1,jdf,jdf1,iyr,iday
     .        , iweek,idow,icrd,julday,id,im,nstart,nstop
     .        , iyb,iyf,idoyb,idoyf,ihrb,ihrf,iminb,iminf
     .        , nsat,itsat(maxsat),nics,nintrs,notl
     .        , numsv,mjdn,nepcht,nepchn,ioerr,iarg,iclarg
c     .       , iudbug

* MOD TAH 210122: Introduced higher sampling rate factor (e.g. for 900-sec sampled
*     t-file to output at 300-seconds, sample_rate = 3
      integer*4 sample_rate   ! Sampling rate multiplier on t-file rate;
                              ! Value must >= 1; default is 1.

      real*8 tb,tstp,te,tn,tf,tn1,tf1,satics(maxorb,maxsat),fmjdn,sdelt
     .     , dt,timdif,docmc(3)
                 
      character*1 ans,gnss
      character*3 orbtyp,spfmt
      character*4 org,icsnam(maxorb),string
      character*5 crdsys,precmod,nutmod,gravmod,frame,srpmod
     .           , eradmod,antradmod
      character*8 otlmod
      character*10 pcvmod
      character*16 spfile,tfile,tfilef,satnam(maxsat)
* MOD TAH 030210: Increased length of clock file name (so that it
*     can be located anywhere).
      character*256 clkfile
      character*16 lowerc,reply
      character*80 version,buff80
      character*256 message
               
c      logical debug
c      data debug/.false./

c Get version number

      call oversn(version)
      write(6,'(a)')' '
      write(message,'(a,a80)') 'Started TTONGS',version
      call report_stat('STATUS','TTONGS','orbits/ttongs',' ',message,0)

c Set default orbtyp, org, and coordinate system

      orbtyp = 'FIT'
      org    = 'SIO'
      crdsys   = 'ITR97'
      spfile = ' '
      clkfile = ' '    
      otlmod = ' '
      pcvmod = ' '                                             
c** temporary set PCV and O-tide models
c**      pcvmod = "IGS05" 
c**      otlmod = "FES2004 "
    

c Units: 

c   ** Be careful, subroutine TROT uses 5,6, 14-17, 20-22, and 30
C

c      T-file (unit applies here to both the input inertial and rotated
c        Earth-fixed T-file--TROT opens and closes the latter as unit 17)
         iut = 17
c      Clock file (set equal zero if not clock file)
         iuclkfile = 18  
c     SP#1 or SP#3 file
         iungs = 19
c      Screen
         iscrn = 6
c      Print output (set zero here for no print)
         iprnt = 25
c      Debug file
c         iudbug = 26                                


c Open the print file 

        open(unit=iprnt,file='ttongs.out',status='unknown'
     .      ,form='formatted', iostat=ioerr)  
        write(iprnt,'(a,a80,/)') 'Started TTONGS',version


c Read the Input 

c     Required (if missing, delineate by ' ' )
c        T-file name    
c        Output format (SP1 or SP3)
c        Organization  (e.g. SIO)       
c        Type of orbit (FIT, EXT, BRD) 
c        Reference frame for Earth-fixed orbit (5 charcters, e.g. ITR97) 
c           Next argument added 061101 (check for confusion with sp3-file name)
c        Ocean loading model (needed to correct CM to CE) (e.g FES2004)
c        PCV model (e.g. IGS05)
c        Output file name (e.g. sio10462.sp3)  
c     Optional ' 
c        Name of clock file (e.g. sio10462.clk)  (delineate by ' ' if time arguments present)
c        Start/Stop times  yyyy doy hh mm    yyyy doy hh mm
         
c Check for command-line arguments (if missing, write the help)

      iarg = iclarg(1,tfile)   
      if( iarg.le.0 ) then              
        write(6,'(1x)') 
        write(6,'(a)') 
     .     'Command-line arguments: '
        write(6,'(a)') 
     .' ttongs [tfile] [fmt] [org] [orbtyp] [frame] [otlmod] [pcvmod]'
     .      ,'[sp-file] [clkfile] [start] [stop] <sample rate>'
        write(6,'(a)') '  where start stop are YYYY DOY HH MM '
        write(6,'(a)') '        <sample rate> to increase rate'
        write(6,'(a)') '  clkfile and start, stop are optional'   
        write(6,'(a)') '  Default sample rate = 1 (optional; int >=1)'
        write(6,'(/,3a)')  'e.g. ',
     .  ' ttongs tpgga9.186 sp3 sio fit itr97 fes2004 IGS05_1402 '
     .       ,'SIO10462.sp3 SIO10462.clk 1999 186 0 0 1999 186 23 45 1'
        write(6,'(1x)') 
        stop
      endif


c Read first the T-file name so that we can get the date from the header

      if( tfile(1:1).ne.'t' ) call report_stat('FATAL'
     .  ,'TTONGS','orbits/ttongs',' ','Invalid entry for T-file name',0)
      call topens(tfile,'old',iut,ioerr)
      if( ioerr.ne.0 ) then
        call report_stat('FATAL','TTONGS','orbits/ttongs',tfile
     .         ,'Error opening T-file, or T-file not found: ',ioerr)
      else
        call report_stat('STATUS','TTONGS','orbits/ttongs',tfile
     .           ,'Opened T-file: ',0)
      endif

c Read the T-file header to get times and satellites
 
      frame = 'UNKWN'
      call thdred ( iut,iscrn,iprnt,nsat,gnss,itsat,satnam
     .            , jdb,tb,jdstp,tstp,sdelt,nepcht,jde,te
     .            , nics,satics,nintrs,icsnam
     .            , precmod,nutmod,gravmod,frame,srpmod
     .            , eradmod,antradmod )     
      close(iut)
c     convert the start time them into gps-week and day of week

c**   Do we want to use the start time or the IC epoch to get the sp3 file name?
c**   For 24-hr files, the T-file will usually start on the previous day.
      call dayjul( jdb, iyr, iday )
      call doygwk(iday,iyr,iweek,idow)   
                             
c Get the organization and file type

      iarg = iclarg(2,spfmt)  
      if( iarg.le.0 ) call report_stat('FATAL','TTONGS'
     .                    ,'orbits/ttongs',' '
     .           ,'Missing command-line argument for output format',0) 
      call uppers(spfmt)
      if( spfmt.ne.'SP1' .and. spfmt.ne.'SP3' ) 
     .        call report_stat('FATAL','TTONGS','orbits/ttongs',' '
     .                         ,'Invalid output format',0)
      call uppers(spfmt)
      iarg = iclarg(3,org)
      if(iarg.le.0) call report_stat('FATAL','TTONGS','orbits/ttongs'
     .      ,' ','Missing command-line argument for organization',0)  
      call uppers(org)

c Create the output file name if not input

      iarg = iclarg(8,spfile)   
      if( iarg.le.0 ) call report_stat('FATAL','TTONGS'
     .                    ,'orbits/ttongs',' '
     .         ,'Missing command-line argument for OTL model',0)  
      if( spfile(1:1).eq.' ') then            
        write(spfile(1:3),'(a3)') org
        write(spfile(4:8),'(i4,i1)') iweek,idow
        if (spfile(4:4).eq." ") write(spfile(4:4),'("0")')
        write(spfile(9:12),'(".",a3)') spfmt
      endif
      spfile=lowerc(spfile) 
      open(unit=iungs,file=spfile,status='unknown',form='formatted'
     .    ,iostat=ioerr)
      if( ioerr.ne.0 ) call report_stat('FATAL','TTONGS','orbits/ttongs'
     .     ,spfile,'Error opening SP file ',ioerr)

                 
c Get the type of orbit, terrestrial frame, and OTL model 

      iarg = iclarg(4,orbtyp) 
      if(iarg.le.0) call report_stat('FATAL','TTONGS','orbits/ttongs'
     .        ,' ','Missing command-line argument for orbit type',0) 
      call uppers(orbtyp)
      if( orbtyp.ne.'FIT' .and. orbtyp.ne.'BRD' .and. orbtyp.ne.'EXT')
     .     call report_stat('FATAL','TTONGS','orbits/ttongs',' '
     .                        ,'Invalid orbit type',0)   
      iarg = iclarg(5,crdsys)
      if(iarg.le.0) call report_stat('FATAL','TTONGS','orbits/ttongs'
     .     ,' ','Missing command-line argument for terrestial frame',0) 
       iarg = iclarg(6,otlmod)    
       if(iarg.le.0) call report_stat('FATAL','TTONGS','orbits/ttongs'
     .     ,' ','Missing command-line argument for OTL model',0) 
      if( otlmod(1:1).eq.' ') call report_stat('WARNING','TTONGS'
     .     ,'orbits/ttongs',' ','OTL model missing--assume NONE',0) 
      call uppers(otlmod)
      iarg = iclarg(7,pcvmod) 
      if(iarg.le.0) call report_stat('FATAL','TTONGS','orbits/ttongs'
     .     ,' ','Missing command-line argument for PCV model',0) 
      if( pcvmod(1:1).eq.' ') call report_stat('WARNING','TTONGS'
     .     ,'orbits/ttongs',' ','PCV model missing--assume NONE',0) 
      call uppers(pcvmod)

      
c Get the clock file if present

      iarg = iclarg(9,clkfile)   
      if( clkfile(1:1).ne.' ' ) then
        open(unit=iuclkfile,file=clkfile,status='old'
     .    ,form='formatted',iostat=ioerr)
        if(ioerr.ne.0) then
          call report_stat('FATAL','TTONGS','orbits/ttongs'
     .              ,clkfile,'Error opening clock file ',ioerr)         
        else
          call report_stat('STATUS','TTONGS','orbits/ttongs'
     .              ,clkfile,'Opened clock file: ',ioerr)  
          read(iuclkfile,'(a)') buff80
          rewind(iuclkfile)
        endif       
      else
         iuclkfile = 0    
         call report_stat('WARNING','TTONGS','orbits/ttongs',' '
     . ,'No clock file available, setting values = 999999.999999',ioerr)
      endif
                
c Get the start/stop times  (write entire file if omitted)

      iyb = 0   
      iarg = iclarg(10,string)     
      read(string,'(i4)',iostat=ioerr) iyb
      if( ioerr.ne.0 )  
     .      call report_stat('STATUS','TTONGS','orbits/ttongs',' '
     .           ,'Error reading year from command-line',ioerr)
* MOD TAH 210122: Set the default sampling rate to match t-file.
      sample_rate = 1
      if( iarg.gt.0 ) then
        iarg = iclarg(11,string) 
        read(string,'(i3)',iostat=ioerr) idoyb
        iarg = iclarg(12,string)   
        read(string,'(i2)',iostat=ioerr) ihrb
        iarg = iclarg(13,string) 
        read(string,'(i2)',iostat=ioerr) iminb    
        if( ioerr.ne.0 )  
     .    call report_stat('STATUS','TTONGS','orbits/ttongs',' '
     .           ,'Error reading start time from command-line',ioerr)
        iarg = iclarg(14,string) 
        read(string,'(i4)',iostat=ioerr) iyf
        iarg = iclarg(15,string)
        read(string,'(i3)',iostat=ioerr) idoyf
        iarg = iclarg(16,string)  
        read(string,'(i2)',iostat=ioerr) ihrf
        iarg = iclarg(17,string)    
        read(string,'(i2)',iostat=ioerr) iminf
* MOD TAH 210122: See if sample rate passed.
        iarg = iclarg(18,string) 
        if( iarg.gt.0 ) then ! Sample_rate passed   
           read(string,*,iostat=ioerr) sample_rate
        else
           sample_rate = 1
        endif

        if( ioerr.ne.0 )  
     .    call report_stat('STATUS','TTONGS','orbits/ttongs',' '
     .         ,'Error reading stop time from command-line',ioerr)
      endif    
       
       
c Determine the start and stop epochs for the SP3 file

c     first shrink the  output file by 5 epochs on each end in order to interpolate accurately
      jdn = jdb
      tn  = tb
      dt = sdelt*5.d0
      call timinc( jdn,tn,dt ) 
      nstart = 6  
      jdf = jdb
      tf = tb                  
c     extra -1 is because start is epoch 1 not epoch 0
      dt = sdelt*(nepcht -5 - 1)    
      jdf = jdb
      tf = tb
      call timinc( jdf,tf,dt )  
      nstop = nepcht - 5
 

c     now shrink it further if an output span was specified

      if( iyb.ne.0 ) then  
         call monday(idoyb,im,id,iyb)
         jdn1 = julday(im,id,iyb)
         tn1 = ihrb*3600.d0 + iminb*60.d0 
         dt = timdif(jdn1,tn1,jdn,tn) 
         if( dt.ge.-1.0d0 ) then
           jdn = jdn1
           tn = tn1  
           nstart = nstart + int(dt/sdelt) 
         else     
           call report_stat('WARNING','TTONGS','orbits/ttongs'
     .                     ,' ','Requested start time too early',0)
         endif   
         call monday(idoyf,im,id,iyf)
         jdf1 = julday(im,id,iyf)
         tf1 = ihrf*3600.d0 + iminf*60.d0    
         dt = timdif( jdf,tf,jdf1,tf1 )
         if( dt.gt.1.d0 ) then
           jdf = jdf1
           tf = tf1
           nstop = nstop - int(dt/sdelt)   
         else
           call report_stat('WARNING','TTONGS','orbits/ttongs'
     .                     ,' ','Requested stop time too late',0)
         endif       
      endif 


C Compute number of epochs and the Modified Julian date for the NGS file header
              
      nepchn = nstop - nstart + 1    
c     ( PEP JD = MJD + 2400000 + 1)
      mjdn = jdn - 2400000 - 1
      fmjdn= tn/86400.d0

c For now, set the number of satellites to be the same as the T-file

      numsv  = nsat

c Write the NGS-file header records

      if (spfmt.eq.'sp3'.or.spfmt.eq.'SP3') then
* MOD TAH 210122: Added sample_rate to call (only sp3).
         call wsp3hd( iungs,iprnt,mjdn,fmjdn,sdelt,nepchn
     .              , numsv,gnss,itsat
     .              , org,orbtyp,crdsys,pcvmod,otlmod, sample_rate )
       else
         read(crdsys(4:5),'(i2)') icrd
         call wsp1hd( iungs,iscrn,mjdn,fmjdn,sdelt,nepchn,numsv,itsat
     .           , org(1:3),orbtyp(1:1),icrd )
      endif
 

c If the T-file is inertial, rotate into the Earth-fixed frame

      if ( frame.ne.'EFIXD' ) then
c       trot will open and then close both files
        close ( iut)
        idir = -1       
c       earth-fixed T-file (opened int TROT)
        tfilef = tfile
        tfilef(6:6) = 'e'
        call trot(tfile,tfilef,idir,frame)
c       Need to open and read the header of the new (Earth-fixed) T-file
        open(unit=iut,file=tfilef,status='old',form='unformatted'
     .      ,iostat=ioerr)
        if (ioerr .ne. 0) then
         call report_stat('FATAL','TTONGS','orbits/ttongs',tfilef,
     .   'Error opening earth-fixed T-file: ',ioerr)
        endif
        frame = 'EFIXD'
        call thdred ( iut,iscrn,iprnt,nsat,gnss,itsat,satnam
     .            , jdb,tb,jdstp,tstp,sdelt,nepcht,jde,te
     .            , nics,satics,nintrs,icsnam
     .            , precmod,nutmod,gravmod,frame,srpmod
     .            , eradmod,antradmod )    
      endif	
		   

c Get the coefficients for ocean-loading CMC corrections
 
      if( otlmod(1:4).ne.'NONE' ) then             
c      coefficients saved in otlcmc; time arguments and docmc dummy this call
c      hard wire the # of components
       notl = 11
       call otlcmc( jde,te,otlmod,notl,1,docmc )
       write(message,'(a,a8)') 
     .   'Converting CM to CE using otlcmc.dat offsets for ',otlmod 
       call report_stat('STATUS','TTONGS','orbits/ttongs',' ',message,0)
      endif


c Read the T-file file and write the NGS-file data records
* MOD TAH 210122: Updated to pass the sample_rate into call.         
      call sdtrit( iut,iungs,iuclkfile,iprnt,spfile,nintrs
     .           , nsat,gnss,itsat,spfmt,otlmod
     .           , jdb,tb,sdelt,nepcht,nstart,nstop, sample_rate )

      call report_stat('STATUS','TTONGS','orbits/ttongs',' ',
     .'Normal end in TTONGS ',0)

      stop
      end
