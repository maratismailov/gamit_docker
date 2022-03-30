Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego. 1997.   All rights reserved.

      Program FIC2RX

c     Translate a single FICA file into RINEX observation and navigation files.  
c     This program is designed to be run interactively or from a shell script.
c     Coordinates may be updated from a GAMIT L-file or GLOBK apr file, and antenna 
c     information from station.info.   This program currently handles TI 4100 data 
c     collected with GESAR or CORE and MACROMETER data translated from ARGO . 
          
c     Command format:
c
c        fic2rx   [FICA-file] [L-file]   [J-file]  
             
c          If there are no arguments, the program will prompt for an input FICA file.
c          If there is only one argument, the coordinates are set to be zero. 
c          The J-file is required for MiniMac data in order to convert the pseudo-range
c            to RINEX convention; should be omitted for other receivers.
c
c     R. King  March, 1997; last modified by King 020807

      
c     A note on conventions.  
c           
c     FICA (and X-) files are in the DOPPLER CONVENTION
c
c             d(phase)        1        d(range)
c               -----   = - ------     -------
c               d(t)        wavelength  d(t)

c     RINEX is defined to be in the PSEUDORANGE CONVENTION
c
c             d(phase)        1        d(range)
c               -----   = + ------     -------
c               d(t)        wavelength  d(t)
c
c     RINEX convention also demands that the initial phase value be zero.


      implicit none              
                  
c     --maximum dimensions for satellites in session
       include '../includes/dimpar.h'  

c     --I/O units
      include '../includes/units.h' 

c     --file names and unit numbers
      include '../includes/makex.h'

c     --station name and coordinates from L-file
      include '../includes/model.h'  

c     --not needed?
c      include '../includes/errflg.h' 
      
c     --unit numbers not contained in makex.h
      integer*4 usestbl

c     --controls for status of data and file
      logical fend,ferr,found_data,weekset,found_header,first_nav,debug
      integer*4 iflag
          
c     --channels allowed by receiver and filled at each epoch
      integer*4 mchan   

c     --epoch counter for observation and navigation records
c      integer*4 nepoch - in model.h
                
c     --FICA/RINEX time tags
      integer*4 kyr,kmo,kdy,kdoy,khr,kmin
      real*8 sec  

c     --date of run
      integer*4 jyr,jmo,jdy

c     --functions in gamit/lib
      integer*4 julday,iarray,iclarg,nblen

c     --initial integers
      real*8 phs10(maxsat),phs20(maxsat)
                      
c     --other reals
      real*8 utcoff,vlight,svdt,dum,decyrs,site_epoch

c     --other integers
      integer*4 iarg,ioerr,isnrf,idum,len,icol,i,j

c     --other characters
      character*1 ay1           
      character*2 ay2
      character*3 rcvrsw 
      character*3 adoy,buff3 
      character*4 site
      character*15 buff15
      character*16 uname 
      character*80 frinexn,wildcard,pickfn
      character*109 version
      character*256 message 

c     other logicals
      logical fcheck,sestbl_reqd,found,first_obs(maxsat),batch,crdfile
 

c     Items from the FICA file
       
c     --current time (GPS week, seconds of week)
      integer*4 igpswk(maxchn)                                     
      real*8 gpssec(maxchn)
c     --start/stop times  (GPS week, seconds of weeks)
      integer*4 iwknstart,iwknstop
      real*8 sowstart,sowstop   
c     --quality factor L1, L2 and tracking mode 
      integer*4 iqvec(maxchn,2),tmode(maxchn)
C     --sat id num PRN
      integer*4     nprn,svid(maxchn)
c     --L1, L2 phase and pseudorange
      real*8 dofl1(maxchn),dofl2(maxchn),prgl1(maxchn),prgl2(maxchn)
c     --SNR
      real*8 denrat(maxchn,2)
c     --week number, from file name or block 9
      integer*4 iwknfile

  
c     --receiver software 
      real*8 swveru

c     --Operator-entered info
      character*8  auser
      character*16 arcvr,aficsite,antht
      character*80 comments(maxlin)
      integer*4 ncomments

c     RINEX defined items

      real*4 rxver
      character*20 rxpgm,rxusr,rxdat
c     --comment
      integer irxcom
      character*60 rxcom(maxlin)
c     --mark name
      character*60 rxmrk
c     --observer
      character*20 rxobs
c     --agency
      character*40 rxagy
c     --receiver serial number,type and SW version
c      character*20 rcvrsn,rcvers (in model.h)
c      character*20 rctype   - in model.h
c     --antenna serial number and type
c      character*20 antsn - in model.h
c      character*20 anttyp - in model.h
c     --antenna offsets
      real*8 anth,ante,antn
c     --cartesian coordinates
      real*8 coords(3)
c     --wavelength factors  
      integer nwave1,nwave2
c     --observation types
      integer nobtyp
      character*2 rxobtp(maxobt)
c     --data interval in seconds
      real*8 sample_int,rxint
c     --data start time
      integer irxyr0,irxmo0,irxdy0,irxhr0,irxmn0
      real*8 rxsec0
c     --data stop time
      integer irxyr1,irxmo1,irxdy1,irxhr1,irxmn1
      real*8 rxsec1   
c     --number of SVs total for header
      integer*4 numsv
c     --satellite numbers for header
      integer*4 totsv(maxsat)
c     --number of obs of each type for each satellite for header
      integer*4 numobs(maxobt,maxsat)
c     --satellite list each epoch
      integer*4 epochsv(maxchn)
c     --data flag
      integer idflag      
c     --observations
      real*8 obs_rx(4)
c     --signal strength indicator
      integer*4 issi(maxchn,4)
C     loss-of-lock indicator
      integer*4 lli
c     --Standard ephemeris values
      real*8 bcephem(16),bcclock(6),subfr1(8)  

c     --speed of light for converting MiniMac ranges
      data vlight/299792.458d0/

c----------------------------------------------------------------------------
   
c  Initialization
c  --------------

c     Remove old versions of the status, warning, and error files
      call report_stat('CLEAR','FIC2RX',' ',' ', ' ',0)
      call report_stat('CLEAR','LIB',' ',' ', ' ',0)   
c     these stored in ..includes/makex.h
      uinfor = 8
      uscren = 6
      usceno = 10
      urinex = 11
      ucoord = 12
      uficaf = 14 
c     this local
      usestbl = 15  
c     J-file for MiniMac only
      usvclk = 21
      debug = .false.
      numsv = 0
      do i=1,maxsat
        totsv(i) = 0 
        phs10(i) = 0.d0
        phs20(i) = 0.d0
      enddo
      antht = ' '      
      do i=1,maxlin
        comments(i) = ' '
        rxcom(i) = ' '
      enddo   
      idflag = 0

c  Identify FIC2RX, user, and agency 
c  ----------------------------------
          
c     open the print file  
      open (unit=uinfor,file='fic2rx.out',form='formatted'
     .     , status='unknown',iostat=ioerr)
      if(ioerr .ne. 0 ) then
        call report_stat('FATAL','FIC2RX','fic2rx','fic2rx.out',
     .  'Error opening FIC2RX print file: ',ioerr)
      endif
      write(uinfor,'(/,a)') 'Program FIC2RX'
c     get the version number and write it to the screen or log file
c     call mversn(uinfor,version)   -- not relevant since code moved to /fica   
c     form is char*40: '9.37 of 97/01/13 07:50:00 (SunOS)'
c      write(uinfor,'(a1)') ' '
c     begin the GAMIT.status file  
c      write(message,'(3a)')  'Started FIC2RX (MAKEX v. ',version,')'
c     call report_stat('STATUS','FIC2RX','fic2rx',' ',message,0) 
      call report_stat('STATUS','FIC2RX','fic2rx',' '
     .  ,'Summary output written to file fic2rx.out',0)


c  Get the name of the FICA file and open it
c  -----------------------------------------
             
c     file name (fficaf) and unit number (uficaf) are stored in makex.h
      batch = .true.
      iarg = iclarg(1,fficaf)
      if( iarg.le.0 ) then
c     if no arguments on the command line, prompt from an 'ls' list  
        batch = .false.
        wildcard = '*.fic' 
        len = 5
        fficaf = pickfn(wildcard,len )
        fficaf = fficaf(1:len)
      else
        call rcpar(1,fficaf)
      endif       
      open( uficaf, file=fficaf, iostat=ioerr, status='old' )
      if( ioerr.eq. 0) then    
         write(message,'(2a)') 'Opened FICA file ',fficaf(1:13)
         call report_stat('STATUS','FIC2RX','fic2rx',' ',message,0)
         write(uinfor,'(a)') message
      else
         call report_stat('FATAL','FIC2RX','fic2rx',fficaf
     .                   ,'Error opening input FICA file',ioerr)
      endif
      site = fficaf(1:4)  
                  

c  Read through the FICA file looking for header information
c  (Blocks 0 and 101 for GESAR) 
c  --------------------------------------------------------------
                               
      found_header = .false.   
      do while ( .not.found_header )
        call rficahd( debug,found_header,rcvrsw,swveru,sample_int
     .              , auser,arcvr,aficsite,antht
     .              , ncomments,comments
     .              , nwave1,nwave2 )                            
        rxint = sample_int
      enddo 
      if( .not.found_header) then  
        write(message,'(a)') 
     .    'No header found on FICA file: RINEX header may be bogus'
        call report_stat('WARNING','FIC2RX','fic2rx',' ',message,0)
        write(uinfor,'(a)') message 
      endif
      rewind (uficaf)                    
      

c  Set the # channels, observables, and wavelenght  according to receiver type
c  --------------------------------------------------------------------------
   
c     --  observation types
      nobtyp = 4
      rxobtp(1) = 'L1'
      rxobtp(2) = 'L2'
      rxobtp(3) = 'P1'
      rxobtp(4) = 'P2'  
  
      if( rcvrsw.eq.'GES'.or.rcvrsw.eq.'COR'.or.rcvrsw.eq.'ROM' ) then 
         write(message,'(2a)') 'Receiver apparently TI, set # channels'
     .        ,' = 4, obs = L1 L2 P1 P2, wavelengths = 1 1' 
         call report_stat('STATUS','FIC2RX','fic2rc',' ',message,0) 
         write(uinfor,'(a)') message    
         nchan = 4                                     
         nobtyp = 4
         rxobtp(1) = 'L1'
         rxobtp(2) = 'L2'
         rxobtp(3) = 'P1'
         rxobtp(4) = 'P2'       
         nwave1 = 1
         nwave2 = 1  

      else if ( rcvrsw.eq.'MIN' ) then    
         write(message,'(2a)') 'Receiver apparently MiniMac, set # '
     .     ,'channels = 8, obs = L1 L2 C1, wavelengths = 1 2'
         call report_stat('STATUS','FIC2RX','fic2rc',' ',message,0) 
         write(uinfor,'(a)') message
         nchan = 8 
         nobtyp = 3    
         rxobtp(1) = 'L1'
         rxobtp(2) = 'L2'
         rxobtp(3) = 'C1'
         nwave1 = 1
         nwave2 = 2  

      else
         nchan = maxchn
         write(message,'(a)')
     .      'Firmware not for TI or MiniMac, add code'
         write(uinfor,'(a)') message
         call report_stat('FATAL','FIC2RX','fic2rx',' ',message,0)
      endif
 

c  Copy the FICA header info into the RINEX comments field
c  -------------------------------------------------------

c     put a standard message into the first two lines
      irxcom = 3
      rxcom(1) = '**FICA files translated into RINEX**' 
      rxcom(2) = '  Contents of FICA header blocks:'
      rxcom(3) = '  ' 
      do i = 1,ncomments
        irxcom = irxcom + 1
        if( irxcom.le.maxlin ) then
           rxcom(irxcom) = comments(i)(1:60)
        else
          write(message,'(a,i3,a)') 
     .      'Number of comment lines exceeds maxlin (',maxlin
     .      ,'); omit remainder'
          call report_stat('WARNING','FIC2RX','fic2rx',' ',message,0) 
          write(uinfor,'(a)') message
          goto 1
        endif   
      enddo
    1 continue
 

c  Read the FICA file again to get the week number, start/stop times and SV list 
c  -----------------------------------------------------------------------------
        
c     We need the week number to interpret the start time since the
c     TI FICA BLK 6 returns only the seconds-of-week.  Set it 
c     intially from the file name; for GESAR files, it will be updated
c     whenever a valid ephemeris record (BLK 9) is read by DOFICA, for
c     TI ROM files, the week number is in the data record (BLK 401).
c     We allow three conventions for FICA file names:
c        ssssyddd.fic
c        ssssy.ddd.fic 
c        ssssyy.ddd
c     See if one of these is used by looking for the period
      icol = index(fficaf,'.')
      if( icol.eq.9 ) then
c        assume the form ssssyddd.fic
         read(fficaf(5:5),'(a1)',iostat=ioerr) ay1 
         if( ioerr.ne.0 ) call report_stat('FATAL','FIC2RX','fic2rx'
     .    ,fficaf,'Error reading year from FICA file name',ioerr)
         read(fficaf(6:8),'(a3)',iostat=ioerr) adoy    
         if( ioerr.ne.0 ) call report_stat('FATAL','FIC2RX','fic2rx'
     .    ,fficaf,'Error reading doy from FICA file name',ioerr)
       elseif (icol.eq.6) then
c        assume the form ssssy.ddd.fic
         read(fficaf(5:5),'(a1)',iostat=ioerr) ay1
         if( ioerr.ne.0 ) call report_stat('FATAL','FIC2RX','fic2rx'
     .    ,fficaf,'Error reading year from FICA file name',ioerr)
         read(fficaf(7:9),'(a3)',iostat=ioerr) adoy    
         if( ioerr.ne.0 ) call report_stat('FATAL','FIC2RX','fic2rx'
     .    ,fficaf,'Error reading doy from FICA file name',ioerr) 
       elseif (icol.eq.7) then
c        assume the form ssssyy.ddd
         read(fficaf(6:6),'(a1)',iostat=ioerr) ay1
         if( ioerr.ne.0 ) call report_stat('FATAL','FICA2RX','fic2rx'
     .     ,fficaf,'Error reading year from FICA file name',ioerr)
         read(fficaf(8:10),'(a3)',iostat=ioerr) adoy
         if( ioerr.ne.0 ) call report_stat('FATAL','FIC2RX','fic2rx'
     .     ,fficaf,'Error reading doy fro FICA file name',ioerr) 
       else
         call report_stat('FATAL','FIC2RX','fic2rx'
     .    ,fficaf,'Unsupported FICA file name',ioerr)     
       endif
c         assume that FICA files fall between 1985 and 1994
      read(ay1,'(i1)',iostat=ioerr) kyr
      if( ioerr.ne.0 ) call report_stat('FATAL','FIC2RX','fic2rx'
     .    ,ay1,'Error reading integer year from character year',ioerr)
      if( kyr.lt.4 ) then
        kyr = 1990 + kyr
      else 
        kyr = 1980 + kyr
      endif
      read(adoy,'(i3)',iostat=ioerr) kdoy  
      if( ioerr.ne.0 ) call report_stat('FATAL','FIC2RX','fic2rx'
     .    ,adoy,'Error reading integer doy from character doy',ioerr)
      khr = 0
      kmin = 0
      sec = 0.d0
      call timcon(-4,iwknfile,sowstart,kyr,kdoy,khr,kmin,sec,utcoff)
      weekset = .true.
      iwknstart = iwknfile  
      write(message,'(a,i4,a)') 
     .     'Initial week number set to ',iwknfile,' from FICA file name'
      write(uinfor,'(a)') message
      call report_stat('STATUS','FIC2RX','fic2rx',' ',message,0)
      found_data  = .false.     
      nepoch = 0      
      fend = .false.
      do while ( .not.fend  )  
        call dofica( debug,iflag,fend,ferr,nprn,svid,tmode
     .             , igpswk,gpssec,dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .             , bcephem,bcclock,subfr1,iwknfile,weekset
     .             , iwknstart,sowstart ) 
c       iflag = 1 data record  ; = 2 ephemeris record    
c*        print *,'DEBUG iflag found_data ',iflag,found_data
        if( iflag.eq. 1 ) then
c           see if first data record
            if( .not.found_data ) then
c               set the start seconds-of-week  (no week # on data record)
                sowstart = gpssec(1)
                found_data = .true. 
                write(uinfor,'(a,i4,f10.2,a)') 'Start time (wk sow = '
     .                    ,iwknstart,sowstart,') set from data record'
            endif
            do i = 1,nchan
               if( svid(i).ne.0 ) then
                 found = .false. 
                 do  j = 1, numsv 
c                    print *,'i j svid ',i,j,svid(i)
                    if( svid(i).eq.totsv(j) ) found = .true.
                 enddo 
c                if this is new SV; if so add it to the list
                 if( .not.found ) then
                   numsv = numsv + 1
                   totsv(numsv) = svid(i) 
                 endif  
                 numobs(1,iarray(svid(i),totsv,numsv)) = 
     .                   numobs(1,iarray(svid(i),totsv,numsv)) + 1
                 idum = iarray(svid(i),totsv,numsv)
c                 print *,'i svid numsv idum numobs',i,svid(i),numsv,idum
c     .                , numobs(1,iarray(svid(i),totsv,numsv))
               endif
            enddo 
        endif    
        if ( iflag.eq.2 .and. .not.found_data .and.
     .       iwknfile.ne.iwknstart ) then 
          iwknstart = iwknfile
          write(uinfor,'(a,i4,a)' ) 'Week number updated from BLK 9 (='
     .         ,iwknfile,')'   
        endif   
        nepoch = nepoch + 1 
c        if( nepoch.gt.20 ) stop
      enddo               
c      print *,'numsv totsv ',numsv,(totsv(i),i=1,numsv)
c      print *,'numobs ',(numobs(1,i),i=1,numsv)
c     finished, so set the stop time
      iwknstop = iwknfile
      sowstop = gpssec(1)
c     now reset 'iwknfile' to the start value for the next read
      iwknfile = iwknstart
      rewind(uficaf)
                     

c  If an L-file is input, get the coordinates from it
c  --------------------------------------------------
      
      if( batch ) then
        iarg = iclarg(2,fcoord)
        if( iarg.le.0 ) then  
         write(message,'(a)') 
     .      'No coordinate file input--coordinates set to zero'
         call report_stat('WARNING','FIC2RX','fic2rx',' ',message,0)
         write(uinfor,'(a)') message
          crdfile = .false.
        else 
          if( .not.fcheck(fcoord) ) then  
             write(message,'(3a)') 'L-file ',fcoord(1:nblen(fcoord))
     .            ,' not found--coordinates set to zero' 
             write(uinfor,'(a)') message
             call report_stat('WARNING','FIC2RX','fic2rx',' ',message,0) 
             write(uinfor,'(a)') message
             crdfile = .false.   
          else
             crdfile = .true.
          endif
        endif
      else
        write(uinfor,'(/,a)') 'Get coordinates from an L-file'
c       if 0 or 1 arguments on the command line, prompt from an 'ls' list
        wildcard = 'l*' 
        len = 2
        fcoord =  pickfn( wildcard,len )  
        fcoord = fcoord(1:len)
        if( fcoord(1:6).eq.'none ') then
           write(message,'(a)') 
     .        'No L-file requested--coordinates set to zero'
           call report_stat('WARNING','FIC2RX','fic2rx',' ',message,0) 
           write(uinfor,'(a)') message
           crdfile = .false.
        elseif( .not.fcheck(fcoord) ) then  
           write(message,'(3a)') 'L-file ',fcoord(1:nblen(fcoord))
     .         ,' not found--coordinates set to zero' 
           write(uinfor,'(a)') message
           call report_stat('WARNING','FIC2RX','fic2rx',' ',message,0) 
           crdfile = .false.
        else
           crdfile = .true.  
        endif
      endif
      if( crdfile ) then
        iul = ucoord
        open ( iul, file = fcoord, iostat=ioerr, status = 'old' )
c       determine if l-file or apr file
        call crd_file_type(fcoord,kfflg)
        if( ioerr.eq.0 ) then
           write(message,'(2a)') 'Opened coordinate file ',fcoord
           write(uinfor,'(a)' ) message     
           call report_stat('STATUS','FIC2RX','fic2rx',' ',message,0)   
           site_epoch = decyrs(kyr,kdoy,0.d0)
           call lread ( site, site_epoch )   
           do i=1,3
            coords(i)= kpos(i) + kvel(i)*(site_epoch-kepoch0)
           enddo
        else    
           call report_stat('FATAL','FIC2RX','fic2rx',fcoord
     .                      ,'Error opening coordinate file: ',ioerr) 
        endif
      endif

         
c  If the receiver is a Mini-Mac, we need to open the J-file to convert PR
c  -----------------------------------------------------------------------
       
      if( rcvrsw.eq.'MIN' ) then  

        if( batch ) then 
          iarg = iclarg(3,fsvclk)
          if( iarg.le.0 ) then       
            call report_stat('FATAL','FIC2RX','fic2rx',' '
     .        ,'Must have J-file as 3rd argument for Mini-Mac',0)
          else 
            if( .not.fcheck(fsvclk) ) then  
               write(message,'(3a)') 'J-file ',fsvclk(1:nblen(fsvclk))
     .              ,' not found' 
               write(uinfor,'(a)') message
               call report_stat('FATAL','FIC2RX','fic2rx',' ',message,0)
            endif
          endif
        else
         write(uinfor,'(/,a)') 'Get SV clock for Mini-Mac from a J-file'
c        if 0 or 1 arguments on the command line, prompt from an 'ls' list
         wildcard = 'j*'
         len = 2
         fsvclk = pickfn( wildcard,len )  
         fsvclk = fsvclk(1:len)
         if( fsvclk(1:6).eq.'none ') then
            call report_stat('FATAL','FIC2RX','fic2rx',' '
     .          ,'Must have J-file for Mini-Mac',0) 
         elseif( .not.fcheck(fsvclk) ) then  
            write(message,'(3a)') 'J-file ',fsvclk(1:nblen(fsvclk))
     .          ,' not found' 
            write(uinfor,'(a)') message
            call report_stat('FATAL','FIC2RX','fic2rx',' ',message,0) 
         endif
        endif
        open ( usvclk, file = fsvclk, iostat=ioerr, status = 'old')
        if( ioerr.eq.0 ) then
          write(message,'(2a)') 'Opened J-file ',fsvclk
          write(uinfor,'(a)' ) message     
          call report_stat('STATUS','FIC2RX','fic2rx',' ',message,0) 
          call readj( usvclk,totsv,numsv,idum,idum,dum,dum,0
     .              , idum,dum,dum,dum,dum,dum )
        else   
          call report_stat('FATAL','FIC2RX','fic2rx',fsvclk
     .                    ,'Error opening J-file: ',ioerr) 
        endif
      endif
               

c  Open the RINEX observation file and write the header
c  ---------------------------------------------------
           
c     get the RINEX file name
      write(ay2,'(i2.2)') mod(kyr,100)
      frinex(1:12) = fficaf(1:4) // adoy // '0.' // ay2 // 'o'
c     open the RINEX file (output)
      open(unit=urinex,
     .        file=frinex,
c*        Allow overwriting, at least while debugging
     .        status='unknown',
c*        Don't allow overwriting of RINEX files--too dangerous
c*     .        status='new',
     .        form='formatted',
     .        iostat=ioerr)
      if (ioerr.ne.0) then    
        write (uinfor,'(a,a20)') 'Error opening file : ', frinex
        call report_stat('FATAL','FIC2RX','fic2rx',frinex
     .                      ,'Error opening file',ioerr)
      else   
         call report_stat('STATUS','FIC2RX','fic2rx',frinex
     .                   ,'Opened file',0)
         write (uinfor,'(a,a20)')  'Opened file: ', frinex
      endif
c     -- RINEX version
      rxver = 2.10
c       identify FIC2RX version in 20 characters
      write(rxpgm,'(a7,a13)') 'FIC2RX ',version(1:13)
c     -- user name, agency, and run date
      call getusr(uname) 
      rxusr = '    '//uname  
      rxagy = '                                       '
      if( fcheck('sestbl.')) then
         open (unit=usestbl,file='sestbl.',status='old',iostat=ioerr)
         if( ioerr.ne.0) call report_stat('FATAL','FIC2RX','fic2rx'
     .                 ,'sestbl.','Error opening file ',ioerr)
         sestbl_reqd = .true.
         call rdsest( 17, 'Processing agency',3,buff3,usestbl
     .              , sestbl_reqd, ioerr )
         rxagy(1:3) = buff3 
         if( ioerr.ne.0 ) then 
             write(message,'(a)') 
     .          'Agency missing from sestbl., set blank in RINEX header'
             call report_stat('WARNING','FIC2RX','fic2rx',' ',message,0) 
             write(uinfor,'(a)') message 
         endif
      else       
         write(message,'(a)') 
     .        'No sestbl. available, set agency blank in RINEX header'
         call report_stat('WARNING','FIC2RX','fic2rx',' ',message,0) 
         write(uinfor,'(a)') message
          rxagy = ' '
      endif
c     -- date of run
      call getdat(jyr,jmo,jdy) 
      write (rxdat,'(i4.4,"/",i2.2,"/",i2.2)') jyr,jmo,jdy 
c     -- monument (station) name
      rxmrk = aficsite // '                                          '
c     -- observer
      rxobs = auser // '            '  
c     -- agency set from read of sestbl. 
c     -- receiver serial number
      rcvrsn(1:16) = arcvr
      rcvrsn(17:20) = '    '
c     -- receiver type  
      if( rcvrsw.eq.'GES'.or.rcvrsw.eq.'COR'.or.rcvrsw.eq.'ROM' ) then
        rctype(1:6) = 'TI4100' 
      elseif (rcvrsw.eq.'MIN' ) then
        rctype(1:6) = 'MIN6AT'
      else   
        write(message,'(a,a3)') 'Receiver type not coded for firmware '
     .         ,rcvrsw  
        write(uinfor,'(a)') message
        call report_stat('WARNING','FIC2RX','fic2rx',' ',message,0)
        rctype(1:6) = '      '
      endif
      rctype(7:20) = '              '
c     -- receiver firmware version
      rcvers(1:3) = rcvrsw  
      write(rcvers(4:9),'(1x,f5.2)') swveru  
      rcvers(10:20) = '           '
c    --  antenna type and serial number
      anttyp = '                    '
      antsn = '                    '   
c    -- antenna height                  
      read(antht,'(f16.0)',iostat=ioerr) anth
      if( ioerr.ne.0 ) then  
        write(message,'(2a)') 'Unreadable anth from FICA: ',anth
        write(uinfor,'(a)') message
        call report_stat('WARNING','FIC2RX','fic2rx',' ',message,0)   
      endif
      ante = 0.d0
      antn = 0.d0               
c     --time of first and last observations
      call timcon( 4,iwknstart,sowstart
     .           , irxyr0,kdoy,irxhr0,irxmn0,rxsec0,utcoff )  
      call monday(kdoy,irxmo0,irxdy0,irxyr0) 
c      print *,'sowstart irxhr0 irxmn0 rxsec0 ',
c     .        sowstart, irxhr0, irxmn0, rxsec0
      call timcon( 4,iwknstop,sowstop
     .           , irxyr1,kdoy,irxhr1,irxmn1,rxsec1,utcoff )  
      call monday(kdoy,irxmo1,irxdy1,irxyr1)  
c     --set the counts for all the observation types the same
      do i=1,numsv
        do j=2,nobtyp
          numobs(j,i) = numobs(1,i)
        enddo
      enddo
c     --write the header
      call wrxhed (urinex,
     .      rxver,rxpgm,rxusr,rxdat,rxcom,irxcom,rxmrk,rxobs,rxagy,
     .      rcvrsn,rctype,rcvers,antsn,anttyp,
     .      coords(1),coords(2),coords(3),
     .      anth,ante,antn,nwave1,nwave2,nobtyp,rxobtp,rxint,
     .      irxyr0,irxmo0,irxdy0,irxhr0,irxmn0,rxsec0,
     .      irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1,
     .      numsv,totsv,numobs )
                                
c     Loop over data records
c     ----------------------
                
      fend = .false.
      nepoch = 0 
      do i=1,maxsat
        first_obs(i) = .true.
      enddo
 
      do while (.not.fend )

        call dofica( debug,iflag,fend,ferr,nprn,svid,tmode
     .             , igpswk,gpssec,dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .             , bcephem,bcclock,subfr1,iwknfile,weekset
     .             , iwknstart,sowstart ) 
c       iflag = 1 data record; (=2 ephemeris record)
        if( iflag.eq.1 ) then
c         convert GPS week,sec to calender date/time
          call timcon( 4,iwknfile,gpssec(1)
     .               , irxyr1,kdoy,irxhr1,irxmn1,rxsec1,utcoff )  
          call monday(kdoy,irxmo1,irxdy1,kyr)  
          irxyr1 = mod(kyr,100)
c         count the number of non-zero channels at this epoch
          mchan = 0
          do i = 1,nchan
            if( svid(i).ne.0 ) then
              mchan = mchan + 1 
              epochsv(i) = svid(i)
            endif
          enddo
c         write the time line
          write(urinex,'(5i3,f11.7,i3,12(i3))')
     .              irxyr1,irxmo1,irxdy1,irxhr1,irxmn1,rxsec1
     .            , idflag,mchan,(epochsv(i),i=1,mchan)     
          do i = 1, nchan
           if( svid(i).ne.0 ) then
c            change phase from Doppler (FICA) to pseudorange (RINEX) convention 
             obs_rx(1) = - dofl1(i) 
             obs_rx(2) = - dofl2(i) 
             obs_rx(3) = prgl1(i)
             obs_rx(4) = prgl2(i) 
c            update the session SV counter and initial phase values
             if ( first_obs(iarray(svid(i),totsv,numsv)) ) then 
c               this is the first occurrence of the SV, set the initial phase
                phs10(iarray(svid(i),totsv,numsv)) = nint(obs_rx(1))
                phs20(iarray(svid(i),totsv,numsv)) = nint(obs_rx(2))
                first_obs(iarray(svid(i),totsv,numsv)) = .false.
             endif  
c            RINEX standard is to remove the integer part of the intial phase;
c            this is controversial, but it's easier to compare final x-files if we do
             obs_rx(1) = obs_rx(1) - phs10(iarray(svid(i),totsv,numsv))
             obs_rx(2) = obs_rx(2) - phs20(iarray(svid(i),totsv,numsv)) 
c            map L1 and L2 signal amplitudes to RINEX ISSI
             call mapamp
     .           (1,rcvrsw,swveru,1,denrat(i,1),isnrf,issi(i,1))
             call mapamp
     .           (1,rcvrsw,swveru,2,denrat(i,2),isnrf,issi(i,2))  
             lli = 0
             issi(i,3) = 0
             issi(i,4) = 0  
c            if a Mini-Mac, convert the pseudorange into RINEX standard
             if( rcvrsw.eq.'MIN' ) then
                call timcon( 1,igpswk(1),gpssec(1)
     .                     , kyr,kdoy,khr,kmin,sec,utcoff )
                call monday(kdoy,kmo,kdy,kyr)
                jdobs= julday(kmo,kdy,kyr)
                tobs = khr*3600.d0 + kmin*60.d0 + sec
                call readj( usvclk,totsv,numsv,svid(i),jdobs
     .                    , tobs,svdt,1,idum,dum,dum,dum,dum,dum )
                prgl1(i) = prgl1(i) - svdt*vlight*1.d3
             endif
c            check for reasonableness to avoid writing asterisks
             if( dabs(obs_rx(1)).gt.9.d9 ) obs_rx(1) = 0.d0
             if( dabs(obs_rx(2)).gt.9.d9 ) obs_rx(2) = 0.d0
             if( dabs(obs_rx(3)).gt.9.d9 ) obs_rx(3) = 0.d0
             if( dabs(obs_rx(4)).gt.9.d9 ) obs_rx(4) = 0.d0 
c            write the data line
             write(urinex,'(4(f14.3,i1.0,i1))') (obs_rx(j),lli,issi(i,j)
     .                                         , j=1,4 )  
            endif
          enddo   
c**debug          if( nepoch.gt.10 ) fend = .true.
          nepoch = nepoch + 1
        endif
    
      enddo                  
      write(message,'(i5,a,a12 )')  nepoch,' epochs written to ',frinex
      write(uinfor,'(a)') message
      call report_stat('STATUS','FIC2RX','fic2rx',' ',message,0) 
      rewind( uficaf )
      close ( urinex )
          

c     Write the navigation file
c     -------------------------
                  
      call report_stat('STATUS','FIC2RX','fic2rx',' '
     .                ,'Looking for navigation records',0) 
      fend = .false.
      ferr = .false.
      first_nav = .true.
      nepoch = 0   
      do while (.not.fend )
                
        call dofica( debug,iflag,fend,ferr,nprn,svid,tmode
     .             , igpswk,gpssec,dofl1,dofl2,prgl1,prgl2,denrat,iqvec
     .             , bcephem,bcclock,subfr1,iwknfile,weekset
     .             , iwknstart,sowstart ) 
c       iflag = 2 ephemeris record; (=1 phase data record)
        if( iflag.eq.2 ) then
          if( first_nav ) then
            frinexn = frinex
            frinexn(12:12) = 'n'  
c           open the RINEX file (output)
            open(unit=urinex,
     .        file=frinexn,
c*            Yes, allow overwriting, at least for debug
     .        status='unknown',
c*            Don't allow overwriting of RINEX files--too dangerous
c*     .        status='new',
     .        form='formatted',
     .        iostat=ioerr)
            if (ioerr.ne.0) then    
              write (uinfor,'(a,a20)') 'Error opening file : ', frinexn
              call report_stat('FATAL','FIC2RX','fic2rx',frinexn
     .                        ,'Error opening file',ioerr)
            else   
              call report_stat('STATUS','FIC2RX','fic2rx',frinexn
     .                        ,'Opened file',0)
              write (uinfor,'(a,a20)')  'Opened file: ', frinexn
            endif
            buff15= 'NAVIGATION DATA'
            write (urinex,'(f9.2,11x,a15,25x,a)',iostat=ioerr)
     .                  rxver,buff15,'RINEX VERSION / TYPE'
            write (urinex,'(3a20,A)',iostat=ioerr)
     .                  rxpgm,rxusr,rxdat,'PGM / RUN BY / DATE' 
            write (urinex,'(a,27x,a7,13x)') 
     .                'Broadcast elements from FICA file','COMMENT'
            write (urinex,'(60x,a20)') 'END OF HEADER       '
            first_nav = .false.
          endif             
          call wrxnav ( urinex,nprn,iwknfile,bcephem,bcclock,subfr1 )
          nepoch = nepoch + 1  
c**debug          if( nepoch.gt.10 ) fend = .true.
        endif
 
      enddo 
      
      if( nepoch.gt.0 ) then 
        write(message,'(i5,a,a12)') nepoch,' epochs written to ',frinexn  
      else
        write(message,'(a)') 'No navigation records found' 
      endif
      write(uinfor,'(a)') message
      call report_stat('STATUS','FIC2RX','fic2rx',' ',message,0) 
      rewind( uficaf )
      close ( urinex ) 
      call report_stat('STATUS','FIC2RX','fic2rx',' ','Normal end',0)
  
      stop
      end


