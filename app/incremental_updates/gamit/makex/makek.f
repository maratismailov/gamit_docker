Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994. All rights reserved.

      Program MAKEK
c
c     make k files from e and x files
c     R.W. King August 1989, from code in MAKEX
c
c     Input files:
c
c       x-file lu 20     X-file
c       e-file lu 21     Orbit file
c       l-file lu 22     L-file of station coordinates
c
c     Output file:
c
c       k-file lu 23      The CLOCK file, assumed created.
c
      implicit none

      include '../includes/dimpar.h' 
      include '../includes/units.h'   
      include '../includes/global.h'
      include '../includes/makex.h'
      include '../includes/model.h'

c
c                        DECLARATIONS
c                        ************
c
      logical         lans,satok,found,warnings,batch,badpr

      character*1      latflag,lonflag
      character*3      rcvrsw,rxobtyp(maxdat)
      character*4      upperc,site
      character*6      adelta
      character*16     uname,satnam(maxsat)
      character*20     rcvnum,antnum
      character*80     wildcard,pickfn,buff80
      character*109    version 
      character*40     version40
      character*16     sitnam16
      character*80     files(100),files2(200)
      character*256    message

      integer*4        iwknx,iwknx0,iwkne,ncount
     .               , irunt(6),ihnsec,ioerr,nsat
     .               , iprnx(maxsat),ierr(maxsat),nblen
     .               , iwknk,imin,id,ih,iy,nfiles,nfiles2,ifile
     .               , ihr,min,idoy,im,latd,latm,lond,lonm
     .               , iserr,icall,itflag,len
     .               , iarg,iclarg,doy,idum,i,j

      real*8           seclat,seclon,sowx,sowx0,coords(3)
     .               , utcoff,delta,secdif
     .               , prgl1(maxsat),bephem(16),xsat(3)
     .               , bclock(6),sowk,sec,sowe 
     .               , subfr1(8),sod
     .               , decyrs,site_epoch,rdum,rdum3(3)

      integer          maxprn
      parameter        (maxprn=31)
c     contrary to usual practice the variable sat_there is indexed by
c     the PRN numbers themselves, not a sequence index, since there is
c     is not a full, sequential listing of satellites in this routine
       
c     K-file headers and format statment (new 160831)
      character*105 kheader1
      character*126 kheader2 
      character*12  ksource
      character*51 afmt
      character*4  kversion

      logical fcheck,lask,sat_there(maxprn)
     .     ,debug/.false./,debuge/.false./
c
c     k-file format     
      data afmt/'(2i4,2i3,f11.7,i6,f15.7,2x,a1,i2.2,3f17.9,i6,f15.7)'/ 

c      kheader1 100  
c        'Ver ' kversn ' Station clock values from '
c          4     4           27 
c      ksource   12 
c      space      1
c      fsvclk    12
c      '  MAKEX ' 8
c      version40 40
                  


 1000 format (1x,a,1x,$)

c                       RUN INITIALIZATION
c                       ******************

C     unit numbers
      uinfor = 0
      uscren = 6
      uxfile = 20
      unav  = 21
c    kin change
      iul = 22
c    end kin change
      uclock = 23

c     get the run time,user name and version number

      call mversn(version)
      version40 = version(1:40)
      call getdat(irunt(1),irunt(2),irunt(3))
      call gettim(irunt(4),irunt(5),irunt(6),ihnsec )
      call getusr(uname)

c     write makek status line
      WRITE(uscren,'(a)')' '
      WRITE(message,5)version
    5 FORMAT('Started MAKEK ',a109)
      call report_stat('STATUS','MAKEK','makek',' ',message,0)

c      write (uscren,10) version,uname(1:nblen(uname)),(irunt(i),i=1,6)
c   10 format (/,1x,'MAKEK ',a,/
c     .   ,1x,'Run by ',a,1x,'on',1x,
c     .   i4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2,//)
 
c                   READ THE INPUT
c                   **************

c     if there are command-line arguments, use them and skip the interactive questions
      iarg = iclarg(1,fnav)   
      if( iarg.gt.0 ) then              
        batch = .true.
cold:    four arguments: nav-file, x-file, l-file, and output interval; the last is optional
crwk 161116: five arguments:
c         nav-file  x-file  l-file gnss out-interval  (last two are optional)
c       overwrite any existing k-files
        iarg = iclarg(2,fxfile)
        if( iarg.le.0 ) call report_stat('FATAL','MAKEK','makek',' ',
     .     'Missing command-line argument for X-file',0)
        iarg = iclarg(3,fcoord)  
        if( iarg.le.0 ) call report_stat('FATAL','MAKEK','makek',' ',
     .     'Missing command-line argument for L-file',0) 
        iarg = iclarg(4,gnss)
        if( iarg.le.0 ) then
          gnss = 'G' 
        endif
        adelta = ' '
        iarg = iclarg(5,adelta)
        if( iarg.gt.0 ) then
          read(adelta,*) delta
        else       
          delta = 600.d0
        endif
      else   
c       no command-line arguments, enter info interactively             
        write(6,*) 
        write(6,*) 
     .     'To run with command-line arguments, abort and enter: '
        write(6,*) 
     .    ' makek [nav-file] [x-file] [l-file] [output interval in sec]'
        write(6,*) ' '
        write(6,*) 
     .          '  All required except interval, which defaults to 600s'
        write(6,*) ' '
        batch = .false.
c       enter the input and output file names
        nfiles2 = 0
        buff80 = 'e*,*'
        call getdir (buff80,100,files,nfiles)
        if(nfiles.gt.0) then
          do i =1,nfiles
             files2(i) = files(i)
          enddo
          nfiles2 = nfiles
        else
          call report_stat('STATUS','MAKEK','makek',' ',
     .       'No e-file found, look for a RINEX nav file',0)
        endif
        buff80 = '*.??n'
        call getdir (buff80,100,files,nfiles)
        if( nfiles.gt.0 ) then
          do  i =1,nfiles
            files2(nfiles2+i) = files(i)
          enddo
          nfiles2 = nfiles2 + nfiles 
        endif
        if( nfiles2.le.0 ) call report_stat('FATAL','MAKEK','makek',' '
     .           ,'No navigation file found',0)        
        write (6,*) 'Available ephemeris files:'
        do  i = 1,nfiles2
           write (*,'(i3,1x,a)') i,files2(i)
        enddo
        write (6,1000) 'Choose a number:'
        read (5,*) ifile
        fnav = files2(ifile) 
  
c       get the X-file name
        write (6,*) 'Select a file containing pseudoranges.'
        wildcard = 'x*.*'
        len = 4
        fxfile  = pickfn ( wildcard,len )
        fxfile = fxfile(1:len)    

c       get the  L-file name
        fcoord='lfile.' 
        if( .not.fcheck(fcoord) ) then
          found = .false.
          do while (.not.found)
            write(uscren,'(a)') ' L-file not found'
            write(uscren,1000) ' Enter L-file name:'
            read(5,'(a)') fcoord 
            if( fcheck(fcoord) ) found =.true.
          enddo
        endif  
       
c       get the K-file interval  
        delta = 600.d0
        write(uscren,1000) 'Enter interval for K-file (sec) [600]:'
        read(5,*) delta

                 
c       wanna debug?
        write(uscren,1000) 'Debug run ?'
        lans = lask()
        if( lans ) then
          write (uscren,'(a)') 'Debug ephemeris ? '
             debuge = lask()
          write(uscren,'(a)') 'Debug x-file also?'
             debug = lask()
        endif           

      endif    

c     DISPLAY INPUT CONROLS AND SITE ID 
c     **********************************

      site=fxfile(2:5)
      site=upperc(site)          
      write(message,'(a,a4)') 'Making K-file for site: ',site
      call report_stat('STATUS','MAKEK','makek',' ',message,0)

                  
c     Create the header lines for the K-file
c     *************************************
                
c     j-file not used with makek 
      fsvclk= ' '      
      kheader1 = ' ' 
      kheader2 = ' '      
      ksource = ' ' 
* rwk 200504: MAKEK hard-wired to use navigation files 
      qsp3 = .false. 
      qnav = .true.
      if( qsp3 ) then
        ksource = fsp3(1:12) 
      elseif( qnav ) then
        ksource = fnav(1:12)
      endif                            
      kversion = ' 2.0'                                                
      kheader1= 'Ver '//kversion//' Station clock values from '//ksource
     .//' '//fsvclk(1:12)//'  MAKEX '//version40
      if( qsp3 ) then
      kheader2 = 'YEAR DOY HR MN  SEC(UTC)   WKNO   SOW(GPST)     PRN
     .OBSERVED PR(sec)    SV CLOCK    SITE CLOCK ' 
      else
      kheader2 = 'YEAR DOY HR MN  SEC(UTC)   WKNO   SOW(GPST)     PRN
     .OBSERVED PR(sec)    SV CLOCK    SITE CLOCK   NAV-FILE WKNO SOW(GPS
     .T)' 
      endif       

c       READ X-FILE HEADER and CALCULATE SITE COORDINATES
c      ************************************************
                 
c     Open the X-file
      open(unit=uxfile,file=fxfile,status='old',iostat=ioerr)
        if (ioerr .ne. 0 ) 
     .    call report_stat('FATAL','MAKEK','makek',fxfile,
     .    'Error opening file: ',ioerr)

c     Read and display the X-file header

      call xhdred ( uxfile,uinfor,uscren,nepoch,inter
     .            , ircint,isessn
     .            , mtime,iy,im,id,ihr,min,sec
     .            , nchan,ischan,satnam
     .            , ndat,dattyp,rxobtyp,lambda
     .            , offarp,sitnam16,rcvrsw,swver,antcod
     .            , rctype,rcvnum,anttyp,antnum
     .            , latflag,latd,latm,seclat
     .            , lonflag,lond,lonm,seclon,height
     .            , ntext,text,gnss )

      if( mtime.eq.1 ) then
         itflag = -2
      elseif( mtime.eq.2 ) then
         itflag = -4
      else             
        write(message,'(a,i3)') 'Bad time flag: ',mtime
        call report_stat('FATAL','MAKEK','makek',' ',message,0)
      endif               
      doy = idoy(iy,im,id)
      sod = ihr*3600.d0 + min*60.d0 + sec 
      call timcon ( itflag,iwknx0,sowx0
     .              , iy,doy,ihr,min,sec,utcoff )

c     Go and get the best estimates of the site from the
c     coordinates file. these coordinates, the best apriori possible
c     are used to compute clock offset parameters.

c     L-file coordinates stored in 'modkin.h' common/kinpar1/
      sitecd = site  
      open(unit=iul,file=fcoord,status='old',iostat=ioerr) 
       if( ioerr.ne.0 ) call report_stat('FATAL','MAKEK','makek'
     .     ,fcoord,'Error opening file: ',ioerr)   
c     determine whether l-file or apr file
      call crd_file_type(fcoord,kfflg)   
      site_epoch = decyrs( iy,doy,sod ) 
      call lread ( site,site_epoch )  
      do i=1,3
        coords(i)= kpos(i) + kvel(i)*(site_epoch-kepoch0)
      enddo
      if(debug) print *,'MAKEK site_epoch kepoch0 kpos kvel '
     .       ,site_epoch,kepoch0,coords,kpos,kvel


c       READ THE NAVIGATION FILE
c      ***********************

c     Open and read the navigation file
      open(unit=unav,file=fnav,status='old',iostat=ioerr)
      if (ioerr .ne. 0 ) then
        call report_stat('FATAL','MAKEK','makek',fnav,
     .  'Error opening file: ',ioerr)
      endif 
c     read in all the E-file records
      icall = 0    
c     prn# not used for icall=0                  
      gnss = satnam(1)(1:1) 
      if(debug) write(*,'(a,i2)') 'MAKEK calling getnav icall ',icall
      call getnav ( debuge,icall,iwknk,sowk,gnss,idum,idum,ischan
     .            , rdum3,rdum,satok,idum,rdum )
      if(debug) then
        print *,'MAKEK after getnav icall iwknk sowk gnss idum '
     .                            , icall,iwknk,sowk,gnss,idum
        print *,'  idum ischan ',idum,ischan
        print *,'  rdum3 rdum satok idum rdum '
     .         ,   rdum3,rdum,satok,idum,rdum
      endif 
        
c       GET THE K-FILE NAME
c       *******************

      call lowers (fxfile)
      fclock = fxfile
      i = index (fclock,'x')
      fclock(i:i) = 'k'
      if ( fcheck(fclock) .and. .not.batch ) then
         write(uscren,1000)
     .   fclock(1:nblen(fclock))//' exists. Overwrite?'
         if (.not. lask()) then
            write (6,1000) (' Enter output K-file name:')
            read(5,'(a)') fclock
         endif
      endif
 
                           
c         OPEN THE K-FILE and WRITE THE HEADERS
c         *************************************
   
      open(unit=uclock,file=fclock,status='unknown',iostat=ioerr)
      if (ioerr .ne. 0 ) then
        call report_stat('FATAL','MAKEK','makek',fclock,
     .  'Error opening file: ',ioerr)
      endif                    
      write(uclock,'(a,/,a,/,a)') kheader1,kheader2,afmt 


c         READ X- and NAV-FILE DATA
c        ***********************

c     Set the initial K-file epoch at the first X-file epoch
      iwknk= iwknx0
      sowk = sowx0
      ncount = 0
      found = .false.
      do i=1,maxprn
        sat_there(i)= .true.
      enddo
c     come here to read one epoch of the X-file
  100 continue
      ncount = ncount+1

      call getpr (debug,iwknx0,sowx0,inter,nchan,ischan,ndat
     .  ,dattyp,iwknx,sowx,nsat,iprnx,ierr,prgl1,ioerr)

c     if ioerr = -1 end of file
      if (ioerr .eq. -1) then
         call report_stat('STATUS','MAKEK','makek',' ',
     .   'End of X-file reached ',ioerr)
         goto 999
      else if (ioerr. gt. 0) then
         call report_stat('FATAL','MAKEK','makek',' ',
     .   'Error reading X-file ',ioerr)
      endif

      if ( debug) then
          call wtime(6,iwknx,sowx,'GPST','Epoch on X-file:')
      endif

c     See if the requested epoch has been found
c     If not, or if no data read another record
c      print *,'iwknx,sowx,iwknk,sowk,nsat: ',iwknx,sowx,iwknk,sowk,nsat
      if (secdif(iwknx,sowx,iwknk,sowk). lt. 0.d0 ) goto 100
      if (nsat.le.0) goto 100

c     as soon as we have found a good epoch, bring the
c     kfile epoch up to this point.
      if (.not.found) then
         iwknk = iwknx
         sowk  = sowx
         found = .true.
      endif

c     The requested epoch has been found.  Get the ephemeris and
c     clock information for the earliest epoch after the requested
c     one, for each satellite.  Then compute a clock correction.

      do 200 i=1,nsat 
c       check for reasonable pseudorange
        if( prgl1(i).eq.0.d0 ) then
           badpr = .true. 
        else
           badpr = .false.
        endif

        icall = 1                       
        gnss = satnam(i)(1:1)
        call getnav ( debuge,icall,iwknk,sowk,gnss,iprnx(i),nchan,ischan
     .              , xsat,svclock(1,i),satok,iwkne,sowe )

        if(debug) then
          print *,'MAKEK after getnav icall iwknk sowk gnss iprnx '
     .                              , icall,iwknk,sowk,gnss,iprnx(i)
          print *,'  nchan ischan ',nchan,ischan
          print *,'  i xsat svclock satok iwkne sowe '
     .         ,     i,xsat,svclock(1,i),satok,iwkne,sowe
        endif 
c       print warning only for first encounter of each satellite
        if( .not.satok ) then
         if( iprnx(i).le.maxprn. and. sat_there(iprnx(i)) ) then   
           write(message,'(a,i3)') 'No ephemeris info for PRN ',iprnx(i)
           call report_stat('WARNING','MAKEK','makek',' ',message,0)             
           sat_there(iprnx(i)) = .false.
         endif
        endif                            
        if(debug) then 
           write(*,'(a,2i3,4d16.8,l1)')  
     .    'MAKEK i iprnx xsat svclock satok '
     .         , i,iprnx(i),xsat,svclock(1,i),satok
         endif 

        if (satok .and. .not.badpr ) then
c         kin change
c         site coordinates are read from l-file and not from x-file, so
c         that makek can be run with updated coordinates without having to update
c         the x-file  


          if( mtime.eq.1 ) then
            itflag = 2
          elseif( mtime.eq.2 ) then
            itflag = 4
          endif   
          call timcon ( itflag,iwknx,sowx
     .                , iy,doy,ihr,min,sec,utcoff )
          sod = ihr*3600.d0 + min*60.d0 + sec 
          site_epoch = decyrs( iy,doy,sod ) 
          call lread ( site,site_epoch ) 
          do j=1,3
           coords(j)= kpos(j) + kvel(j)*(site_epoch-kepoch0)
          enddo      
          if(debug) print *,'MAKEK coords ',coords
          call stnclk ( debug,gnss,iwknx,sowx,iprnx(i),prgl1(i)
     .                , coords,xsat,svclock(1,i)
     .                , rclock,iserr,warnings ) 
          if(debug) then
             write(*,'(a,2d16.8,i3)') 
     .          'MAKEK after stnclk svclock rclock iserr '
     .          ,              svclock(1,i),rclock,iserr
          endif 
          if (iserr .eq. 0) then   
             if( qsp3 ) iwkne = 0   
             call wclock ( uclock,afmt,gnss,iprnx(i),iwknx,sowx,
     .       prgl1(i),svclock(1,i),rclock,iwkne,sowe  )
          endif
        endif
  200 continue

c     Update the K-file epoch
c     Use the next "even" epoch to have the best chance of matching the
c     (usually) hourly epochs of the broadcast ephemeris
      call secsum(iwknk,sowk,delta,iwknk,sowk)
      sowk= delta*dint( sowk/delta )
      call timcon(1,iwknk,sowk,iy,id,ih,imin,sec,utcoff )

c     now read the next record
      goto 100

  999 if( warnings ) then
        call report_stat('WARNING','MAKEK','makek',fclock,
     .  '**Warnings issued for bad clocks',0)
      endif

      call report_stat('STATUS','MAKEK','makek',fclock,
     .'Normal end to MAKEK',0)
      stop

      end
