Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995. All rights reserved.
     
      Program BCTOT

C    Written by Yehuda Bock, November 1987
C    Modified by R. King  May 1989; October 1995; January 2001   
c    Modified by R. King December 2015, August 2016  for GNSS


c    Read a Broadcast Ephemeris (E-) file in FICA Block 9 or RINEX nav-file format
c    and create a T-file in an earth-fixed system as well as an inertial system

      implicit none

      include '../includes/dimpar.h'     
      include '../includes/orbits.h'

      logical cmdline,found,debug,goodrec
                        
      character*1 gnss,yawbias

      integer*4 jds,jde,jdf,iye,idne,ihe,imine,iys,idns,ihs,mins
     .        , iwkne,iwkn,iewkn,itwkn
     .        , iesat,ixsat,itsat,nsat,nprn,iflag,nintrs
     .        , iexpnt,nbrd,ntbrd,isat,nesat,icall
     .        , ibcerr,idir,i,j,k,nics,icbrd
     .        , iepoch,inter,nepoch,closest_epoch
     .        , mon,day,yr,hr,min
     .        , i00(3),i01(3),ioerr,jyr,jdoy,iper
     .        , isvn,frqchn,svnstart(5),svnstop(5)
     .        , ifyr,isessn   
     .        ,  check_flg ,nssat,ssats(maxsat)
     .        , iarg,iclarg,nblen  

      real*8 ts,te,tf,ephem,bclock,bephem,tephem,secs,trans_sow
     .     , utcoff,satics,sece,satprm,delt,sow,sowe,subfr1
     .     , acheck,tsow,offset_hrs,x,t00(3),t01(3),stepsize
     .     , sbmass,yawrate
     .     , sowclk,tf3(3),sec,span
      real*8 antpwr   ! Antenna transmit power (W)

 

c      Units
      integer*4 iscrn,iterm,iprnt,iut1,inut,ipole,iutin,iutout
     .        , iux,iubc,iusess  
  
c  Units:
c    Assigned in OPENB:
c      iterm =   5
c      iscrn =   6
c      iux   =  14
c      iubc  =  15
c      iutin =  16
c      iutout=  17
c      inut  =  20
c      iut1  =  21
c      ipole =  22
c      iprnt =  30
c    Assigned in BCTOT:
c      iscrn  =  6
c      iusess = 40 (session.info)
c    Assigned in GET_MODELS:
c      iusest = 41  (sestbl.)

 
c     this eventually to be moved to includes/orbits.h
      integer maxhed                                                  
      parameter(maxhed=20)

      character*1  lowerc,upperc,ans,aflag
      character*3  cdoy
      character*4  icsnam(maxorb),cyr
      character*5  frame,precmod,nutmod,gravmod,srpmod,eframe
     .             ,eradmod,antradmod 
      character*6  prog
      character*16 bcfile,xfile,tfile,tfilef,gfile,gfilef,tfin,printfile
     .            ,tmpnam,satnam(maxsat)
      character*20 antbody 
      character*120 version
      character*256 message                 

      dimension bephem(16,maxbrd,maxsat), ephem(16), satprm(6)
     .        , bclock(6), satics(maxorb,maxsat),ixsat(maxsat)
     .        , iesat(maxsat),itsat(maxsat)
     .        , itwkn(maxbrd,maxsat), iewkn(maxbrd,maxsat)
     .        , tsow(maxbrd,maxsat)
     .        , iexpnt(maxsat),nbrd(maxsat),ntbrd(maxsat)
     .        , tephem(16,maxbrd,maxsat),subfr1(8),x(3,maxsat)
     .        , sowclk(maxbrd,maxsat) 

c      function
       logical fcheck
       integer*4 julday 

       data debug/.false./
        
      iscrn = 6
   
c        Get version number
      call oversn(version)

      WRITE(6,'(a)')' '
      WRITE(message,5)version
    5 FORMAT('Started BCTOT ',a120)
      call report_stat('STATUS','BCTOT','orbits/bctot',' ',message,0)


c----------- Read the input----------------------------------
c
c     bctot  [year] [day-of-year]  [navigation file]  [t-file]  [gnss]   [x-file]
c               (required)          (required)     (required  (optional) (optional)
c
c     With command-line arguments, start/stop times must be available from 
c     session.info or an x-file
c
c     No arguments implies interactive, times can be entered manually
c     
 
      cmdline = .false. 
      iarg = iclarg(1,cyr)  
      if( iarg.gt.0 ) then
        cmdline = .true. 
        frame = ' '
        read(cyr,*,iostat=ioerr) yr
        if( ioerr.ne.0 ) call report_stat
     .      ('FATAL','BCTOT','orbits/bctot',' '
     .      ,'Error reading command-line year',ioerr)
        iarg = iclarg(2,cdoy)
        read(cdoy,*,iostat=ioerr) jdoy  
        if( ioerr.ne.0 ) call report_stat
     .      ('FATAL','BCTOT','orbits/bctot',' '
     .      ,'Error reading command-line day',ioerr)
        iarg = iclarg(3,bcfile)   
        iarg = iclarg(4,tfile)            
        if( iarg.le.0 ) 
     .    call report_stat('FATAL','BCTOT','orbits/bctot',' '
     .                    ,'Too-few command-line arguments',0)
c       BC file name must be standard format, without full path for logic in BCTOT
        if( nblen(bcfile).gt.12 ) call report_stat('FATAL','BCTOT'
     .       , 'orbits/bctot',' '
     .       , 'nav file name too long --cannot use full path',0)    
c       construct efixed t-file name
        tfilef=tfile
        iper=index(tfile,'.')
        tfilef(iper-1:iper-1)='e'

        iarg = iclarg(5,gnss)
        if( iarg.le.0 ) then  
          gnss = 'G'
        endif

        iarg = iclarg(6,tmpnam)
        if( iarg.gt.0 ) then  
             xfile = tmpnam
        else
            xfile = " "
        endif                                              

      else 
        write(*,'(a,/,2a,/,2a )') 
     .     'Run-string: '
     .    ,'bctot  [year] [day-of-year]  [navigation file]  [t-file]'  
     .    ,'[gnss]  [x-file]  '
     .    ,'        (required)             (required)       (required)'
     .    ,'   (optional)  (optional)'

         stop
       endif

c------ Get the model and integration inputs from defaults or sestbl. and session.info,
c         (or read 5th and 6th lines of input file)
           
      call get_models( frame, precmod, nutmod, gravmod, srpmod, eradmod
     .               , antradmod, stepsize, delt, nics, icsnam ) 
 
c------ Assign unit numbers and open files used by multiple /orbits routines

      tfin = '                '
      printfile = 'bctot.out'          
      call openb( iterm,iscrn,iprnt,iutin,iutout,iubc,iux
     .          , inut,iut1,ipole,tfin,tfilef,bcfile,xfile,printfile )


c--------Get the start, stop, and IC epochs for the T-file 
c          (requires session.info, created by sh_bcfit)
                                                 
      if( fcheck('session.info') ) then   
        iusess = 40 
        open(unit=iusess,file='session.info',status='old'
     .          , iostat=ioerr )
        if( ioerr.ne.0 ) then
          call report_stat('FATAL','BCTOT','bctot'
     .         ,' ','Error opening session.info',ioerr)
        else
         call report_stat('STATIS','BCTOT','bctot'
     .         ,' ','Opened session.info',ioerr)
        endif
        check_flg = 2
        isessn = 1
        call rsesfo ( iusess,.false.,check_flg,yr,jdoy,isessn
     .              , hr,min,inter,nepoch,nssat,ssats  )
        close(unit=iusess)      
c       get month, day from day of yr
        call monday(jdoy,mon,day,yr)  
c       session.info info returns only minutes, assume seconds = 0
        sec = 0.d0
      else  
        call report_stat('FATAL','BCTOT','bctot',' '
     .          ,'No session.info, cannot get start times',0 )  
      endif
c     get PEP JD and sec-of-day for start, stop, and IC epoch
      jds = julday(mon,day,yr) 
      ts = 3600.d0*hr + 60.d0*min + sec
      jdf = jds               
      tf = ts 
      span = dble(inter*(nepoch)) 
      call timinc(jdf,tf,span)  
      jde = jds 
      te = ts
      call timinc(jde,te,span/2.d0)
      if(debug) print *,'start epoch stop from session.info '
     .                 ,jds,ts,jde,te,jdf,tf
c     convert Julian day of ephemeris epoch to year and day number
      call dayjul(jde,iye,idne)
      call ds2hms( iye,idne,te,ihe,imine,sece )
c     convert the ICs from GPS yr/day/hms to GPS week number and seconds
      call timcon( -4,iwkne,sowe,iye,idne,ihe,imine,sece,utcoff )
      write(iprnt,'(/,1x,a,i4,1x,i3,2x,2i3,f4.0,2x,i4,1x,f10.2)')
     .        'Midpoint of span for IC epoch: '
     .       , iye,idne,ihe,imine,sece,iwkne,sowe
               

c------ Read all ephemeris values into storage-------------------------
                       
      nesat = 0
      do i=1,maxsat
        iesat(i) = 0
        nbrd(i)  = 0
      enddo                            

c     read one record of the navigation file
      icall = 0
c     begin loop on navigation records
   20 call reade( iubc,icall,gnss
     .          , iflag,trans_sow,nprn,iwkn,ephem,bclock,subfr1 )
c      --reade returns GPS wk,sow in GPST, having converted Glonass UTC time
      icall = 1 
      if(debug)  print *,'BCTOT nprn,iwkn sow ephem iflag '
     .         ,nprn,iwkn,bclock(1),ephem,iflag
      if( iflag.eq.-1 ) goto 50
                     
c     If a PRN is not in the existing list add it unless the record is bad
      goodrec = .true.
      if( nprn.le.0.or.nprn.gt.32.or.iflag.eq.1 ) goodrec = .false.
      if( gnss.eq.'G' ) then 
        acheck = ephem(4)*ephem(4)
        if (dabs(acheck) .lt. 1.d7 .or. dabs(acheck) .gt. 1.d8) then
           write (message,'(a,i3,a,f7.0,a,d12.4)')
     .        'Bad BC ephemeris record   prn=',nprn,' sow='
     .       ,ephem(1),' a=',acheck
           call report_stat('WARNING','BCTOT','orbits/bctot',' '
     .                      ,message,0)
           goodrec = .false.
         endif 
       endif
       if( goodrec ) then 
c       if the PRN is not in the existing list, add it
        found = .false.
        if( nesat.gt.0 ) then
          do  isat= 1, nesat
            if( nprn.eq.iesat(isat) ) found = .true.
          enddo 
        endif
        if ( .not.found ) then
          nesat = nesat + 1   
cd          print *,'BCTOT nprn nesat ',nprn,nesat
          if( nesat.gt.maxsat) then
             write(message,30)nesat,maxsat
30           format('Number of satellites NESAT: ',i4,
     .       ' Exceeds program dimensions, MAXSAT:',i4)
             call report_stat('FATAL','BCTOT','orbits/bctot',' '
     .                        ,message,0)
          endif
          iesat(nesat) = nprn
        endif

c       Store the ephemeris message by satellite;    
        if( debug ) print *,'Storing messages, nesat = ',nesat
        do  isat = 1, nesat
          if( nprn.eq.iesat(isat) ) then
            nbrd(isat) = nbrd(isat) + 1  
            if( nbrd(isat).gt.maxbrd ) then
               write(message,40)nbrd(isat),iesat(isat),maxbrd
40            format('Number of broadcast values (',i3,
     .         ') for satellite PRN ',i3,' exceeds MAXBRD (',i4,')')
              call report_stat('FATAL','BCTOT','orbits/bctot',' '
     .                         ,message,0)
            endif 
            if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .          gnss.eq.'J'.or.gnss.eq.'I' ) then
              if(debug) print *
     .            ,'BCTOT isat nbrd iwkn ',isat,nbrd(isat),iwkn
              iewkn(nbrd(isat),isat) = iwkn
              do  i=1,16
                bephem(i,nbrd(isat),isat) = ephem(i)
              enddo 
            elseif( gnss.eq.'R') then
              do i=1,9
                bephem(i,nbrd(isat),isat) = ephem(i)
              enddo
              sowclk(nbrd(isat),isat) = bclock(1)
            endif  
          endif
        enddo

c     end if on valid navigation record; go read another
      endif
      goto 20

   50 continue
      call report_stat('STATUS','BCTOT','orbits/bctot',bcfile
     .,'Successfully read broadcast navigation file: ',0)
             
C------ Check the correspondence of E-file and X-file satellite --------
C       (messages displayed to screen in subroutine chksvs)

      if(upperc(xfile(1:1)).ne.'X' ) then
         call report_stat('WARNING','BCTOT','orbits/bctot',' ',
     .    'No X-file, using all satellites on the nav-file',0)
         nsat=nesat
         do i=1,nsat
            iexpnt(i) = i
         enddo
      else if (upperc(xfile(1:1)).eq.'X' ) then
         call chksvs( nesat,iesat,nsat,ixsat,iexpnt)
      endif 
      if( debug ) print *,'nsat iexpnt ',nsat,(iexpnt(i),i=1,nsat)

c       Move storage to the X-file slots

      do k=1,nsat 
        if( iexpnt(k).ne.0 ) then
          itsat(k)= iesat(iexpnt(k))
          ntbrd(k) = nbrd(iexpnt(k))
          if(debug) print *, 'isat iexpnt iesat nbrd '
     .       ,k,iexpnt(k),iesat(iexpnt(k)),nbrd(iexpnt(k))
          do j=1,ntbrd(k)
            itwkn(j,k)=iewkn(j,iexpnt(k))
            if(debug) print *,'isat ibrd iexpnt iewkn '
     .        , k,j,iexpnt(k),iewkn(j,iexpnt(k))
            do  i=1,16
              tephem(i,j,k)= bephem(i,j,iexpnt(k))
            enddo
            if( gnss.eq.'G'.or.gnss.eq.'C'.or.gnss.eq.'E'.or.
     .          gnss.eq.'I' ) then
              tsow(j,k) = bephem(1,j,iexpnt(k)) 
            elseif( gnss.eq.'R') then
              tsow(j,k) = sowclk(j,iexpnt(k))
            else
              write(message,'(a,a1,a)') 'GNSS ',gnss,'not yet supported'
              call report_stat('WARNING','BCTOT','orbits/bctot',' '
     .                    ,message,0) 
            endif
          enddo
        endif
      enddo
      write(iprnt,'(a)') 'Number of records for each SV' 
      do isat = 1,nsat
        write(iprnt,'(i3,i5)') itsat(isat),ntbrd(isat)      
        if(debug ) then 
           print *,'itsat nbrd ',itsat(isat),ntbrd(isat)
           do j=1,ntbrd(i)
             print *,bephem(1,j,isat)
           enddo
        endif 
      enddo

c------ Construct the satellite name from the GNSS and PRN 
                           
      do i=1,nsat
* MOD TAH 190702: Added antpwr to snav_read call
         call svnav_read( -1,iye,idne,ihe,imine,gnss,itsat(i),
     .        isvn,frqchn,antbody,sbmass,yawbias,yawrate, antpwr, 
     .        svnstart,svnstop )
cd        print *,'i iyr idoy ihr imin gnss itsat isvn '
cd     .         , i,iyr,idoy,ihr,imin,gnss(i),itsat(i),isvn
cd        print *,'antbody',antbody
        if( gnss.eq.'G' ) then
c         GPS:  e.g, G23   53 IIR     
* MOD TAH 200606: Made I2.2 and flexiable a format to allow
*         possible use of all antbody length (applied to all
*         writes below).  
          write(satnam(i),'(a1,i2.2,3x,i2,1x,a)')
     .              'G',itsat(i),isvn,trim(antbody(7:))
        elseif( gnss.eq.'R' ) then 
c         GLONASS: e.g.  R23  701 M
          write(satnam(i),'(a1,i2.2,2x,i3,1x,a)')
     .              'R',itsat(i),isvn,trim(antbody(9:)) 
        elseif( gnss.eq.'E' ) then 
c         GALILEO: e.g.  E23 0201 0A
          write(satnam(i),'(a1,i2.2,1x,i4,1x,a)')
     .              'E',itsat(i),isvn,trim(antbody(9:)) 
        elseif( gnss.eq.'C' ) then 
c         BEIDOU: e.g.  C23      2I
          write(satnam(i),'(a1,i2.2,1x,i4,1x,a)')
     .              'C',itsat(i),isvn,trim(antbody(8:)) 
        elseif( gnss.eq.'I' ) then
c         IRNSS: e.g.   I07    7 IGESO
          write(satnam(i),'(a1,i2.2,3x,i2,1x,a)')
     .              'I',itsat(i),isvn,trim(antbody(7:)) 
        else
          write(message,'(a,a1,i2.2 )')  
     .        'Satellite name not defined for ',gnss,itsat(i)
          call report_stat('WARNING','BCTOT','orbits/bctot',' '
     .                    ,message,0) 
        endif
      enddo                          


c------ Display the satellites selected and their reference times at the IC epoch-
c       (convert Glonass times to GPST)
                                              
      if( debug ) 
     .   write(*,'(/,2a)')
     .    ' PRN        Message reference epoch at IC epoch '
     .   ,'     Offset from IC epoch (hrs) '
      write(iprnt,'(/,a,/,2a)')
     .   '------------------------------------------------------------'
     .   ,' PRN        Message reference epoch at IC epoch '
     .   ,'     Offset from IC epoch (hrs) '   
      if( debug ) print *,'Getting record closest to IC epoch '
      do isat=1,nsat                         
         if(debug) print *,'isat ntbrd iwkne sowe  '
     .                      ,isat,ntbrd(isat),iwkne,sowe
         icbrd = closest_epoch ( ntbrd(isat),itwkn(1,isat),tsow(1,isat)
     .                        , iwkne,sowe,offset_hrs )
         call timcon( 4,itwkn(icbrd,isat),tsow(icbrd,isat)
     .              , iye,idne,ihe,imine,sece,utcoff )
         if(debug) print *,'icbrd itwkn tsow iye idne ihe imine sece '
     .        ,icbrd,itwkn(icbrd,isat),tsow(icbrd,isat)
     .        , iye,idne,ihe,imine,sece 
         if( debug ) 
     .   write(*,'(i3,4x,i4,2x,f8.0,3x,2i4,2x,2i3,f6.2,6x,f8.2 )')
     .             itsat(isat),itwkn(icbrd,isat),tsow(icbrd,isat)
     .           , iye,idne,ihe,imine,sece
     .           , offset_hrs    
         write(iprnt,'(1x,i3,4x,i4,2x,f8.0,3x,2i4,2x,2i3,f6.2,6x,f8.2)')
     .             itsat(isat),itwkn(icbrd,isat),tsow(icbrd,isat)
     .           , iye,idne,ihe,imine,sece
     .           , offset_hrs
      enddo
        

C------- Get the T-file header information and write the header--------------------------

C        Compute a set of Cartesian 'initial conditions' for the T- and G-files

      do isat = 1,nsat
        icbrd = closest_epoch ( ntbrd(isat),itwkn(1,isat),tsow(1,isat)
     .                        , iwkne,sowe,offset_hrs )
        if( debug ) then
           print *,'IC epoch index for PRN ',itsat(isat),icbrd
           print *,'itwkn tsow ',
     .      itwkn(icbrd,isat),tephem(1,icbrd,isat)
          print *,'tephem: ',(tephem(i,icbrd,isat),i=1,16)
        endif 
        call brdxyz ( iwkne,sowe,itwkn(icbrd,isat),tephem(1,icbrd,isat)
     .              , satprm,ibcerr,itsat(isat))  
        do i=1,6
           satics(i,isat)=satprm(i)
        enddo
c**     velocites are not dynamically consistent: error in brdxyz or rotcrd?
c**     Interpolate to get values for G-file
C       Default radiation pressure parameters
        if(nics.gt.6) then
          satics(7,isat) = 1.d0
          do j = 8,15
            satics(j,isat)= 0.d0
          enddo
        endif
      enddo
                   

c--------Write the  Earth-fixed and inertial T-files

      eframe = 'EFIXD'
      nintrs = 3
      nepoch= ( (jdf-jds)*86400.d0 + (tf-ts) ) / delt + .000001d0 + 1
      write(message,'(a,a16,a)')
     .              ' Earth-fixed T-file ',tfilef,' being created'   
      call report_stat('STATUS','BCTOT','orbits/bctot',' ',message,0)
      if( debug ) then
        print *,'BCTOT models ',precmod,nutmod,gravmod,srpmod
     .               ,eradmod,antradmod 
        print *,'     delt nepoch ',delt,nepoch     
      endif
      call thdrit( iutout,jde,te,jds,ts,jdf,tf,delt,nepoch,nintrs
     .           , nsat,gnss,itsat,satnam,satics,nics,bcfile,icsnam
     .           , precmod,nutmod,gravmod,eframe,srpmod,eradmod
     .           , antradmod )
c     convert julian day of initial epoch to year and day number
      call dayjul(jds,iys,idns)
c     write the header for BC ephemeris offsets
      write(iprnt,'(/,a)') '-------------------------------------------'
      write(iprnt,'(/,a,/)')' ** BC ephemeris epoch offsets in hours **'

c     loop over all epochs
      do iepoch=1,nepoch

         if(iepoch.eq.1) then
c           convert start (GPST) epoch to GPST week number and seconds
            call ds2hms( iys,idns,ts,ihs,mins,secs )
            call timcon( -4,iwkn,sow,iys,idns,ihs,mins,secs,utcoff )
         else
c           add delta seconds
            call secsum(iwkn,sow,delt,iwkn,sow)
         end if
c        write the epoch number and time to the print file
         call timcon(4,iwkn,sow,iys,idns,ihs,mins,secs,utcoff )
         write(iprnt,'(1x,a,3x,i4,3x,2i4,2i3,f4.0,i6,f10.2)') 
     .     ' PRN   Offset     Epoch: '
     .        ,iepoch,iys,idns,ihs,mins,secs,iwkn,sow

c       loop over satellites
        do isat = 1,nsat
          icbrd = closest_epoch ( ntbrd(isat),itwkn(1,isat),tsow(1,isat)
     .                          , iwkn,sow,offset_hrs )
          call brdxyz ( iwkn,sow,itwkn(icbrd,isat),tephem(1,icbrd,isat)
     .                , satprm,ibcerr,ntbrd(isat) )  
c         if (isat.eq.1 ) then
c          print *,'iwkne sowe itwkn icbrd ',iwkne,sowe,itwkn(icbrd,isat)
c     .                ,icbrd
c          print *,'tephem ',(tephem(i,icbrd,isat),i=1,16) 
c          print *,'satprm ',satprm
c         endif
c         Store satellite position
          x(1,isat)=satprm(1)
          x(2,isat)=satprm(2)
          x(3,isat)=satprm(3)
c         write offsets to the print file
          aflag = ' '
          if( dabs(offset_hrs).gt.2.d0) aflag = '*'
          write(iprnt,'(2x,i2,f8.2,1x,a1)') itsat(isat),offset_hrs,aflag
c       end loop on satellites
        enddo

c       write the T-file record for this epoch
        write(iutout) ((x(i,isat),i=1,nintrs),isat=1,nsat)

c     end loop on epochs
      enddo

      call report_stat('STATUS','BCTOT','orbits/bctot',tfilef
     .,'Successfully wrote earth-fixed T-file: ',0)

c     trot and gmake will open and close the files
      call closeb(iutin,iutout,inut,iut1,ipole)
     
c     Now write out an inertial tfile if required
      if( frame.ne.'EFIXD' ) then
        idir = 1  
        call trot(tfilef,tfile,idir,frame)
c       skip this since reported in trot:
c        call report_stat('STATUS','BCTOT','orbits/bctot',tfile
c     .  ,'Successfully wrote inertial T-file: ',0)
      endif
     

c------- Generate a G-file at the central epoch-------------------------------------

      prog = 'BCTOT '
      if( frame.eq.'EFIXD' )  then
         call report_stat('WARNING','BCTOT','orbits/bctot',' ',
     .   'Earth fixed G-file is being created ',0)
         call gmake( prog,tfilef,gfilef,itsat )
       else  
         call gmake( prog,tfile,gfile,itsat )
         call arcinp( tfile,gfile,jds,ts,jdf,tf,stepsize )
      endif

      call report_stat('STATUS','BCTOT','orbits/bctot',' ',
     .    'Normal end in BCTOT',0)

      stop
      end                                                                    

c############################################################################

      Subroutine get_models( frame, precmod, nutmod, gravmod, srpmod
     .                     , eradmod,antradmod
     .                     , stepsize, delt, nics, icsnam )


      implicit none
                                
      include '../includes/dimpar.h'
      include '../includes/orbits.h'

      logical reqd,fcheck

      character*4  icsnam(maxorb)
      character*4  ics_sphrc(9),ics_srxyz(9),ics_berne(15),ics_bern2(12)
      character*5  frame,precmod,nutmod,gravmod,srpmod,eradmod,antradmod
     .             ,buf5
      character*8  buf8  
      character*80 nut_header
      character*256 message

      integer*4 mchkey,iusest,nics,ioerr,ill,i
                                 
      real*8 delt,stepsize
                   

c      data icsflg/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
c     .             ,'RAD1','RAD2','RAD3','    ','    ','    '
c     .             ,'    ','    ','    '/
      data ics_sphrc/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     .              ,'DRAD','YRAD','ZRAD'/
      data ics_srxyz/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     .              ,'XRAD','YRAD','ZRAD'/
      data ics_berne/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     .              ,'DRAD','YRAD','BRAD','DCOS','DSIN','YCOS'
     .              ,'YSIN','BCOS','BSIN'/ 
      data ics_bern2/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     .              ,'DRAD','YRAD','BRAD','XCN1','XCN3','ZCN1'/
c

          
c     Unit number for the sestbl.
      iusest = 41
         

c     If a sestbl. is available, read it to get the values

      if (fcheck('sestbl.')) then
        open(iusest,file='sestbl.',status='old',iostat=ioerr)
        if ( ioerr .ne. 0) then
          call report_stat('FATAL','BCTOT','orbits/bctot',' '
     .                    ,'Error opening sestbl. ',ioerr)
        endif    
        reqd = .false.
        call rdsest(14,'Inertial frame',5,buf5,iusest,reqd,ill) 
        if( ill .ne. 0 ) then   
            write(message,'(a,i3)') 'Error reading sestbl. : ',ill
            call report_stat('FATAL','BCTOT','orbits/bctot',' '
     .                       ,message,0) 
        endif
        if ( buf5.ne.'     ' ) then
          frame = buf5
        else
          frame = 'J2000'
        endif
        if(frame.ne.'B1950'.and.frame.ne.'J2000') 
     .      call report_stat('FATAL','BCTOT','orbits/bctot',frame,
     .        'Invalid inertial frame in sestbl.: ',0)
        reqd = .false.
        call rdsest( 25, 'Inertial Reference System', 5, buf5, 21,
     .               reqd, ill)
        if ( ill.ne.0 ) call report_stat('WARNING','BCTOT',
     .          'orbits/bctot', 'sestbl.',
     .          'Using default Inertial Referecne Frame: IAU76',0)
        if (buf5.ne.'     ') then 
            precmod = buf5 
            nutmod = buf5
        else
            if(frame.eq.'B1950') precmod = 'IAU68'
            if(frame.eq.'J2000') precmod = 'IAU76' 
        endif    
        reqd = .false.
        call rdsest(24,'Reference System for ARC',5,buf5,iusest,reqd
     .              ,ill)
        if( ill .ne. 0 ) then
           call report_stat('FATAL','BCTOT','orbits/bctot',' '
     .                      ,'Error reading sestbl. : ',ill)
        endif
        if ( buf5.ne.'     ' ) then
           gravmod = buf5
        else 
           gravmod = 'EGM08'
        endif
        call rdsest(23,'Radiation Model for ARC',5,buf5,iusest,reqd,ill)
        if( ill .ne. 0 ) then   
            write(message,'(a,i3)') 'Error reading sestbl. : ',ill
            call report_stat('FATAL','BCTOT','orbits/bctot',' '
     .                       ,message,0)
        endif    
        if ( buf5.ne.'     ' ) then
          srpmod = buf5       
        else
          srpmod = 'BERNE'
        endif
        call rdsest(21,'Earth radiation model',5,buf5,iusest,reqd,ill)
        if( ill .ne. 0 ) then   
            write(message,'(a,i3)') 'Error reading sestbl. : ',ill
            call report_stat('FATAL','BCTOT','orbits/bctot',' '
     .                       ,message,0)
        endif
        if ( buf5.ne.'     ' ) then
          eradmod = buf5     
        else
          eradmod = 'NONE '
        endif        
        call rdsest(20,'Antenna thrust model',5,buf5,iusest,reqd,ill)
        if( ill .ne. 0 ) then   
            write(message,'(a,i3)') 'Error reading sestbl. : ',ill
            call report_stat('FATAL','BCTOT','orbits/bctot',' '
     .                       ,message,0)
        endif
        if ( buf5.ne.'     ' )  then
          antradmod = buf5   
        else
          antradmod = 'NONE '
        endif
        call rdsest(24,'Tabular interval for ARC',8,buf8,iusest,reqd
     .             ,ill)
        if( ill .ne. 0 ) then    
          write(message,'(a,i3)') 'Error reading sestbl; ill=',ill
          call report_stat('FATAL','BCTOT','orbits/bctot',' ',message,0)
        endif                
        if ( buf8.ne.'        ' )  then
           read(buf8,'(f8.0)',iostat=ioerr) delt
           if( ioerr.ne.0 ) call report_stat('FATAL','BCTOT'
     .           ,'orbits/bctot',' ','Error reading delt',ioerr)
        else
           delt = 900.d0
        endif
        call rdsest(8,'Stepsize',8,buf8,iusest,reqd,ill)
        if( ill .ne. 0 ) then   
          write(message,'(a,i3)') 'Error reading sestbl. : ',ill
          call report_stat('FATAL','BCTOT','orbits/bctot',' ',message,0)
        endif
        if ( buf8.ne.'        ' ) then
           read(buf8,'(f8.2)',iostat=ioerr) stepsize
           if( ioerr.ne.0 ) call report_stat('FATAL','BCTOT'
     .           ,'orbits/bctot',' ','Error reading stepsize',ioerr)
        else
          stepsize = 75.d0 
        endif         
      else
c     No sestbl. available
        call report_stat('WARNING','BCTOT','orbits/bctot',' ',
     .    'No sestbl. available, set default models',ioerr)
        frame = 'J2000'
        if( frame(1:1).eq.' ' ) frame = 'J2000'  
        gravmod = 'IGS92'
        srpmod  = 'SPHRC'
        eradmod = 'NONE '
        antradmod = 'NONE '
        stepsize = 75.d0
        delt = 900.d0
      endif    
c     assign a precession model: J2000 -> IAU76, B1950 -> IAU68
c precmod already set appropriately above in Inertial Reference Frame section of sestbl. read
c      if(frame.eq.'J2000') precmod = 'IAU76'
c      if(frame.eq.'B1950') precmod = 'IAU68' 
c     get the nutation model from the nutabl file header if nutabl used
      if( .not.fcheck('nbody') ) then
        open( unit=20,file='nutabl.',status='old',iostat=ioerr)
        if( ioerr.ne.0 ) 
     .    call report_stat('FATAL','BCTOT','orbits/bctot','nutabl.',
     .      'Error opening nutation table: ',ioerr)
        read(20,'(a)',iostat=ioerr) nut_header
        if( ioerr.ne.0) call report_stat('FATAL','BCTOT','orbits/bctot'
     .     ,' ','Error reading first line of nutation file to get model'
     .     ,ioerr)
        close(20)  
* MOD TAH 200303: Changed default IAU0A which is fixdrv default. (Override 
*        sestbl. if needed
        if(mchkey(nut_header,'IAU20',80,5).gt.0) then
          nutmod='IAU0A'
        else
          nutmod='IAU80'
        endif  
      else
        if ( precmod .eq. 'IAU68' ) then
          nutmod = 'IAU80'
        else if  ( precmod .eq. 'IAU76' ) then
c       IAU00 will get nutations from MHB_2000; IAU0A uses the SOFA routines.
* MOD TAH 200303: Changed default IAU0A which is fixdrv default. (Override 
*        sestbl. if needed
          nutmod = 'IAU0A'
        endif
      endif

c     Load IC names according to Radiation Model for ARC sestbl. entry

      if (srpmod .eq. 'SPHRC' .or. srpmod .eq. 'SRDYZ') then
        nics = 9
        do i=1,nics
          icsnam(i) = ics_sphrc(i)
        enddo
      elseif (srpmod .eq. 'SRDYB') then
        nics = 9
        do i=1,nics
          icsnam(i) = ics_berne(i)
        enddo
      elseif (srpmod .eq. 'SRXYZ') then
        nics = 9
        do i=1,nics
          icsnam(i) = ics_srxyz(i)
        enddo
      elseif (srpmod .eq. 'BERNE'.or. srpmod .eq. 'BERN1'.or.
     .        srpmod .eq. 'UCLR1'  ) then
        nics = 15
        do i=1,nics
          icsnam(i) = ics_berne(i)
        enddo         
      elseif (srpmod .eq. 'BERN2') then  
        nics = 12
        do i=1,nics
          icsnam(i) = ics_bern2(i)
        enddo                               
      endif

      return
      end


 


