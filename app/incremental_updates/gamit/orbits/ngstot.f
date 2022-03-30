Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1994. All rights reserved.
      Program NGSTOT
C
C Written by R.W. King & Yehuda Bock  March 1988
C
C Read an NGS Standard Product emphemeris file and create a T-file
C in an earth-fixed system or an inertial system
C
c  mod mar 1990 by pch - remove call to ngsprg (program header bs)
c  mod jul 1990 by ss  - remove argument FIRST from THDRED
c  mod oct 1990 by ss
c  mod nov 1991 by rwk to make similar to TTONGS, particularly as
c      pertains to inputing start, stop times
c  mod mar 1993 by peng - add output SP#3 format
c  mod apr 1995 by pt - hardwire models for tfile headers as follows:
c                       precession : IAU76
c                       nutation   : IAU80
c                       gravity    : (blank)
c  mod may 1995 by pt - user decides on B1950 or J2000 inertial frame. This
c      assigns the precession model IAU68 for B1950, IAU76 for J2000
c  mod rwk 2006  correct from a center-of-solid-Earth (CE) to center-of-mass 
c                (CM) system, considering ocean tides
c  mod rwk 2014  remove the option to input an x-file, and replace the
c                command-line argument with a gnss selection code
c  mod lei 2015  read frame and gravity from sestbl. 
c  mod rwk 2015  allow reading of any (single) GNSS set of SVs 
c  mod rwk 2018  don't open the nutation file if not needed 
c  mod rwk Apr 2019 Change the initial t-/g- file to use ECOM1 rather than SPHRC 

      implicit none

      include '../includes/dimpar.h'

      integer*4 jds,jdf,jde,mjds,iyear,imonth,iday,ihr,imin,iyr,nblen
     .        ,  numsp3sv,numsat,itsat,msat,nepchs
     .        ,  nintrs,iungs,iut,iux,idir,icall,nics,iscrn,ierr
     .        , iy,iy1,idoy,idoy1,iweek,idow,julday,frqchn,svnstart(5)
     .        , svnstop(5),iarg,iclarg,isvn,notl,ioerr,i,j,k
      integer*4 ill   ! added by lei
      logical reqd    ! added by lei
      character*5 buf5 ! added by lei

c     function
      integer*4 mchkey

      real*8 sec, delts, frcts, dt, te, tf, ts, a, pjd
      real*8 satics(maxorb,maxsat),accsat(maxsat),pvsigb,clksigb
     .     , docmc(3),sbmass,yawrate
      real*8 antpwr   ! Antenna transmit power (W)
      

      character*1 lowerc,upperc,spver,gnss_sel   
     .          , yawbias
      character*2 linesym
      character*4 icstfl(6),srpnam(13),icsnam(maxorb)
      character*5 precmod,nutmod,gravmod,frame,srpmod,eframe
     .          , eradmod,antradmod
      character*6 prog
      character*8 otlmod
      character*16 spfile,tfile,gfile,tfilef,satnam(maxsat)
      character*20 antbody 
      character*80 nut_header
      character*120 version
      character*256 message

      logical sat_ok(maxsat),fcheck
 
      logical debug/.false./

      dimension itsat(maxsat)

      data icstfl/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'/
      data nutmod/'     '/,gravmod/'     '/,srpmod/'     '/
     .   , eradmod/'     '/,antradmod/'     '/

c Initialization
c     Screen
      iscrn = 6
c     X-file
      iux = 14
c     SP-file
      iungs = 15
c     T-file
      iut = 16  
      spfile = " "
      tfile  = " "

c Get version number

      call oversn(version)
c*      write(iscrn,'(a)')' '
      write(message,5) version
    5 format('Started NGSTOT ',a120)
      call report_stat('STATUS','NGSTOT','orbits/ngstot',' ',message,0)
        
c  Read the input

c  if there are no command-line arguments, echo the help 
c
  
      iarg = iclarg(1,spfile) 
      if( iarg.le.0 ) then
         write(*,'(/a,/,a,/)')
     .       'ngstot  [sp3file]  [t-file]  [gnss-code] '
     .      ,'        (required) (required)  (optional)'
* MOD TAH 200606: Stop rather than fatal generated below when no
*     arguments are given
         stop 'ngstot: No arguments given'
      else
c       file name must be standard format, without full path for logic in NGSTOT
        if( nblen(spfile).gt.16 ) call report_stat('FATAL','NGSTOT'
     .       , 'orbits/ngstot',' '
     .       , 'SP3 file name too long --cannot use full path',0)
        iarg = iclarg(2,tfile)                                  
        if( iarg.le.0 )  
     .           call report_stat('FATAL','NGSTOT','orbits/ngstot',' '
     .              ,'Missing command-line argument for T-file name',0)
        iarg = iclarg(3,gnss_sel)
        if( iarg.le.0 ) then 
          gnss_sel = 'G'
        else
          gnss_sel = upperc(gnss_sel) 
        endif
      endif
 

c  Report the command line
 
      write(message,'(a,a16,a,a16,a,a1)') 'Reading ',spfile
     .     ,' to write ',tfile,' for GNSS = ',gnss_sel
      call report_stat('STATUS','NGSTOT','orbits/ngstot',' ',message,0)
 

c Open the SP3 and T- files

      call ngsopn( spfile,tfile,tfilef,iungs,iux,iut )
         
c Check for consistency between t-file day-of-year and sp3-file day-of-week

c     do only if file names have the dates in the usual format: orgwwwwd.sp3 
      read(spfile(4:8),'(i4,i1)',iostat=ierr,err=10) iweek,idow
      read(tfile(6:6),'(i1)',iostat=ierr,err=10) iy
      read(tfile(8:10),'(i3)',iostat=ierr,err=10) idoy  
c     branch to 10 with no error issued if the file format is non-standard
      if( (iweek.gt.0.and.iweek.lt.2000) .and. 
     .    (idow.ge.0.and.idow.le.6)  .and.
     .    (iy.ge.0.and.iy.le.9) .and.
     .    (idoy.gt.0.and.idoy.le.366) ) then
        idoy1 = 0
        call doygwk(idoy1,iy1,iweek,idow) 
        iy1 = mod(iy1,10)
        if( iy1.ne.iy .or. idoy1.ne.idoy ) then
          write(message,'(a,1x,i1,i4,a,1x,i1,i4)') 
     .        'SP3 filename implies yr,doy = ',iy1,idoy1
     .       ,' but T-file name implies yr,doy = ',iy,idoy
          call report_stat('WARNING','NGSTOT','orbits/ngstot',' '
     .                     ,message,0)
        endif
      endif

c Determine the file format using the first two characters of the first line

c   SP1    ' # '
c   SP3-a  '# '  or '# a'
c   SP3-b  '#b'
c   SP3-c  '#c'

   10  read(iungs,'(a2)') linesym
      if ( linesym.eq." #") then  
        spver = '1'
      elseif ( linesym.eq.'# '.or.linesym.eq.'#a' ) then
        spver = 'a'
      elseif ( linesym.eq.'#b' ) then
        spver = 'b'
      elseif ( linesym.eq.'#c' ) then
        spver = 'c'         
      elseif ( linesym.eq.'#d' ) then 
        spver = 'd' 
      else
        call report_stat('FATAL','ORBITS','ngstot',' '
     .                  ,'Unrecognized SP orbit format',0 )
      endif
      rewind(iungs)

c Clean satics array
	do i=1,maxorb
		do j=1,maxsat
			satics(i,j)=0.0
		enddo
	enddo

c Read the SP file header to get times and satellites

      if (spver.eq."1") then
        call rsp1hd( iungs,iyear,imonth,iday,ihr,imin,sec,delts
     .             , mjds,frcts,nepchs,numsat,itsat )

      else
                       
        if( debug ) print *,'NGSTOT calling rsp3hd '
        call rsp3hd( iungs,gnss_sel
     .             , iyear,imonth,iday,ihr,imin,sec
     .             , delts,mjds,frcts,nepchs
     .             , numsp3sv,numsat,itsat,accsat
     .             , pvsigb,clksigb,otlmod )

      endif 
      if( debug ) then 
        print *,'NGSTOT gnss_sel numsp3sv numsat clksigb otlmod '
     .              ,   gnss_sel,numsp3sv,numsat,clksigb,otlmod
        print *,' itsat',(itsat(i),i=1,numsat)
        print *,'       accsat  ',(accsat(i),i=1,numsat)
      endif

c Convert SP times to PEP JD
                
      ts = dble(ihr)*3600.d0 + dble(imin)*60.d0 + sec  
      jds = julday(imonth,iday,iyear)
c**    no longer convert T-file times to UTC - they are now GPST whenever written
c      gpsutc= taiutc (jds) - 19.d0
c      call timinc( jds,ts,-gpsutc )
c      avoid roundoff by assuming that the epochs are even seconds
      ts = anint(ts)
       
c Compute the end time for the T-file header
      dt = delts*(nepchs-1)
      jdf= jds
      tf= ts
      call timinc( jdf,tf,dt )
c       Avoid roundoff by assuming that the epochs are even seconds
      tf = anint(tf)

c Compute the initial condition epoch as the middle of the observation span
      dt = delts*(nepchs/2)
      jde= jds
      te= ts
      call timinc( jde,te,dt ) 

c If initial condition time is within 1 hour of midday make it exactly 12:00
      if ( dabs(te-43200.d0) .le. 3600.d0 ) then
        te = 43200.d0
      endif 

c Avoid roundoff by assuming that the epochs are even seconds
      te = anint(te)
             
c If ocean-tidal loading correction used, read in the coefficients
c  (saved in /otlcmc; docmc not calculated during this call)

      if( otlmod(1:1).ne.' '.and.otlmod(1:4).ne.'NONE') then                         
c      hard-wire # components
       notl = 11                             
       call otlcmc( jds,ts,otlmod,notl,1,docmc )
       write(message,'(a,a8)') 
     .   'Converting CE to CM using otlcmc.dat offsets for ',otlmod 
       call report_stat('STATUS','NGSTOT','orbits/ngstot',' ',message,0)
      endif

c Loop thru the SP file to get initial conditions near the center of the span
          
cd      print *,'calling GETICS jde te ',jde,te 
      call getics( iungs,spver,jde,te,gnss_sel
     .           , numsp3sv,numsat,itsat,satics )              
      if( debug ) then
        write(*,'(a,i3,100i3)') 
     .   'In NGSTOT numsp3sv ', numsp3sv
        write(*,'(a,i3,100i3)') 
     .   'In NGSTOT numsat itsat ', numsat,(itsat(i),i=1,numsat)
        write(*,'(6(/,1x,3f11.3,3f9.4))') 
     .     ((satics(i,j),i=1,6),j=1,numsat)
      endif
  

c It's possible that an SV has been added to the sp3 file with bogus positions
c (usually 0.0); check here the ICs and the eliminate bad SVs from the T-file
c RWK 150922: Allow for high-altitude SVs 
                       
      do i = 1, numsat 
        sat_ok(i) = .true.
        a = dsqrt(satics(1,i)**2 + satics(2,i)**2 + satics(3,i)**2)  
c       if( a.le.2.d4 .or. a.ge.3.d4 ) then 
        if( a.le.2.d4 .or. a.ge.5.d4 ) then
          write(message,'(a,i2,a,i2,a)') 'Bad ICs for sat ',i
     .                   ,' (PRN',itsat(i),')--omit from g/t files'
          call report_stat('WARNING','NGSTOT','orbits/ngstot',' '
     .                    ,message,0) 
          sat_ok(i) = .false.   
        endif
      enddo  
      msat = numsat
      i = 1         
      if( debug ) then
        print *,'orig IC list numsat ',numsat
        do j=1,numsat
          print *,'j itsat sat_ok satics(1) '
     .            ,j,itsat(j),sat_ok(j),satics(1,j)
        enddo
      endif
      do while ( i.le.msat ) 
         if(.not.sat_ok(i) ) then 
            msat = msat - 1
            do j = i,numsat-1
               itsat(j) = itsat(j+1) 
               sat_ok(j)=sat_ok(j+1)
               do k = 1,6
                 satics(k,j) = satics(k,j+1)
               enddo
            enddo   
         endif  
      i = i + 1
      enddo
      numsat = msat   
      if( debug ) then
        print *,'after ok check numsat ',numsat
        do j=1,numsat      
          print *,'j itsat sat_ok satics(1) '
     .            ,j,itsat(j),sat_ok(j),satics(1,j)
        enddo
      endif

          
c Construct the 16-character satellite names

                  
      pjd = jde + te/86400.d0
      call pjdhms( pjd,iyr,idoy,imonth,iday,ihr,imin )
      do i=1,numsat                              
cd        print *,'calling SVNAV_READ i gnss itsat iyr idoy ihr imin '
cd     .                        ,i,gnss_sel,itsat(i),iyr,idoy,ihr,imin
* MOD TAH 190702: Added place holder for antpwr to snav_read call
        call svnav_read( -1,iyr,idoy,ihr,imin,gnss_sel,itsat(i),
     .       isvn,frqchn,antbody,sbmass,yawbias,yawrate,antpwr, 
     .       svnstart,svnstop )
cd        print *,'i iyr idoy ihr imin gnss itsat isvn '
cd     .         , i,iyr,idoy,ihr,imin,gnss_sel,itsat(i),isvn
cd        print *,'antbody',antbody
        if( gnss_sel.eq.'G' ) then
c         GPS:  e.g, G23   53 IIR 
* MOD TAH 200606: Made I2;2 and flexiable a format to allow
*         possible use of all antbody length (applied to all
*         writes below).   
* MOD TAH 201031: Added explict end to antbody string to avoid
*         runtimes errors on some compilers. 
          write(satnam(i),'(a1,i2.2,3x,i2,1x,a)')
     .              'G',itsat(i),isvn,trim(antbody(7:13))
        elseif( gnss_sel.eq.'R' ) then 
c         GLONASS: e.g.  R23  701 M
          write(satnam(i),'(a1,i2.2,2x,i3,1x,a)')
     .              'R',itsat(i),isvn,trim(antbody(9:15)) 
        elseif( gnss_sel.eq.'E' ) then 
c         GALILEO: e.g.  E23 0201 0A
          write(satnam(i),'(a1,i2.2,1x,i4,1x,a)')
     .              'E',itsat(i),isvn,trim(antbody(9:15))
        elseif( gnss_sel.eq.'C' ) then 
c         BEIDOU: e.g.  C23      2I
          write(satnam(i),'(a1,i2.2,1x,i4,1x,a)')
     .              'C',itsat(i),isvn,trim(antbody(8:14)) 
        else
          write(message,'(a,a1,i2.2 )')  
     .        'Satellite name not defined for ',gnss_sel,itsat(i)
          call report_stat('WARNING','NGSTOT','orbits/ngstot',' '
     .                    ,message,0) 
        endif
      enddo

c Set the scratch T-file to be earth-fixed

      eframe = 'EFIXD'

c Write the T-file header records
                           
c     positions only - no velocities or partials
      nintrs = 3  
      do i=1,6
        icsnam(i) = icstfl(i)
      enddo         
      srpmod = 'ECOM1' 
      call assign_srpnames(srpmod,nics,srpnam)
      do i=7,15
        icsnam(i) = srpnam(i-6)
      enddo 
 
c       write(*,'('using ',a5,' and ',a5)')

c     put the ICs into a CM system if ocean-tide CMC correction used in sp3
      if( otlmod(1:1).ne.' ') then
        call otlcmc( jde,te,otlmod,notl,2,docmc )
        if( debug ) then  
          write(6,'(a,9x,3f8.1)') 'Otide CM-CE (m) at IC epoch:  '
     .          , (docmc(i),i=1,3)
c         this corrects from CE to CM, but the sign that works (+)  seems
c         to be the opposite of what Scherneck says on his web page. 
c         Gerd Gendt (GFZ) confirms this.
          do j=1,numsat
            do i=1,3
              satics(i,j) = satics(i,j) + docmc(i)/1.d3
            enddo
          enddo
        endif          
      endif

c====================== commented out by L. Wang ===================================
cc     set the reference frame and gravity model 
c      frame = 'J2000'
c      precmod = 'IAU76'       
c      gravmod = 'EGM96'  
c======================= end-of-comment ============================================
cxxxxxxxxxxxxxxxxxxxxxx  added by lei xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
c     set the reference frame and gravity model by reading sestbl.
c     -- set the default values
* MOD TAH 200304: Changed default values to EGR08 and IAU0A from EGM08 abd IAU76.
      frame = 'J2000'
      gravmod = 'EGR08'                
      precmod = 'IAU0A'
c     -- open sesstion table
      open( unit=21, file='sestbl.',status='old',iostat=ioerr)
      if( ioerr.ne.0 ) then
        call report_stat('WARNING','NGSTOT','orbits/ngstot','sestbl.',
     .    'Errors in opening session table, assuming defaults',ioerr)

      else

c     -- read frame and set precmod
        reqd = .false.
        call  rdsest( 14, 'Inertial frame', 5, buf5, 21,
     .                reqd, ill )
        if ( ill.ne.0 ) call report_stat('WARNING','NGSTOT',
     .          'orbits/ngstot', 'sestbl.','Errors in sestbl entries',0)
        if(buf5.ne.'     ') frame = buf5
        if(frame.ne.'B1950'.and.frame.ne.'J2000') then  
          call report_stat('WARNING','NGSTOT','orbits/ngstot','sestbl.'
     .       ,'Invalid inertial frame, assuming J2000',0)  
          frame = 'J2000';  
        endif
        frame(1:1) = upperc(frame(1:1))
        if(frame.eq.'B1950') precmod = 'IAU68'
* MOD TAH 200304: Default J2000 to IAU0A precmod.
*       if(frame.eq.'J2000') precmod = 'IAU76'
        if(frame.eq.'J2000') precmod = 'IAU0A'
c     -- read gravity
        reqd = .false.
        call rdsest( 24, 'Reference System for ARC', 5, buf5, 21,
     .               reqd, ill)
        if ( ill.ne.0 ) call report_stat('WARNING','NGSTOT',
     .          'orbits/ngstot', 'sestbl.','Errors in sestbl entries',0)
        if(buf5.ne.'     ') gravmod = buf5
        if( gravmod.ne.'WGS84' .and. gravmod.ne.'WGS72' .and.
     .      gravmod.ne.'MERIT' .and. gravmod.ne.'IGS92' .and.
     .      gravmod.ne.'EGM96' .and. gravmod.ne.'EGM08' .and.
     .      gravmod.ne.'EGR08' )  then
           call report_stat('WARNING','NGSTOT','orbits/ngstot','sestbl.'
     .      ,'Invalid Reference System for ARC, assuming EGM08',0)    
            gravmod = 'EGM08'
        endif
        reqd = .false.
        call rdsest( 25, 'Inertial Reference System', 5, buf5, 21,
     .               reqd, ill)
        if ( ill.ne.0 ) call report_stat('WARNING','NGSTOT',
     .          'orbits/ngstot', 'sestbl.',
     .          'Using default Inertial Referecne Frame: IAU0A',0)
        if (buf5.ne.'     ') then 
            precmod = buf5 
            nutmod = buf5
        else
            if(frame.eq.'B1950') precmod = 'IAU68'
* MOD TAH 200304: Set default is buf5 is blank to IAU0A
            if(frame.eq.'J2000') precmod = 'IAU0A' 
        endif    
        close( 21 )

      endif
      
cxxxxxxxxxxxxxxxxxxxxx end of modification xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                          
c     open the nutation file if old-style
      if( .not.fcheck('nbody') ) then
        open( unit=20,file='nutabl.',status='old',iostat=ioerr)
        if( ioerr.ne.0 ) 
     .     call report_stat('FATAL','NGSTOT','orbits/ngstot','nutabl.',
     .    'Error opening nutation table: ',ioerr)
        read(20,'(a)',iostat=ioerr) nut_header
        if(ioerr.ne.0) call report_stat('FATAL','NGSTOT','orbits/ngstot'
     .   ,' ','Error reading first line of nutation file to get model'
     .   ,ioerr)
         close(20)
         if(mchkey(nut_header,'IAU20',80,5).gt.0) then
* MOD TAH 200303: Changed default IAU0A which is fixdrv default. (Override 
*        sestbl. if needed
* MOD TAH 200304: If nutabl. used, the IAU00 is correct nutation model.
           nutmod='IAU00' 
         else
           nutmod='IAU80'
         endif      
      else
* MOD TAH 200303: 
c       Old IAU00 will call MHB_2000 for nutations
c       New 10.71 IAU0A will call the SOFA routines.
        if ( precmod .eq. 'IAU76' ) then
* MOD TAH 200303: Changed default to IAU0A which is fixdrv default. (Override 
*        sestbl. if needed
* MOD TAH 200304: If still IAU76, then selected by usere so keep IAU00 nutation
*        series.
          nutmod='IAU00'
        else if (precmod .eq. 'IAU68' ) then
           nutmod='IAU80'
        else        
          nutmod = precmod
        endif
      endif
      if(debug) print*,'precmod, nutmod: ',precmod,' ',nutmod
      
c     set the non-gravitational force parameters to nominal values
      do j=1,numsat
         satics(7,j)= 1.D0
         satics(8,j)= 0.D0
         satics(9,j)= 0.D0
      enddo

c     write the earth fixed T-file header

      if( debug ) 
     .   print *,'calling THDRIT jde te jds ts jdf tf '
     .       ,               jde,te,jds,ts,jdf,tf 

      call thdrit( iut,jde,te,jds,ts,jdf,tf,delts,nepchs,nintrs
     .         , numsat,gnss_sel,itsat,satnam,satics,nics,spfile,icsnam
     .         , precmod,nutmod,gravmod,eframe,srpmod,eradmod,antradmod)
!     print*,'precmod, nutmod 1: ',precmod,' ',nutmod       

c Read the SP file and write the earth-fixed T-file data records

      call tdtrit( iungs,iut,spver,jds,ts,jdf,tf,delts,nepchs,nintrs
     .           , gnss_sel,numsp3sv,numsat,itsat
     .           , pvsigb,clksigb,otlmod )


c Now write out the inertial T-file

      close( iut )
      idir = 1     
      call report_stat('STATUS','NGSTOT','orbits/ngstot',tfile
     .                  ,'Writing inertial T-file ',0)
      call trot(tfilef,tfile,idir,frame)

c Now write out the G-file
      icall = 1
      prog = 'NGSTOT'
      call gmake( prog,tfile,gfile,itsat,icall )

      call report_stat('STATUS','NGSTOT','orbits/ngstot',' ',
     .'Normal end to NGSTOT ',0)

      stop
      end
