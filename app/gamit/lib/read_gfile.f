Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995/2012.  All rights reserved.

      Subroutine read_gfile( gfname,iug,jde,te,frame,precmod,nutmod
     .                     , gravmod,srpmod,eradmod,antradmod,time_type 
     .                     , nsat,nics,icsnam,satnam,satics )

c       Read all values from the g-file (based on old /lib/ghdred.f)
c       R. King May 2012

c Input                        
c  gfname  c*16 : g-file name
c  iug     i*4  : logical unit for g-file (opened by calling program)

c Output
c   jde                    i*4   PEP Julian Day of IC epoch
c   te                     r*8   Seconds-of-day of IC epoch 
c   frame                  c*5   Inertial frame  (B1950 or J2000)
c   precmod                c*5   Precession model  (IAU69 or IAU76)
c   nutmod                 c*5   Nutation model for (IAU80 or IAU00)
c   gravmod                c*5   Gravity model   
c   srpmod                 c*5   Model for direct solar radiation pressure
c   eradmod                c*5   Model for Earth radiation pressure
c   antradmod              c*5   Model for antenna thrust 
c   nsat                   i*4   Number of satellites
c   nics                   i*4   Number of parameters (6 ICs + radiation pressure)
c   icsnam(maxorb)         c*4  Names of parameters
c   satnam(maxsat)         c*16  SV names (currently use only 6 characters, e.g. PRN 13)
c   satics(maxorb,maxsat)  r*8  Values of Initial conditions and parameter values
c                                           
      implicit none

      include '../includes/dimpar.h'
                        
      character*1 gnss,yawbias 
      character*4 icsnam(maxorb),srpnam(9),time_type
      character*5  precmod,frame,srpmod,nutmod,gravmod,eradmod,antradmod
      character*16 gfname,satnam(maxsat)
      character*20 antbody 
      character*81 buf81
      character*80 prog_name
      character*256 message

      integer*4 iug,nsat,nics,id,ih,im,idoy,min,iy,jde,iprn,isvn,frqchn
     .        , svnstart(5),svnstop(5),len,ioerr,i,j

      real*8 sec,te,utcoff,satics(maxorb,maxsat),sbmass,yawrate
      real*8 antpwr  ! Transmission power (W)

      logical eoh,eos,firstcall/.true./

c       functions
      integer*4 julday,rcpar
      real*8 taiutc


c  Get calling program name for report_stat

      len =  rcpar(0,prog_name)
           
c  Open the g-file
          
      open (unit=iug,file=gfname,status='old',iostat=ioerr)
      if (ioerr .ne. 0) then
        call report_stat('FATAL','ARC','filopn',gfname,
     .  'Error could not open G-file',ioerr)
      endif

c  Read the first line of the G-file for IC epoch, time-type, and models

      buf81 = ' '
      read(iug,'(a81)',iostat=ioerr) buf81 
      if( buf81(3:3).eq.' ') then                           
c       old format, 2-digit year
        read(buf81,'(i2,1x,i3,1x,i2,1x,i2,f3.0,20x,a4,7(1x,a5))'
     .        ,iostat=ioerr) iy,idoy,ih,min,sec,time_type
     .                     , frame,precmod,srpmod,nutmod,gravmod
     .                     , eradmod,antradmod
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .    ,'lib/read_gfile','Error reading pre-1995 g-file',gfname,0)
        call fix_y2k(iy)
      else
c       new format, 4-digit year
        read(buf81,'(i4,1x,i3,1x,i2,1x,i2,f3.0,18x,a4,7(1x,a5))'
     .        ,iostat=ioerr) iy,idoy,ih,min,sec,time_type
     .                    , frame,precmod,srpmod,nutmod,gravmod
     .                    , eradmod,antradmod    
        if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .    ,'lib/read_gfile','Error reading 1st line of g-file',gfname,0)
      endif
    

c   Read the second line to get the parameter number and names
                        
      read(iug,'(a81)',iostat=ioerr) buf81     
      read(buf81(1:2),'(i2)',iostat=ioerr) nics
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .    ,'lib/read_gfile','Error reading 2nd line of g-file',gfname,0)
      do i = 1,maxorb
        icsnam(i) = ' '
      enddo
      if( nics.gt.maxorb ) then
       write(message,'(a,i2,a,i2)') 'Parameter number on G-file ('
     .         ,nics,') exceeds MAXORB (',maxorb,')'
       call report_stat('FATAL',prog_name,'lib/read_gfile',gfname
     .         ,message,0)
      endif           
      if( nics.gt.19 ) call report_stat('FATAL',prog_name
     .    ,'lib/read_gfile'
     .    ,'No more than 15 parameters allowed on g-file',gfname,0)
      read(buf81(3:77),'(19(1x,a4))',iostat=ioerr) 
     .    (icsnam(i),i=1,nics)  
      if( ioerr.ne.0 ) call report_stat('FATAL',prog_name,
     .    'lib/read_gfile','Error reading parameter names from g-file'
     .    ,gfname,0)
           

c  Set the default model names (mostly for pre-1995 g-files) if blank
   
      if( time_type.eq.'     ')   time_type = 'UTC '
      if( frame.eq.'     ' )      frame = 'B1950'
      if( precmod.eq.'     ')     precmod = 'IAU68'
      if( nutmod.eq.'     ')      nutmod =  'IAU80'
      if( gravmod.eq.'     ')     gravmod = 'IGS92'
      if( eradmod.eq.'     ')     eradmod = 'NONE '
      if( antradmod.eq.'     ')   antradmod = 'NONE'
c     No longer allow a blank (pre-1996) for the radiation pressure model since
c     the spherical model is no longer supported.
c       1  ECOM1 (formerly BERNE) constant + once-per-rev direct, Y, B
c       2  ECOM2 same as ECOM1 plus direct 2-per-rev and 4-per-rev
c       7  UCLR1  Univ College London harmonic model (same parameters as ECOM1)
c       8  UCLR2  Univ College London grid model (same parameters as ECOM1)
      if( srpmod.ne.'ECOM1'.and.srpmod.ne.'BERNE'.and.srpmod.ne.'ECOM2'
     .    .and.srpmod.ne.'ECOMC'.and.
     .         srpmod.ne.'UCLR1'.and.srpmod.ne.'UCLR2' ) then
        write(message,'(a,a5,a)') 'Radiation-pressure model (',srpmod
     .    ,') not supported'
        call report_stat('FATAL',prog_name,'lib/read_gfile',' '
     .                  ,message,0)
      endif

c  Compute JD, seconds of day from the g-file times
          
      call monday(idoy,im,id,iy)
      jde = julday( im,id,iy)
      te = ih*3600.d0 + min*60.d0 + sec
c     UTC no longer allowed by FIXDRV
      if( time_type.eq.'UTC') then
        utcoff = taiutc(jde) - 19.d0
        call timinc(jde,te,utcoff)  
        write(message,'(a)') 'G-file epoch is UTC, changed to GPST' 
        call report_stat('WARNING',prog_name,'lib/read_gfile',gfname
     .                  , message,0)
      endif


c  Read through comments to the end of the header
        
      eoh = .false.
      do while( .not.eoh ) 
        read(iug,'(a)',iostat=ioerr) buf81
        if( ioerr.ne.0 )  call report_stat('FATAL',prog_name
     .      ,'lib/read_gfile',gfname,'Error finding END on g-file',0)
        if( buf81(1:3).eq.'END' ) eoh = .true.
      enddo

c  Read the ICs for each satellite
                           
      eos = .false. 
      nsat = 0 
      do while( .not.eos ) 
        read(iug,'(a)',iostat=ioerr) buf81
        if( buf81(1:3).eq.'END'.or.ioerr.eq.-1 ) then
          eos = .true.                    
        elseif( ioerr.ne.0 ) then
          if( ioerr.ne.0 )  call report_stat('FATAL',prog_name
     .      ,'lib/read_gfile',gfname,'Error finding END on g-file',0)
        else 
          nsat = nsat + 1
          satnam(nsat) = buf81(1:16) 
c         if pre-gnss names, translate them here
          if( satnam(nsat)(1:3).eq.'PRN' ) then
            if( firstcall ) then 
              call report_stat('WARNING',prog_name,'lib/read_gfile'
     .          ,gfname
     .          ,'Old-style satnam on g-file, translating to GNSS style'
     .         ,0)
              firstcall = .false.
            endif
            gnss = 'G'
             read(satnam(nsat)(5:6),'(i2)') iprn 
cd             print *,'READ_GFILE read nsat, satnam iprn ',nsat
cd     .           , satnam(nsat),iprn
* MOD TAH 190702: Added antpwr to snav_read call
             call svnav_read( -1,iy,idoy,ih,min,gnss,iprn,isvn,frqchn
     .                      , antbody,sbmass,yawbias,yawrate, antpwr
     .                      , svnstart,svnstop )
            write(satnam(nsat),'(a1,i2,2x,i2,1x,a8)')
     .          'G',iprn,isvn,antbody(7:14)
cd             print *,'nsat satnam ',nsat,satnam(nsat)
          endif    
          do i=1,nics
            read(iug,'(d20.0)',iostat=ioerr) satics(i,nsat)  
            if( ioerr.ne.0 ) then
              write(message,'(a,i2,a,i2)') 
     .           'Error reading IC ',i,' for satellite ',nsat
              call report_stat('FATAL',prog_name,'lib/read_gfile'
     .              ,gfname,message,0) 
            endif
          enddo  
cd          print *,'READ_GFILE iprn isvn antbody nics satics '
cd     .           , iprn,isvn,antbody,nics,(satics(i,nsat),i=1,nics)
        endif
      enddo

cDEBUUG  
cd      print *,'DEBUG in READ_GFLE:'
cd      print *,'gfname,iug ',gfname,iug
cd      print *,'times ',iy,idoy,ih,min,sec,time_type,jde,te
cd      print *,'models ',frame,precmod,nutmod,gravmod,srpmod
cd     .          ,eradmod,antradmod
cd      write(*,'(a,20(1x,a4))') 'icsnam ',(icsnam(i),i=1,nics)
cd      print *,'nsat nics ',nsat,nics
cd      do j=1,nsat
cd        write(*,'(a,20d20.13)') satnam(j),(satics(i,j),i=1,nics)
cd      enddo
cd      print *,'debug stop in read_gfile '
cd      stop

      return
      end
