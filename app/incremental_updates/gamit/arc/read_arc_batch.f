      Subroutine read_arc_batch

c     Read the ARC batch file, storing the values in arc.h
c     Code moved from arc, get_sats, arcmrg, and lib/ghdred
c     All variables stored in common in arc.h
c     R. King 15 May 2012

      implicit none

      include '../includes/dimpar.h'
      include '../includes/global.h'
      include '../includes/arc.h'
                                   
      integer*4 nkpsat,ioerr,iy,im,id,ih,idoy,min,isvn,frqchn
     .        , svnstart(5),svnstop(5),linelength,i
                                    
      real*8 sec
                            
      character*1 yawbias 
      character*16 string     
      character*128 line 
      character*256 message 
 
      logical eos,firstsat/.true./,debug/.false./

c       Functions           
      integer*4 julday,nblen
      character*3 upperc       
      character*6 lowerc 
* MOD TAH 200328: Added string for read idbprt (optional on line).
      character*128 sarcout   ! Line with arcout name.  The optional
                              ! idbprt at end of line causes issues when
                              ! read. 
      character*16 sidbprt    ! String used to read idbprt from line 
                              ! in case it is not there.     
      logical kbit


c Read the list of satellites to be integrated

      do i=1,maxsat
        satnam(i) = ' '
      enddo
      nsats = 0
      eos = .false.
      do while (.not.eos ) 
        string = ' ' 
        read(5,'(a16)',iostat=ioerr) string
        if ( ioerr.ne. 0 ) 
     .    call report_stat('FATAL','ARC','read_arc_batch','string',
     .   'Error reading satellites from the ARC input batch file',ioerr)
        if( upperc(string(1:3)) .eq.'END' ) then
          eos=.true.
        else   
          nsats = nsats + 1
          satnam(nsats) = string  
        endif
      enddo
          
c Read the models and integration intervals
              
      read(5,'(a)',iostat=ioerr)  line 
      if( ioerr.ne.0 ) call report_stat('FATAL','ARC','read_arc_batch'
     . ,' ','Error reading batch file model line',ioerr)
cd      print *,'BATCH line:'
cd      print *,line 
      linelength = nblen(line)  
cd      print *,'linelength ',linelength 
      read(line(1:72)
     .     ,'(a5,1x,a5,1x,f6.1,1x,f8.4,2x,a4,2x,a20,a5,1x,a5,1x,a5)'
     .     ,iostat=ioerr)  gravmod,srpmod,delt,diint
     .     ,time_type,frame_name,precmod,eradmod,antradmod 
       if(ioerr.ne.0) call report_stat('FATAL','ARC','read_arc_batch'
     .    ,line,'Error reading 1st 72 columns of ARC  batch file',ioerr)

* MOD TAH 190622: Removed code below and replaced with more flexible version

* MOD TAH 171224: Made input more flexible.
* MOD RWK 170723: Restore the original (TAH mod didn't work for some reason)
C     if(linelength.gt.72) then
C       read(line(73:81),'(3i3)',iostat=ioerr) gravdeg,etidedeg,otidedeg
* TAH        read(line(73:),*,iostat=ioerr) gravdeg,etidedeg,otidedeg
cd        print *,'READ_ARC_BATCH gravdeg etidedeg otidedeg '
cd     .     ,gravdeg,etidedeg,otidedeg     
C       if(ioerr.ne.0) call report_stat('FATAL','ARC','read_arc_batch'
C    .    ,line,'Error reading grav/tide degree from ARC  batch file'
C    .    ,ioerr)
C       write(message,'("Degree options Grav ",i3,"   Etide ",I3
C    .        ,"   Otide ",I3)')  gravdeg,etidedeg,otidedeg
C       call report_stat('STATUS','ARC','read_arc_batch',' ',message,0)
C     else                                           
C       gravdeg = 12
C       etidedeg = 4
C       otidedeg = 12 
C       write(message,'(2a)') 'Old-style batch file, set degree of '
C    .     ,'gravity, etide, and otide harmonics to defaults (12 4 12)'
C       call report_stat('WARNING','ARC','read_arc_batch',' ',message,0)
C     endif
**RWK 181205: Allow reading of 'lbody' to control inclusion of Venus and Jupiter
C     if(linelength.gt.81) then
C       read(line(82:83),'(i2)',iostat=ioerr) lbody
C       if(ioerr.ne.0) call report_stat('FATAL','ARC','read_arc_batch'
C    .    ,line,'Error reading lbody from ARC  batch file',ioerr)
**      this next line to catch case of the '1' being read I2 as '10'
C       if(lbody.ne.0) then 
C         lbody = 1
C         call report_stat('STATUS','ARC','read_arc_batch',' '
C    .      ,'Venus and Jupiter used in the integration',0)
C       else
C         call report_stat('STATUS','ARC','read_arc_batch',' '
C    .      ,'Venus and Jupiter not used in the integration',0)
C       endif 
C     else
C       lbody = 0
C        call report_stat('STATUS','ARC','read_arc_batch',' '
C    .      ,'Venus and Jupiter not used in the integration',0)
C     endif 
* END REMOVED CODE.

      if(linelength.gt.72) then
        read(line(73:),*,iostat=ioerr) gravdeg,etidedeg,otidedeg,
     .        lbody
* MOD TAH 190622: Set lbody to zero if not passed.  
        if( ioerr.eq.-1 ) then  ! lbody not on line so set to zero
           lbody = 0
           ioerr = 0
        end if 
        if(ioerr.ne.0) call report_stat('FATAL','ARC','read_arc_batch'
     .     ,line,'Error reading grav/tide degree from ARC  batch file'
     .     ,ioerr)
      else                                           
        gravdeg = 12
        etidedeg = 4
        otidedeg = 12 
        lbody = 0
        write(message,'(2a)') 'Old-style batch file, set degree of '
     .     ,'gravity, etide, and otide harmonics to defaults (12 4 12)'
        call report_stat('WARNING','ARC','read_arc_batch',' ',message,0)
      endif
      if( lbody.eq.0 ) then
         call report_stat('STATUS','ARC','read_arc_batch',' '
     .      ,'Venus and Jupiter not used in the integration',0)
      endif 

      write(message,'("Degree options Grav ",i3,"   Etide ",I3
     .        ,"   Otide ",I3," LBODY ",I1)')  
     .            gravdeg,etidedeg,otidedeg, lbody

      call report_stat('STATUS','ARC','read_arc_batch',' ',message,0)
   
c Assign values where blank in old batch files

      if( time_type.eq.'    ' ) time_type = 'UTC '
      if( frame_name.eq.'                    ') 
     .     frame_name = 'INERTIAL     1950.0 '
c       determine the inertial frame
      if( index(frame_name,'5').ne.0)  then
        frame='B1950'
      elseif( index(frame_name,'2').ne.0)  then
        frame='J2000'
      else
        call report_stat('FATAL','ARC','arc',frame,
     .  'Unrecognised frame in ARC input batch file',0)
      endif        
c     if no precession model specified, assign IAU68 to B1950, IAU76 to J2000
      if( precmod.eq.'     '.and.frame.eq.'B1950') precmod = 'IAU68'
      if( precmod.eq.'     '.and.frame.eq.'J2000') precmod = 'IAU76'
c       check precession model is valid
      if( precmod .ne. 'IAU76' .and. precmod .ne.'IAU68' .and. 
     .    precmod .ne. 'IAU0A' .and. precmod .ne.'IAU0C' .and.
     .    precmod .ne. 'IAU06' .and. precmod .ne.'IAU6A') then
        call report_stat('FATAL','ARC','arc',precmod,
     .  'Unrecognised precession model in ARC input batch file',0)
      elseif((frame.eq.'J2000'.and.precmod.eq.'IAU68').or.
     .       (frame.eq.'B1950'.and.precmod.eq.'IAU76')) then
        write(message,'(a,a5,a,a5,a)') 'Frame ',frame,' and precession '
     .       ,precmod,' are incompatible'
        call report_stat('FATAL','ARC','read_arc_batch',' ',message,0)
      endif       
c** RWK 120512: Removed from the model line the hidden flags to turn off
c               solid-earth tides and lunar eclipses, now controlled by
c               the idbprt binary-code number in the print-file line
      if( eradmod.eq.'     ') eradmod = 'NONE '
      if( antradmod.eq.'     ') antradmod = 'NONE'


c Check consistency of integration and tabular intervales

      if((delt.ne.1350.0.and.delt.ne.900.0).or.(diint.ne.75.0.and.
     .    diint.ne.168.75)) then
        write(message,'(3a,f8.4,a,f6.1)') 'Integration error may occur'
     .     ,' since standard stepsize or tabular interval not selected:'
     .     ,' stepsize ',diint,' tabular interval ',delt
        call report_stat('WARNING','ARC','read_arc_batch',' ',message,0)
      endif
      if ( dmod(delt,diint).ne.0.0 ) then
        write(message,'(2a,f6.1,f8.4)') 'Tabular interval not integer '
     .    ,' divisible by the integration stepsize: ',delt,diint
        call report_stat('FATAL','ARC','read_arc_batch',' ',message,0)
      endif

 
c Read the name of the print file 
* MOD TAH 200328: Modified read to be more general in reading the 
*     (non-documented) idprt variable.  Read line free format and
*     check for EOF on line
*     Reomoved all lines with C comment below.
C     read(5,'(a)',iostat=ioerr) afname
C     if(ioerr.ne.0) call report_stat('FATAL','ARC','read_arc_batch'
C    .  ,' ','Error reading print-file name from ARC batch file',ioerr)
c       Although the arc print-file name is declared 16-characters,
c       in practice we assume that it is of the form 'arcout.DDD' or,
c       in any case, less than 14 characters so that we can hide a
c       3-digit number in the last three places to invoke debug.
C     if( afname(14:16).ne.'   ') then
C        read(afname(13:16),'(i4)',iostat=ioerr) idbprt
C        if( ioerr.ne.0 ) call report_stat('FATAL','ARC','arc',afname
C    .     ,'Non-blank char at end of afname must be integers for debug'
C    .     ,ioerr)                 
C        afname(13:16) = '   '                   
C     else 
C       idbprt = 0
C     endif
*     Get the arcout name and debug option and see what is on line
*     Line must be read in two parts of esle the missing idbprt value
*     will be left handing in the read if read from file and not present. 
      read(*,'(a)',iostat=ioerr) sarcout
      read(sarcout,*,iostat=ioerr) afname, sidbprt
      if( len_trim(sidbprt).eq.0 .or. ioerr.eq.-1 ) then
         idbprt = 0
      else
         read(sidbprt,*,iostat=ioerr) idbprt
         if( ioerr.ne.0 ) call report_stat('FATAL','ARC',
     .         'read_arc_batch', ' ',
     .        'Error IDBPRT from ARC batch file',ioerr)
      endif

c     warn if tides or lunar eclipses turned off
      if( kbit(idbprt,11)) then 
        write(message,'(a)') 'Solid_Earth tide effect turned off'
        call report_stat('WARNING','ARC','arc',' ',message,0)
      endif  
c      check hidden flag for turning off lunar eclipses
      if( kbit(idbprt,12))  then 
        write(message,'(a)') 'Lunar eclipses turned off'
        call report_stat('WARNING','ARC','arc',' ',message,0)
      endif     

 
c Read the name of the g-file 

      read(5,'(a)',iostat=ioerr) gfname
      if(ioerr.ne.0) call report_stat('FATAL','ARC','read_arc_batch'
     .  ,' ','Error reading g-file name from ARC batch file',ioerr)
                           
c Skip the blank line formerly used for the name of the x-file

      read(5,'(a)',iostat=ioerr) 

c Read the start and stop times for the integration

      read(5,'(i4,1x,i3,1x,2(i2,1x),f8.5)',iostat=ioerr) 
     .     iy,idoy,ih,min,sec 
      if(ioerr.ne.0) call report_stat('FATAL','ARC','read_arc_batch'
     .  ,' ','Error reading start time from ARC batch file',ioerr)
      call monday(idoy,im,id,iy)
      call check_y2k(iy)
      jdb= julday( im,id,iy )
      tb= ih*3600.d0 + min*60.d0 + sec
      read(5,'(i4,1x,i3,1x,2(i2,1x),f8.5)',iostat=ioerr) 
     .     iy,idoy,ih,min,sec
      if(ioerr.ne.0) call report_stat('FATAL','ARC','read_arc_batch'
     .  ,' ','Error reading stop time from ARC batch file',ioerr)
      call monday(idoy,im,id,iy)
      call check_y2k(iy)
      jdf= julday( im,id,iy )
      tf= ih*3600.d0 + min*60.d0 + sec

c Translate old-style PRN names to new-style antbody names
               
cd      print *,'READ_ARC_BATCH nsats satnam ',nsats,satnam(nsats)
      do i=1,nsats
        if( satnam(i)(1:3).eq.'PRN' ) then  
          if(firstsat ) then   
            call report_stat('WARNING','ARC','lib/read_gfile',gfname
     .  ,'Old-style satnam on ARC batch file, translating to GNSS style'
     .        ,0)
            firstsat = .false.
          endif
          gnss = 'G'
          read(satnam(i)(5:6),'(i2)') iprn 
cd        print *,'READ_GFILE read nsats, satnam iprn ',nsats
cd .             , satnam(nsats),iprn
* MOD TAH 190702: Added antpwr to snav_read call
          call svnav_read( -1,iy,idoy,ih,min,gnss,iprn,isvn,frqchn
     .                   , antbody,sbmass,yawbias,yawrate, antpwr
     .                   , svnstart,svnstop )
* MOD TAH 200606: Made I2.2 and flexiable a format to allow
*         possible use of all antbody length.  This seems to be only
*         for GPS (PRN name) and antbody(7:) should be fine (see cases
*         for writing antbody(x:) in ngstot.f 
          write(satnam(i),'(a1,i2.2,2x,i2,1x,a)')
     .        'G',iprn,isvn,trim(antbody(7:))
cd          print *,'After translation satnam ',satnam(i)
        endif
      enddo
    
c Read the flag indicating whether partials are to be integrated ('Y'/'N')

      read(5,'(a)') apar

                                                                         
c Read the name of the t-file

      read(5,'(a)',iostat=ioerr) tfname
      if(ioerr.ne.0) call report_stat('FATAL','ARC','read_arc_batch'
     .  ,' ','Error reading t-file name from ARC batch file',ioerr)

      if(debug) then                                  
        print *,'debug in read_arc_batch nsats ',nsats 
        write(*,'(20(a16,1x))') (satnam(i),i=1,nsats) 
        print *,'delt diint ',delt,diint
        print *,'models ',frame,precmod,nutmod,gravmod,srpmod
     .          ,eradmod,antradmod
        print *,'afname ',afname
        print *,'gfname ',gfname
        print *,'start ',jdb,tb
        print *,'end   ',jdf,tf
        print *,'partial flag ',apar
        print *,'tfname ',tfname
        print *,'debug stop in read_arc_batch '
        stop
      endif

      return
      end




                                 



























