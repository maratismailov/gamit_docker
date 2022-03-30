      Subroutine svnav_read( idir,iyr,iday,ihr,imin,gnss,isat
     .           , jsat,frqchn,antbody,sbmass,yawbias,yawrate, antpwr,
     .             start,stop ) 
c
c PURPOSE: Determine the GNSS SVN from PRN or vica versa, and the antenna/body
c          -type, frequency (for GLONASS), satellite mass, maximum yaw rate, 
c           and sun sensor bias status.
c
c PARAMETERS:                        
c         IN: idir    : idir = 1 SVN to PRN, idir = -1 PRN to SVN  I*4
c             iyr     : 4-digit year                               I*4
c             iday    : 3-digit day of year                        I*4
c             ihr     : 2-digit hr                                 I*4
c             imin    : 2-digit minute                             I*4
c             gnss    : 1-character GNSS code (G R E C J I)        C*1
c             isat    : input SVN or PRN number                    I*4
c
c        OUT: jsat    : output SVN or PRN number                   I*4
c             frqchn  : frquency channel (GLONASS only)            I*4
c             antbody : antenna/body-type (rcvant_tab/ANTEX std)   C*20
c             sbmass  : S/C mass in grams (if old-style will be 0) R*8
c             yawbias : Yaw bias                                   C*1
c             yawrate : maximum yaw rate of the input prn or svn   R*8 
c             antpwr  : Transmission power from igs_metadata       R*8
c
c CREATED  27 August 2014 by R. King from earlier version Mar 1995 by S McClusky
c
c There are three formats possible:
c   Old-style GAMIT table with header begining '  NSN/PRN', no longer supported
c   New-style GAMIT table with header 'svnav.dat  Version  2.0'
c   IGS Satellite Metadata SINEX with header beginning '%=SNX' 
c Since there are no new-style GAMIT tables beyond version 2.0, we'll use 
c the version variable for subsequent release of the SINEX file 

c The new-style GAMIT format is:

c svnav.dat  Version  2.0         
c SYS SVN  PRN CHAN ANT/BODY               MASS(G) YAW BIAS  YAW RATE  START           STOP             COMMENTS
c  G     1   4   0  BLOCK I                 453800.     U      0.1999  1978  53  0  0  2100   1  0  0                                           
c  R   701   6  99  GLONASS-M              1450000.     U      0.2500  2003 344  0  0  2010 118  0  0   #channel 99 indicates unknown  
c  R   731  22  -3  GLONASS-M              1450000.     U      0.2500  2010  60  0  0  2100   1  0  0
c  D     1  30   0  BEIDOU-2M              1000000.     Y      0.3000  2007 103  0  0  2100   1  0  0  #mass unknown
c  E   101  11   0  GALILEO-1               700000.     Y      0.4000  2011 294  0  0  2100   1  0  0
 
c Non-blank first column is a header or comment.   
                                                                                 
c NOTE 1: Since ARC integrations using cross a day boundary to provide interpolation,
c         there is a logical problem if the SVN/PRN assignment changes on the day
c         boundary. Since we'll want to keep the assignment that corresponds to 
c         the observation day (IC epoch), if a subsequent entry is found, we'll ignore 
c         the change.

* MOD TAH 190702: Added return of transmission power from igs_metadata.snx.  If
*     svnav.dat linked, 0.0 is returned since value is not known.

      implicit none

                  
      character*1  gnss,gnss_tmp,yawbias         
      character*10 svid 
      character*20 antbody 
      character*80 prog_name   
      character*128 record
      character*256 message

      integer*4 lun,idir,isat,jsat,frqchn
     .         ,iyr,iday,ihr,imin,time(5),start(5),stop(5),power
     .         ,iprn_tmp,isvn_tmp,len,rcpar,itemp,ioerr,i
c
      real*8 sbmass,yawrate,version
      real*8 antpwr    ! Tranmission power (W).
                    
      logical eof,found   

* MOD TAH 190627: Cleaned up code indenation and made use of svnav.dat
*     or igs_metadata_WWWW.snx more explicit.
      logical svnav_file  ! Set true is svnav.dat linked to svnav.dat
                          ! as opposed to igs_metadata_WWWW.snx file.  

c     function
      integer*8 itimdif

c Get calling program name for report_stat
      len = rcpar(0,prog_name)
                  
cd      print *,'SVNAV_READ idir iyr iday ihr imin isat '
cd     .      ,  idir,gnss,iyr,iday,ihr,imin,isat  

c Put the requested time into an array for comparing with file entires
      time(1) = iyr
      time(2) = iday
      time(3) = ihr
      time(4) = imin
      time(5) = 0                            
c     times read from file have no seconds
      start(5) = 0
      stop(5) = 0 

c Open the svnav.dat file
      lun = 69 
      open(unit=lun,file='svnav.dat',form='formatted',status='old'
     .    ,iostat=ioerr ) 
      If( ioerr.ne.0 )  call report_stat('FATAL',prog_name
     .   ,'lib/svnav_read','svnav.dat','File not found ',ioerr)
                                       
c Get the version number from the first head line
      read(lun,'(a)') record               
cd      print *,'1st record ',record 
      if( record(3:5).eq.'NSN') then
        call report_stat( 'FATAL',prog_name,'lib/svnav_read'
     .    ,'svnav.dat','Old version of file not supported',0)
      elseif( record(12:18).eq.'Version') then
        read( record(20:23),'(f4.0)') version
c       this must be 2.0 (Make explicit test)
* MOD TAH 190627: Make explicit this is svnav.dat
        if( version.ne.2.0 ) then
           call report_stat('FATAL',prog_name,'lib/snvav_read',
     .          'svnav.dat','Wrong version of snvav.dat ',
     .          int(version*10))
        else
          svnav_file = .true.
        endif 
      elseif( record(3:5).eq.'SNX') then
        read( record(7:10),'(f4.0)')  version 
* MOD TAH 190627: Make explicit this is not svnav.dat
        svnav_file = .false.
      else
        call report_stat( 'FATAL',prog_name, 'lib/svnav_read',
     .               'svnav.dat','Unknown version of file ',0)
      endif
          
c See if a GAMIT-style svnav.dat file
* MOD TAH 190627: Explicitly test for svnav.dat (old 
C     if( version.eq.2.0 ) then 
      if( svnav_file ) then 

c        Skip the column head line
         read(lun,'(a)') record
         if (ioerr .ne. 0) 
     .   call report_stat('FATAL',prog_name,'lib/snvav_read',
     .             'svnav.dat','Error reading svnav.dat file',ioerr)
                       

c        Loop over svnav.dat file to find the requested entry
         eof = .false.    
         found = .false.
         do while (.not.found .and. .not.eof )
           read(lun,'(a)',iostat=ioerr) record   
           if( ioerr.eq.-1 ) then
             eof = .true.
           elseif( ioerr.ne.0 ) then
              call report_stat('FATAL',prog_name,'lib/snvav_read'
     .         ,'svnav.dat','Error reading svnav.dat file',ioerr)
           elseif (record(1:1).ne.' ') then
c            comment 
             continue 
           else
             read(record,'(1x,a1,2x,i4,2x,i2,2x,i2,2x,a20,f11.0,5x
     .              , a1,f12.4,2(i6,i4,2i3))')  
     .         gnss_tmp,isvn_tmp,iprn_tmp,frqchn,antbody,sbmass
     .       , yawbias, yawrate
     .       , (start(i),i=1,4),(stop(i),i=1,4)

c            check if prn or svn is being requested
             if ( idir.lt.0 ) then  
c              SVN from PRN   
               if ( gnss_tmp.eq.gnss .and. iprn_tmp.eq.isat ) then
                 if( itimdif(start,time).le.0 .and. 
     .            itimdif(stop,time).ge.0 ) then           
                   found = .true.
cd                   print *,'found ',found 
                   jsat = isvn_tmp
                 endif
               endif
             elseif (idir.gt.0 ) then 
c              PRN from SVN
c              Beidou does not have official SVNS but svnav.dat has 
c              MGEX-assigned values
c              if( gnss.eq.'C' )  call report_stat('FATAL',prog_name
c    .   ,'lib/svnav_read','svnav.dat','Beidou has no SVN defined',0)
               if ( gnss_tmp.eq.gnss .and. isvn_tmp.eq.isat ) then
                 if( itimdif(start,time).ge.0 .and. 
     .            itimdif(stop,time).le.0 ) then
                   found = .true.
                   jsat = iprn_tmp
                 endif
               endif
             else
               call report_stat('FATAL',prog_name,'lib/svnav_read',' '
     .                      ,'idir=0, something wrong',0)
             endif
           endif
         enddo
                                   
c Stop if SV not found in svnav.dat file
         if( .not.found ) then 
           if( idir.lt.0 ) then 
             write(message,'(a,a1,i2,a)') 'PRN ',gnss,isat
     .           ,' not found on svnav.dat file'
           elseif( idir.gt.0 ) then
             write(message,'(a,a1,i3,a)') 'SVN  ',gnss,isat
     .           ,' not found on svnav.dat file'
           endif             
cd           print *,'prog_name(1:7)',prog_name(1:7)
           if( prog_name(1:7).eq.'dcbtab2' .or. 
     .         prog_name(1:3).eq.'arc') then
               call report_stat('WARNING',prog_name,'lib/svnav_read',
     .                     ' ', message,0)
           else
               call report_stat('FATAL',prog_name,'lib/svnav_read',
     .                    ' ', message,0)
           endif
         endif  

* MOD TAH 1900702: Set antpwr to 0 since not in svnav.dat
         antpwr = 0.d0

c If not a GAMIT-style file, read a SINEX files
  
      ELSE 
        call read_svsinex( lun,idir,iyr,iday,ihr,imin,gnss,isat
     .                   , jsat,frqchn,antbody,svid,sbmass
     .                   , yawbias,yawrate,power ) 
*       Convert power I4 to R8 antpwr
        antpwr = power
      endif
      close(lun)    

c     DEBUG
c     stop

      end




