      Subroutine read_svsinex( lun,idir,iyr,iday,ihr,imin,gnss,isat
     .                       , jsat,frqchn,antbody,svid,sbmass
     .                       , yawbias,yawrate,power)

c PURPOSE: Read an IGS Satellite metadata SINEX file, returning the values
c          needed for the SV at the input epoch

c PARAMETERS:
c       IN:   lun     : logical unit number for read               I*4
c             idir    : iidir = 1 SVN to PRN, idir = -1 PRN to SVN I*4
c             iyr     : 4-digit year                               I*4
c             iday    : 3-digit day of year                        I*4
c             ihr     : 2-digit hr                                 I*4
c             imin    : 2-digit minute                             I*4
c             gnss    : 1-character GNSS code (G R E C J I)        C*1
c             isat    : input SVN or PRN number                    I*4
c
c        OUT: jsat    : output SFN or PRN number                   I*4 
c             frqchn  : frquency channel (GLONASS only)            I*4
c             antbody : antenna/body-type (rcvant_tab/ANTEX std)   C*20
c             svid    : SINEX name for satellite                   C*10 - not yet used
c             sbmass  : S/C mass in grams (if old-style will be 0) R*8
c             yawbias : Yaw bias                                   C*1
c             yawrate : maximum yaw rate of the input prn or svn   R*8 
c             power   : transmitted power in watts                 I*4

c CREATED 14 November 2017 by R. King

      implicit none

                  
      character*1  gnss,gnss_tmp,yawbias         
      character*3  prn,          ! PRN string read from sinex file
     .             prn_in        ! PRN string for satellite we are looking for
                         ! or the correct one when svn found (depends on idir)                        
      character*4  svn,          ! SVN string read from sinex file
     .             svn_in        ! SVN string for satellite we are looking for
                         ! or the correct one when prn found.

      character*10 cospar_id   ! ID not used in GAMIT
      character*6  SatCat      ! Catalog number; not used in GAMIT.
      character*20 antbody_in  ! Body type read from metadata snx

      integer*4 idir     ! Direction of PRN->SVN or visa-vers
      integer*4 yr1,  yr2   ! Year read from sinex file
      integer*4 doy1, doy2  ! DOY read from sinex file
      real*8    sod1, sod2  ! Seconds of day read from sinex file
      real*8    secs        ! Needed for call to ds2hms: Added TAH 190701

      integer*8 itimdif     ! GAMIT function to return time difference in seconds?
                                                    
      character*10 svid 
      character*20 antbody 
      character*80 prog_name   
      character*128 record
      character*256 message

      integer*4 lun,isat,jsat,frqchn,power,iyr,iday,ihr,imin,time(5)
     .        , start(5),stop(5),len,rcpar,ioerr,i
c        
      real*4 version
      real*8 sbmass,yawrate

      logical found  ! Used to indicate that blocks and SVNs have been found.
      logical debug / .false. /   ! Turn on ouput.
                    
c Get calling program name for report_stat
      len = rcpar(0,prog_name)
                                                                    
c Put the requested time into an array for comparing with file entires
      time(1) = iyr
      time(2) = iday
      time(3) = ihr
      time(4) = imin
      time(5) = 0                            
c     times read from file have no seconds
      start(5) = 0
      stop(5) = 0  
      frqchn = 0   ! Initialize value for when GLONASS no used.

c Get the SVN/PRN correspondences  ( SATELLITE/PRN  block ) 
* Start read the sinex file looking for the start of SATELLITE/PRN  block
*     Once found, look for our specific prn or svn. 
      rewind(lun)
      found = .false.
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record
        if(ioerr.ne.0 ) then 
            call report_stat('FATAL',prog_name,'lib/read_svsinex', 
     .         'svnav.dat',
     .         'Failed to find SATELLITE/PRN block',ioerr)
        else
            if( record(12:14).eq.'PRN' ) then
               found = .true.
            endif 
        endif 
      enddo  

* MOD TAH 190627: Don't check direction yet we will do that below as
*     the PRN or SVN is found.  We need to find the SVN number because
*     most the file is based on SVN not (non-unique) PRN. 
C     if( idir.eq.-1 ) then
c       PRN to SVN
*     Depending on direction either assign PRN straing (G04) or the 
*     svn string (G074) so we can search the records.
*     Assuming here isat in input PRN or SVN and jsat is output PRN or SCN
      if( idir.eq.-1 ) then  ! PRN to SVN
         write(prn_in,'(a1,I2.2)') gnss, isat
      else                   ! SVN to PRN 
         write(svn_in,'(a1,I3.3)') gnss, isat
      endif

*     Now read the file to see if can find the entry 
      found = .false.
      do while ( .not.found ) 
        read(lun,'(a)',iostat=ioerr) record
        if( record(1:4).eq.'-SAT'.or.ioerr.ne.0 ) then 
           call report_stat('FATAL',prog_name, 
     .           'lib/read_svsinex','svnav.dat',
     .           'Failed to find SV in PRN block of SINEX file',ioerr)    
        else
*          Only try to decode record if it is not a comment
           if( record(1:1).eq.' ' ) then
              read(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,a3)'
     .            ,iostat=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,prn 

              If( ioerr.ne.0 )  call report_stat('FATAL',prog_name,
     .             'lib/read_svsinex','svnav.dat',
     .              'Error decoding PRN block of SINEX file',ioerr)
              if( yr2.eq.0000 ) then 
                  yr2 = 2100
                  doy2 = 1
              endif
              start(1) = yr1
              start(2) = doy1
              call ds2hms(yr1,doy1,sod1,start(3),start(4),secs) 
              start(5) = nint(secs)
              stop(1) = yr2
              stop(2) = doy2
              call ds2hms(yr2,doy2,sod2,stop(3),stop(4),secs) 
              stop(5) = nint(secs)
*             See if have a match to the prn_in or svn_in we are
*             looking for
              if( idir.eq.-1 ) then     ! PRN_in passed, see if match
c                PRN to SVN: Compare prn string from file with prn_in
                 if( prn.eq.prn_in .and.
     .              itimdif(start,time).le.0 .and.
     .              itimdif(stop,time).ge.0 ) then
                   found = .true. 
                   read(svn(2:4),'(i3)') jsat
                   svn_in = svn   ! Save so that we use to find other 
                                  ! entries such as power, mass and frequencies.
                 endif  
              elseif( idir.eq.1 ) then
c               SVN to PRN: Compare svn string from file with svn_in
                if( svn.eq.svn_in .and.
     .              itimdif(start,time).le.0 .and.
     .              itimdif(stop,time).ge.0 ) then
                   found = .true. 
                   read(prn(2:3),'(i2)') jsat
                   prn_in = prn    ! Save in case we need.
                endif     
              else 
                If( ioerr.ne.0 )  call report_stat('FATAL',prog_name, 
     .              'lib/read_svsinex','svnav.dat', 
     .              'idir not -1 or 1 in call to read_svsinex',ioerr)
              endif
           end if
        endif
      enddo

c Get the frequency channel for Glonass ( SATELLITE/FREQUENCY_CHANNEL block)
*     This only needed if gnss == R for Glonass.
      if( gnss(1:1).eq.'R' ) then 
         rewind(lun)
         found = .false.
         do while(.not.found) 
           read(lun,'(a)',iostat=ioerr) record              
           if( ioerr.ne.0 ) then
              call report_stat('FATAL',prog_name,'lib/read_svsinex', 
     .           'svnav.dat',
     .           'Failed to find FREQUENCY_CHANNEL SINEX block',ioerr)
           elseif( record(12:15).eq.'FREQ' ) then
             found = .true.
           endif          
         enddo
         found = .false.                                 
         do while(.not.found) 
            read(lun,'(a)',iostat=ioerr) record
            if( record(1:4).eq.'-SAT'.or.ioerr.ne.0 ) then 
                write(message,220) svn_in, prn_in
 220            format('Failed to find SV ',a,' PRN ',a,1x,
     .                 'in FREQUENCY_CHANNEL SINEX block')
                call report_stat('FATAL',prog_name, 
     .             'lib/read_svsinex','svnav.dat', message, ioerr)
            else             
*           Only try to decode record if it is not a comment
               if( record(1:1).eq.' ' ) then
                 read(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,i4)',
     .              iostat=ioerr) svn,yr1,doy1,sod1,
     .                                yr2,doy2,sod2, frqchn
                 if(ioerr.ne.0) call report_stat('FATAL',prog_name,
     .               'lib/read_svsinex','svnav.dat', 
     .               'Reading GLONASS frequency',ioerr)                   
                 if( yr2.eq.0000 ) then 
                    yr2 = 2100
                    doy2 = 1
                 endif 
                 start(1) = yr1
                 start(2) = doy1
                 call ds2hms(yr1,doy1,sod1,start(3),start(4), secs) 
                 start(5) = nint(secs)
                 stop(1) = yr2
                 stop(2) = doy2
                 call ds2hms(yr2,doy2,sod2,stop(3),stop(4), secs) 
                 stop(5) = nint(secs)
                 if( svn.eq.svn_in .and.
     .               itimdif(start,time).le.0 .and.
     .               itimdif(stop,time).ge.0 ) then
                     found = .true. 
                  endif  
               endif
            endif
         enddo
      endif         ! gnss = R

c Get the SATELLITE/IDENTIFIER block
c     +SATELLITE/IDENTIFIER                                                           
c     *                                                                               
c     *SVN_ COSPAR ID SatCat Block__________ Comment__________________________________
c     *                                                                               
c      G001 1978-020A  10684 GPS-I           Launched 1978-02-22; NAVSTAR 1
c      G002 1978-047A  10893 GPS-I           Launched 1978-05-13; NAVSTAR 2
c      G003 1978-093A  11054 GPS-I           Launched 1978-10-07; NAVSTAR 3

      rewind(lun)
      found = .false.
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record   
        if(ioerr.ne.0 ) then
           call report_stat('FATAL',prog_name,'lib/read_svsinex'
     .       ,'svnav.dat'
     .       ,'Failed to find IDENTIFIER SINEX block',ioerr)
        elseif( record(12:21).eq.'IDENTIFIER' ) then
          found = .true.
        endif          
      enddo
      found = .false.                                 
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record
        if( record(1:4).eq.'-SAT'.or.ioerr.ne.0 ) then 
           call report_stat('FATAL',prog_name, 
     .        'lib/read_svsinex','svnav.dat',
     .        'Failed to find SV in IDENTIFIER SINEX block',ioerr)
        else             
*         Only try to decode record if it is not a comment
          if( record(1:1).eq.' ' ) then
             read(record,'(1x,a4,1x,a9,1x,A6,1x,a15)'
     .           ,iostat=ioerr) svn, cospar_id,SatCat, antbody_in
             if(ioerr.ne.0) call report_stat('FATAL',prog_name,
     .                'lib/read_svsinex','svnav.dat', 
     .                'Reading satellite id ',ioerr) 
*            Only one line per satellite (unique ID)  
             if( svn.eq.svn_in ) then
                 found = .true. 
*                Replace names to make consistent with GAMIT codes
*                Need to do GPS replacement in two steps due to 
*                trailing space
                 call sub_char( antbody_in,'GPS-','BLOCK_')
                 call sub_char( antbody_in,'_',' ')
*                Due to sub_char allowing multiple replacements
*                we can't directly replace'GLO' with 'GLONASS'
*                (sub_char knows and won't replace, so we need
*                to do in two steps
                 call sub_char( antbody_in,'GLO','XXX')
                 call sub_char( antbody_in,'XXX','GLONASS')     
                 call sub_char( antbody_in,'GAL','XXX')
                 call sub_char( antbody_in,'XXX','GALILEO')
                 call sub_char( antbody_in,'BDS','BEIDOU')
                 antbody = antbody_in
             endif  
          endif
        endif 
      enddo
             

c Get the satellite mass ( SATELLITE/MASS block )

      rewind(lun)
      found = .false.
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record   
        if(ioerr.ne.0 ) then
           call report_stat('FATAL',prog_name,'lib/read_svsinex'
     .       ,'svnav.dat'
     .       ,'Failed to find MASS SINEX block',ioerr)
        elseif( record(12:15).eq.'MASS' ) then
          found = .true.
        endif          
      enddo
      found = .false.                                 
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record
        if( record(1:4).eq.'-SAT'.or.ioerr.ne.0 ) then 
           call report_stat('FATAL',prog_name, 
     .        'lib/read_svsinex','svnav.dat',
     .        'Failed to find SV in MASS SINEX block',ioerr)
        else             
*         Only try to decode record if it is not a comment
          if( record(1:1).eq.' ' ) then
             read(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,f9.3)'
     .           ,iostat=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,sbmass 
             if(ioerr.ne.0) call report_stat('FATAL',prog_name,
     .                'lib/read_svsinex','svnav.dat', 
     .                'Reading satellite mass ',ioerr)                   
             if( yr2.eq.0000 ) then 
                 yr2 = 2100
                 doy2 = 1
             endif
             start(1) = yr1
             start(2) = doy1
             call ds2hms(yr1,doy1,sod1,start(3),start(4),secs) 
             start(5) = nint(secs)
             stop(1) = yr2
             stop(2) = doy2
             call ds2hms(yr2,doy2,sod2,stop(3),stop(4), secs) 
             stop(5) = nint(secs) 
             if( svn.eq.svn_in .and.
     .           itimdif(start,time).le.0 .and.
     .           itimdif(stop,time).ge.0 ) then
                 sbmass = sbmass*1.d3   ! Convert to grams from kg 
                                        ! (grams expected by GAMIT) 
                 found = .true. 
             endif  
          endif
        endif 
      enddo

c Get the yaw bias and rate.
* MOD TAH 190701: Implemented with new SATELLITE/YAW_BIAS_RATE block added
*     but TAH to igs_metadata.snx (using script sh_svnav_yaw_to_igsmeta)

      rewind(lun)
      found = .false.
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record   
        if( record(1:4).eq.'-SAT'.or.ioerr.ne.0 ) then 
           if(ioerr.ne.0 ) then
              call report_stat('FATAL',prog_name,'lib/read_svsinex'
     .          ,'svnav.dat'
     .          ,'Failed to find YAW_BIAS_RATE SINEX block',ioerr)
           endif 
        elseif( record(12:19).eq.'YAW_BIAS' ) then
          found = .true.
        endif          
      enddo
      found = .false.                                 
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record
        if( record(1:4).eq.'-SAT'.or.ioerr.ne.0 ) then 
          call report_stat('FATAL',prog_name, 
     .        'lib/read_svsinex','svnav.dat',
     .        'Failed to find SV in YAW block',ioerr)
        else             
*         Only try to decode record if it is not a comment
          if( record(1:1).eq.' ' ) then
             read(record,'(1x,a4,2(1x,i4,1x,i3,1x,F5.0),4x,a1,1x,f8.4)'
     .           ,iostat=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,
     .                          yawbias, yawrate
             if(ioerr.ne.0) call report_stat('FATAL',prog_name,
     .                'lib/read_svsinex','svnav.dat', 
     .                'Reading yaw rate ',ioerr)                   
             if( yr2.eq.0000 ) then 
                 yr2 = 2100
                 doy2 = 1
             endif
             start(1) = yr1
             start(2) = doy1
             call ds2hms(yr1,doy1,sod1,start(3),start(4), secs) 
             start(5) = nint(secs)
             stop(1) = yr2
             stop(2) = doy2
             call ds2hms(yr2,doy2,sod2,stop(3),stop(4), secs) 
             stop(5) = nint(secs) 
             if( svn.eq.svn_in .and.
     .           itimdif(start,time).le.0 .and.
     .           itimdif(stop,time).ge.0 ) then
                  found = .true. 
             endif  
           endif
         end if
      enddo


c Get the transmitter power ( SATELLITE/TX_POWER block)

      rewind(lun)
      found = .false.
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record 
        if(ioerr.ne.0 ) then
           call report_stat('FATAL',prog_name,'lib/read_svsinex',
     .         'svnav.dat','Failed to find TX_POWER SINEX block',ioerr)
        elseif( record(12:15).eq.'TX_P' ) then
          found = .true.
        endif 
      enddo   
      found = .false.                                 
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record 
        if( record(1:4).eq.'-SAT'.or.ioerr.ne.0 ) then 
            call report_stat('WARNIUG',prog_name,
     .         'lib/read_svsinex','svnav.dat',
     .         'Failed to find ' // svn_in // ' in TX_POWER block',
     .         ioerr)
           power = 185   ! Set nominal power
           found = .true.
        else             
*          Only try to decode record if it is not a comment
           if( record(1:1).eq.' ' ) then
              read(record,'(1x,a4,2(1x,i4,1x,i3,1x,F5.0),1x,i4)'
     .           ,iostat=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,power  
              if( yr2.eq.0000 ) then 
                  yr2 = 2100
                  doy2 = 1
              endif
              if(ioerr.ne.0) then
                  write(*,'(a,1x,a)') 'TX_POWER block: ',
     .               trim(record)
                  call report_stat('FATAL',prog_name,
     .              'lib/read_svsinex','svnav.dat',
     .               'Failed to Reading in TX_POWER block',ioerr)
                  power = 0.d0
              endif

              start(1) = yr1
              start(2) = doy1
              call ds2hms(yr1,doy1,sod1,start(3),start(4), secs) 
              start(5) = nint(secs)
              stop(1) = yr2
              stop(2) = doy2
              call ds2hms(yr2,doy2,sod2,stop(3),stop(4), secs) 
              stop(5) = nint(secs)
              if( svn.eq.svn_in .and.
     .             itimdif(start,time).le.0 .and.
     .             itimdif(stop,time).ge.0 ) then
                   found = .true. 
              endif  
           endif
         endif 
      enddo

      if( debug ) then
         write(*,310) svn_in, prn_in, antbody, start, stop, frqchn,
     .                sbmass,  power, yawbias, yawrate 
 310     format('META ',a4,1x,a3,1x,a20,1x,
     .          ' Dates ',2(I4,1x,I3,3(1x,I2.2),1x),
     .          ' Meta ',i3,1x,F9.1,' g',I6,' W ',1x,a1,1x,F8.4)
      endif

      return
      end

  



      





            



      
       
               


