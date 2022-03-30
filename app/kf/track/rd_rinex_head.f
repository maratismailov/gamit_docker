CTITLE READ_RINEX_HEAD
 
      subroutine read_rinex_head(obs_file, debug, rxv  )

      implicit none
 
*     Routine to read the header from the data files.
*     Returns will be:
*     Approximate site position
*     station name
 
      include '../includes/xfile_def.h'

* obs_file -- Full name of rinex file to open and read header from

      character*(*) obs_file
      
* Local variables

*    swver  - Rinex software version
      character*20 swver
            
*   ierr    - IOSTAT error
*  trimlen  - Length of string
*   indx    - Position of substring in a string
 
      integer*4 ierr, jerr, trimlen, indx, i, k
      integer*4 nline, lent  ! additional lines need to read
                             ! obervable types and last index

      real*4 rxv             ! Rinex version 

*   eoh     - End of header flags
*   rinex   - True if rinex file
       logical eoh, rinex

       integer*4 debug  ! debug start value to see if debug output
 
*   line    - Line read from file
 
 
      character*256 line
      character*1 cr
      character*2 datatypes(xf_maxdat)

 
****  Open the data file
      cr = char(13)
      call openfile(obs_lu, obs_file, 'old',' ','ASCII', ierr)   
      call report_error('IOSTAT',ierr,'open',obs_file, 1,
     .        'read_rinex_file')

*     Start looping over the header
      eoh = .false.
      rinex = .false.
      xf_ntext = 0
      jerr = 0
      do while ( .not.eoh )
          read(obs_lu,'(a)', iostat=ierr) line 
          call sub_char(line,cr,' ')
          call casefold(line)
          call report_error('IOSTAT',ierr,'read', obs_file,1,
     .                'End of hedader not found')
          if( trimlen(line).lt.10) eoh = .true.
          if( index(line,'END OF HEAD').gt.0 ) eoh = .true.

*         Check rinex verions
          if( index(line,'RINEX VERSION').gt.0 ) then
              read(line,*) rxv
              write(*,120) trim(obs_file), rxv
 120          format('Data File ',a,' Rinex Version ',F5.2)
              if( rxv.ge.3.0 ) then
                 call report_stat('FATAL','SVPOS',
     .                'RINEX 3 not suppported yet',obs_file,'',0)
              end if
          end if

*         See if we find Rinex file
          indx = index(line,'OBSERVATION')
          if( indx.gt.12 ) then
               rinex = .true.
C              read(line,*,iostat=jerr) xf_isessn 
C              write(*,*) jerr, xf_isessn
               xf_isessn = 1
               call report_error('IOSTAT',jerr,'decod',line,0,'RXHEAD')
          end if
          
*         For version 2
          indx = index(line,'O                   G')
          if( indx.gt.12 ) then
               rinex = .true.
               read(line, *,iostat=jerr) xf_isessn
               call report_error('IOSTAT',jerr,'decod',line,0,'RXHEAD')
          end if
          
*         See if marker name
          indx = index(line,'MARKER NAME         ')
          if( indx.gt.12 ) xf_makernam = line(1:8)

*         See if operator name
          indx = index(line,'OBSERVER / AGENCY   ')
          if( indx.gt.12 ) xf_openam = line(1:40)
                    
*         See if receiver type and serial # , software version
          indx = index(line,'REC # / TYPE / VERS')
          if( indx.gt.12 ) then
              read(line, 100,iostat=jerr ) xf_rcvnum, xf_rctype, swver
 100          format(a20,a20,a20)
               call report_error('IOSTAT',jerr,'decod',line,0,'RXHEAD')
          end if
          
*         See if antenna type and serial #.
          indx = index(line,'ANT # / TYPE       ')
          if( indx.gt.12 ) then
            read(line, 100, iostat=jerr ) xf_antnum, xf_anttyp
               call report_error('IOSTAT',jerr,'decod',line,0,'RXHEAD')
          end if
                     
*         See if position
          indx = index(line,'APPROX POSITION XYZ')
          if( indx.gt.12) then
               read(line,*,iostat=jerr) xf_pos  
               call report_error('IOSTAT',jerr,'decod',line,0,'RXHEAD')
               xf_posflg = 'XYZ'
          end if

*         See if antenna offset to reference point.
          indx = index(line,'ANTENNA: DELTA H/E/N')
          if( indx.gt.12 )  read(line, *,iostat=jerr ) xf_offarp
               call report_error('IOSTAT',jerr,'decod',line,0,'RXHEAD')

  
*         See if data types
          indx = index(line,'# / TYPES OF OBSERV')
          if( indx .gt.0 ) then
* MOD TAH 090114: Allow for more than 9 entries on a line
              read(line,*,iostat=jerr) xf_ndat,
     .                    (datatypes(i),i=1, min(xf_ndat,9)) 
              xf_msat = 0
              call report_error('IOSTAT',jerr,'decod',line,0,'RXHEAD')
              if( xf_ndat.gt.xf_maxdat ) then
                  write(*,180) xf_ndat, xf_maxdat 
 180              format('**DISASTER** Too many observables. ',
     .                   I3,' in current RX file, ',i3,
     .                   ' set in xf_maxdat')
                  stop 'TRACK DISASTER: Too many observables'
              end if

* MOD TAH 150411: More general code that can extended when more 
*             observables are needed.
****          See if we need to read next line (if xf_ndat<=9 only one
*             line is needed and so no aditional lines are read.
              nline = int((xf_ndat-1)/9)
              do k = 1, nline 
                  read(obs_lu,'(a)',iostat=ierr) line  
                  call report_error('IOSTAT',ierr,'read',obs_file,
     .                               0,'RXHEAD')
                  lent = 9
                  if( k.eq.nline ) lent = xf_ndat-nline*9
                  read(line,*,iostat=jerr)
     .                    (datatypes(i+k*9),i=1, lent)
                  call report_error('IOSTAT',jerr,'decod',line,0,
     .                    'Read Datatypes')

              enddo 

* MOD TAH 100717: Changed C1 to type 5 and test later and use only if
*             P1 not available.
* MOD TAH 130425: Added C2 as type 8, used if P2 is not available
* MOD TAH 180322: Moved this cocde because assignment will depend on
*             GNSS system being processed.  Assign_dattp is new subroutine
C             xf_mdat = 0
C             do i = 1, xf_ndat 
C                xf_dattyp(i) = 0
C                if (datatypes(i).eq.'L1') xf_dattyp(i) = 1
C                if (datatypes(i).eq.'L2') xf_dattyp(i) = 2
C                if (datatypes(i).eq.'P1') xf_dattyp(i) = 3
C                if (datatypes(i).eq.'P2') xf_dattyp(i) = 4 
C                if (datatypes(i).eq.'C1') xf_dattyp(i) = 5               
C                if (datatypes(i).eq.'D1') xf_dattyp(i) = 6
C                if (datatypes(i).eq.'D2') xf_dattyp(i) = 7
C                if (datatypes(i).eq.'C2') xf_dattyp(i) = 8
C                if (xf_mdat.lt.xf_dattyp(i)) 
C    .                               xf_mdat = xf_dattyp(i)
C             end do
              xf_mdat = 6   ! Expect observations through to C2 
              call assign_dattyp( datatypes )
              if( debug.ge.-2  ) then 
                 print *,'OBS: # xf_ndat ', xf_ndat ,' Types ',
     .                                 datatypes(1: xf_ndat) 
C                do i = 1,4
C                  print *,'XF_DATTYP SYS ',i,' MAP ', xf_dattyp(1:6,i)
C                end do
              end if
          end if

*         See if first obs. epoch.
          indx = index(line,'TIME OF FIRST OBS   ')
          if( indx.gt.12 ) then
c            read(line, *) xf_srttime, xf_srtsec
             read(line, *, iostat=jerr) xf_srttime, xf_srtsec
               call report_error('IOSTAT',jerr,'decod',line,0,'RXHEAD')
             call ymdhms_to_mjd(xf_srttime, xf_srtsec, xf_start_mjd)  
          end if
          
*         See if  last obs. epoch.
          indx = index(line,'TIME OF LAST OBS   ')
          if( indx.gt.12 ) then
              read(line, *, iostat=jerr) xf_endtime, xf_endsec
               call report_error('IOSTAT',jerr,'decod',line,0,'RXHEAD')
              call ymdhms_to_mjd(xf_endtime, xf_endsec, xf_end_mjd)
          endif

*         See if  # of satellites.
          indx = index(line,'# OF SATELLITES    ')
          if( indx.gt.12 )  then
              read(line, *, iostat=jerr) xf_nsat
               call report_error('IOSTAT',jerr,'decod',line,0,'RXHEAD')
          end if
          

*         See if  collection inteval.
          indx = index(line,'INTERVAL            ')
          if( indx.gt.12 )    read(line, *,iostat=jerr) xf_ircint
               call report_error('IOSTAT',jerr,'decod',line,0,'RXHEAD')

      end do

 
*     See if we found rinex file
      if( .not. rinex ) then
          write(*,170) obs_file(1:trimlen(obs_file))
 170      format(a,' Does not appear to be RINEX file')
          stop 'read_rinex_head: Wrong type data file'
      end if
 
****  Thats all
      return
      end
 
