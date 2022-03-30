      program trackRTr

      implicit none

*     Program to read rinex files and feed into the trackRT processing
*     engine.  Here we have another command line reader and code to
*     read rinex files.

      include 'trackRT.h'         ! Common block
      include 'trackRTObs.h'      ! Real-time data structures 

      integer*4 ep   ! Counter for number epochs processed.
      integer*4 pass_debug(10)   ! Debug array passed from read_bat

****  OK Decode the runstring

      call trtr_coml

      call read_batrt( 0 , pass_debug )
     
      num_rtobs = 1
      rx_not_open = .true.
      ep = 0

      do while ( num_rtobs.gt. 0 )

         call rx_to_obsInt
         call procobsRT( num_rtobs, ep )

      end do

      end

CTITLE TRTR_COML

      subroutine trtr_coml

      implicit none

      include 'trackRT.h'         ! Common block

****  Local

      integer*4 i   ! Loop counter
     .,   ierr      ! IOSTAT error
     .,   n, ns         ! Counter for finding sites names
 
      character*256 coml   ! Command line argument
      character*4 insite   ! Name of site from rxfile name

      logical done  ! Set true when command line finished
     .,       rxend ! Set true at end of list of rx files


****  Loop over the runstring arguments
      done = .false.
      i = 0
      num_site = 0
      bat_file = ' '
      prt_root = ' '
      ref_code = ' '
      rt_machine = '127.0.0.1'
      rt_port = 0
      do while ( .not. done )
          i = i + 1
          call rcpar(i,coml)
          if( coml(1:2).eq.'-m' ) then
             i = i + 1
             call rcpar(i,rt_machine)
          else if ( coml(1:2).eq.'-p' ) then
             i = i + 1
             call rcpar(i,coml)
             read(coml,*,iostat=ierr) rt_port
          else if( coml(1:2).eq.'-f' ) then 
             i = i + 1
             call rcpar(i,bat_file)
          else if( coml(1:2).eq.'-n' ) then 
             i = i + 1
             call rcpar(i,prt_root)
          else if( coml(1:2).eq.'-r' ) then 
             i = i + 1
             call rcpar(i,ref_code)
             call casefold(ref_code)
*            See if matches an existing site
             ns = 0
             do n = 1, num_site
                if( ref_code.eq.site_names(n) ) then
                    ns = n
                    exit
                endif
             end do
             if( ns.eq.0 ) then
                ! New site, so add
                num_site = num_site + 1
                if( num_site.gt.max_site ) then
                   call report_stat('fatal','trackrt','trtr_coml',
     .                ' ','Too many sites: Max allowed', max_site)
                endif
                site_names(num_site) = ref_code
             end if

          else if ( coml(1:2).eq.'-d' ) then
*            Start scanning for rinex files 
             rxend = .false.
             do while ( .not. rxend )
                 i = i + 1
                 call rcpar(i,coml)
                 if( coml(1:1).eq.'-' ) then
                     rxend = .true.
                     i = i - 1
                 else if( len_trim(coml).gt.0 ) then  ! Extract file name

*                    Extract the site name from the RXfile to see if
*                    already given (only one file per site).
                     n = len_trim(coml)
                     do while ( n.gt.1 )
                         n = n - 1
                         if( coml(n:n).eq.'/' ) then
                             insite = coml(n+1:n+4)
                             n = -1
                         endif
                     enddo
                     if( n.eq.1 ) insite = coml(1:4)
                     call casefold(insite)
*                    See if we have this site already
                     ns = 0
                     do n = 1, num_site
                        if( insite.eq.site_names(n) ) then
                            ns = n
                            exit
                        endif
                     end do
                     if( ns.eq.0 ) then
                        ! New site, so add
                        num_site = num_site + 1
                        ns = num_site
                     end if
                     if( num_site.gt.max_site ) then
                        call report_stat('fatal','trackrt','trtr_coml',
     .                      ' ','Too many sites: Max allowed', max_site)
                     endif
                     site_names(ns) = insite
                     if( site_names(ns).eq.ref_code ) then
                         site_type(ns) = 1   ! Leave as kinematic
                     else   ! Mark as kinematic
                         site_type(ns) = 1
                     end if
*                    Save the rinex file name
                     rx_file(ns) = coml
 
                endif
                if( len_trim(coml).eq.0 ) rxend = .true.
            end do
          end if
          if( len_trim(coml).eq.0 ) done = .true.
      end do
*
****  Make sure we have a reference site
      if( ref_code(1:1).eq.' ' ) then
         ref_code = site_names(1)
         site_type(1) = 1
      end if

      num_prc = num_site

*     Now see if have all we need.
      if ( len_trim(bat_file).eq.0 ) then
         write(*,120) 
 120     format('No command file; use -f option')
         call proper_runstring('trackRT.hlp','trackRT',1)
      endif

*     Report setup 
      write(*,210) trim(rt_machine), rt_port
 210  format('TrackRTr: Data steam from ',a,' Port ',i5)
      write(*,230) trim(bat_file)
 230  format('TrackRTr: Command file   : ',a)
      if( len_trim(prt_root).gt.0 ) then
         write(*,240) trim(prt_root)
 240     format('TrackRTr: Print root     : ',a)
      else
         write(*,240) 'stdout'
      endif

      if( len_trim(ref_code).gt.0 ) then
         write(*,250) ref_code
 250     format('TrackRTr: Reference Site : ',a,/,
     .          'TrackRTr: Sites to be processed')
      else
         ref_code = ' '
      endif 

      do i = 1, num_prc
         write(*,270) i,site_names(i)
 270     format('TrackRT: ',i2,' Site ',a)
      end do

      end
 
 
CTITLE RX_TO_OBSINT

      subroutine rx_to_obsInt

      implicit none

*     Routine to read rinex file records and return the data
*     in the TrackRT data structure.  Normally this is done
*     with trackRT_saveobs routine in the realtime comm module.

      include 'trackRT.h'         ! Common block
      include 'trackRTObs.h'      ! Real-time data structures 

      integer*4 obs_lu  ! Unit numbers for the rinex file
     .,   i,j,k,l       ! Loop counters
     .,   ierr, jerr    ! IOSTAT errors
     .,   trimlen       ! Length of string
     .,   indx          ! Index in string
     .,   date(5)       ! YMDHM current
     .,   id            ! Record ID read from rx file
     .,   rx_msat, tm_msat  ! number of satellite records at
                        ! current epic
     .,   rx_iprn(max_sat)  ! PRN number
     .,   itmp            ! Temp counter for multiple lines
     .,   ilck(max_sat)   ! Lock flag
     .,   flags(max_rxtype, max_sat)  ! Flags by satellite
     .,   num, n, nl        ! tmp integers
     .,   gpsweek       ! GPS week

      real*8 sectag     ! Seconds tag
     .,    mjd_curr     ! MJD of current obs
     .,    vals(max_rxtype,max_sat)  ! Values read from RX file
     .,    gpssec       ! GPS seconds of week

      real*8 rx_getObs  ! Function to return selected data type

     
 
      character*80 line  ! Line from rinex file
     .,            line2 ! Second line for extended lines
     .,            save_line1(max_site), save_line2(max_site)
                         ! Saved version of PRN lines

      character*1 cr     ! Cariage return
     .,       rx_svtype(max_sat) ! Satellite types

      character*2 datatypes(max_rxtype)  ! L1,L2 etc data types

      logical eoh   ! End of header
     .,       sync  ! Set true when the times in current RX file 
                    ! are within 10 msec of reference site
     .,       read_already(max_site)  ! Set true when header line already
                    ! read.  Happens is reference site is
                    ! missing an epoch

      save save_line1, save_line2, read_already


      cr = char(13)

****  See if the rinex files have been open yet
      if( rx_not_open ) then
          rx_not_open = .false.
          do i = 1, num_site
             obs_lu = 100 + i
             open(obs_lu,file=rx_file(i),status='old',iostat=ierr)
             call report_error('IOSTAT',ierr,'open',rx_file(i),
     .            1,'rx_to_obsInt')

*            Now skip over the header.  Here we do not use the header
*            because this information is supplied externally.  We do
*            need to extract the types of data in the files
             eoh = .false. 
             do while ( .not.eoh )
                read(obs_lu,'(a)', iostat=ierr) line 
                call sub_char(line,cr,' ')
                call casefold(line)
                call report_error('IOSTAT',ierr,'read', rx_file(i),1,
     .                      'End of header not found')
                if( trimlen(line).lt.10) eoh = .true.
                if( index(line,'END OF HEAD').gt.0 ) eoh = .true.
   
*               See if data types
                indx = index(line,'# / TYPES OF OBSERV')
                if( indx .gt.0 ) then
* MOD TAH 090114: Allow for more than 9 entries on a line
                    read(line,*,iostat=jerr) rx_ndat(i),
     .                          (datatypes(j),j=1, max(rx_ndat(i),9)) 
                    rx_msat = 0
                    call report_error('IOSTAT',jerr,'decod',line,
     .                    0,'RXHEAD')
                    if( rx_ndat(i).gt.32 ) then
                        write(*,180) rx_ndat, 32 
 180                    format('**DISASTER** Too many observables. ',
     .                          I3,' in current RX file, ',i3,
     .                          ' set in rx_maxdat')
                        stop 'TRACK DISASTER: Too many observables'
                    end if
****                See if we need to read next line
                    if( rx_ndat(i).gt. 9 ) then
                        read(obs_lu,'(a)',iostat=ierr) line  
                        call report_error('IOSTAT',ierr,'read',
     .                                   rx_file(i),  0,'rx_to_obsInt')
                        read(line,*,iostat=jerr) rx_ndat(i),
     .                          (datatypes(j),j=10, rx_ndat(i))
                    endif 

                    do j = 1, rx_ndat(i) 
                       rx_dt(j,i) = 0
                       if (datatypes(j).eq.'L1') rx_dt(j,i) = 1
                       if (datatypes(j).eq.'L2') rx_dt(j,i) = 2
                       if (datatypes(j).eq.'P1') rx_dt(j,i) = 3
                       if (datatypes(j).eq.'P2') rx_dt(j,i) = 4 
                       if (datatypes(j).eq.'C1') rx_dt(j,i) = 5               
                       if (datatypes(j).eq.'C2') rx_dt(j,i) = 6               
                       if (datatypes(j).eq.'S1') rx_dt(j,i) = 7
                       if (datatypes(j).eq.'S2') rx_dt(j,i) = 8
                       if (datatypes(j).eq.'D1') rx_dt(j,i) = 9
                       if (datatypes(j).eq.'D2') rx_dt(j,i) = 10
                    end do
                end if
             end do
             read_already(i) = .false.
          end do   ! Looping over stations
      end if       ! Rinex files not open

****  Files are open; loop over stations returning data current time (reference
*     site site 1 must have data)
      num_rtobs = 0
      do i = 1, num_site
c     do i = num_site, 1, -1
c          print *,'Site ',i, rx_ndat(i), mjd_curr, num_rtobs, 
c     .            read_already(i), sync
          obs_lu = 100 + i

*         Loop to allow syncing of observations with the first one
          sync = .false.
          do while ( .not.sync )
              if( .not. read_already(i) ) then 
                  read(obs_lu,'(a80)', iostat=ierr) line
                  if( ierr.ne.0 ) RETURN
              else
                  line = save_line1(i)
              end if
*             Check to see if rinex header line (from concatinated
*             file).  If so then skip new header
              call casefold(line)
              indx = index(line,'OBSERVATION DATA')
              if( indx.gt.0 ) then
                  write(*,'(a,1x,a)') 'New RINEX header found in',
     .                  trim(rx_file(i))
*                 Skip to end of header
                  eoh = .false.
                  do while ( .not.eoh .and. ierr.eq.0 )
                      read(obs_lu,'(a)', iostat=ierr) line 
                      call sub_char(line,cr,' ')
                      call casefold(line)
                      call report_error('IOSTAT',ierr,'read', 
     .                      rx_file(i),1,'End of header not found')
                      if( trimlen(line).lt.10) eoh = .true.
                      if( index(line,'END OF HEAD').gt.0 ) eoh = .true.
                  end do
*                 Now read next line
                  read(obs_lu,'(a80)', iostat=ierr) line
              end if

****          Continue reading file.
              read(line,210,iostat=ierr) date, sectag, id, rx_msat,
     .               (rx_svtype(j),rx_iprn(j),j=1,min(rx_msat,12))
 210          format(5i3,f11.7,i3,i3,24(a1,i2))
c              write(*,210) date, sectag, id, rx_msat,
c     .               (rx_svtype(j),rx_iprn(j),j=1,min(rx_msat,12))
* MOD TAH 061220: Only read the extended line if the ID is not
*     between 2 and 4 (comment and stop/go records)
              if( rx_msat.gt.12 .and. (id.lt.2 .or. id.gt.4)) then
                  tm_msat = rx_msat - 12
                  if( .not.read_already(i) ) then
                     read(obs_lu,'(a)',iostat=ierr) line2
                  else
                     line2 = save_line2(i)
                  endif
                  read(line2,215,iostat=ierr) 
     .                (rx_svtype(j+12),rx_iprn(j+12),j=1,tm_msat)
 215              format(32x,12(a1,i2))  
              end if
    
*             Convert the times to MJD and make sure matches current value
              call ymdhms_to_mjd(date, sectag, mjd_curr)
              gpsweek = int((mjd_curr - 44244)/7)
              gpssec  = (mjd_curr - (gpsweek*7 + 44244))*86400.d0

*             See if times match
              if( num_rtobs.gt.0 ) then
*                 OK: we have data already, see if this obs is the same
*                 or before or after
                  if( mjd_curr-RT_MJD_obs(1).lt.
     .                                    -0.01d0/86400.d0 ) then
*                    This observation is before current one so we need
*                    skip over data and get to next line
                     do j = 1,rx_msat*(int((rx_ndat(i)-1)/5)+1)
                        read(obs_lu,'(a80)', iostat=ierr) line
                        if( ierr.ne.0 ) then
                            write(*,250) ierr, site_names(i), mjd_curr
 250                        format('Unexpected IOSTAT error ',i5,
     .                         ' read ',a,' rinex file, MJD ',F14.6)
                            RETURN
                        end if
                     end do
                     read_already(i) = .false.
                  elseif( mjd_curr-RT_MJD_obs(1).ge.
     .                                          +0.01d0/86400.d0 ) then
*                    See if data is latter in which case we save lines and
*                    wait
                     sync = .true.
                     read_already(i) = .true.
                     save_line1(i) = line
                     save_line2(i) = line2
                  else
                     read_already(i) = .false.
                     sync = .true.
                  endif
              else
*                  This is reference first site so say we are sync but that
*                  we have not read header already
                   read_already(i) = .false.
                   sync = .true.
              endif

              if( sync .and. .not.read_already(i) ) then
*                 Read and save the data
*                 Now loop over the data records.
                  itmp = min(rx_ndat(i),5)
                  nl = (rx_ndat(i)-1)/5+1
                  do j = 1, rx_msat
                     do k = 1, nl 
                        read(obs_lu,'(a)',iostat=ierr) line
                        call report_error('IOSTAT',ierr,'read',
     .                         rx_file,0, 'rx_to_obsInt')

                        call sub_char(line,cr," ") 
                        if( k.lt.nl ) then
                            itmp = min(5*k,rx_ndat(i))
                        else
                            itmp = rx_ndat(i)
                        endif
                        if( ierr.eq.0 )
     .                  read(line,120,iostat=jerr) (vals(l,j), 
     .                      ilck(l), flags(l,j), l =(k-1)*5+1,itmp)
 120                    format(5(f14.3,i1,i1))
                        call report_error('IOSTAT',jerr,'read',line,0,
     .                         'rx_to_obsInt')
                     end do
                  end do

****              OK: Now save the data into the track_RT obsInternal
*                 arrays
                  num = rx_ndat(i)
                  do j = 1, rx_msat
                      num_rtobs = num_rtobs + 1
                      n = num_rtobs
                      RT_flags(n) = flags(1,j)
                      RT_StatID(n) = site_names(i)
                      RT_satSys(n) = rx_svtype(j)  ! Satellite System ('G' or 'R')
                      RT_satNum(n) = rx_iprn(j)   ! Satellite Number (PRN for GPS NAVSTAR)
                      RT_slot(n)   = 0             ! Slot Number (for Glonass)
                      RT_GPSWeek(n) = gpsweek        ! Week of GPS-Time
                      RT_GPSweeks(n) = gpssec        ! Second of Week (GPS-Time)
                      RT_C1(n) = rx_getObs(num,vals(1,j),rx_dt(1,i),5) ! CA-code pseudorange (meters)
                      RT_C2(n) = rx_getObs(num,vals(1,j),rx_dt(1,i),6) ! CA-code pseudorange (meters)
                      RT_P1(n) = rx_getObs(num,vals(1,j),rx_dt(1,i),3) ! P1-code pseudorange (meters)
                      RT_P2(n) = rx_getObs(num,vals(1,j),rx_dt(1,i),4) ! P2-code pseudorange (meters)
                      RT_L1(n) = rx_getObs(num,vals(1,j),rx_dt(1,i),1) ! L1 carrier phase (cycles)
                      RT_L2(n) = rx_getObs(num,vals(1,j),rx_dt(1,i),2) ! L2 carrier phase (cycles)
                      RT_slip_cnt_L1(n) = 0     ! L1 cumulative loss of continuity indicator (negative value = undefined)
                      RT_slip_cnt_L2(n) = 0     ! L2 cumulative loss of continuity indicator (negative value = undefined)
                      RT_lock_timei_L1(n) = 0   ! L1 last lock time indicator                (negative value = undefined)
                      RT_lock_timei_L2(n) = 0   ! L2 last lock time indicator                (negative value = undefined)
                      RT_S1(n) = rx_getObs(num,vals(1,j),rx_dt(1,i),7) ! L1 signal-to noise ratio
                      RT_S2(n) = rx_getObs(num,vals(1,j),rx_dt(1,i),8) ! L2 signal-to noise ratio
                      RT_SNR1(n) = 0  ! L1 signal-to noise ratio (mapped to integer)
                      RT_SNR2(n) = 0  ! L2 signal-to noise ratio (mapped to integer)

                      RT_MJD_obs(n) = mjd_curr   ! MLD of epopch

*                     Set the initial error flags (Maybe RT_flag could be used too?)
                      RT_errflag(n) = 0 

*                     Now make sure we have 'P1' and 'P2' ranges
                      if( RT_P1(n).eq.0 ) then
                          RT_P1(n) = RT_C1(n)
                          call sbit(RT_flags(n),16,1)
                      endif
                      if( RT_P2(n).eq.0 ) then
                          RT_P2(n) = RT_C2(n)
                          call sbit(RT_flags(n),17,1)
                      endif

                  end do
              end if 
          end do
      end do

****  Thats all
      return
      end

CTITLE RX_GETOBS

      real*8 function rx_getObs( num, vals, dattyp, code )
      implicit none

      integer*4 num            ! Number of data types available      
      real*8 vals(num)         ! Values read from rinex file
      integer*4 dattyp(num)    ! Code for data type 1-9
      integer*4 code           ! Code type to return ( 0 if not available)

      integer*4 j  ! Loop counter 

      rx_getObs = 0.d0
      do j = 1, num
         if( dattyp(j).eq.code ) rx_getObs = vals(j) 
      end do

***** Thats all
      return
      end


