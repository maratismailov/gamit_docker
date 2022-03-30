      subroutine trt_cpcoml( argc, port, machine, cpp_cmdfile, 
     .          cpp_prtroot, cpp_refsite, cpp_nums, cpp_prcsite ) 

      implicit none

*     Routine to copy C-read command line arguments to fortran common

      include 'trackRTObs.h'
      include 'trackRT.h'

* PASSED Variables

      integer*4 argc   ! Number of arguments
     .,         port   ! Port for comms
     .,         cpp_nums   ! Number of reference sites
     .,         trimlen    ! Length of sting
      logical endfnd    ! Set true when \0 found in cmdfile name

      character*80 machine  ! Comms machine
      character*256 cpp_cmdfile  ! Track batch file
      character*256 cpp_prtroot  ! root string for print file

      character*4 cpp_refsite    ! Reference site
     .,           cpp_prcsite(max_site) ! Site to process

* LOCAL 
      integer*4 i, j  ! Loop counter

****  Copy over values
      write(*,100) trackRT_version
 100  format(/,'++ TrackRT Version ',a,' ++',/)
      if( argc.le.1 ) then
*        No command line, so print help
         call proper_runstring('trackRT.hlp','trackRT',1)
      end if

*     Save the informaiton passed 
      rt_port = port
*     Find end of machine name
      endfnd = .false.
      rt_machine = machine
      i = 0
      do while( i.lt.80 .and. .not. endfnd )
         i = i + 1
         if( ichar(rt_machine(i:i)).eq.0 ) then
             endfnd = .true.
             rt_machine(i:) = ' '
         endif
      enddo 
*     Now find end of command file name
      endfnd = .false.
      i = 0
      bat_file = cpp_cmdfile
      do while( i.lt.256 .and. .not. endfnd )
         i = i + 1
         if( ichar(bat_file(i:i)).eq.0 ) then
             endfnd = .true.
             bat_file(i:) = ' '
         endif
      enddo 

*     Now find end of prt root string
      endfnd = .false.
      i = 0
      prt_root = cpp_prtroot
      do while( i.lt.256 .and. .not. endfnd )
         i = i + 1
         if( ichar(prt_root(i:i)).eq.0 ) then
             endfnd = .true.
             prt_root(i:) = ' '
         endif
      enddo 

*     Copy over the ref site name
      ref_code = cpp_refsite
      call casefold(ref_code)
      if( ichar(ref_code(1:1)).eq.0 ) then
          ref_code = ' '
          j = 0
      else
          site_names(1) = ref_code
          site_type(1) = 1   ! Leave as float for moment
          j = 1
      end if
      num_prc = cpp_nums
      do i = 1, num_prc
         if( cpp_prcsite(i).ne.ref_code ) then
             j = j + 1
             site_names(j) = cpp_prcsite(i)
             call casefold(site_names(j))
             site_type(j) = 1  ! Mark as kinematic
         endif
      end do
*     Save number of stations
      num_prc = j
      num_site = num_prc

*     Now see if have all we need.
      if ( len_trim(bat_file).eq.0 ) then
         write(*,120) 
 120     format('No command file; use -f option')
         call proper_runstring('trackRT.hlp','trackRT',1)
      endif

*     Report setup 
      write(*,210) trim(rt_machine), rt_port
 210  format('TrackRT: Data steam from ',a,' Port ',i5)
      write(*,230) trim(bat_file)
 230  format('TrackRT: Command file   : ',a)
      if( len_trim(ref_code).gt.0 ) then
         write(*,250) ref_code
 250     format('TrackRT: Reference Site : ',a,/,
     .          'TrackRT: Sites to be processed')
      else
         ref_code = ' '
      endif 

      do i = 1, num_prc
         write(*,270) i,site_names(i)
 270     format('TrackRT: ',i2,' Site ',a)
      end do
 
      return
      end

      subroutine saveobs( nobs, rtobs )

      implicit none

*     F90 routine that interfaces to the c++ routine trackRTComm.cpp
*     to save one data record from the real-time stream.  

      include 'trackRTObs.h'
!      use trackRTMod

* PASSED variables
      integer*4 nobs  ! Number of observations so far

      type ( obsInternal ) rtobs

****  The number of observations is passed from the C++ Comm
*     routine
!      rtobsArray(nobs) = rtobs
      RT_flags(nobs) =        rtobs%flags         
      RT_StatID(nobs) =       rtobs%StatID        ! Station ID
      RT_satSys(nobs) =       rtobs%satSys        ! Satellite System ('G' or 'R')
      RT_satNum(nobs) =       rtobs%satNum        ! Satellite Number (PRN for GPS NAVSTAR)
      RT_slot(nobs) =         rtobs%slot          ! Slot Number (for Glonass)
      RT_GPSWeek(nobs) =      rtobs%GPSWeek       ! Week of GPS-Time
      RT_GPSweeks(nobs) =     rtobs%GPSweeks      ! Second of Week (GPS-Time)
      RT_C1(nobs) =           rtobs%C1            ! CA-code pseudorange (meters)
      RT_C2(nobs) =           rtobs%C2            ! CA-code pseudorange (meters)
      RT_P1(nobs) =           rtobs%P1            ! P1-code pseudorange (meters)
      RT_P2(nobs) =           rtobs%P2            ! P2-code pseudorange (meters)
      RT_L1(nobs) =           rtobs%L1            ! L1 carrier phase (cycles)
      RT_L2 (nobs) =          rtobs%L2            ! L2 carrier phase (cycles)
      RT_slip_cnt_L1(nobs) =  rtobs%slip_cnt_L1   ! L1 cumulative loss of continuity indicator (negative value = undefined)
      RT_slip_cnt_L2(nobs) =  rtobs%slip_cnt_L2   ! L2 cumulative loss of continuity indicator (negative value = undefined)
      RT_lock_timei_L1(nobs) =rtobs%lock_timei_L1 ! L1 last lock time indicator                (negative value = undefined)
      RT_lock_timei_L2(nobs) =rtobs%lock_timei_L2 ! L2 last lock time indicator                (negative value = undefined)
      RT_S1(nobs) =           rtobs%S1            ! L1 signal-to noise ratio
      RT_S2(nobs) =           rtobs%S2            ! L2 signal-to noise ratio
      RT_SNR1(nobs) =         rtobs%SNR1          ! L1 signal-to noise ratio (mapped to integer)
      RT_SNR2(nobs) =         rtobs%SNR2          ! L2 signal-to noise ratio (mapped to integer)

!     Compute the MLD of this obs
      RT_MJD_obs(nobs) =  44244 + RT_GPSWeek(nobs)*7 + 
     .                    RT_GPSweeks(nobs)/86400 

!     Set the initial error flags (Maybe RT_flag could be used too?)
      RT_errflag(nobs) = 0 

*     Now make sure we have 'P1' and 'P2' ranges
      if( RT_P1(nobs).eq.0 ) then
          RT_P1(nobs) = RT_C1(nobs)
          call sbit(RT_flags(nobs),16,1)
      endif
      if( RT_P2(nobs).eq.0 ) then
          RT_P2(nobs) = RT_C2(nobs)
          call sbit(RT_flags(nobs),17,1)
      endif

     
      num_rtobs = nobs

c     write(*,120) nobs, rtobs%satSys, rtobs%satNum, rtobs%StatID,
c    .             rtobs%GPSWeek, rtobs%GPSWeeks,
c    .             rtobs%P1, rtobs%C1
c 120  format('Save ',i4,1x,a1,I2.2,1x,a9,1x,I4,1x,F12.4,1x,
c     .       'P1 C1',2F15.3)

!      print *,'Num Data ', nobs

!      print *,'SAT ',rtobs%satNum,' ',rtobs%flags, ' ', rtobs%satSys
!      print *,'Data L1 ',rtobs%P1, rtobs%P2, rtobs%L1, rtobs%L2

      return 
      end

      subroutine procobs( nobs )

      implicit none

*     F90 routine that will process the realtime data saved at the
*     current epoch  

      include 'trackRTObs.h'
!      use trackRTMod

* PASSED
      integer*4 nobs  ! Number of observations from trackRTComm; should
                      ! match value saved in common

* LOCAL
      integer*4 i, j   ! Loop counters

!     Make sure that the number of obs match
      write(*,120) nobs, num_rtobs
 120  format('ProcObs Number of obs ',2I6)

!      Just list the data we have
      do i = 1, num_rtobs, 20
!         write(*,220) i, rtobsArray(i)%satSys, 
!     .             rtobsArray(i)%satNum,  rtobsArray(i)%StatID,
!     .             rtobsArray(i)%GPSWeek, rtobsArray(i)%GPSWeeks,
!     .             rtobsArray(i)%P1
         write(*,220) i, RT_satSys(i), 
     .             RT_satNum(i),  RT_StatID(i),
     .             RT_GPSWeek(i), RT_GPSWeeks(i),
     .             RT_P1(i)
 220     format('Proc ',i4,1x,a1,I2.2,1x,a9,1x,I4,1x,F13.5,1x,
     .          ' P1 ',F15.3)
      end do
      write(*,320)
 320  format('End of Process',/)
c      write(*,'(a)') 'Type number to continue'
c      read(*,*) j

      return
      end 
