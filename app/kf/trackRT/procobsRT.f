      subroutine procobsRT( num_rtpass, ep )

      implicit none

*     TrackRT and trackRTr main processing engine.  Here we start with the
*     time of this set of measurements and update the filter state.

      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 num_rtpass    ! Number of data passed (should
                              ! equal num_rtobs)
      integer*4 ep            ! Counter for number of epochs of data processed.

* LOCAL
      integer*4 i, j, blk   ! Loop counter
      integer*4 iter
      integer*4 ierr     ! IOSTAT error on update file
      integer*4 pass_debug(10)   ! Debug array passed from read_bat
      logical  OK
      character*256 line  

      if( ep.ge.usr_nepochs .and. usr_nepochs.gt. 0 ) 
     .       stop 'USR_NEPOCHS Reached'

****  Make sure values match
      if( num_rtpass.ne.num_rtobs ) then
         write(line,110) num_rtpass, num_rtobs
 110     format('NUM_RTPASS ',i4,' and NUM_RTOBS ',i4,
     .          ' do not match; return')
         call report_stat('WARNING','TrackRT','procobsRT',' ',
     .        line, num_rtpass*10000+num_rtobs)
         return 
      endif

****  Make sure we have data
      if( num_rtobs.eq.0 ) RETURN

****  See if before requested start time (default use_start=0)
      if ( RT_MJD_obs(1).lt. usr_start ) RETURN 

****  If this is epoch zero, then initialize the estimation system
      if( ep.eq.0 ) then
          call init_estRT
          kf_prev_mjd = RT_MJD_obs(1)
          kf_step = 0
          call read_dcbsRT(dcb_file, kf_prev_mjd)
      end if

****  Test the epochs in this group.  If there has been a delay in processing
*     a couple of epochs may be here together
      call split_epoch(ep) 

****  Main stages of the processing
      do blk = 1, num_blk

*        Increment counter for number of epochs
         ep = ep + 1

*        Set the start and stop records for this epoch
         if( blk.eq.1 ) then
            sblk = 1
         else
            sblk = blks(blk-1)+1
         endif

         eblk = blks(blk)

*        Sort the data so that sites are in acending order
* MOD 100505: Now done in split_epoch
         call get_blk( blk, ep) 

*****    Make sure we have reference station data
         if( ssblk(1).eq.0 ) then
              call report_stat('WARNING','TrackRT','procobsRT',' ',
     .            'No Reference station data,  Epoch ',ep)
              exit
         endif
 
****     Compute the time step taken
         kf_step = (RT_MJD_obs(sblk)-kf_prev_mjd)*86400

****     If the time step is zero; then there a problem.  Dump
*        what is happening and skip this zero step
         if( abs(kf_step).lt.1.d-2 .and. ep.gt. 1 ) then
*            Step seems to be zero.  Dump what is happening.
             write(*,120) kf_step, blk, num_blk, sblk, eblk, 
     .           blks(1:num_blk)
 120         format('ERROR: ZERO Time step ',E13.3,' sec in block ',
     .              I2.2,'/',I2.2,' Start/Stop ',2i4,' Blocks ',
     .              100I5)
             write(*,140) (j, RT_sitenum(j), site_names(RT_sitenum(j)),
     .              RT_satNum(j),RT_GPSWeek(j), RT_GPSweeks(j),
     .              RT_MJD_obs(j),j=1,num_rtobs)
 140         format('RT OBS Times',/,
     .          1000('Ent ',i4,' Site ',i3,1x,a,1x,' G ',i2.2,
     .               ' Week/Sec ',I4,1x,F16.6,1x,' MJD ',F14.6,/))
             exit   ! Skip and try next block
         end if

         kf_prev_mjd = RT_MJD_obs(sblk)
         if( (ep.ge.debug(1) .or. ep.ge.debug(3)) .and.
     .       (ep.le.debug(2) .or. ep.le.debug(4)) )
     .   write(*,220) ep, blk, num_blk, sblk, eblk, num_rtobs, kf_step
 220     format(/,140('_'),/,
     .          'Processing Epoch ',i5,' Blk/Num ',I2.2,'/',I2.2,
     .          ' Start/End ',2i5,' Obs ',i4,
     .          ' Time Step ',F6.2,' seconds')

*        See if we need to re-read a command file
         if( needupd ) then
*           See if we have already read the file 
*           See if file exists
            open(100,file=upd_file,status='old',iostat=ierr)
            
            if( ierr.eq.0 .and. .not.updread ) then
                call read_batrt(1, pass_debug)
                updread = .true.

****            See if reset are needed
                call reset_state(ep, kf_prev_mjd)

            elseif( ierr.ne.0 ) then ! File no longer exits, reset
                    ! so that it will be read later if it comes back
                 updread = .false.
            end if
            close(100)
         end if

*        See if we need to update the sp3 file information.  This will
*        open and read new sp3 files as needed.
         call update_sp3

*        See if we need to update or open the output files based on the
*        rollover interval
         call update_outfs(ep)

*****    Based on current coordinates get the clock offsets
*        at each station
         call get_clkRT(ep)  


*        For the sites we are processing; compute the raw omc.  Note here
*        that the direct station clock errors are not added here
         call comp_rawomc(ep)

*        DEBUG Ouput
         if( ep.ge.debug(5) .and. ep.lt.debug(6)) then 
            do i = sblk, eblk, 1
               write(*,320) ep, i, rt_ord(i), RT_satSys(i), 
     .            RT_satNum(i),  RT_StatID(i),
     .            RT_GPSWeek(i), RT_GPSWeeks(i),
     .            RT_omcraw(:,i), RT_azel(2,i), RT_errflag(i)
 320           format('RT_OMCRAW ',i5,1x,i3,1x,i3,1x,a1,1x,I2.2,1x,
     .            a9,1x,I4,1x,F13.5,1x,' OMC ',4F15.3,1x,F7.2,1x, O3)
            end do
         end if

****     Now add and test the estimates of phase ambquities.
         call get_bflags(ep)
        
****     Now we loop until this epoch looks good
         OK = .false.
         iter = 0
         do while ( .not. OK )

*           Set OK true and then any problem or update will set it false.
            iter = iter + 1
            if( iter.gt.10 ) OK = .true.

****        Now update solution estimates.  The OK logical returns true if the
*           update is OK no process noise or data editing were needed.
            call update_posrt(ep, OK, iter )

*****       If the solution is OK (ie., no more ambiquities and editing) then
*           output values.  (In future versions, the fixed-lag filter would
*           added here).
            if( OK ) then
               call output_outfs(ep) 
               call flush
            end if

         end do

****     See if need to report status
         call report_statusRT(lus, ep )

      end do  ! Loop over blocks

****  Thats all
      return
      end

CTITLE SPLIT_EPOCH
  
      subroutine split_epoch(ep)
 
      implicit none

*     Routine to split current block of data into epochs.  Normally
*     one epoch at a time is expected but when there are delays
*     multiple epochs are possible. 

*     Routine also assigns site numbers to each entry.  If site number
*     can't found then we can't process data (could be extra sites coming
*     to port and so is not fatal)

*     RT_ord pointer list so that sites are ordered also generated here.

      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

      integer*4 ep    ! Current epoch number

      integer*4 i, j, k, n   ! Loop counter
      integer*4 ns, ne  ! Start and stop for station
      integer*4 indx      ! Index in string
      integer*4 finind(max_rtobs) ! Index for where observations need to go
                          ! so that RT_ arrays contain only sites being processed
      integer*4 sortvar(max_rtobs)  ! Variable that combines time and site number so
                          ! that we can sort by time and then site.  Should handle
                          ! 99 sites and 10-Hz data
     .,         sortind(max_rtobs)  ! Index orginally 1-num_rtobs which once sorted
                          ! tells us where values need to go.
      integer*4 fin_rtobs ! Final number of observations with sites being 
                          ! processed,

      real*8 curr_GPSs     ! Current seconds of week
      integer*4 curr_GPSw  ! current GSP week

      character*128 line

***** Remove the sites that are not being processed.
      fin_rtobs = 0

****  Assign site numbes and get list of data to be processed.
      do i = 1, num_rtobs
         indx = 1
         call get_cmd(RT_StatID(i), site_names, num_site, j, indx)
         if( j.gt.0 ) then
             RT_sitenum(i) = j
             fin_rtobs = fin_rtobs + 1
             finind(i) = fin_rtobs
         else
             RT_sitenum(i) = 0
             finind(i) = 0
         end if
      end do

***** Now move the observation arrays to remove un-used sites.
      if( num_rtobs.ne.fin_rtobs ) then
         do i = 1, num_rtobs
            if( finind(i).ne.i .and. finind(i).gt.0 ) then
*               Move the observations
                j = finind(i)
                call move_RT(i,j)
            end if
         end do
      end if

****  Now reset number if data
      num_rtobs = fin_rtobs

****  Scan through the times to make that all the data
*     are new 
      curr_GPSs = RT_GPSweeks(1)
      curr_GPSw = RT_GPSweek(1)
      if( ep.gt.0 ) then 
*         Now try to find time of new data and remove any old
*         data
          fin_rtobs = 0
          do i = 1, num_rtobs
****         See if new
             if( RT_MJD_obs(i).gt.kf_prev_mjd+5.d-3/86400.d0 ) then
                fin_rtobs = fin_rtobs + 1
                finind(i) = fin_rtobs
                if( RT_GPSweek(i)*7.d0*86400.d0+
     .              RT_GPSweeks(i).lt. 
     .              curr_GPSw*7.d0*86400.d0+curr_GPSs ) then
                    curr_GPSs = RT_GPSweeks(i)
                    curr_GPSw = RT_GPSweek(i)
                end if
             else
                write(line,110) ep, i, num_rtobs, RT_MJD_obs(i), 
     .                kf_prev_mjd  
 110            format('Time earlier or equal previous time Ep ',i6,
     .              ' Obs ',i4,' of ',i4,
     .              ' Current MJD ',F13.6,' Previous MJD ',F13.6)
                call report_stat('WARNING','TrackRT','split_epoch',' ',
     .              line, ep)
                write(lus,'(a)') trim(line)
*
                finind(i) = 0
             end if
          end do

*****     Now move the observation arrays to remove un-used sites.
          if( num_rtobs.ne.fin_rtobs ) then
             do i = 1, num_rtobs
                if( finind(i).ne.i .and. finind(i).gt.0 ) then
*                  Move the observations
                   j = finind(i)
                   call move_RT(i,j)
                end if
             end do
****         Now reset number if data
             num_rtobs = fin_rtobs
          end if
      end if

****  All data is now at a later time, now scan to sort data into
*     time order.  We do this with RT_ord array which gives the ordered 
*     list of data.
*     Form the sort index
      do i = 1,num_rtobs
          sortvar(i) = nint( (RT_GPSweek(i)-curr_GPSw)*7.d0*86400.d0 +
     .                (RT_GPSweeks(i)-curr_GPSs)*100 )*100+
     .                 RT_sitenum(i)
          sortind(i) = i
      end do

*     Now sort on sortvar(i) using a bubble sort
      do i = 1, num_rtobs - 1 
         do j = 1, num_rtobs - 1
            if( sortvar(j).gt.sortvar(j+1) ) then ! Next value is less
*               so swap values
                call switch_I4(sortvar(j), sortvar(j+1), k)
                call switch_I4(sortind(j), sortind(j+1), k)
            endif
         end do
      end do

****  Now re-arrange data into sorted form  
      call sort_blk( ep, sortind, sortvar) 

****  Final check but probably not neeed now.
      curr_GPSs = RT_GPSweeks(1)
      num_blk = 1
      do i = 2, num_rtobs
*         Allow up to 5 ms of "slack" in times.
          if( RT_GPSweeks(i).gt. curr_GPSs+5.d-3 ) then
              blks(num_blk) = i -1
              num_blk = num_blk + 1
              curr_GPSs = RT_GPSweeks(i)
          elseif( RT_GPSweeks(i).lt. curr_GPSs-5.d-3 ) then
              write(line,120) ep, i, RT_GPSweeks(i), curr_GPSs 
 120          format('Time runs backwards Ep ',i6,' Obs ',i4,
     .            ' GPS Week Seconds  ',2F13.6)
              call report_stat('WARNING','TrackRT','split_epoch',' ',
     .            line, ep)
          end if
      end do
      blks(num_blk) = num_rtobs

****  Now get the block boundaries.
      n = 0
      do k = 1, num_blk

*        Set the start and stop records for this epoch
         if( k.eq.1 ) then
            sblk = 1
         else
            sblk = blks(k-1)
         endif
         eblk = blks(k)

****     Now get site order in this block
         do j = 1, num_site
            do i = sblk,eblk
               if( RT_sitenum(i).eq.j ) then
                  n = n + 1
                  RT_ord(n) = i
               endif
            end do
         end do
      end do

*     Debug
      if( ep.ge.debug(3) .and. ep.le. debug(4)  )
     .write(*,220) num_site, num_blk,  blks(1:num_blk)
 220  format('+++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     .       'RT_ORD ',i2,' Sites; num_blk ',i3,' Blocks ',10(i4,1x))

***   Thats all
      return
      end 

CTITLE MOVE_RT

      subroutine move_RT(s,d)

      implicit none

*     Routine to move data from index s to index d.  s must be greater than d.

      include 'trackRTObs.h'      ! Real-time data structures 

* PASSED
      integer*4 s   ! Sourece index
     .,         d   ! Destination index

      if( s.lt.d ) then
          write(*,120) s,d
 120      format('Fatal error: Trying to move data up. Source ',i4
     .,          ' Destination ',i4)
          stop 'BAD DATA ORDERING'
      end if

****  Start move
  
      RT_flags(d)  =  RT_flags(s)    
      RT_StatID(d) =  RT_StatID(s)         
      RT_satSys(d) =  RT_satSys(s)         
      RT_satNum(d) =  RT_satNum(s)         
      RT_slot(d)   =  RT_slot(s)           
      RT_GPSWeek(d)    = RT_GPSWeek(s)        
      RT_GPSweeks(d)   = RT_GPSweeks(s)       
      RT_C1(d)         = RT_C1(s)             
      RT_C2(d)         = RT_C2(s)             
      RT_P1(d)         = RT_P1(s)             
      RT_P2(d)         = RT_P2(s)             
      RT_L1(d)         = RT_L1(s)             
      RT_L2(d)         = RT_L2 (s)            
      RT_slip_cnt_L1(d)   = RT_slip_cnt_L1(s)    
      RT_slip_cnt_L2(d)   = RT_slip_cnt_L2(s)    
      RT_lock_timei_L1(d) = RT_lock_timei_L1(s)  
      RT_lock_timei_L2(d) = RT_lock_timei_L2(s)  
      RT_S1(d)            = RT_S1(s)             
      RT_S2(d)            = RT_S2(s)             
      RT_SNR1(d)          = RT_SNR1(s)           
      RT_SNR2(d)          = RT_SNR2(s)           

      RT_sitenum(d)       = RT_sitenum(s) 
      if( RT_sitenum(d).le.0 ) then
         write(*,220) s,d,  RT_StatID(s)
  220    format('ERROR: Zero station number with move from ',i4,
     .          ' to ',i4,' StatID ',a)
         stop 'BAD SITE number'
      end if      

      RT_MJD_obs(d)       = RT_MJD_obs(s)        
      RT_errflag(d)       = RT_errflag(s)        

****  Thats all
      return
      end

CTITLE SORT_BLK

      subroutine sort_blk( ep, sortind, sortvar ) 

      implicit none

*     Routine to sort data into ascending site order.  This routine
*     takes the raw OW data between sblk and eblk and sort the data
*     so that the stations are in the order of the sites. The following
*     arrays are set and re-set:
*     RT_ord: Index in it points to where the observation needs to go
*             (at the end it is reset to run from sblk-eblk sequentially
*     sslk(ns), eblk(ns) -- these arrays give the start and end of the 
*             obs number at each site.  At the end these run in sequence.

      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep    ! Current epoch number
     .,   sortind(max_rtobs)   ! sorted list 
     .,   sortvar(max_rtobs)   ! variable sorted on (for debug)

* LOCAL


!     Copy of the data arrays 
      integer*4 CP_flags(max_rtobs)  
      character*21 CP_StatID(max_rtobs)       ! Station ID
      character*1 CP_satSys(max_rtobs)        ! Satellite System ('G' or 'R')
      integer*4 CP_satNum(max_rtobs)          ! Satellite Number (PRN for GPS NAVSTAR)
      integer*4 CP_slot(max_rtobs)            ! Slot Number (for Glonass)
      integer*4 CP_GPSWeek(max_rtobs)         ! Week of GPS-Time
      real*8    CP_GPSweeks(max_rtobs)        ! Second of Week (GPS-Time)
      real*8 CP_C1(max_rtobs)                 ! CA-code pseudorange (meters)
      real*8 CP_C2(max_rtobs)                 ! CA-code pseudorange (meters)
      real*8 CP_P1(max_rtobs)                 ! P1-code pseudorange (meters)
      real*8 CP_P2(max_rtobs)                 ! P2-code pseudorange (meters)
      real*8 CP_L1(max_rtobs)                 ! L1 carrier phase (cycles)
      real*8 CP_L2 (max_rtobs)                ! L2 carrier phase (cycles)
      integer*4 CP_slip_cnt_L1(max_rtobs)     ! L1 cumulative loss of continuity indicator (negative value = undefined)
      integer*4 CP_slip_cnt_L2(max_rtobs)     ! L2 cumulative loss of continuity indicator (negative value = undefined)
      integer*4 CP_lock_timei_L1(max_rtobs)   ! L1 last lock time indicator                (negative value = undefined)
      integer*4 CP_lock_timei_L2(max_rtobs)   ! L2 last lock time indicator                (negative value = undefined)
      real*8 CP_S1(max_rtobs)                 ! L1 signal-to noise ratio
      real*8 CP_S2(max_rtobs)                 ! L2 signal-to noise ratio
      integer*4 CP_SNR1(max_rtobs)            ! L1 signal-to noise ratio (mapped to integer)
      integer*4 CP_SNR2(max_rtobs)            ! L2 signal-to noise ratio (mapped to integer)
      integer*4 CP_sitenum(max_rtobs)         ! Site number in list of sites (zero if not
                                              ! in list and hence not used)
      real*8 CP_MJD_obs(max_rtobs)            ! MJD of epopch
      integer*4 CP_errflag(max_rtobs)         ! Error flags

* LOCAL
      integer*4 i,j, k 
      integer*4 ns, ne  ! Start and stop site boundaries 
      integer*4 isblk(max_site), ieblk(max_site)  ! Initial start and stop of each
                                              ! sites data in a block

****  First move the current block of data to its copy.  To be safe, we move all 
*     data in all blocks.

****  Start move
      do j = 1, num_rtobs
        CP_flags(j)  =  RT_flags(j)    
        CP_StatID(j) =  RT_StatID(j)         
        CP_satSys(j) =  RT_satSys(j)         
        CP_satNum(j) =  RT_satNum(j)         
        CP_slot(j)   =  RT_slot(j)           
        CP_GPSWeek(j)    = RT_GPSWeek(j)        
        CP_GPSweeks(j)   = RT_GPSweeks(j)       
        CP_C1(j)         = RT_C1(j)             
        CP_C2(j)         = RT_C2(j)             
        CP_P1(j)         = RT_P1(j)             
        CP_P2(j)         = RT_P2(j)             
        CP_L1(j)         = RT_L1(j)             
        CP_L2(j)         = RT_L2 (j)            
        CP_slip_cnt_L1(j)   = RT_slip_cnt_L1(j)    
        CP_slip_cnt_L2(j)   = RT_slip_cnt_L2(j)    
        CP_lock_timei_L1(j) = RT_lock_timei_L1(j)  
        CP_lock_timei_L2(j) = RT_lock_timei_L2(j)  
        CP_S1(j)            = RT_S1(j)             
        CP_S2(j)            = RT_S2(j)             
        CP_SNR1(j)          = RT_SNR1(j)           
        CP_SNR2(j)          = RT_SNR2(j)           
        CP_sitenum(j)       = RT_sitenum(j) 
        CP_MJD_obs(j)       = RT_MJD_obs(j)        
        CP_errflag(j)       = RT_errflag(j)
      end do 

****  Now find out where the boundaies are in each block
C     do i = 1, num_site
C        ns = 10000
C        ne = 0
C        isblk(i) = 0  ! Initialize
C        ieblk(i) = 0
C        do j = sblk,eblk
C           if( RT_sitenum(j).eq.i .and. j.lt.ns ) then
C               ns = j
C           end if
C           if( RT_sitenum(j).eq.i .and. j.ge.ne ) then
C               ne = j
C           end if
C        end do
C        if( ne.gt.0 .and. ns.lt.1000) then 
C           if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
C    .      write(*,240) blk,i,(j,RT_sitenum(j),RT_satNum(j),j=ns,ne)
C240        format('RT_INT Blk ',i2,' Site ',i2,' Site/Sats: ',
C    .           100(i3,1x,i2,1x,i2.2,' |'))
C           isblk(i) = ns
C           ieblk(i) = ne
C        end if     
C     end do

****  OK now copy data back in the correct order
C     k = 0 
C     do i = 1, num_site
C        ssblk(i) = 0
C        seblk(i) = 0
C        if( isblk(i).ne.0 ) then 
C           ssblk(i) = k+1
C           do j = isblk(i), ieblk(i) 
C              k = k + 1
C              seblk(i) = k 
C              RT_flags(k)  =  CP_flags(j)    
C              RT_StatID(k) =  CP_StatID(j)         
C              RT_satSys(k) =  CP_satSys(j)         
C              RT_satNum(k) =  CP_satNum(j)         
C              RT_slot(k)   =  CP_slot(j)           
C              RT_GPSWeek(k)    = CP_GPSWeek(j)        
C              RT_GPSweeks(k)   = CP_GPSweeks(j)       
C              RT_C1(k)         = CP_C1(j)             
C              RT_C2(k)         = CP_C2(j)             
C              RT_P1(k)         = CP_P1(j)             
C              RT_P2(k)         = CP_P2(j)             
C              RT_L1(k)         = CP_L1(j)             
C              RT_L2(k)         = CP_L2 (j)            
C              RT_slip_cnt_L1(k)   = CP_slip_cnt_L1(j)    
C              RT_slip_cnt_L2(k)   = CP_slip_cnt_L2(j)    
C              RT_lock_timei_L1(k) = CP_lock_timei_L1(j)  
C              RT_lock_timei_L2(k) = CP_lock_timei_L2(j)  
C              RT_S1(k)            = CP_S1(j)             
C              RT_S2(k)            = CP_S2(j)             
C              RT_SNR1(k)          = CP_SNR1(j)           
C              RT_SNR2(k)          = CP_SNR2(j)           
C              RT_sitenum(k)       = CP_sitenum(j) 
C              RT_MJD_obs(k)       = CP_MJD_obs(j)        
C              RT_errflag(k)       = CP_errflag(j)
C           end do
C        end if
C     end do

****  OK now copy data back in the correct order
      do i = 1, num_RTobs
         k = sortind(i)
         RT_flags(i)  =  CP_flags(k)    
         RT_StatID(i) =  CP_StatID(k)         
         RT_satSys(i) =  CP_satSys(k)         
         RT_satNum(i) =  CP_satNum(k)         
         RT_slot(i)   =  CP_slot(k)           
         RT_GPSWeek(i)    = CP_GPSWeek(k)        
         RT_GPSweeks(i)   = CP_GPSweeks(k)       
         RT_C1(i)         = CP_C1(k)             
         RT_C2(i)         = CP_C2(k)             
         RT_P1(i)         = CP_P1(k)             
         RT_P2(i)         = CP_P2(k)             
         RT_L1(i)         = CP_L1(k)             
         RT_L2(i)         = CP_L2 (k)            
         RT_slip_cnt_L1(i)   = CP_slip_cnt_L1(k)    
         RT_slip_cnt_L2(i)   = CP_slip_cnt_L2(k)    
         RT_lock_timei_L1(i) = CP_lock_timei_L1(k)  
         RT_lock_timei_L2(i) = CP_lock_timei_L2(k)  
         RT_S1(i)            = CP_S1(k)             
         RT_S2(i)            = CP_S2(k)             
         RT_SNR1(i)          = CP_SNR1(k)           
         RT_SNR2(i)          = CP_SNR2(k)           
         RT_sitenum(i)       = CP_sitenum(k) 
         RT_MJD_obs(i)       = CP_MJD_obs(k)        
         RT_errflag(i)       = CP_errflag(k)
      end do

****  OK, now reset RT_ord
      do j = 1, num_RTobs
         RT_ord(j) = j
      end do

      if( ep.ge.debug(3) .and. ep.le.debug(4) )  then
         do i = 1, num_RTobs
            write(*,340) i, RT_sitenum(i),site_names(RT_sitenum(i)),
     .           RT_satNum(i), sortind(i), sortvar(i),
     .           RT_MJD_obs(i)
 340        format('SORT ',i4,' Site ',i2,1x,a,' G ',I2.2,
     .             ' Indices ',i4,1x,i8,' MJD ',F13.6)
         end do
      end if

***** Thats all
      return
      end

      subroutine get_blk( blk, ep )


      implicit none

*     Routine to get the boundaries between each site in the current
*     block.
*     RT_ord: Index in it points to where the observation needs to go
*             (at the end it is reset to run from sblk-eblk sequentially
*     sslk(ns), eblk(ns) -- these arrays give the start and end of the 
*             obs number at each site.  At the end these run in sequence.

      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep    ! Current epoch number
     .,   blk         ! Current block number

* LOCAL
      integer*4 i,j, k 
      integer*4 ns, ne  ! Start and stop site boundaries 
      integer*4 isblk(max_site), ieblk(max_site)  ! Initial start and stop of each
                                              ! sites data in a block

****  Now find out where the boundaies are in each block
      do i = 1, num_site
         ns = 10000
         ne = 0
         ssblk(i) = 0  ! Initialize
         seblk(i) = 0
         do j = sblk,eblk
            if( RT_sitenum(j).eq.i .and. j.lt.ns ) then
                ns = j
            end if
            if( RT_sitenum(j).eq.i .and. j.ge.ne ) then
                ne = j
            end if
         end do
         if( ne.gt.0 .and. ns.lt.1000) then 
            if( ep.ge.debug(3) .and. ep.le.debug(4) ) 
     .      write(*,240) blk,i,(j,RT_sitenum(j),RT_satNum(j),j=ns,ne)
 240        format('RT_INT Blk ',i2,' Site ',i2,' Site/Sats: ',
     .           100(i3,1x,i2,1x,i2.2,' |'))
            ssblk(i) = ns
            seblk(i) = ne
         end if     
      end do

****  Thats all
      return
      end 

      

