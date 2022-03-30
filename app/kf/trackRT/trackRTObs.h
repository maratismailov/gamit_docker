*     TrackRT common declaration
* MOD TAH 050211: Post-2.5 version of BNC which casts ASCII data does not need
*     the structure defined below.
* MOD TAH 050211: Added L5 for GPS
* MOD TAH 032917: Updated from BNC Ver 2.9 and above (added BNC_Version) 

      include 'obsInternal.h'   ! Data structure definition: Not USED trackRT V 1.10
                                ! but retained so that trackRTB will still run.

!     Set maximum number of observations all
      integer*4 max_rtobs   ! max size of obsInternal array
      integer*4 max_blks    ! maximum number of epochs that can be in one block of data
      integer*4 maxbytes    ! Largest ASCII buffer size allowed (must matck trackRTComm.h 
                            ! value
            
      parameter ( max_rtobs = 1000 )
      parameter ( max_blks  = 64   )
      parameter ( maxbytes  = 363000 ) 

! Declarations for common block

      integer*4 num_rtobs   ! Number of data records at this epoch
      integer*4 BNC_Vers    ! BNC_Version (if not set, determined from data)

!     Array elements to be saved  
      integer*4 RT_flags(max_rtobs)  
      character*21 RT_StatID(max_rtobs)       ! Station ID
      character*1 RT_satSys(max_rtobs)        ! Satellite System ('G' or 'R')
      integer*4 RT_satNum(max_rtobs)          ! Satellite Number (PRN for GPS NAVSTAR)
      integer*4 RT_slot(max_rtobs)            ! Slot Number (for Glonass)
      integer*4 RT_GPSWeek(max_rtobs)         ! Week of GPS-Time
      real*8    RT_GPSweeks(max_rtobs)        ! Second of Week (GPS-Time)
      real*8 RT_C1(max_rtobs)                 ! CA-code L1 pseudorange (meters)
      real*8 RT_C2(max_rtobs)                 ! CA-code L2 pseudorange (meters)
      real*8 RT_C5(max_rtobs)                 ! CA-code L5 pseudorange (meters)
      real*8 RT_P1(max_rtobs)                 ! P1-code pseudorange (meters)
      real*8 RT_P2(max_rtobs)                 ! P2-code pseudorange (meters)
      real*8 RT_L1(max_rtobs)                 ! L1 carrier phase (cycles)
      real*8 RT_L2 (max_rtobs)                ! L2 carrier phase (cycles)
      real*8 RT_L5 (max_rtobs)                ! L5 carrier phase (cycles)
      integer*4 RT_slip_cnt_L1(max_rtobs)     ! L1 cumulative loss of continuity indicator (negative value = undefined)
      integer*4 RT_slip_cnt_L2(max_rtobs)     ! L2 cumulative loss of continuity indicator (negative value = undefined)
      integer*4 RT_slip_cnt_L5(max_rtobs)     ! L5 cumulative loss of continuity indicator (negative value = undefined)
      integer*4 RT_lock_timei_L1(max_rtobs)   ! L1 last lock time indicator                (negative value = undefined)
      integer*4 RT_lock_timei_L2(max_rtobs)   ! L2 last lock time indicator                (negative value = undefined)
      real*8 RT_S1(max_rtobs)                 ! L1 signal-to noise ratio
      real*8 RT_S2(max_rtobs)                 ! L2 signal-to noise ratio
      real*8 RT_S5(max_rtobs)                 ! L5 signal-to noise ratio
      integer*4 RT_SNR1(max_rtobs)            ! L1 signal-to noise ratio (mapped to integer)
      integer*4 RT_SNR2(max_rtobs)            ! L2 signal-to noise ratio (mapped to integer)
      integer*4 RT_SNR5(max_rtobs)            ! L5 signal-to noise ratio (mapped to integer)

      integer*4 RT_sitenum(max_rtobs)         ! Site number in list of sites (zero if not
                                              ! in list and hence not used)

      real*8 RT_MJD_obs(max_rtobs)            ! MJD of epopch

      integer*4 RT_errflag(max_rtobs)         ! Error flags
!     Error flags (but bit)
!     1  -- Flagged bad in clock estimate
!     2  -- P1-P2 too large 
!     3  -- No SD difference to this site from earlier sites
!     4  --
!     5  -- Data below elevation cutoff
!     6  -- PRN has been excluded (either not in SP3 or user) 

!     RT_flags definitions.  (Lower order bit not know)
!     Bit   Meaning
!     1   -- P1 range copied from C1 slot
!     2   -- P2 range copied from C2 slot 

!     Section for multiple blocks
      integer*4 num_blk         ! Number of blocks of epochs
      integer*4 sblk, eblk      ! Current start and end of block
      integer*4 blks(max_blks)  ! End counter for each block 

! Common declaration

      common / trackRT_I4 / num_rtobs, BNC_Vers 
     .,        RT_flags, RT_satNum, RT_slot, RT_GPSWeek 
     .,        RT_slip_cnt_L1, RT_slip_cnt_L2, RT_slip_cnt_L5
     .,        RT_lock_timei_L1, RT_lock_timei_L2
     .,        RT_SNR1, RT_SNR2, RT_SNR5
     .,        num_blk, sblk, eblk, blks 
     .,        RT_sitenum, RT_errflag

      common / trackRT_ch / RT_StatID, RT_satSys

      common / trackRT_R8 / RT_GPSweeks, RT_C1, RT_C2,  RT_C5 
     .,        RT_P1, RT_P2, RT_L1, RT_L2, RT_L5
     .,        RT_S1, RT_S2, RT_S5, RT_MJD_obs

*******************************************************************************
*     TrackRT processing variables
      integer*4 RT_ambpnt(max_rtobs) ! Pointer to entry on bf_ent that
              ! applies to each one-way
      integer*4 RT_ord(max_rtobs) ! Pointers to have the data listed
              ! in site order. Important when double differences are formed
              ! (Also site order changes in RT data stream).

      real*8 RT_omcraw(4,max_rtobs) ! Obs-minus-computed for 
              ! L1 phase (cycles), L2 phase (cycles), P1 range (m) and
              ! P2 range (m)
      real*8 RT_azel(2,max_rtobs)     ! Azimuth/Elev angle of data (deg)

      common / trackRTpr_I4 / RT_ambpnt, RT_ord
    
      common / trackRTpr_R8 / RT_omcraw, RT_azel      
