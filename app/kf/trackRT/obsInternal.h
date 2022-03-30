!  Include definition from BNC/RTCM/GPSDecoder.h
! class t_obsInternal {
!  public:

!   t_obsInternal() : 
!     flags(0),
!     satSys(' '),
!     satNum(0),
!     slot(0),
!     GPSWeek(0),
!     GPSWeeks(0.0),
!     C1(0.0),
!     C2(0.0),
!     P1(0.0),
!     P2(0.0),
!     L1(0.0),
!     L2(0.0),
!     slip_cnt_L1(-1),
!     slip_cnt_L2(-1),
!     lock_timei_L1(-1),
!     lock_timei_L2(-1),
!     S1(0.0),
!     S2(0.0),
!     SNR1(0),
!     SNR2(0) {
!     StatID[0] = '\x0';
!   }
!   int    flags;
!   char   StatID[20+1];  // Station ID
!   char   satSys;        // Satellite System ('G' or 'R')
!   int    satNum;        // Satellite Number (PRN for GPS NAVSTAR)
!   int    slot;          // Slot Number (for Glonass)
!   int    GPSWeek;       // Week of GPS-Time
!   double GPSWeeks;      // Second of Week (GPS-Time)
!   double C1;            // CA-code pseudorange (meters)
!   double C2;            // CA-code pseudorange (meters)
!   double P1;            // P1-code pseudorange (meters)
!   double P2;            // P2-code pseudorange (meters)
!   double L1;            // L1 carrier phase (cycles)
!   double L2;            // L2 carrier phase (cycles)
!   int    slip_cnt_L1;   // L1 cumulative loss of continuity indicator (negative value = undefined)
!   int    slip_cnt_L2;   // L2 cumulative loss of continuity indicator (negative value = undefined)
!   int    lock_timei_L1; // L1 last lock time indicator                (negative value = undefined)
!   int    lock_timei_L2; // L2 last lock time indicator                (negative value = undefined)
!   double S1;            // L1 signal-to noise ratio
!   double S2;            // L2 signal-to noise ratio
!   int    SNR1;          // L1 signal-to noise ratio (mapped to integer)
!   int    SNR2;          // L2 signal-to noise ratio (mapped to integer)
! };

! F90 structure definition of the t_obsInternal object
!
      type :: obsInternal 
          integer*4 flags 
          character*21 StatID       ! Station ID
          character*1 satSys        ! Satellite System ('G' or 'R')
          integer*4 satNum          ! Satellite Number (PRN for GPS NAVSTAR)
          integer*4 slot            ! Slot Number (for Glonass)
          integer*4 GPSWeek         ! Week of GPS-Time
          real*8    GPSweeks        ! Second of Week (GPS-Time)
          real*8 C1                 ! CA-code pseudorange (meters)
          real*8 C2                 ! CA-code pseudorange (meters)
          real*8 P1                 ! P1-code pseudorange (meters)
          real*8 P2                 ! P2-code pseudorange (meters)
          real*8 L1                 ! L1 carrier phase (cycles)
          real*8 L2                 ! L2 carrier phase (cycles)
          integer*4 slip_cnt_L1     ! L1 cumulative loss of continuity indicator (negative value = undefined)
          integer*4 slip_cnt_L2     ! L2 cumulative loss of continuity indicator (negative value = undefined)
          integer*4 lock_timei_L1   ! L1 last lock time indicator                (negative value = undefined)
          integer*4 lock_timei_L2   ! L2 last lock time indicator                (negative value = undefined)
          real*8 S1                 ! L1 signal-to noise ratio
          real*8 S2                 ! L2 signal-to noise ratio
          integer*4 SNR1            ! L1 signal-to noise ratio (mapped to integer)
          integer*4 SNR2            ! L2 signal-to noise ratio (mapped to integer)
      end type

