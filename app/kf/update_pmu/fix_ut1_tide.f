 
      program fix_ut1_tide

*     This program will fix the UT1 tidal signals in old KalObs files
*     which were run with calc ver 7.0 to <7.3 and had the tidal UT1
*     flag set to +1.  After these experiments have been run through
*     update_pmu, they do not have the short period UT1 corrections applied
*     anymore.  This program will add these corrections.
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
 
 
*   i,j,k       - Loop counters
*   ierr        - Fmp file error flag
*   len_runstring   - Length of runstring
*   rcpar       - HP runstring reading utility
*   trimlen     - Length of string
 
 
      integer*4 i, ierr, len_runstring, rcpar, trimlen
 
*   kbit        - Bit checking function
 
      logical kbit

*   dut1       - Tide contribution to UT1 (tsec)
*   delta_tau  - Change in delay (ps)
*   delta_tau_dot - Change in rate (fs/s)

      real*8 dut1, delta_tau, delta_tau_dot
 
*   KO_name     - Name of the Kalobs file to be updated
 
 
      character*132 KO_name
 
***** Get the Kalobs file name
 
      len_runstring = rcpar(1, KO_name )
 
*                                     ! Name not given, STOP
      if( len_runstring.eq.0 ) then
          stop ' FIX_UT1_TIDE: Incorrect Runstring, add KalObs name'
      end if
 
***** Now open the KalObs file
 
      call open_KalObs( KO_name, ierr )
      call report_error('FmpOpen',ierr,'open',KO_name, 1, 
     .                  'FIX_UT1_TIDE')
 
***** Read header
 
      call RW_KalObs_header('R', ierr)
      call report_error('Fmpread',ierr,'read',KO_name, 1, 
     .                  'FIX_UT1_TIDE')

***** See if file needs fixing
      if( kbit(data_notes,11) ) then
*         Already fixed, so get out
          write(*,100) KO_name(1:trimlen(KO_name))
 100      format(' KalObs file ',a,' already has been fixed',/,
     .           ' FIX_UT1_TIDE terminating with no change to data')
          stop ' FIX_UT1_TIDE: ALready fixed'
      end if
 
***** Now loop over data file, updating as we go
 
      do i = 1, num_obs
 
          call RW_KalObs_block('R','DATA',site, ierr, i )
          call report_error('FmpRead',ierr,'read','data block',1,
*                                     ! Kill if error
     .                      'FIX_UT1_TIDE')

          call short_period_ut1( epoch, dut1 )

*         Update the value of the pmu_calc for UT1 by adding short period terms
          pmu_calc(3) = pmu_calc(3) + dut1*15.d3
*                                                 ! UT1 delay 
          delta_tau = dut1*15.d3 * pmu_part(3,1)

*****     Change theoretical value of group delay, phase delay, and
*         single-band delay.
*                                                              ! Group
          db_theoretical(1) = db_theoretical(1) + delta_tau
*                                                              ! Phase
          db_theoretical(2) = db_theoretical(2) + delta_tau
*                                                              ! SB delay
          db_theoretical(3) = db_theoretical(3) + delta_tau

*****     Calculate change to phase delay rate, assuming that the
*         contribution due to the PM-UT1 rate is zero
          delta_tau_dot = dut1*15.d3 * pmu_part(3,2)

*****     Add this contribution to the theoretical phase_delay rate
          db_theoretical(4) = db_theoretical(4) + delta_tau_dot

*         write out data record
          call RW_KalObs_block('W','DATA',site, ierr,i)
          call report_error('FmpWrite',ierr,'writ','data block',1,
*                                     ! kill of error
     .                      'FIX_UT1_TIDE')
 
      end do
 
***** Update data notes and write the header back out
*     Now update the header so that we know it is fixed
      call sbit( data_notes, 11, 1)
      ut1_tide = -1
 
      call increment_Kalver
      call RW_KalObs_header('W', ierr)
      call close_KalObs( ierr )
 
***** Thats all
      end
