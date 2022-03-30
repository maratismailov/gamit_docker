CTITLE APPLY_EDITS
 
      subroutine apply_edits

      implicit none 
 
 
*     Routine to apply the editing conditions to the data as they are
*     read from the file.  This routine sets the bits in DATA_FLAG to
*     reflect any "problems" with the data.   This value of DATA_FLAG
*     is later masked with DATA_MASK to see if the observation will be
*     used in the solution.
*     Some of the bits in DATA_FLAG are either set when the data is read
*     in or during the application of the contributions and medium
*     calibrations.
*                                 10:24 AM  FRI., 13  MAR., 1987
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   i,j,k       - Loop counters
 
      integer*4 i
 
*   kbit        - Bit checking function
 
 
      logical kbit
 
****  See if we do not meet FRNGE quality

* MOD TAH 901227 See if quality code is XX if it is then it is not
*     available so we should ignore.
      if( FRNGE_qual(1:2).ne.'XX' ) then
*                                                   ! Check for letter qual
*                                                   ! code also.
*                                                   ! Obs. no good
          if( FRNGE_qual(2:2).lt.min_FRNGE_quality .or.
     .        FRNGE_qual(2:2).ge.'A'   ) then
              call sbit(data_flag,1,1)
          else
              call sbit(data_flag,1,0)
          end if
      else
          call sbit(data_flag,1,0 )
      end if
 
***** Check the number of channels for this observation
*                                          ! Usual number of channels have
      if( kbit(data_notes,9) ) then
*                                          ! been set, so check
          if( num_channels.ne. usual_num_channels ) then
*                                          ! Set SOLVK interactive active
              call sbit(data_flag, 7, 1)
*                                          ! for the moment since we are
*                                          ! not using
          end if
      end if
 
***** See if we meet elevation angle limits
 
      do i = 1,2
*                                                       ! Below cutt_off
          if( elev(i,1).lt.elev_cutoff(site(i)) ) then
              call sbit( data_flag,14,1)
          end if
      end do
 
***** See if source deleted
      if( kbit(down_source,source) ) then
          call sbit(data_flag,12,1)
      end if
 
***** See if inside the down time region for this baseline
*     Loop over the number of down time entries
      do i = 1, down_num
 
*         Are we inside the time interval for this period
          if( epoch.ge.(down_times(1,i)+start_epoch) .and.
*                                                            ! Inside
     .        epoch.le.(down_times(2,i)+start_epoch) ) then
*                                                 ! time interval
 
*             See if this baseline is the correct one for this time
              IF( (site(1).eq.down_sites(1,i) .or.
     .             down_sites(1,i).eq.999999 )         .AND.
     .            (site(2).eq.down_sites(2,i) .or.
*                                                           ! donot use
     .             down_sites(2,i).eq.999999 )         )  THEN
*                                                      ! this baseline
 
                  call sbit( data_flag,13,1)
              END IF
 
*             Check for reverse baseline
              IF( (site(2).eq.down_sites(1,i) .or.
     .             down_sites(1,i).eq.999999 )         .AND.
     .            (site(1).eq.down_sites(2,i) .or.
*                                                           ! donot use
     .             down_sites(2,i).eq.999999 )         )  THEN
*                                                      ! this baseline
 
                  call sbit( data_flag,13,1)
              END IF
 
*                     ! in time region
          end if
*                     ! Looping over the number of down times
      end do
 
***** Thats all
      return
      end
 
