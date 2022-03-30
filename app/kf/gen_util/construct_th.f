CTITLE CONSTRUCT_THEORETICAL
 
      subroutine construct_theoretical

      implicit none 
 
 
*     Routine to compute the theorectical delay for the observation
*     being processed.  The theoretical stored in the data files
*     has no contributions added to it all.  In this rotuine we
*     add those contributions which the user has requested.  We
*     also put the ambiguities into the theoretical as well.
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   i           - loop counter
*   idum        - Dummy set to make phase_cal_cont applied always
 
      integer*4 i, idum
 
*   kbit        - SOLVE function to see if bit is set.
 
      logical kbit
 
*   temp_cont(4)- Temporary contribution to fix bug in Calc Ver 6.0
*               - for general relativity contribution.
 
      real*4 temp_cont(4)
 
***** Copy data base theorecticals into user values, also copy data
*     sigmas into user space
      do i = 1,4
          theoretical(i) = db_theoretical(i)

* MOD TAH 900129: Check size of sigma to avoid overflow in calcs.
          if( db_sigma(i).gt.1.d6 ) then
              db_sigma(i) = 1.d6
          end if
          sigma(i)       = db_sigma(i)
      end do
 
*     Reset the user specific parts of the data_flag
 
      call reset_data_flag( data_flag )
 
***** Add ambiquities to the delays for group and phase by subtracting
*     from theorectical
 
      theoretical(1) = theoretical(1) - num_grp_amb*grp_amb
      theoretical(2) = theoretical(2) - num_phs_amb*phs_amb
 
***** Now add the baseline contributions (if requested by the user)
 
      call add_cont( cont_baseline,avail_baseline,1,theoretical,
     .               etd_cont    )
      call add_cont( cont_baseline,avail_baseline,2,theoretical,
     .               ptd_cont    )
      call add_cont( cont_baseline,avail_baseline,3,theoretical,
     .               ocean_cont  )
      call add_cont( cont_baseline,avail_baseline,4,theoretical,
     .               axo_cont    )
 
C     if( Calc_ver.eq.6.0 ) then
C         temp_cont(1) = -gen_rel_cont(1)
C         temp_cont(2) = -gen_rel_cont(2)
C     else
C         temp_cont(1) =  gen_rel_cont(1)
C         temp_cont(2) =  gen_rel_cont(2)
C     end if
 
C     call add_cont( cont_baseline,avail_baseline,5,theoretical,
C    .               gen_rel_cont)
 
      call get_gen_rel_cont( temp_cont )
 
      call add_cont( cont_baseline,avail_baseline,5,theoretical,
     .               temp_cont )
 
***** Now add source structure
      call add_struct( cont_structure,avail_structure, source,
     .                 theoretical,source_struc_cont)
 
***** Now add cable calibration, feed rotation, and phase calibration
 
      do i = 1,2
*                                                      ! Cable cal
          call add_site_cont( cont_site(site(i)), avail_site(site(i)),
     .        1,theoretical, cable_cont,     i, 1,1,1)
*                                                      ! Feed rotation corr.
          call add_site_cont( cont_site(site(i)), avail_site(site(i)),
     .        2,theoretical, feed_rot_cont,  i, 0,1,0)
          idum = -1
*                                                         ! Make sure applied
*                                                         ! Phase cal cont.
          call add_site_cont( cont_site(site(i)), idum,
     .        3,theoretical, phase_cal_cont,
     .        i, -(2*i-3), 0,0)
*                ! NOTE: This causes cont to be subtracted from the theo.
*                !       at both sites 1 and 2 (same sign convention as GSFC)

*                                 ! Add radial ocean load
          call add_site_cont( cont_site(site(i)), avail_site(site(i)),
     .        5,theoretical, ocean_rad    ,  i, 1,1,1)
*                                 ! Add horizontal ocean load
          call add_site_cont( cont_site(site(i)), avail_site(site(i)),
     .        6,theoretical, ocean_horz   ,  i, 1,1,1)
      end do
 
***** Now do the ionospheric delay
 
      call add_ion
 
***** Now do the site dependent terms
 
      do i = 1,2
          call add_medium (i)
      end do

***** Now see if we should add clock jumps.  Note: we do not include change
*     in rate.
      if( kbit(data_notes,10) ) then
          do i = 1,3
             theoretical(i) = theoretical(i) + clock_jump(2,1) 
     .                                       - clock_jump(1,1)
          end do
      end if
      theoretical(4) = theoretical(4) + clock_jump(2,2) 
     .                                - clock_jump(1,2)
 
***** Now include extended earth tide partials
* MOD TAH 870730: Moved to calculation of the etd partial to KALUPD to
*     decrease size of KALGN.
C     call etd_partial
 
 
***** Thats all
      return
      end
 
