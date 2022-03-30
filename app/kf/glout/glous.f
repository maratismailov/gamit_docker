CTITLE GLOUS
 
      subroutine glous

      implicit none  
 
*     This segment stablizes the solution.  At the moment it just
*     fixes the translation, and the RA origins
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   ema_data(max_vma_space) - Allocation of memory for
*               - covariance matrix and solution vector
*   i,j         - Loop counter
*   iout        - Output LU
*   options     - options passed through runstring
*   used(8)     - Temporary variable for saying that sites or sources
*               - are used.  Later will be replaced by user option
 
      integer*4 ema_data(1), i,j, iout, options, used(8)
 
*   kbit        - Bit testing function
 
      logical kbit
 
*   dummy           - Dummy value for getting from ema to main
*   origin_gain(max_glb_parn)   - Kalman filter gain for
*               - origin fixing
*   pmu_parts(3,3,max_glb_sites)    - PMU partial derivatives
*   rate_trans(6)   - Rate translation origin (set to minimize the
*                   - horizontal motions)
*   translation(6)  - XYZ origin fix.  Set to mininize the horizontal
*                   - displacements
*   ut1_main        - Main memory value of UT1-AT
 
      real*8 dummy, origin_gain(max_glb_parn),
     .    pmu_parts(3,3,max_glb_sites), rate_trans(6), translation(6),
     .    ut1_main
 
 
      common / progcon / iout, options
 
 
      common / globc_ema / ema_data
 
***** Stabilize the translation origin (position and velocity)
 
*     Get the origin first
 
      do i = 1,8
*                         ! Tell software to use everything
          used(i) = -1
      end do

      do i = 1,6
         rate_trans(i) = 0.d0
         translation(i) = 0.d0
      end do
 
*     Get the pmu partials
      ut1_main = gut1_apr(1)
      do i = 1, gnum_sites
*                         ! Loop over XYZ
          do j = 1,3
              call pmu_part(i,pmu_parts(1,j,i), apr_val_site(1,1,i),
     .            j,ut1_main)
          end do
      end do
 
      IF( KBIT(OPTIONS,6) ) THEN
 
      call Fix_system( 0.d0, apr_val_site, 2, 'EMA', translation,
     .     ema_data(isol_parm), parn_site, 2, used, gnum_sites,
     .     ema_data(isol_parm+2*num_glb_parn), -1, pmu_parts)
 
 
      call apply_trans( parn_site(1,1,1),2, ema_data(isol_parm),
     .                  gnum_sites, translation)
 
 
      call apply_rot( parn_site(1,1,1), 2, ema_data(isol_parm),
     .                pmu_parts, gnum_sites, translation(4))
 
 
*                     ! Loop over XYZ
      do i = 1, 3
 
          translation(i) = 0
          do j = 1, gnum_sites
              call get_from_ema(dummy, ema_data(isol_parm),
     .             parn_site(i,1,j) )
              translation(i) = translation(i) + dummy
          end do
*                                                       ! Should cause no
          translation(i) = -translation(i)/gnum_sites
*                                         ! change to the answres
          call glb_stabilize( parn_site(i,1,1), gnum_sites, 6,
     .         ema_data(icov_parm), num_glb_parn, num_glb_parn,
     .         ema_data(isol_parm), origin_gain, used, translation(i))
      end do
 
      do i = 1,6
          rate_trans(i) = translation(i)
      end do
      END IF
 
*     Stabilize the rates
 
      IF( KBIT(OPTIONS,7) ) THEN
      call Fix_system( 0.d0, apr_val_site, 2, 'EMA', rate_trans,
     .     ema_data(isol_parm), parn_site(1,2,1), 2, used, gnum_sites,
     .     ema_data(isol_parm+2*num_glb_parn), -1, pmu_parts)
 
*                     ! Loop over XYZ dot
      do i = 1, 3
          call glb_stabilize( parn_site(i,2,1), gnum_sites, 6,
     .         ema_data(icov_parm), num_glb_parn, num_glb_parn,
     .         ema_data(isol_parm), origin_gain, used, rate_trans(i))
      end do
 
      call apply_rot( parn_site(1,2,1), 2, ema_data(isol_parm),
     .                pmu_parts, gnum_sites, rate_trans(4))
 
 
      END IF
 
*     Now see if we should fix final orientation (with all rates estimated)
      if( (rate_trans(1) .ne.0 .and. rate_trans(2) .ne.0 .and.
     .     rate_trans(3) .ne.0) .OR.
     .    (translation(1).ne.0 .and. translation(2).ne.0 .and.
*                                          ! Fix orienation because all site
     .     translation(3).ne.0)    ) then
*                                          ! position or rates have been
*                                          ! estimated.
*                         ! x and y wobble
          do i = 1,2
*                                         ! Stabilize the value
              if( parn_wob(i).ne.0 ) then
                  call get_from_ema( rate_trans(i), ema_data(isol_parm),
     .                               parn_wob(i) )
                  rate_trans(i) = -rate_trans(i)
                  call glb_stabilize( parn_wob(i), 1, 1,
     .                ema_data(icov_parm), num_glb_parn, num_glb_parn,
     .                ema_data(isol_parm), origin_gain, -1,
     .                rate_trans(i) )
                  call get_from_ema( dummy, ema_data(isol_parm),
     .                parn_wob(i) )
                  rate_trans(i) = dummy + rate_trans(i+3)
                  call put_to_ema( ema_data(isol_parm), rate_trans(i),
     .                parn_wob(i) )
              end if
          end do
*         UT1 value
*                                     ! Stabilize the value
          if( parn_ut1(1).ne.0 ) then
              call get_from_ema( rate_trans(3), ema_data(isol_parm),
     .                           parn_ut1(1) )
              rate_trans(3) = -rate_trans(3)
              call glb_stabilize( parn_ut1(1), 1, 1,
     .            ema_data(icov_parm), num_glb_parn, num_glb_parn,
     .            ema_data(isol_parm), origin_gain, -1,
     .            rate_trans(3) )
              call get_from_ema( dummy, ema_data(isol_parm),
     .             parn_ut1(1) )
              rate_trans(3) = dummy + rate_trans(6)
              call put_to_ema( ema_data(isol_parm), rate_trans(3),
     .             parn_ut1(1) )
          end if
 
      end if
 
*     Stabilize the RA origin
      IF( KBIT(OPTIONS,11) ) then
 
          call glb_stabilize( parn_source(1,1,1), gnum_sources, 4,
     .         ema_data(icov_parm), num_glb_parn, num_glb_parn,
     .         ema_data(isol_parm), origin_gain, used, 0.d0 )
 
      END IF
 
***** Now check that we covariance matrix is OK
 
      call check_covar(ema_data(icov_parm), num_glb_parn, num_glb_parn)
 
***** Thats all
      return
      end
 
