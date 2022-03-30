CTITLE GLB_BAK_WRIT1
 
      subroutine glb_bak_writ1( cov_sav_parm, sol_sav_parm,
     .                          cov_obs, sol_obs, obs_corr )

      implicit none  
 
*     Routine to write out the back solution results
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glorg_common.h'
      include '../includes/glb_hdr_def.h'
 
*   i       - Loop counter
*   iel     - row number for markov element
*   ierr    - IOSTAT error
*   iout    - Output unit number
*   loc_site_parn(3,max_glb_sites)  - Parameter numbers to be
*           - used when sites output
*   loc_source_parn(2,max_glb_sources) - Parameter numbers to
*           - be used when source output
*   loc_svs_parn(max_svs_elem, max_glb_svs) - Parameter numbers to be used
*             for satellite orbits
*   loc_parn(max_glb_sites) - General local parameter number
*   loc_ind_mar(3,max_glb_sites) - Local indication of sites to
*             to be printed to bak file.
 
      integer*4 ierr, iout, loc_site_parn(3,max_glb_sites),
     .    loc_source_parn(2,max_glb_sources), loc_parn(max_glb_sites),
     .    loc_svs_parn(max_svs_elem,max_glb_svs), i, j, 
     .    loc_ind_mar(3,max_glb_sites)
 
*   cov_sav_parm(num_glb_parn,num_glb_parn) - Covariance matrix
*           - of final results
*   cov_obs(cnum_parn,cnum_parn)    - Contains postfit observation
*                               - covariance matrix
*   obs_corr(cnum_parn)         - Corrections to the observed values
*                               - to remove translations and
*                               - rotations.
*   sol_sav_parm(num_glb_parn)  - solution vector
*   sol_obs(cnum_parn)          - post fit residuals
 
      real*8 cov_sav_parm(num_glb_parn,num_glb_parn),
     .    cov_obs(cnum_parn,cnum_parn), obs_corr(cnum_parn),
     .    sol_sav_parm(num_glb_parn), sol_obs(cnum_parn) 

      logical kbit

      character*14 st_param_label(6)
 
      data st_param_label / 'X coordinate  '
     .,                    'Y coordinate  '
     .,                    'Z coordinate  '
     .,                    'X-axis        '
     .,                    'Y-axis        '
     .,                    'Z-axis        ' /
  
 
***** Write out the header first
      iout = 202
 
      call write_bak_header( iout, 0 )
 
***** Check the covariance matrix for negative diagonal
 
      call check_covar ( cov_sav_parm, num_glb_parn, num_glb_parn )

* MOD TAH 190526: See if we will write position summary
      IF( KBIT(bak_opts,14) ) THEN
          call write_glb_pos( iout, 0, use_pos , 
     .          cov_sav_parm, sol_sav_parm, glb_inp_file )
      END IF
      
 
***** Write all of the markov elements
 
      write(iout,125)
  125 format(/' MARKOV PARAMETER ESTIMATES')
      write(iout,130)
 
  130 format('  #',t10,'PARAMETER',t44,'Estimate',t60,'Adjustment',
     .                t72,'Sigma')
 
***** Start with sites, First check if markov ! Positions

*     Allow the markov only output to be over-riden by the user
*     selected sites to be printed.
      call set_mar_prt(loc_ind_mar) 
 
      call check_est( parn_site(1,1,1), 6, 3, gnum_sites,
     .                loc_site_parn, 3, loc_ind_mar, 1, 3*gnum_sites  )

      call keep_loc( loc_site_parn, 3, 3, ltog_sites, cnum_sites,
     .               gnum_sites)
 
C     call write_loc_par( iout, bak_opts, 7, loc_site_parn, 3, -1,
      call write_loc_par( iout,        0, 7, loc_site_parn, 3, -1,
     .        gnum_sites, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
 
*     Now site rates
      call check_est( parn_site(1,2,1), 6, 3, gnum_sites,
     .                loc_site_parn, 3, ind_mar, 1, num_glb_mar )

      call keep_loc( loc_site_parn, 3, 3, ltog_sites, cnum_sites,
     .               gnum_sites)
      call write_loc_par( iout,       0, 28, loc_site_parn, 3, -1,
     .        gnum_sites, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)

****  Now do the translation parameters
      do i = 1, 3
         if( mar_tran(i,1).ne.0 .or. mar_tran(i,3).ne.0 ) then 
            call out_glbl(iout, parn_tran(i,1), 0.0d0,
     .           sol_sav_parm, cov_sav_parm, num_glb_parn, 1.d0,
     .           'TRANSLTN', st_param_label(i), '(m) ', 0)
         end if
      end do
      do i = 1,3
         if( mar_tran(i,2).ne.0 ) then 
            call out_glbl(iout, parn_tran(i,2), 0.0d0,
     .           sol_sav_parm, cov_sav_parm, num_glb_parn, 1.d0,
     .           'TRANRATE', st_param_label(i), '(m/yr) ', 0)
         end if
      end do 

****  Now do the rotation parameters
      do i = 1, 3
         if( mar_rot(i,1).ne.0 .or. mar_rot(i,3).ne.0 ) then 
            call out_glbl(iout, parn_rot(i,1), 0.0d0,
     .           sol_sav_parm, cov_sav_parm, num_glb_parn, 1.d0,
     .           'ROTATION', st_param_label(i+3), '(mas)', 0)
         end if
      end do
      do i = 1,3
         if( mar_rot(i,2).ne.0 ) then 
            call out_glbl(iout, parn_rot(i,2), 0.0d0,
     .           sol_sav_parm, cov_sav_parm, num_glb_parn, 1.d0,
     .           'ROT_RATE', st_param_label(i+3), '(ms/y)', 0)
         end if
      end do 

****  Now do the scale parameters
* MOD TAH 130716: Test for white noise scale 
      if( mar_scale(1).ne.0 .or. mar_scale(3).ne.0 ) then
          call out_glbl(iout, parn_scale(1), 0.0d0,
     .         sol_sav_parm, cov_sav_parm, num_glb_parn, 1.d0,
     .           'SCALE   ', '(ppb) ',' ', 0)
      end if

      if( mar_scale(2).ne.0 ) then
           call out_glbl(iout, parn_scale(2), 0.0d0,
     .         sol_sav_parm, cov_sav_parm, num_glb_parn, 1.d0,
     .           'SCALRATE', '(ppb/yr) ',' ', 0)
      end if

* MOD TAH 040703: Check for atmospheric delay
      call check_est( parn_atm, 1, 1, gnum_sites,
     .                loc_parn, 1, ind_mar, 1, num_glb_mar )
      call keep_loc( loc_site_parn, 3, 3, ltog_sites, cnum_sites,
     .               gnum_sites)
 
      call write_loc_par( iout,       0, 61, loc_parn, 1, -1,
     .        gnum_sites, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
       

***** Now do axis offset
 
      call check_est( parn_axo(1,1), 2, 2, gnum_sites,
     .                loc_parn, 1, ind_mar, 1, num_glb_mar )
      call keep_loc( loc_site_parn, 3, 3, ltog_sites, cnum_sites,
     .               gnum_sites)
 
      call write_loc_par( iout,       0, 10, loc_parn, 1, -1,
     .        gnum_sites, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
 
 
***** Now do source positions
 
      call check_est( parn_source(1,1,1), 4, 2, gnum_sources,
     .                loc_source_parn, 2, ind_mar,1, num_glb_mar )
      call keep_loc( loc_source_parn, 2, 2, ltog_sources, cnum_sources,
     .               gnum_sources)
 
      call write_loc_par( iout,       0, 11, loc_source_parn, 2, -1,
     .        gnum_sources, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
 
***** Now do source position rates
 
      call check_est( parn_source(1,2,1), 4,2, gnum_sources,
     .                loc_source_parn, 2, ind_mar, 1,num_glb_mar )
      call keep_loc( loc_source_parn, 2, 2, ltog_sources, cnum_sources,
     .               gnum_sources)
 
      call write_loc_par( iout,       0, 31, loc_source_parn, 2, -1,
     .        gnum_sources, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
 
***** Now do satellite ephemerides position rates.  Here we output the
*     Satellite parameters even if they are not markov since they can 
*     change when ephemeris epoch change.

      do i = 1, gnum_svs
         do j = 1,max_svs_elem
            loc_svs_parn(j,i) = parn_svs(j,i)
         end do
      end do 
      call keep_loc( loc_svs_parn, max_svs_elem, max_svs_elem, 
     .               ltog_svs, cnum_svs, gnum_svs)

      call write_loc_par( iout,       0, 51, loc_svs_parn, 
     .        max_svs_elem, -1,
     .        gnum_svs, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)

 
***** Wobble and UT1
 
      call check_est( parn_wob, 1, 1, 8,
     .                loc_parn, 1, ind_mar, 1,num_glb_mar )
 
      call write_loc_par( iout,       0, 13, loc_parn, 1, -1,
     .        8, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
 
***** UT1
 
      call check_est( parn_ut1,1, 1, 6,
     .                loc_parn, 1, ind_mar, 1,num_glb_mar )
 
      call write_loc_par( iout,       0, 14, loc_parn, 1, -1,
     .        6, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)

****  Write out the correlations
      call write_pmu_corr( iout, cov_sav_parm, num_glb_parn,
     .                     parn_wob(1), parn_wob(2), parn_ut1(1))

****  Check multi-day PMU
      if( kbit(mul_pmu_opt,32) ) then
          call write_loc_par( iout,       0, 56, loc_parn, 1, -1,
     .        6, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
          call write_loc_par( iout,       0, 57, loc_parn, 1, -1,
     .        6, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
          call write_loc_par( iout,       0, 58, loc_parn, 1, -1,
     .        6, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
      end if
      call write_iers(iout,cov_sav_parm, sol_sav_parm)
 
***** Nutation angles
 
      call check_est( parn_nut_ang, 1,1, 8,
     .                loc_parn, 1, ind_mar, 1, num_glb_mar )
 
      call write_loc_par( iout,       0, 15, loc_parn, 1, -1,
     .        8, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
 
 
***** NOW DO THE STATE TRANSISSION PARAMETERS
 
      write(iout,175, iostat=ierr)
  175 format(/' TRANSISSION PARAMETER ESTIMATES')
      write( iout, 130)
 
***** Start with sites, First check if markov ! Positions
 
C     call check_est( parn_site(1,1,1), 6, 3,  gnum_sites,
C    .                loc_site_parn, 3, trans_row, 3,num_trans_row )
 
C     call write_loc_par( iout, bak_opts, 7, loc_site_parn, 3, -1,
C     call write_loc_par( iout,        0, 7, loc_site_parn, 3, -1,
C    .        gnum_sites, cov_sav_parm, sol_sav_parm, num_glb_parn,
C    .        obs_corr, 0)
 
***** Now do axis offset
 
      call check_est( parn_axo(1,1), 2, 1, gnum_sites,
     .                loc_parn, 1, trans_row, 3,num_trans_row )
 
      call write_loc_par( iout,       0, 10, loc_parn, 1, -1,
     .        gnum_sites, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
 
 
***** Now do source positions
 
      call check_est( parn_source(1,1,1), 4, 2,  gnum_sources,
     .                loc_source_parn, 2, trans_row,3, num_trans_row )
 
      call write_loc_par( iout,       0, 11, loc_source_parn, 2, -1,
     .        gnum_sources, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
 
 
***** Wobble and UT1
 
      call check_est( parn_wob, 1, 1,  8,
     .                loc_parn, 1, trans_row, 3,num_trans_row )
 
      call write_loc_par( iout,       0, 13, loc_parn, 1, -1,
     .        8, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
 
***** UT1
 
      call check_est( parn_ut1, 1, 1, 6,
     .                loc_parn, 1, trans_row, 3, num_trans_row )
 
      call write_loc_par( iout,       0, 14, loc_parn, 1, -1,
     .        6, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
 
***** Nutation angles
 
      call check_est( parn_nut_ang, 1, 1, 8,
     .                loc_parn, 1, trans_row, 3, num_trans_row )
 
      call write_loc_par( iout,       0, 15, loc_parn, 1, -1,
     .        8, cov_sav_parm, sol_sav_parm, num_glb_parn,
     .        obs_corr, 0)
 
***** Thats all
      return
      end
 
CTITLE GLB_BAK_WRIT2
 
      subroutine glb_bak_writ2( cov_sav_parm, sol_sav_parm,
     .                          cov_obs, sol_obs, obs_corr )
 
      implicit none 
 
*     Routine to write out the back solution results
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   i       - Loop counter
*   iel     - row number for markov element
*   ierr    - IOSTAT error
*   iout    - Output unit number
*   loc_site_parn(3,max_glb_sites)  - Parameter numbers to be
*           - used when sites output
*   loc_source_parn(2,max_glb_sources) - Parameter numbers to
*           - be used when source output
*   loc_parn(max_glb_sites) - General local parameter number
 
 
      integer*4 i, iout, loc_site_parn(3,max_glb_sites),
     .    loc_source_parn(2,max_glb_sources), loc_parn(max_glb_sites)
 
*   cov_sav_parm(num_glb_parn,num_glb_parn) - Covariance matrix
*           - of final results
*   cov_obs(cnum_parn,cnum_parn)    - Contains postfit observation
*                               - covariance matrix
*   obs_corr(cnum_parn)         - Corrections to the observed values
*                               - to remove translations and
*                               - rotations.
*   sol_sav_parm(num_glb_parn)  - solution vector
*   sol_obs(cnum_parn)          - post fit residuals
 
      real*8 cov_sav_parm(num_glb_parn,num_glb_parn),
     .    cov_obs(cnum_parn,cnum_parn), obs_corr(cnum_parn),
     .    sol_sav_parm(num_glb_parn), sol_obs(cnum_parn)
 
 
 
***** Write out the header first
      iout = 202
 
***** NOW DO THE POSTFIT RESIDUALS
 
      if( compute_glb_res ) then
 
 
*         Check the covariance matrix for negative elements
          call check_covar( cov_obs, cnum_parn, cnum_used )
 
          write(iout,400)
  400     format(/' LOCAL PARAMETER POSTFIT RESIDUALS')
          write(iout,410)
  410     format('  #',t10,'PARAMETER',t46,'Observed',t60,
     .           'Residual',t72,'Sigma')
 
*****     Start with SITE positions
          call clear_loc_parn( loc_site_parn, 3*cnum_sites)
 
*                         ! Loop over XYZ
          do i = 1,3
 
              call codes_to_parn(6+i, gpar_codes, cnum_used,
     .                           loc_site_parn(i,1), 3)
          end do
 
          call write_loc_par( iout, 0, 7, loc_site_parn, 3,
     .            ltog_sites, cnum_sites, cov_obs , sol_obs,
     .            cnum_parn, obs_corr, cnum_parn)
 
 
*****     Check SITE position rates
          call clear_loc_parn( loc_site_parn, 3*cnum_sites)
 
*                         ! Loop over XYZ
          do i = 1,3
 
              call codes_to_parn(27+i, gpar_codes, cnum_used,
     .                           loc_site_parn(i,1), 3)
          end do
 
          call write_loc_par( iout, 0, 28, loc_site_parn, 3,
     .            ltog_sites, cnum_sites, cov_obs , sol_obs,
     .            cnum_parn, obs_corr, cnum_parn)
 
 
*****     Axis offset residuals
 
          call clear_loc_parn( loc_parn, cnum_sites)
 
          call codes_to_parn( 10, gpar_codes, cnum_parn,
     .                       loc_parn, 1)
 
          call write_loc_par( iout, 0,10, loc_parn, 1,
     .            ltog_sites, cnum_sites, cov_obs , sol_obs,
     .            cnum_parn, obs_corr, cnum_parn)
 
 
*****     Do the source positions
          call clear_loc_parn( loc_source_parn, 2*cnum_sources)
 
*                         ! Loop over RA and Dec
          do i = 1,2
 
              call codes_to_parn(10+i, gpar_codes, cnum_used,
     .                           loc_source_parn(i,1), 2)
          end do
 
          call write_loc_par( iout, 0,11, loc_source_parn, 2,
     .            ltog_sources, cnum_sources, cov_obs , sol_obs,
     .            cnum_parn, obs_corr, cnum_parn)
 
*****     Do the source positions rates
          call clear_loc_parn( loc_source_parn, 2*cnum_sources)
 
*                         ! Loop over RA and dec
          do i = 1,2
 
              call codes_to_parn(30+i, gpar_codes, cnum_used,
     .                           loc_source_parn(i,1), 2)
          end do
 
          call write_loc_par( iout, 0,31, loc_source_parn, 2,
     .            ltog_sources, cnum_sources, cov_obs , sol_obs,
     .            cnum_parn, obs_corr, cnum_parn)
 
 
*****     Do the wobble components
          call clear_loc_parn( loc_parn, 8)
 
          call codes_to_parn(13, gpar_codes, cnum_used,
     .                       loc_parn, 1)
 
          call write_loc_par( iout, 0,13, loc_parn, 1,
     .            -1, 8, cov_obs , sol_obs,
     .            cnum_parn, obs_corr, cnum_parn)
 
*****     Do UT1
          call clear_loc_parn( loc_parn, 6)
 
          call codes_to_parn(14, gpar_codes, cnum_used,
     .                       loc_parn, 1)
 
          call write_loc_par( iout, 0,14, loc_parn, 1,
     .            -1, 6, cov_obs , sol_obs,
     .            cnum_parn, obs_corr, cnum_parn)
 
*****     Do the nutation angles
          call clear_loc_parn( loc_parn, 8)
 
          call codes_to_parn(15, gpar_codes, cnum_used,
     .                       loc_parn, 1)
 
          call write_loc_par( iout, 0,15, loc_parn, 1,
     .            -1, 8, cov_obs , sol_obs,
     .            cnum_parn, obs_corr, cnum_parn)
 
*****     Do the Gamma residual
          call clear_loc_parn( loc_parn, 1)
 
          call codes_to_parn(27, gpar_codes, cnum_used,
     .                       loc_parn, 1)
 
          call write_loc_par( iout, 0,27, loc_parn, 1,
     .            -1, 1, cov_obs , sol_obs,
     .            cnum_parn, obs_corr, cnum_parn)
 
*                 ! residual computed
      END IF
 
***** Thats all
      return
      end
 
