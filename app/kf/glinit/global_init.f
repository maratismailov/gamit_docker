CTITLE GLOBAL_INIT
 
      subroutine global_init

      implicit none 
 
 
*     This routine initalizes some of the variable used by the global
*     Kalman filter, and clears the ema area.
 
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
 
*   EmaSize   - Function to return the size of an ema area.
 
 
      integer*4 EmaSize, i,j
 
***** Clear the ema area first
 
      call clear_ema( glb_ema, EmaSize(glb_ema, last_glb_ema) )
 
***** Set up the seasonal frequencies and dampings
 
*                             ! X Wobble annual (cycles/year)
      sw_freq(1) = 1.0
*                             ! Y Wobble annual (cycles/year)
      sw_freq(2) = 1.0
 
*                             ! X damping time ~ 1000 days (years)
      sw_damp(1) = 0.365
*                             ! X damping time ~ 1000 days (years)
      sw_damp(2) = 0.365
 
*     UT1 values
 
*                             ! annual (cycles/year)
      su_freq(1) = 1.0
*                             ! semi-annual (cycles/year)
      su_freq(2) = 2.0
 
*                             ! annual damping time ~ 1000 days (years)
      su_damp(1) = 0.365
*                             ! semiannual damping time ~ 1000 days (years)
      su_damp(2) = 0.365
 
***** Free mode data
 
*                             ! Wobble period in days
      wob_period = 432
*                             ! Q factor for wobble
      wob_Q      = 100
 
*                             ! FCN period in days
      nut_period = 435
*                             ! Q factor for FCN
      nut_Q      = 100
 
***** Some logicals for the solution
 
      glb_bak_soln    = .false.
      glb_glb_tides   = .false.
      compute_glb_res = .false.
      clear_bo_site   = .false.
      no_direct_copy  = .false.
      rad_reset       = .true.
      old_irw         = .false.
      uni_wght        = .false.
      del_scratch     = .false.
 
***** Initialize some of the file names (so they wont be used)

      glb_sol_file = glb_sol_default  
      pmu_table_file = ' '
      nut_table_file = ' '
      apr_table_file = ' '
      glb_out_file   = ' '
      glb_bak_file   = ' '
      glr_cmd_file   = ' '
      glb_org_file   = ' '
      nut_inp_file   = ' '
      plan_inp_file  = ' '
      pmu_inp_file   = ' '
      sd_inp_file    = ' '
      glr_sol_file   = ' '

      if( .not.make_svs_file ) glb_svs_file   = ' '
      svs_mar_file   = ' '
      gdescription   = 'Globk Analysis'

      num_eq = 0
      num_renames = 0
      eq_reset = .false.
      num_apr_files = 0
      type_svs_elem = 0
      num_mul_pmu = 0

*     Set epoch tolerances
      tol_mul_pmu = 0.01d0
      tol_svs_eph = 0.10d0

*     Initialize the satellite antenna apriori
      do i = 1, max_glb_svs
         do j = 1,3
           apr_val_ant(j,i) = -999.d0
         end do
      end do

*     Clear the times_used arrays
      do i = 1, max_glb_sites
          times_used(i) = 0
      end do

*     Set the max_chi**2 increment allowed
      max_chi_inc = 100.0

*     Set the max prefit residual (meters) 
      max_prefit_diff = 10000.0 

*     Set max rotation allowed (make default large) (mas)
      max_eop_rot = 10000.0

*     Clear the number of non-secular terms and downweights
      num_nonsec = 0
      num_est_nons = 0
      num_ss = 0

* MOD TAH 120714: Set appload_mod to show not set by user yet
      appload_mod = 0

* NOD TAH 200220: See IERS10 as default mean-pole model (may
*     change after ITRF2020 is apdoted.
      mean_pole_def = 'IERS10'

***** Thats about all
      return
      end
 
