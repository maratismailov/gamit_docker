CTITLE GLB_INC_POSTFIT
 
      subroutine glb_inc_postfit( cov_obs, sol_obs, cov_copy)

      implicit none 
  
*     Routine to increment the postfit residual Chi**2 for the
*     global solution, during the back solution.
*
*     The chi**2 is computed rigously using:
*        2          T
*     Chi  = sol_obs  cov_obs sol_obs
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
      include '../includes/glb_hdr_def.h'
 
*   i,j     - Loop counters
*   ipivot(max_glb_parn)    - Pivot rows for invert_vis
 
      integer*4 i, ipivot(max_glb_parn)
 
*   cov_obs(cnum_parn,cnum_parn)    - Covariance matrix
*                   - of data
*   cov_copy(cnum_parn,cnum_parn) - Copy of the obervation
*                   - covariance matrix. (So that we can invert to
*                   - to increment chi**sq)
 
*   scale(max_glb_parn)     - Scale array for invert_vis
 
*   sol_obs(cnum_parn)              - Prefit residuals after
*                   - update for current estimate of parameters
*   temp_sum        - Temporary summation variable
 
      real*8 cov_obs(cnum_parn,cnum_parn),
     .    cov_copy(cnum_parn,cnum_parn), scale(max_glb_parn),
     .    sol_obs(cnum_parn), temp_sum
 
 
      logical first_call
 
      save first_call
 
 
 
      common ipivot, scale
 
      data  first_call  / .true. /
 
***** If this is first call clear the summation variables
 
      if( first_call ) then
          sum_post_chi     = 0.d0
          sum_post_chi_num = 0.d0
          first_call = .false.
      end if
 
***** Save a copy of COV_OBS
 
      do i = 1, cnum_used
          call DWMOV( cov_obs(1,i),1, cov_copy(1,i),1, cnum_used)
      end do
 
***** Now scale the diagonal of cov_copy so that the matrix will not
*     singular
 
      call dwsmy(1.001d0,cov_copy(1,1), cnum_parn+1,
     .                   cov_copy(1,1), cnum_parn+1, cnum_used)
 
 
***** Invert to the copy of cov_obs
 
      call invert_vis( cov_copy, sol_obs, scale, ipivot, cnum_used,
*                                     ! No solution vector computed
     .                 cnum_parn, 0 )
 
****  Now compute increment to chi**2
 
      sum_loc_chi = 0.d0
 
      do i = 1, cnum_used
          call DWDOT( temp_sum, cov_copy(1,i),1, sol_obs,1,
     .                cnum_used )
          sum_loc_chi = sum_loc_chi + sol_obs(i)*temp_sum
      end do
 
      if( sum_loc_chi.gt.0 ) then
          sum_post_chi     = sum_post_chi + sum_loc_chi
          sum_post_chi_num = sum_post_chi_num + cnum_used
 
*         Should print a warning message if sum_loc_chi is < than zero
      end if
 
****  That's all
      return
      end
 
