CTITLE GLB_RESIDUALS
 
      subroutine glb_residuals( cov_sav_parm, sol_sav_parm,
     .                          a_part, part_pnt, row_copy,
     .                          cov_obs, sol_obs, temp_gain)

      implicit none  
 
*     Routine to compute the postfit residuals to the input
*     set of parameters
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
      include '../includes/glb_hdr_def.h'
 
 
*   i,j             - Loop counters
*   ic,ir           - Row and col counters
*   part_pnt(2,max_glb_deriv,cnum_parn) - Pointers to the
*                   - non-zero partial derivatives
 
      integer*4 i, ic,ir, part_pnt(2,max_glb_deriv,cnum_parn)
 
*   a_part(max_glb_deriv,cnum_parn) - Non_zero partial
*                   - derivatives
*   cov_obs(cnum_parn, cnum_parn)   - data covariance matrix
*                   - will be replaced with postfit cov. matrix
*   cov_sav_parm(num_glb_parn,num_glb_parn) - Final covariance
*                   - matrix of the parameters
 
*   row_copy(num_glb_parn)  - Copy of a row of cov_sav_parm
*   scr_real(max_glb_parn)  - Scratch area for calculations
*   sol_obs(cnum_parn)  - Prefit residuals, which will be updated
*   sol_sav_parm(num_glb_parn)  - Final corrections to parmaters
*   temp_gain(cnum_parn,num_glb_parn) - Temporary copy of gain
*                   - matrix
 
      real*8 a_part(max_glb_deriv,cnum_parn),
     .    cov_obs(cnum_parn, cnum_parn),
     .    cov_sav_parm(num_glb_parn,num_glb_parn),
     .    row_copy(num_glb_parn), scr_real(max_glb_parn),
     .    sol_obs(cnum_parn), sol_sav_parm(num_glb_parn),
     .    temp_gain(cnum_parn,num_glb_parn)
 
      common scr_real
 
***** First dot the partials into the parameters so that we get
*     the corrections to the original data.
 
      do i = 1, cnum_used
          call glb_dot_part( scr_real, sol_sav_parm,1, a_part(1,i),
     .         part_pnt(1,1,i), indx_pnt(i) )
 
          sol_obs(i) = sol_obs(i) - scr_real(1)
      end do
 
***** Now do the postfit covariance matrix
 
*           m+1 T
***** Form C   A  by rows
*           m+1
 
*                                 ! Loop over rows
      do ir = 1, num_glb_parn
*                                 ! Loop over columns
          do ic = 1, cnum_used
 
              call glb_dot_part( scr_real(ic), cov_sav_parm(1,ir),1,
     .             a_part(1,ic), part_pnt(1,1,ic), indx_pnt(ic) )
 
          end do
 
*         Now move the column into temp_gain
          call dwmov( scr_real,1, temp_gain(1,ir),1, cnum_used)
      end do
 
*                          m+1 T
*     Now finish forming AC   A   and add to Q
*                          m+1
 
      do ic = 1, cnum_used
*                                 ! Only do lower diagonal and copy
          do ir = ic, cnum_used
 
              call glb_dot_part( scr_real(ir),
     .             temp_gain(ic,1),cnum_parn, a_part(1,ir),
     .             part_pnt(1,1,ir), indx_pnt(ir) )
          end do
 
*         Subrract this contribution to cov_obs and then move values to
*         upper diagonal
 
          do ir = ic, cnum_used
              cov_obs(ir,ic) = cov_obs(ir,ic) - scr_real(ir)
          end do
 
*         Move to lower diagonal
          call DWMOV(cov_obs(ic+1,ic),1, cov_obs(ic,ic+1),cnum_parn,
     .               cnum_used-ic)
      end do
 
***** Thats all
      return
      end
 
