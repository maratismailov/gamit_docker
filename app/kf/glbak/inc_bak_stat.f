CTITLE INC_BAK_STATS

      subroutine inc_bak_stats ( cov_obs, sol_obs )

      implicit none  
 
*     Routine increment the statistics for each class of parameter
*     which we have based on the parameter type.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   i       - Loop counter
*   indx    - Index for decoding code
*   type    - Type of parameter from code
 
      integer*4 i, indx, type
 
*   cov_obs( cnum_parn, cnum_parn)  - Data postfit covariance matrx
*   sol_obs( cnum_parn )            - Postfit residuals
*   wgh                             - Weight for current observation
 
      real*8 cov_obs( cnum_parn, cnum_parn), sol_obs( cnum_parn ), wgh
 
*   first_call  - Indicates that we have not yet called
 
      logical first_call
 
 
      save first_call
 
      data first_call / .true. /
 
 
***** Initialize values if first call
      if( first_call ) then
          do i = 1, max_chi_types
              sum_type_wgh(i) = 0
              sum_type_res(i) = 0
              sum_type_num(i) = 0
          end do
          first_call = .false.
      end if
 
***** Now loop over the parameter codes
 
      do i = 1, cnum_used
 
          call decode_code( gpar_codes(i), type, indx )
 
*                                                ! Last three are reserved
          if( type-6.le.max_chi_types-3 ) then
*                                                ! for NEU
              wgh = 1.d0/cov_obs(i,i)
              sum_type_wgh(type-6) = sum_type_wgh(type-6) + wgh
              sum_type_res(type-6) = sum_type_res(type-6) +
     .                               sol_obs(i)**2*wgh
              sum_type_num(type-6) = sum_type_num(type-6) + 1
          end if
      end do
 
***** Thats all
      return
      end

