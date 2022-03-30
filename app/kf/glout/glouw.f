CTITLE GLOUW
 
      subroutine glouw

      implicit none 
 
*     This segment writes out the final solution to the output
*     device with the options passed through the runstring
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   ema_data(max_vma_space) - Allocation of memory for
*               - covariance matrix and solution vector
*   i           - Loop counter
*   iout        - Output LU
*   options     - options passed through runstring
 
      integer*4 ema_data(1), iout, options, i
 
*   kbit        - Function for testing if bit is on
 
      logical kbit

*   use_sites   - Dummy array needed to indicate sites used in
*                 origin for call to write_glb_sum.  Here we just
*                 set to zero,
 
      integer*4 use_site(max_glb_site_wrds) 

      common / progcon / iout, options
 
 
      common / globc_ema / ema_data
 
***** Start,
 
      call write_glb_header( iout, options,
     .                       ema_data(icov_parm), ema_data(isol_parm) )
  
*     Now the parameter estimates

      IF( KBIT(OPTIONS,5) ) THEN
*        clear the dummy use_site array
         do i = 1, max_glb_site_wrds
            use_site(i) = 0
         end do
         call write_glb_sum( iout, options, use_site, 
     .                       ema_data(icov_parm), ema_data(isol_parm),
     .                       list_file )

      END IF

*     See if positions summary to output
      IF( KBIT(OPTIONS,14) ) THEN
*        clear the dummy use_site array
         do i = 1, max_glb_site_wrds
            use_site(i) = 0
         end do
         call write_glb_pos( iout, options, use_site, 
     .                       ema_data(icov_parm), ema_data(isol_parm),
     .                       list_file )

      END IF

****  Write log estimates
      call write_log_sum(iout, options, use_site, 
     .                   ema_data(icov_parm), ema_data(isol_parm) )

*     Write out the standard output solution

      call write_glb_params( iout, options, ema_data(icov_parm),
     .                                      ema_data(isol_parm) )

*     Now write the baseline value information
*     Temporary removal
      if( kbit(options,2) ) then
          call write_glb_basel( iout, options, 1, ema_data(icov_parm),
     .                                        ema_data(isol_parm) )
          call write_glb_bcomp( iout, options, 1, ema_data(icov_parm),
     .                                        ema_data(isol_parm) )
      end if
 
*     Write out baseline rate information
      if( kbit(options,3) ) then
C         call write_glb_basel( iout, options, 2, ema_data(icov_parm),
C    .                                            ema_data(isol_parm) )
          call write_glb_bcomp( iout, options, 2, ema_data(icov_parm),
     .                                            ema_data(isol_parm) )
      end if
 
*     Now write out correlations
      call write_glb_corel( iout, options,  ema_data(icov_parm),
     .                                      ema_data(isol_parm) )

*     Write a blank line at the end
      write(iout,'(1x)')
 
***** Thats all
      return
      end
 
