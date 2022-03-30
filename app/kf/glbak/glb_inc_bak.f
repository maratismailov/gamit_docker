CTITLE glb_inc_bak
 
      subroutine glb_inc_bak( cov_obs, sol_obs )

      implicit none  
 
*     Routine to increment the prefit residual Chi**2 for the
*     global solution during the back solution.
*
*     The chi**2 is computed rigously using:
*        2          T
*     Chi  = sol_obs  cov_obs sol_obs
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
      include '../includes/globk_markov.h'
      include '../includes/glb_hdr_def.h'
 
*   i,j     - Loop counters
 
 
      integer*4 i
 
*   cov_obs(cnum_parn,cnum_parn)    - Covariance matrix
*                   - of data
*   sol_obs(cnum_parn)              - Prefit residuals after
*                   - update for current estimate of parameters
*   loc_sum         - Sum of prefits for this experiment
*   temp_sum        - Temporary summation variable
 
      real*8 cov_obs(cnum_parn,cnum_parn), sol_obs(cnum_parn), loc_sum,
     .    temp_sum
     
*   outline   - Line to be output by report_stat
*   progname  - Program name

      character*256 outline, progname

*   lenprog  - Length of program name
*   rcpar    - Returns runstring

      integer*4 lenprog, rcpar    
     
      logical first_call
 
      save first_call, progname
 
 
      data  first_call  / .true. /

***** If first call get the program name for reporting Chi**2
*     increments.
 
      if( first_call ) then
          first_call = .false.
          lenprog = rcpar(0, progname )
      end if
      
****  Now compute increment to chi**2
 
      loc_sum = 0.0d0
      do i = 1, cnum_used
          call DWDOT( temp_sum, cov_obs(1,i),1, sol_obs,1,
     .                cnum_used )
          loc_sum = loc_sum + sol_obs(i)*temp_sum
      end do
 
* MOD TAH 980519: Save the chi**2 increment so that it can be
*     written to srt_file.
      dchi_save = loc_sum/cnum_used

      write( log_unit,100 ) glb_inp_file, cnum_used, loc_sum/cnum_used
 100  format(' For ',a40,' Chi**2 for ',i4,
     .       ' parameters is ',f10.3)
      if( log_unit.ne.6 ) then 
           write( *,100 ) glb_inp_file, cnum_used, loc_sum/cnum_used
      end if

      if( abs(loc_sum/cnum_used).lt.1000.d0 ) then
          write(outline,110) cnum_used, loc_sum/cnum_used
 110      format('Chi**2 Increment for ',i5,' dof ',F7.3)
      else
          write(outline,120) cnum_used, loc_sum/cnum_used
 120      format('Chi**2 Increment for ',i5,' dof ',d11.3)
          call report_stat('warning',progname,'glbak',glb_inp_file,
     .                 outline,0)
      end if
      call report_stat('status',progname,'glbak',glb_inp_file,
     .                 outline,0)
 
*     Save the local increment to chi**2
      sum_loc_chi = loc_sum
 
****  That's all
      return
      end
 
