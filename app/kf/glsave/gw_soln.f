CTITLE GW_SOLN
      subroutine gw_soln ( cov_parm )

      implicit none  
 
*     Routine to write the covariance matrix and solution vector
*     to the disk file.  To do this the VWRIT routine is used
*     to write directly from ema to disk (via the scratch common area)
 
      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
 
      include '../includes/globk_common.h'
 
*   ierr            - Fmp Error
*   i,j,k           - Loop counters
*   FmpClose        - Closes files
 
      integer*4 ierr
 
*   cov_parm(*)     - Covariance matrix with solution at the
*                   - end (Use a large value to force compiler into
*                   - using I*4 calcualtion of address)
      real*8 cov_parm(*)

* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

      data I8 / 1 /

 
***** Add apriori values back into solution
 
*                                                ! DEBUG, For use in
      ierr = 0
*                                                ! DEBUG to get solution
      if( ierr.ne.0 ) then
*                                                ! listed.
          call list_soln( cnum_parn, cov_parm, 1 )
      end if
 
      call glb_tot_solution( cov_parm((I8*cnum_parn)*cnum_parn+1) )
 
*     Now write out the matrices
* MOD TAH 190603: Change to writd8 to allow for >32767 parameters. 
      if( cnum_parn.gt.32767 ) then 
         call writd8(cglb_dcb, ierr, cov_parm, (I8*128)*cnum_par_vals,
     .              crec_par_vals)
         call report_error('WRITD8',ierr,'writ','COV_PARM',0,'GW_SOLN')
      else    ! Used original "small" code
         call writd(cglb_dcb, ierr, cov_parm, 128*cnum_par_vals,
     .              crec_par_vals)
         call report_error('WRITD',ierr,'writ','COV_PARM',0,'GW_SOLN')
      endif
  
* DEBUG
 
C     write(1,50) (crun_time(i), i=1,5), glb_file
C 50  format(/' Run ',5i3,' Global ',a20,/,
C    .        ' Parameter values ')
 
      if ( ierr.ne.0 ) then
          call list_soln( cnum_parn, cov_parm, 1 )
      end if
 
***** Thats all
      return
      end
 
