CTITLE GLB_BAK_FILTER
 
      subroutine glb_bak_filter( cov_sav_parm, sol_sav_parm, jmat,
     .                           cov_parm, sol_parm, row_copy,
     .                           part_pnt, a_part,
     .                           temp_gain, kgain, cov_obs, sol_obs,
     .                           glb_used, first_soln )

      implicit none   
 
*     Routine to do the back Kalman filtering.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   part_pnt(2, max_glb_deriv, cnum_parn) - pointers to the
 
      integer*4 part_pnt(2, max_glb_deriv, cnum_parn)
 
*   a_part(max_glb_deriv,cnum_parn) - compressed form of the
*           - partials matrix
 
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix for
*           - the global parameters
*   cov_obs(cnum_parn, cnum_parn)   - covariance matrxo for the
*           - parameters read from disk
 
*   cov_sav_parm(num_glb_parn,num_glb_parn) - Covariance matrix for
*           - the global parameters read from forward solution
*   jmat(num_glb_parn,num_glb_parn) - Gain matrux used for averaging
*           - the forward and bak solutions
 
*   kgain(cnum_parn,num_glb_parn)   - Kalman filter Gain matrix
 
*   row_copy(num_glb_parn)          - Copy of one row of the
*                           - covariance matrix.  Used when column
*                           - in state transsion matrix is less than
*                           - row
*   sol_parm(num_glb_parn)  - The solution vector for the global
*           - parameters
 
*   sol_obs(cnum_parn)  - Solution vector read from disk
 
*   sol_sav_parm(num_glb_parn)  - The solution vector for the global
*           - parameters read from the forward solution
 
*   temp_gain(cnum_parn,num_glb_parn) - Working space for computing
*           - the Kalman gain.
 
      real*8 a_part(max_glb_deriv,cnum_parn),
     .    cov_parm(num_glb_parn,num_glb_parn),
     .    cov_obs(cnum_parn, cnum_parn),
     .    cov_sav_parm(num_glb_parn,num_glb_parn),
     .    jmat(num_glb_parn,num_glb_parn),
     .    kgain(cnum_parn,num_glb_parn), row_copy(num_glb_parn),
     .    sol_parm(num_glb_parn), sol_obs(cnum_parn),
     .    sol_sav_parm(num_glb_parn),
     .    temp_gain(cnum_parn,num_glb_parn)

* glb_used - Logical to indicate glb file used
* first_soln - Set true for first day so no averaging to get cov_sav_parm 

      logical glb_used, first_soln 

      integer*4 i,j
 
 
***** Increment the filter for the next experiment being added.
 
*     Predict the values at the next step
 
      call glb_predict( cov_parm, sol_parm, row_copy )
 
****  Now average the foward and back solutions to get the smoothed
*     results
* MOD TAH 190528: Only form average if not first solution (first solution
*     is just the end of the forward solution).
!     do i = 1,num_glb_parn,100
!         write(*,120) 'BEF',i,sol_sav_parm(i), 
!    .                      (cov_sav_parm(i,j),j=i,i+10)
120       format(a,' SOL ',I5,F10.6,' COV ',11E18.6)
!         write(*,120) 'INP',i,sol_parm(i), 
!    .                      (cov_parm(i,j),j=i,i+10)
!     enddo 
      if( .not.first_soln )
     .call glb_average( cov_sav_parm, sol_sav_parm, jmat,
     .                  cov_parm, sol_parm )
!    call glb_average( cov_sav_parm, sol_sav_parm, jmat,
!    .                  cov_parm, sol_parm )

!      do i = 1,num_glb_parn,100
!          write(*,120) 'AFT',I,sol_sav_parm(i), 
!     .                      (cov_sav_parm(i,j),j=i,i+10)
!     enddo 

***** Now continue as usual (Same as forward filter)
*     Now compute the filter gain for this new data
 
      call glb_kalman_gain( cov_parm, a_part, part_pnt, cov_obs,
     .                      temp_gain, kgain )
 
*     Now add the information from this experiment
 
      call glb_update( cov_parm, sol_parm, a_part, part_pnt,
     .                 temp_gain, kgain, sol_obs, cov_obs,'N',
     .                 glb_used )

* MOD TAH 930208: Increment prefit residual chi**2
* MOD TAH 960905: Chi**2 increments moved into glb_update
c     call glb_inc_bak( cov_obs, sol_obs )  
 
***** Thats all
      return
      end
 
