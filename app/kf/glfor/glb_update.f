CTITLE GLB_UPDATE
 
      subroutine glb_update( cov_parm, sol_parm, a_part, part_pnt,
     .           temp_gain, kgain, sol_obs, cov_obs, save_chi,
     .           glb_used )
 
      implicit none 

 
*     The routine adds the contribution of the current observations
*     into the covariance matrix and the solution vector.
*
*     This is done using:
*      m+1    m           m
*     x    = x    + K(y-Ax   )
*      m+1    m+1         m+1
*
*     and
*
*      m+1    m        m
*     C    = C    - KAC
*      m+1    m+1      m+1
*
*            m+1
*     where x    is sol_parm  (the solution vector)
*            m+1
*
*           y    is the data vector (sol_obs)
*           A    is the partial matrix (a_part)
*           K    is the Kalman filter gain (kgain)
*            m
*     and   C    is the parameter covariance matrix (cov_parm)
*            m+1
*
*
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
      include '../includes/glb_hdr_def.h'
 
*   ir,ic       - Row and column counters
*   part_pnt(2,max_glb_deriv,cnum_parn) - Pointer to which
*                   - parameters as used in the partials matrix
 
 
      integer*4 ir,ic, part_pnt(2,max_glb_deriv,cnum_parn)
 
*   a_part(max_glb_deriv,cnum_parn)     - Compressed partials
*                   - matrix
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix of the
*                   - parameters being estimated
*   cov_obs(cnum_parn,cnum_parn)        - Covariance matrix of the
*                   - parameters from the SOLVK solution
*   kgain(cnum_parn,num_glb_parn)       - Kalman filter gain
*                   - matrix
*   scr_real(max_glb_parn)              - Temporary storage for
*                   - rows of matrix during calculations
*   sol_parm(num_glb_parn)              - Solution vector
*   sol_obs(cnum_parn)                  - Data vector
*   temp_gain(cnum_parn,num_glb_parn)   - Temporay matrix as
*                   - decribed above
 
      real*8 a_part(max_glb_deriv,cnum_parn),
     .    cov_parm(num_glb_parn,num_glb_parn),
     .    cov_obs(cnum_parn,cnum_parn), kgain(cnum_parn,num_glb_parn),
     .    scr_real(max_glb_parn), sol_parm(num_glb_parn),
     .    sol_obs(cnum_parn), temp_gain(cnum_parn,num_glb_parn)

* dchi  - Change in chi**2

      real*8 dchi
      
* save_chi  - Y if values should be saved.

      character*(*) save_chi

* glb_used  - Logival to denote if globla file used in solution

      logical glb_used 
 
 
*     Compute the prefit residuals using our best estimates of the
*     parameters

*     Make sure cov_obs is symmetric
      call sym_mat( cov_obs, cnum_parn)
 
      do ir = 1, cnum_used
          call glb_dot_part( scr_real(1), sol_parm(1),1,
     .         a_part(1,ir), part_pnt(1,1,ir), indx_pnt(ir) )
 
*         Subtract from observed value
          sol_obs(ir) = sol_obs(ir) - scr_real(1)
      end do

*     Increment the chi**2 value/

      call glb_inc_prefit( cov_obs, sol_obs, dchi, save_chi ) 

*     Now add the information from this experiment
      if( dchi.ge.0.d0 .and. dchi.lt.max_chi_inc ) then

 
*         Now multiple by Kalman gain to get the update to the parameter
*         values
 
          do ir = 1, num_glb_parn
              call DWDOT(scr_real(1), kgain(1,ir),1, sol_obs(1),1,
     .                   cnum_used)
 
*             Update parameter value
              sol_parm(ir) = sol_parm(ir) + scr_real(1)
          end do
 
****      Now update the covariance matrix
          do ic = 1, num_glb_parn
              do ir = 1, num_glb_parn
                  call DWDOT(scr_real(ir), temp_gain(1,ic),1,
     .                       kgain(1,ir),1, cnum_used)
              end do
 
*             Now subrtact from column of cov_parm
              do ir = ic, num_glb_parn
                  cov_parm(ir,ic) = cov_parm(ir,ic) - scr_real(ir)
              end do
 
*             Move values to upper diagonal
              call DWMOV(cov_parm(ic+1,ic),1,
     .                   cov_parm(ic,ic+1),num_glb_parn,
     .                   num_glb_parn-ic)
          end do
          glb_used = .true. 
      else
          glb_used = .false. 
          write(*,400) dchi
 400      format('** Chi**2 change is ',d16.8,' Not including this',
     .           ' data')
          if( log_unit.ne.6 ) write(log_unit,400) dchi
      end if

*     Make sure cov_parm is symetric
      call sym_mat( cov_parm, num_glb_parn )
 
***** Thats it, solution updated
C     call out_cov( cov_parm, num_glb_parn, 'FINAL POINT')
      return
      end
 
CTITLE SYM_MAT

      subroutine sym_mat ( mat, dim )

      implicit none 

*     Routine to make matrix symmetric

* dim  - Dimension of matrix
* mat  - Matrix to be checked.  The abs() of smallest off-diagonal is set into
*        matrix.

      integer*4 dim
      real*8 mat(dim,dim)

* i, j  - Loop counters

      integer*4 i, j

****  Loop over matrix setting off-diagonal terms
      do i = 1, dim-1
         do j = 1, i-1
            if( abs(mat(i,j)).lt.abs(mat(j,i)) ) then
                mat(j,i) = mat(i,j)
            else
                mat(i,j) = mat(j,i)
            end if
         end do
      end do

***** Thats all
      return
      end


