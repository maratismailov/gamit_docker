
CTITLE UPDATE_SOL

      subroutine update_sol(na, ne, dNL1, dNL2)

      implicit none

*     Routine to update the sol_vec, cov_param and amgiq values
*     when a term has been resolved.

      include '../includes/const_param.h'   
      include 'track_com.h'

* PASSED VARIABLES
      integer*4 na   ! Ambiguity number that was just resolved
     .,         ne   ! Number of equations (1 == LC, 2 = L1+L2)
     .,         dNL1, dNL2 ! Number of cycles change at L1 and L2

* LOCAL

      real*8 lambda(2)  ! Wavelengths

      integer*4 i,j   ! Loop counters
     .,         ipivot(2)  ! Used in inversion
     .,         np(2) ! Parameters to be forced.

*   cov_col(max_parm, 2)    - The columns of the
*               -  covariance matrix for the parameters
*               - being forced.
*   sol_col(2)      - The change in the parameter
*               - estimates needed to get the forced values.
*   var_force(2) - Variances of focced parameters
*   dsol(2)     - Variance with which value should be forced.
*   avat(2,2)   - The multiplication of
*               - the partials matrix of (00001000,000001000)
*               - and the covaraiance matrix.
*   equ_gn(max_parm)  - Kalman gain matrix for
*               - getting solution.
*   scale(2)        - Scaling vector for solution. (passed
*               - to invert_vis)
*   dchi        - Change in Chi**2
*   sol_save(2) - Saved values of solution vector
 
      real*8 cov_col(max_parm, 2), sol_col(2),
     .    dsol(2), var_force(2),
     .    avat(2,2),
     .    equ_gn(max_parm,2), scale(2), dchi, sol_save(2)
 

****  Compute value to which sol_vec will need to be forced
      np(1) = amb_parn(1,na)
      if( np(1).le.0 ) then
          print *,'ERROR in update_sol Amb ',na,np
          RETURN
      endif

      lambda(1) = (vel_light/fR1)
      lambda(2) = (vel_light/fR2) 

****  OK get change value
      if( ne.eq. 1 ) then
          dsol(1) = (lcf1*dNL1+lcf2*dNL2)*lambda(1)/
     .                  ((lcf1+lcf2)*lambda(1))
          var_force(1) = 0.d0
          sol_save(1) = sol_vec(np(1))

          call force_track( cov_parm, sol_vec, max_parm, num_parm,
     .         cov_col, sol_col, np, 1,  dsol,
     .         var_force, dchi, avat, equ_gn, ipivot, scale)
          print *,'Force1 ',na, np(1),sol_save(1), sol_vec(np(1)),
     .                     dsol(1),dchi

*         Now resave ambiq values 
          do i = 1,num_ambs
             if( amb_parn(1,i).gt.0 ) then
                ambiq_float(1,i) = sol_vec(amb_parn(1,i))
                ambiq_flsig(1,i) = sqrt(cov_parm(amb_parn(1,i), 
     .                                           amb_parn(1,i))) 
             endif
          end do
*     else do L1 an dL2
      else if( amb_parn(2,na).gt.0 ) then
          np(2) = amb_parn(2,na)
          dsol(1) = dNl1/lambda(1)
          dsol(2) = dNL2/lambda(2)
          var_force(1) = 0.d0
          var_force(2) = 0.d0
          sol_save(1) = sol_vec(np(1))
          sol_save(2) = sol_vec(np(2))

          call force_track( cov_parm, sol_vec, max_parm, num_parm,
     .         cov_col, sol_col, np, 2,  dsol,
     .         var_force, dchi, avat, equ_gn, ipivot, scale)
          print *,'Force2 ',na, np,sol_save, sol_vec(np(1)),
     .                     sol_vec(np(2)),
     .                     dsol,dchi

*         Now resave ambiq values 
          do i = 1,num_ambs
             do j = 1,ne
                if( amb_parn(j,i).gt.0 ) then
                   ambiq_float(j,i) = sol_vec(amb_parn(j,i))
                   ambiq_flsig(j,i) = sqrt(cov_parm(amb_parn(j,i), 
     .                                              amb_parn(j,i))) 
                endif
             end do
          end do
      else
          print *,'ERROR in update_sol amb ',na,' Expected L2'
      endif

****  Thats all
      return
      end


      subroutine force_track( cov_parm, sol_parm, dim_parm, num_parm,
     .    cov_col, sol_col, nc, num_force, force_values,
     .    force_var, dchi, avat, equ_gn, ipivot, scale )

      implicit none
 
*     This routine will force parameters in a covariance matrix to have
*     specific values.  The parameters to be forced are given in the
*     NC array.  There are num_force of them and the values to be forced
*     to are given in the force_values array.  If local array mapping
*     is to be then map_force should be called before this routine
*     setup the mapping of the arrays needed.
 
*   dum_parm    - Dimensions of cov_param
*   num_parm        - Number of parameters in cov_parm.  Cov_parm
*               - is assumed to be dimensioned to this size.
*   num_force   - Number of parameters to be focesd to specific
*               - values
*   nc(num_force)   - The list of parameters to be made equal.
*   ipivot(num_force)   - The pivot elements for the matrix
*               - inversion
 
      integer*4 dim_parm, num_parm, num_force, 
     .          nc(num_force), ipivot(num_force)
 
*   cov_parm(num_parm,num_parm) - Full covariance matrix to
*               - be modified.
*   sol_parm(num_parm)      - The solution vector.  The nc
*               - elements of this vector will be set equal.
*   cov_col(num_parm, num_force)    - The columns of the
*               -  covariance matrix for the parameters
*               - being forced.
*   sol_col(num_force)      - The change in the parameter
*               - estimates needed to get the forced values.
*   force_values(num_force) - The values the forced parameters
*               - should take on.
*   force_var(num_force) - Variance with which value should be forced.
*   avat(num_force,num_force)   - The multiplication of
*               - the partials matrix of (00001000,000001000)
*               - and the covaraiance matrix.
*   equ_gn(num_parm,num_force)  - Kalman gain matrix for
*               - getting solution.
*   scale(num_force)        - Scaling vector for solution. (passed
*               - to invert_vis)
*   dchi        - Change in Chi**2
 
      real*8 cov_parm(dim_parm,dim_parm), sol_parm(dim_parm),
     .    cov_col(dim_parm, num_force), sol_col(num_force),
     .    force_values(num_force), force_var(num_force),
     .    avat(num_force,num_force),
     .    equ_gn(dim_parm,num_force), scale(num_force), dchi
 
* LOCAL PARAMETERS
 
*   i,j,k       - Loop counters
 
      integer*4 i,j,k
 
*   dsol, dcov  - Summation variables for computing corrections
*               - to the solution vector and covariance matrix.
 
 
      real*8 dsol, dcov
 
*                        T
***** Start, form the AVA  matrix, where A is of the form:
*     y = Ax (y is vector of zero observations)
*     and
*     A =  0 0 0 1 0 0 0 0 0 0....... to num_parm
*          0 0 0 0 0 0 0 1 0 0
*          0 0 0 0 0 0 0 0 1 0
*     where the above form would set parameters 4,8, and 9.
*     For num_force parameters being equated, there are num_force
*     rows in A.
*
*     The above form is much simpler to compute is we save the
*     columns of the covariance matrix for thhose parameters to
*     be forced first.
 
      do i = 1, num_parm
          do j = 1, num_force
              cov_col(i,j) = cov_parm(i,nc(j))
          end do
      end do
 
*     Save the forced parts of the solution vector as well
      do j = 1, num_force
          sol_col(j) = sol_parm(nc(j))
      end do
 
*               T
****  Now do AVA
*
      do i = 1, num_force
          do j = 1, num_force
              avat(i,j) = cov_col(nc(i),j)
* MOD TAH 000302: Added variance to diagonal for constrained 
*             forcing.              
              if( i.eq.j ) avat(i,i) = avat(i,i) + force_var(i) 
          end do
      end do
 
****  Now invert this matrix.  (If "sort of equal" was desired we could
*     add value to diagonal now representing variance of y above). Pass
*     zero as number in solution vector, we dont want to multiply.
*     kgain below is dummy argument.
 
      call invert_vis(avat, equ_gn, scale, ipivot, num_force,
     .               num_force, 0 )

**** Before continuing compute the change in Chi**2 due to condition
      dchi = 0.d0
      do i = 1, num_force  
         do j = 1, num_force 
            dchi = dchi + sol_col(i)*avat(i,j)*sol_col(j)
         end do
      end do
      dchi = dchi/num_force  

 
*     Now form the Kalman gain, equ_gn given by
*                T     T -1
*     equ_gn = VA  (AVA )
*
      do i = 1, num_parm
          do j = 1, num_force
 
*             Do the multiply (could use VIS but stick to straight)
*             call dwmul(equ_gn(i,j), col_col(i,1), num_force,
*    .                    avat(1,j),1, num_force)
 
              equ_gn(i,j) = 0.d0
              do k = 1, num_force
                  equ_gn(i,j) = equ_gn(i,j) + cov_col(i,k)*
     .                                     avat(k,j)
              end do
          end do
      end do
 
****  Now get the change to the solution vector
*
*     x  = x  -  equ_gn*(force_values-sol_col)
*      n    o
*
      do i = 1,num_parm
          dsol = 0.d0
          do j = 1, num_force
              dsol = dsol + equ_gn(i,j)*(force_values(j)-sol_col(j))
          end do
          sol_parm(i) = sol_parm(i) + dsol
      end do
 
*     Now update the covariance matrix
*
*     V  = V - equ_gn*cov_col
*      n    o
 
      do i = 1, num_parm
          do j = 1, num_parm
 
*             Do summation loop
              dcov = 0.d0
              do k = 1, num_force
                  dcov = dcov + equ_gn(i,k)*cov_col(j,k)
              end do
              cov_parm(i,j) = cov_parm(i,j) - dcov
          end do
      end do
 
****  Thats all
      return
      end
