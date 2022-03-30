      subroutine force_parms( cov_parm, sol_parm, num_parm,
     .    cov_col, sol_col, nc, num_force, force_values,
     .    force_var, dchi, avat, equ_gn, ipivot, scale )
 
      implicit none 

*     This routine will force parameters in a covariance matrix to have
*     specific values.  The parameters to be forced are given in the
*     NC array.  There are num_force of them and the values to be forced
*     to are given in the force_values array.  If local array mapping
*     is to be then map_force should be called before this routine
*     setup the mapping of the arrays needed.
 
*   num_parm        - Number of parameters in cov_parm.  Cov_parm
*               - is assumed to be dimensioned to this size.
*   num_force   - Number of parameters to be focesd to specific
*               - values
*   nc(num_force)   - The list of parameters to be made equal.
*   ipivot(num_force)   - The pivot elements for the matrix
*               - inversion
 
      integer*4 num_parm, num_force, nc(num_force), ipivot(num_force)
 
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
 
      real*8 cov_parm(num_parm,num_parm), sol_parm(num_parm),
     .    cov_col(num_parm, num_force), sol_col(num_force),
     .    force_values(num_force), force_var(num_force),
     .    avat(num_force,num_force),
     .    equ_gn(num_parm,num_force), scale(num_force), dchi
 
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
 
CTITLE MAP_FORCE
 
      subroutine map_force(first, icov_col, isol_col, iavat,
     .            iequ_gn, iscale, iipivot, max_space, num_parm,
     .            num_force)

      implicit none 
 
*     This routine sets up the mapping for an Interger*4 array
*     to contain the variables we need for equating parameters.
*
*   first       - location in array of first available word
*   icov_col    - Start of the cov_col matrix
*   isol_col    - Start of the sol_col vector
*   iavat       - Start of the avat matrix
*   iequ_gn     - Start of the equ_gn matrix
*   iscale      - Start of the scale vector
*   iipivot     - Start of the ipivot vector
*   max_space   - maximum number of words allowed to be used.
 
*   num_parm        - NUmber of parameters in solution
*   num_force       - number of conditions
 
      integer*4 first, icov_col, isol_col, iavat, iequ_gn, iscale,
     .    iipivot, max_space, num_parm, num_force
 
* LOCAL VARIABLES
 
*   final       - Last word used.
 
      integer*4 final

* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

      data I8 / 1 /
 
***** Start the mapping.
      icov_col = first
      isol_col = icov_col   + 2*(I8*num_parm)*num_parm
      iavat     = isol_col  + 2*num_parm
      iequ_gn   = iavat     + 2*num_force*num_force
      iscale    = iequ_gn   + 2*num_parm*num_force
      iipivot   = iscale    + 2*num_force
      final     = iipivot   + num_force
 
***** See if have enough memory
*                                         ! Not enough memory.
      if( final.gt.max_space ) then
          write(*,100) final, max_space
 100      format(' MAP_EQUATE ERROR: Not enough memory: ',i9,
     .        ' Integer*4 words required,',/,
     .        19x,'and only ',i8,' words available.')
          stop '  MAP_EQUATE ERROR: Not enough memory'
      end if
 
****  Thats all
      return
      end
 
