      subroutine equate_parms( cov_parm, sol_parm, num_parm,
     .    cov_dcol, sol_dcol, nc, num_cond, dchi, avat, equ_gn,
     .    ipivot, scale, var )

      implicit none 
 
*     This routine will force parameters in a covariance matrix to have
*     the same value.  The parameters to be equated are given the
*     NC array.  There are num_cond of them.  If local array mapping
*     is to be then map_equate should be called before this routine
*     setup the mapping of the arrays needed.
 
*   num_parm        - Number of parameters in cov_parm.  Cov_parm
*               - is assumed to be dimensioned to this size.
*   num_cond        - Number of parameters to be set equal to each
*               - other.
*   nc(num_cond)    - The list of parameters to be made equal.
*   ipivot(num_cond-1)  - The pivot elements for the matrix
*               - inversion
 
      integer*4 num_parm, num_cond, nc(num_cond), ipivot(num_cond-1)
 
*   cov_parm(num_parm,num_parm) - Full covariance matrix to
*               - be modified.
*   sol_parm(num_parm)      - The solution vector.  The nc
*               - elements of this vector will be set equal.
*   cov_dcol(num_parm, num_cond-1)  - Difference in the columns
*               - of the covariance matrix between the first
*               - condition parameter and the reset.
*   sol_dcol(num_cond-1)        - difference in the solution
*               - vector for the first and remaining parameters.
*   dchi        - Change in Chi**2 due to condition applied
*   avat(num_cond-1,num_cond-1) - The multiplication of
*               - the partials matrix of ( ... 1...-1....) and
*               - the covaraiance matrix.
*   equ_gn(num_parm,num_cond-1) - Kalman gain matrix for
*               - getting solution.
*   scale(num_cond-1)       - Scaling vector for solution. (passed
*               - to invert_vis)
*   var         - Variance to be appear to the constraint
 
      real*8 cov_parm(num_parm,num_parm), sol_parm(num_parm),
     .    cov_dcol(num_parm, num_cond-1), sol_dcol(num_cond-1),
     .    dchi, 
     .    avat(num_cond-1,num_cond-1), equ_gn(num_parm,num_cond-1),
     .    scale(num_cond-1), var
 
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
*     A =  0 0 0 1 0 0 -1  0  0  0....... to num_parm
*          0 0 0 1 0 0  0 -1  0  0
*          0 0 0 1 0 0  0  0  0 -1
*     where the above form would equate parameters 4, 7, 8, and 10.
*     For num_cond parameters being equated, there are num_cond-1
*     rows in A.
*
*     The above form is much simpler to compute is we difference the
*     columns of the covariance matrix first.  We choose to difference
*     from the first paramter but any could be used.
* MOD TAH 150223: Make sure none of the equate entries is zero:
      do i = 1, num_cond
         if( nc(i).le.0 .or. nc(i).gt. num_parm) then
             write(*,110) num_cond, nc(1:num_cond)
 110         format('ERROR in condition parameter number: NumCond ',
     .              I5,' List ',100I5,/,'Not completing')
             return
         end if
      end do

      do i = 1, num_parm
          do j = 1, num_cond-1
              cov_dcol(i,j) = cov_parm(i,nc(1)) - cov_parm(i,nc(j+1))
          end do
      end do
*     Form difference of solution vector as well
      do j = 1, num_cond-1
          sol_dcol(j) = sol_parm(nc(1)) - sol_parm(nc(j+1))
      end do
 
*               T
****  Now do AVA
*
      do i = 1, num_cond-1
          do j = 1, num_cond-1
              avat(i,j) = cov_dcol(nc(1),j) - cov_dcol(nc(i+1),j)
              if( i.eq.j ) then
                  avat(i,i) = avat(i,i) + var
              end if
          end do
      end do
 
****  Now invert this matrix.  (If "sort of equal" was desired we could
*     add value to diagonal now representing variance of y above). Pass
*     zero as number in solution vector, we dont want to multiply.
*     kgain below is dummy argument.

*     Check avat to make sure that there are no zeros
      do i = 1,  num_cond-1
         if( avat(i,i).lt.1.d-16) then
* MOD TAH 030116: Moved message upto output equate.f routine
c             write(*,100) i, num_cond
c100         format(' ** WARNING ** Diagonal for element ',i3,
c    .              ' (of ',i3,' conditions) is zero',/,
c    .              '    NOT COMPLETING THIS EQUATE')
             dchi = -1
             return
        end if
      end do
 
      call invert_vis(avat, equ_gn, scale, ipivot, num_cond-1,
     .               num_cond-1, 0 )

**** BEfore continuing compute the change in Chi**2 due to condition
      dchi = 0.d0
      do i = 1, num_cond-1
         do j = 1, num_cond-1
            dchi = dchi + sol_dcol(i)*avat(i,j)*sol_dcol(j)
         end do
      end do
      dchi = dchi/(num_cond-1)

*     Now finish up solution 
*     Now form the Kalman gain, equ_gn given by
*                T     T -1
*     equ_gn = VA  (AVA )
*
      do i = 1, num_parm
          do j = 1, num_cond-1
 
*             Do the multiply (could use VIS but stick to straight)
*             call dwmul(equ_gn(i,j), col_dcol(i,1), num_cond-1,
*    .                    avat(1,j),1, num_cond-1)
 
              equ_gn(i,j) = 0.d0
              do k = 1, num_cond-1
                  equ_gn(i,j) = equ_gn(i,j) + cov_dcol(i,k)*
     .                                     avat(k,j)
              end do
          end do
      end do
 
****  Now get the change to the solution vector
*
*     x  = x  -  equ_gn*sol_dcol
*      n    o
*
      do i = 1,num_parm
          dsol = 0.d0
          do j = 1, num_cond-1
              dsol = dsol + equ_gn(i,j)*sol_dcol(j)
          end do
          sol_parm(i) = sol_parm(i) -dsol
      end do
 
*     Now update the covariance matrix
*
*     V  = V - equ_gn*cov_dcol
*      n    o
 
      do i = 1, num_parm
          do j = 1, num_parm
 
*             Do summation loop
              dcov = 0.d0
              do k = 1, num_cond-1
                  dcov = dcov + equ_gn(i,k)*cov_dcol(j,k)
              end do
              cov_parm(i,j) = cov_parm(i,j) - dcov
          end do
      end do
****  Thats all
      return
      end
 
CTITLE MAP_EQUATE
 
      subroutine map_equate(first, icov_dcol, isol_dcol, iavat,
     .            iequ_gn, iscale, iipivot, max_space, num_parm,
     .            num_cond)

      implicit none 
 
*     This routine sets up the mapping for an Interger*4 array
*     to contain the variables we need for equating parameters.
*
*   first       - location in array of first available word
*   icov_dcol   - Start of the cov_dcol matrix
*   isol_dcol   - Start of the sol_dcol vector
*   iavat       - Start of the avat matrix
*   iequ_gn     - Start of the equ_gn matrix
*   iscale      - Start of the scale vector
*   iipivot     - Start of the ipivot vector
*   max_space   - maximum number of words allowed to be used.
 
*   num_parm        - NUmber of parameters in solution
*   num_cond        - number of conditions
 
      integer*4 first, icov_dcol, isol_dcol, iavat, iequ_gn, iscale,
     .    iipivot, max_space, num_parm, num_cond
 
* LOCAL VARIABLES
 
*   final       - Last word used.
 
      integer*4 final

* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

      data I8 / 1 /
 
***** Start the mapping.
      icov_dcol = first
      isol_dcol = icov_dcol + 2*(I8*num_parm)*num_parm
      iavat     = isol_dcol + 2*num_parm
      iequ_gn   = iavat     + 2*(num_cond-1)*(num_cond-1)
      iscale    = iequ_gn   + 2*num_parm*(num_cond-1)
      iipivot   = iscale    + 2*(num_cond-1)
      final     = iipivot   + num_cond-1
 
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
 
