CTITLE 'INVERT_VIS'
 
      subroutine invert_vis( norm_eq, sol_vec, scale, ipivot, nsize,
     .                       ndim, nsol)
c
c
c     Matrix inversion taken from Kalman_lib.ftn on HP1000.
* MOD TAH 980611: Changed to not scale the matrix if scale is passed
*     with -1 for all values from 1 to nsize.  This is to handle 
*     non-symetric, non-postive definate matrices.
c
*            i                     - Loop counters
*   ndim                  - dimension of matrix
*   ipivot(ndim)          - Pivot elements for inversion
*   nsol	                - Number of solutions to be computed
*   nsize                 - Size of matrix to be inverted
 
      integer*4 i, ndim, ipivot(ndim), nsol	, nsize
*
*      norm_eq(ndim,ndim)    - Symetric matrix to be inverted
*   sol_vec(ndim,nsol)    - Solution vector
*   scale(ndim)           - Scale needed to normalize the matrix
      real*8 norm_eq(ndim,ndim), sol_vec(ndim,nsol), scale(ndim)
*

*   no_scale -- Logical set true if we should not scale the matrix
*               (indicated by the scale factors being preset to -1)
      logical no_scale 
c

*     Check to see if we need to scale
      no_scale = .true. 
      do i = 1, nsize
         if( scale(i).ne. -1.d0 ) no_scale = .false.
      end do
c
***** Scale the matrx before inversion
      if( .not.no_scale ) then 
          call dwvmv( norm_eq(1,1), ndim+1, scale, 1, nsize)
c
*         Check the scale factors
          do i = 1, nsize
             if( scale(i).ne.0 ) then
                 scale(i) = 1.d0/sqrt(abs(scale(i)))
             else
                 scale(i) = 1.d0
                 norm_eq(i,i) = 1.d0
             end if
          end do
c
*****     Now scale the matric
          call scale_matrix( norm_eq, sol_vec, scale, nsize, nsol, ndim)
      end if 
c
***** Now do the Gauss Elimination
      call Gauss_elim( norm_eq, sol_vec, ipivot, nsize, nsol, ndim)
c
***** Descale the matrix
      if( .not. no_scale ) then
          call scale_matrix( norm_eq, sol_vec, scale, nsize, nsol, ndim)
      end if
c
***** Thats all
      return
c
      END
 
