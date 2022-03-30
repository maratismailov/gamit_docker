CTITLE 'SCALE_MATRIX'
 
      subroutine scale_matrix( norm_eq, sol_vec, scale, nsize, nsol,
     .                         ndim)
c
c
*     Routine to scale the matrix before inversion
c
*   i,j                   - Loop counters
*   ndim                  - dimension of matrix
*   nsol	                - Number of solutions to be computed
*   nsize                 - Size of matrix to be inverted
 
      integer*4 i,j, ndim, nsol	, nsize
*
*   norm_eq(ndim,ndim)    - Symetric matrix to be inverted
*   sol_vec(ndim,nsol)    - Solution vector
*   scale(ndim)           - Scale needed to normalize the matrix
      real*8 norm_eq(ndim,ndim), sol_vec(ndim,nsol), scale(ndim)
*
c
c
****  Start by scaling the rowa
      do i = 1, nsize
         call dwsmy( scale(i), norm_eq(i,1), ndim, norm_eq(i,1), ndim,
     .               nsize)
         do j = 1, nsol
            sol_vec(i,j) = sol_vec(i,j)*scale(i)
         end do
      end do
c
***** Scale the columns
      do j = 1, nsize
         call dwsmy( scale(j), norm_eq(1,j),1, norm_eq(1,j),1, nsize)
      end do
c
***** Thats all
      return
      END
 
 
