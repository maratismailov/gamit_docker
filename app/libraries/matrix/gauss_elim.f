CTITLE 'GAUSS_ELIM'
 
      subroutine gauss_elim( norm_eq, sol_vec, ipivot, nsize, nsol,
     .                       ndim)

      implicit none
c
c
*     Gauss ellimination routine taken from Kalman_lib.ftn on HP1000
c
*          i              - Loop counters
*   ndim                  - dimension of matrix
*   nsize                 - Size of matrix to be inverted
*   icol, irow            - Row and column counters
*   jcol,jrow             - other row and solumns counters
*   ipivot(ndim)          - Pivot elements for inversion
*   imab                  - Index to max pivot
*   imax, jmax            - Row and colmum with max values
*   nsol                  - Number of solutions to be computed
      integer*4 i, ndim, nsize, icol, irow, jrow, ipivot(ndim),
     .    imab, imax, jmax, nsol
*
c
*       norm_eq(ndim,ndim)    - Symetric matrix to be inverted
*   pivot                 - Value of pivot element
*   sol_vec(ndim,nsol)    - Solution vector
*   temp                  - Temporary value
      real*8 norm_eq(ndim,ndim), pivot, sol_vec(ndim,nsol), temp
*
c
c
****   Initialize the pivot counters
      do i = 1, nsize
         ipivot(i) = i
      end do
c
***** Start looping on diagonal for searching for pivot
      do irow = 1, nsize
         icol = ipivot(irow)
c
*        Find maximum value in this row
         call dwmab( imab, norm_eq(icol, irow), ndim, nsize-irow+1)
         imax = imab + irow - 1
c
*        Get the pivot
         pivot = norm_eq(icol, imax)
c
*        See if pivot is too small
         if( abs(pivot).lt.1.d-12) then
             write(*,100) irow, pivot, 1.d-12
  100        format(/' Pivot element too small for row ',i3,/,
     .               ' Value is ',d15.6,', tolerancc is ',d15.6)
         end if
c
*****    Interchange the columns for the pivot
         jmax = ipivot(imax)
*                                  ! Interchange
         if( imax.ne.irow ) then
             ipivot(irow) = ipivot(imax)
             ipivot(imax) = icol
c
*            Swap matrix elements
             call dwswp( norm_eq(1,irow),1, norm_eq(1,imax), 1, nsize)
         end if
c
*****    Zero out the elements of this column and form col of inverse
         norm_eq(icol,irow) = norm_eq(jmax,irow)
         norm_eq(jmax,irow) = 1.d0
c
*        Form the column
         call dwsmy(1.d0/pivot, norm_eq(1,irow),1, norm_eq(1,irow),1,
     .              nsize)
c
         do jrow = 1, nsize
            if( jrow.ne.irow) then
                temp = - norm_eq(icol,jrow)
                norm_eq(icol,jrow) = norm_eq(jmax, jrow)
                call dwpiv(temp, norm_eq(1,irow),1, norm_eq(1,jrow),1,
     .                     norm_eq(1,jrow),1, nsize)
                norm_eq(jmax,jrow) = temp/pivot
            end if
         end do
c
*        Do solution
         do jrow = 1, nsol
            temp = - sol_vec(icol,jrow)
            sol_vec(icol,jrow) = sol_vec(jmax, jrow)
            call dwpiv( temp, norm_eq(1,irow),1, sol_vec(1,jrow),1,
     .                  sol_vec(1,jrow),1, nsize)
            sol_vec(jmax,jrow) = temp/pivot
         end do
c
*     end do loop over rows
      end do
c
****  Now change sign of solution vectoe
      if( nsol.gt.0 ) then
* MOD TAH 980609: Only change the sign of the used portion of the
*         solution vector (used to be single call for ndim*nsol)
          do jrow = 1, nsol
              call dwsmy( -1.d0, sol_vec(1,jrow),1, 
     .                           sol_vec(1,jrow),1, nsize)
          end do
      end if
c
***** Thats all
      Return
      END
 
