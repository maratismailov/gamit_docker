CTITLE CHECK_COVAR
 
      subroutine check_covar( cov, dim, num )

      implicit none 
 
 
*     Routine to check the diagonal elements of a covariance
*     matrix for elements less than zero.  If the values if less than
*     abs(1.d-6) then the value is set to zero.  Any other negative
*     value will cause a warning to be printed and the value will
*     be set to its abs value
 
*   dim     - Dimension of the covariance matrix
*   i       - Loop counter
*   num     - Number of elements to be checked
 
      integer*4 dim, i, num
 
*   cov(dim,dim)  - EMA covariance matrix
 
      real*8 cov(dim,dim)
 
 
***** Scan down the diagonal
 
      do i = 1, num
          if( cov(i,i).lt.0 ) then
 
*             See if small (probably due to rounding error for a parameter
*             with zero apriori variance)
*                                            ! Just set to zero
              if( cov(i,i).gt. -1.d-5 ) then
                  cov(i,i) = 0.d0
              else
                  write(*,100) i,cov(i,i)
  100             format(' **WARNING** Row ',i5,' has negative',
     .                   ' diagonal ',d15.6,/,
     .                   '   Value set to ABS in CHECK_COVAR')
                  cov(i,i) = abs(cov(i,i))
              end if
          end if
      end do
 
***** Thats all
      return
      end
 
