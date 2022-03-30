CTITLE OUT_LOC_CORREL
 
      subroutine out_loc_correl( iout, covar, dimc, num, work, title)

      implicit none 
 
 
*     Routine to compute and output the correlations between the parameters
*     in COVAR.
 
 
*   dimc        - Dimension of the covariance matrix
*   i,j,k       - Loop counters
*   iout        - the output device LU number
*   num         - Number of parameters from COVAR whose correlations
*               - should be output
*   trimlen     - HP string length function
 
      integer*4 dimc, i,j, iout, num, trimlen
 
*   covar(dimc,dimc)    - Covariance matrix
*   work(dimc,dimc)     - Working area in which the correlations
*                       - will be stored (may by the same array as
*                       - COVAR)
 
      real*8 covar(dimc,dimc), work(dimc,dimc)
 
*   title               - Title to be output at the beginning of the
*                       - line
 
      character*(*) title
 
****  Compute the correlations
 
      do i = 1, num - 1
          do j = i+1, num
             if( covar(i,i).gt.0 .and. covar(j,j).gt.0 ) then
                 work(i,j) = covar(i,j)/sqrt(covar(i,i)*covar(j,j))
             else
                 work(i,j) = 0.d0
             end if
          end do
      end do
 
*     Now write out the correlations
 
      write(iout,100) title(1:max(1,trimlen(title))),
     .                ((work(i,j),j=i+1,num),i=1,num-1)
  100 format(5x,a,t39,10(3f13.4,:/,t39) )
 
***** Thats all
      return
      end
 
