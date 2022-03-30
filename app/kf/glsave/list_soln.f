 
CTITLE LIST_SOLN
 
      subroutine list_soln( num, cov, inc )
 
      implicit none 

* MOD TAH 190603: Introduced inc to reduce output 
*     List covariance matrix and solution
 
*   i,j,k - Loop counters
*   num   - Dimension of cov
 
      integer*4 i, num, inc
 
 
      real*8 cov(num,num)
 
 
****  Write results
      if( inc.le.0 ) inc = 1
      print *,'LIST_SOLN ',num,inc, cov(1,1),cov(1,2)
 
      do i = 1, num, inc
          write(*,100) i, cov(i,num+1), cov(i,i)
  100     format('GLSAVE LIST_SOLN',1x,i5,'.',f22.8,' VAR ',e22.4)
      end do
 
****  Thats all
      return
      end
