CTITLE OUT_COV
 
      subroutine out_cov ( cov_parm, dim, title )

      implicit none  
 
*     Routine to write covariance matrix and solution vector
 
*   dim         - Size of matrix
*   i,j,k       - Loop counters
 
      integer*4 dim, i
 
*   cov_parm(dim, dim) - Covariance matrix
 
      real*8 cov_parm(dim, dim)
 
*   title       - description of stage we are at
 
      character*(*) title
 
 
*     Write title and solution vector
 
      write(*,100) title, (i, cov_parm(i,dim+1),i=1,dim)
  100 format(' At ',a,1000(:/,5(1x,i4,1x,e16.6)))
      write(*,'(a)') 'SQRT Diagonal'
      write(*,150) (i,sqrt(cov_parm(i,i)),i=1,dim)
  150 format(5(1x,i5,1x,d16.6))
C     write(1,160) (i,cov_parm(i,i),i=1,dim)
C 160 format(5(1x,i2,1x,f12.9))
 
 
C     do i = 1, dim
C         write(1,200) i,(cov_parm(j,i),j=1,dim)
C 200     format(i3,1x,10(5(f12.6,1x),:/,4x))
C     end do
 
***** Thats all
      return
      end
 
