      subroutine transp(A,B,IROW,ICOL)
C
C Transpose a Matrix

      implicit none

      integer*4 icol,irow,i,j
      real*8 a,b
C
      dimension a(irow,icol),b(icol,irow)
C
      do 10 i=1,irow
      do 10 j=1,icol
      b(j,i)=a(i,j)
   10 continue
C
      return
      end
