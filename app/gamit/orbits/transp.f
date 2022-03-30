      subroutine transp(A,B,IROW,ICOL)
C
C Transpose a Matrix
C
      implicit real*8(a-h,o-z)

      integer*4 irow,icol,i,j

      real*8 a,b

      dimension a(irow,icol),b(icol,irow)

      do 10 i=1,irow
      do 10 j=1,icol
      b(j,i)=a(i,j)
   10 continue

      return
      end
