Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE MATTRN(A)
C
C     TRANSPOSE MATRIX A, WHERE A IS A 3X3 MATRIX
C     RICK ABBOT - DECEMBER 1984
C
C     PARAMETERS
C          A       I  INPUT MATRIX
C      
      implicit none
      real*8 a,b
      integer i,j
      DIMENSION A(3,3),B(3,3)
C
      DO 100 I=1,3
      DO 100 J=1,3
  100 B(J,I)=A(I,J)
      DO 200 I=1,3
      DO 200 J=1,3
  200 A(I,J)=B(I,J)
C
      RETURN
      END
