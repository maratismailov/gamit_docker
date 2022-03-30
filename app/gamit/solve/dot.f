Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE DOT(A,B,DOTP)
C
C     DOT PRODUCT : DOT=A*B
C
      implicit none

      real*8 A(9),B(9),dotp
C
      DOTP= A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
C
      RETURN
      END
