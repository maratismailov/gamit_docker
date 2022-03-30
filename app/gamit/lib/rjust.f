      SUBROUTINE RJUST(N,S)
C     SUBROUTINE TO RIGHT JUSTIFY STRING OF N CHARS

      implicit none

      integer*4 i,j,k,n
      CHARACTER*1 S(N)
C
      J=0
      DO 100 I=N,1,-1
         IF(S(I).NE.' ') GO TO 101
         J=J+1
100   CONTINUE
C
101   CONTINUE
      IF(J.EQ.0) RETURN
      DO 200 K=N-J,1,-1
         S(K+J)=S(K)
200   CONTINUE
C
      DO 300 K=1,J
300   S(K)=' '
C
      RETURN
      END
