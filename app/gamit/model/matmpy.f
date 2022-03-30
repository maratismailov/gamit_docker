      SUBROUTINE MATMPY(A,B,C,NROW,NCOLA,NCOLB)
C
C     PREMULTIPLY MATRIX B BY MATRIX A WITH RESULTS IN MATRIX C
C     RICK ABBOT - JANUARY 1983
C
C     PARAMETERS
C          A       I  INPUT MATRIX
C          B       I  INPUT MATRIX
C          C       O  OUTPUT MATRIX
C          NROW    I  NUMBER OF ROWS OF MATRIX A
C          NCOLA   I  NUMBER OF COLUMNS OF MATRIX A
C          NCOLB   I  NUMBER OF COLUMNS OF MATRIX B
C
      implicit none
      REAL*8 A(NROW,NCOLA),B(NROW,NCOLB),C(NROW,NCOLB)
c**   The above should read B(NCOLA,NCOLB),but this version is now
c**   replaced by the /orbits version, moved to the library.  The
c**   error probably had no effect since all matrices used would be
c**   3 x 3 (but I didn't check this) -- rwk 92/1/30
C
      DO 300 I=1,NROW
      DO 200 J=1,NCOLB
      SUM = 0.D0
      DO 100 N=1,NCOLA
      SUM = SUM + A(I,N)*B(N,J)
  100 CONTINUE
      C(I,J) = SUM
  200 CONTINUE
  300 CONTINUE
C
      RETURN
      END
