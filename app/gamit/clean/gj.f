      SUBROUTINE GJ(A,B,N,EP)
C
C     ----     DND  890102
C
C        Inverse a symmetric positive matrix using
C        Gauss-Jordan approach
C
C  INPUT:
C     N     :  dimension of matrix A
C     A     :  input A; output inverse(A)
C     B     :  working vector
C     EP    :  error index
C              EP = 1  :  OK
C              EP = -1 :  matrix is not positive
C
      implicit none
C
      integer*4 index,index1,index2,ep,i,j,k,kk,n
      real*8 a,b,g,w,small

      DIMENSION A(N*(N+1)/2),B(N)

      SMALL=1.0D-15
C
      EP=1
      DO 10 K=1,N
         KK=N-K+1
         W=A(1)
         IF (W.LE.SMALL) GOTO 30
         DO 20 I=2,N
            INDEX=I*(I-1)/2+1
            G=A(INDEX)
            IF (I-KK) 2,2,3
  3         B(I)=G/W
            GOTO 4
  2         B(I)=-G/W
  4         DO 20 J=2,I
            INDEX1=(I-1)*(I-2)/2+J-1
            INDEX2=I*(I-1)/2+J
            A(INDEX1)=A(INDEX2)+G*B(J)
 20      CONTINUE
         INDEX=N*(N+1)/2
         A(INDEX)=1.0D0/W
         DO 10 I=2,N
         INDEX=N*(N-1)/2+I-1
         A(INDEX)=B(I)
 10   CONTINUE
C      DO 40 I=2,N
C      DO 40 J=1,I-1
C      A(J,I)=A(I,J)
C 40   CONTINUE
      RETURN
 30   EP=-EP
CD     PRINT*,'Matrix is not positive!'
      RETURN
      END
