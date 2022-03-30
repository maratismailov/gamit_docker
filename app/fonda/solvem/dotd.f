      SUBROUTINE DOTD(A,B,C,MA,NA,MB,NB,MC,NC,L,M,N,IFLAG)
C
C     DIMENSIONS OF ARRAYS IN CALLING PROGRAM ARE :
C     A(MA,NA) , B(MB,NB) , C(MC,NC)
C     IE. IF B IS A VECTOR THEN IT IS DIMENSIONED B(MB,1)
C
C     THE SIZES OF THE ARRAYS TO BE MULTIPLIED ARE:
C     A(L,M)  *  B(M,N)  =  C(L,N)
C
C     IFLAG  =  1  A * B = C
C            =  2  A**TR * B = C
C            =  3  A * B**TR = C
C
      real*8 A(MA,NA) , B(MB,NB) , C(MC,NC)
      integer i,j,k,l,m,n,ma,na,mb,nb,mc,nc,iflag
C
C     C=0
C
      DO 1 K=1,N
      DO 1 I=1,L
1     C(I,K)=0.0D0
      GO TO (10,20,30) , IFLAG
C
C     A * B = C
C
   10 DO 11 K=1,N
      DO 11 J=1,M
      DO 11 I=1,L
   11 C(I,K)=C(I,K) + A(I,J)*B(J,K)
      RETURN
C
C     A**TR * B = C
C
   20 DO 21 K=1,N
      DO 21 I=1,L
      DO 21 J=1,M
   21 C(I,K)=C(I,K) + A(J,I)*B(J,K)
      RETURN
C
C     A * B**TR = C
C
   30 DO 31 J=1,M
      DO 31 K=1,N
      DO 31 I=1,L
   31 C(I,K)=C(I,K) + A(I,J)*B(K,J)
      RETURN
      END



