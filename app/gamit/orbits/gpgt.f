Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE GPGT(G,P,NRG,NCG,GPGTM,WORK1,WORK2)
C
C     MATRIX MULTIPLICATION - ERROR PROPAGATION FORM
C     G - GENERAL STORAGE MATRIX
C     P - SYMMETRIC STORAGE MATRIX
C     GPGTM - PRODUCT (SYMMETRIC MODE)
C     NRG - NUMBER OF ROWS IN G
C     NCG - NUMBER OF COLUMNS IN G
C     WORK1 - WORK VECTOR,DIMENSION NRG
C     WORK2 - WORK VECTOR,DIMENSION NCG

      implicit none

      INTEGER*4 NRG,NCG,I,J,K,L,M,IND,J2,IU,IV,IY,IQ,KL
      REAL*8 G(*),P(*),GPGTM(*),WORK1(*),WORK2(*),W1
C
      IND=(NRG*(NRG+1))/2
      DO 1 I=1,IND
      GPGTM(I)=0.D0
    1 CONTINUE

      DO 2 I=1,NCG
      WORK2(I)=0.D0
2     CONTINUE

      DO 3 I=1,NRG
      WORK1(I)=0.D0
3     CONTINUE

      DO 10 I=1,NRG
C
C     LOOP THRU COLUMNS OF P
      DO 20 J=1,NCG
         J2=(J*J-J)/2
C     LOOP THRU ELEMENTS OF P BY ROW
         W1=0.0D0
         DO 15 K=1,NCG
            IV=I+NRG*(K-1)
            IF(K-J) 16,17,17
   16       IU=K+J2
            GO TO 18
   17       IU=J+(K*K-K)/2
   18       W1=W1+G(IV)*P(IU)
   15    CONTINUE
         WORK2(J)=W1
 20   CONTINUE
C
      DO 35 L=I,NRG
         W1=0.D0
         DO 30 M=1,NCG
            IY=NRG*(M-1)+L
            W1=W1+WORK2(M)*G(IY)
 30      CONTINUE
         WORK1(L)=W1
 35   CONTINUE
C
C     ADD CONTRIBUTION TO ROW OF GPGT
      DO 40 KL=I,NRG
      IQ=I+(KL*KL-KL)/2
      GPGTM(IQ)=GPGTM(IQ)+WORK1(KL)
   40 CONTINUE
C
   10 CONTINUE
C
      RETURN
      END
