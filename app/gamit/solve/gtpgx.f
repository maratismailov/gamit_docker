Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE GTPGX(G,P,NRG,NCG,GTPG,WORK1,WORK2)
C
C     MATRIX MULTIPLICATION - QUADRATIC FORM
C     SEQUENTIAL STORAGE
C     G - GENERAL STORAGE MATRIX
C     P - SYMMETRIC STORAGE MATRIX
C     NRG - NUMBER OF ROWS IN G
C     NCG - NUMBER OF COLUMNS IN G
C     WORK1 - WORK VECTOR,DIMENSION NRG
C     WORK2 - WORK VECTOR,DIMENSION NCG

      implicit none

      integer*4 ncg,nrg,ir,k,iv,j2,iu,iy,iz,kl,iq,i,j,l,m

      REAL*8 G(*),GTPG(*),WORK1(*),WORK2(*),P(*),w1
C
      DO 10 I=1,NCG
C
      IR=NRG*(I-1)
C     LOOP THRU COLUMNS OF P
      DO 20 J=1,NRG
         J2=(J*J-J)/2
C        LOOP THRU ELEMENTS OF P BY ROW
         W1=0.0D0
         DO 15 K=1,NRG
            IV=IR+K
            IF(K-J) 16,17,17
   16       IU=K+J2
            GO TO 18
   17       IU=J+(K*K-K)/2
   18       W1=W1+G(IV)*P(IU)
   15    CONTINUE
         WORK1(J)=W1
 20   CONTINUE
C
      DO 35 L=I,NCG
         W1=0.D0
         IY=NRG*(L-1)
         DO 30 M=1,NRG
            IZ=IY+M
            W1=W1+WORK1(M)*G(IZ)
   30    CONTINUE
         WORK2(L)=W1
 35   CONTINUE
C
C     ADD CONTRIBUTION TO ROW OF GTPG
      DO 40 KL=I,NCG
      IQ=I+(KL*KL-KL)/2
      GTPG(IQ)=GTPG(IQ)+WORK2(KL)
   40 CONTINUE
C
   10 CONTINUE
C
      RETURN
      END
