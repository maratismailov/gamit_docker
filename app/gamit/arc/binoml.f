Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE BINOML(nn)
C         I        I!
C     COMPUTE BINOMIAL COEFFICIENTS OF THE FORM       =  ----------
C                                                  K      K!(I-K)!
C
C     UP TO A MAXIMUM VALUE OF I AND K
C     RICK ABBOT -  NOVEMBER 1984
C
      implicit none

      include '../includes/dimpar.h'   
      include '../includes/arc.h'

      integer nn,i,k,loop

      BINOM(1,1)=1.D0
      BINOM(2,1)=1.D0
      BINOM(2,2)=1.D0
      DO 25 I=3,nn
      BINOM(I,1)=1.D0
      LOOP=I-1
      DO 20 K=2,LOOP
   20 BINOM(I,K)=BINOM(I-1,K-1)+BINOM(I-1,K)
   25 BINOM(I,I)=1.D0
C
      RETURN
      END
