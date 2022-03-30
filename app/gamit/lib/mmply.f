Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
       SUBROUTINE MMPLY(A,B,R,L,M,N)
C
C   FORM THE MATRIX PRODUCT R=AB
C   THE MATRICES A AND B ARE RETURNED UNCHANGED

      implicit none 

      integer*4 ind1,ind2,ind3,i,j,k,l,m,n

      REAL*8 A(*),B(*),R(*)
C
C       DIMENSION A(L,M),B(M,N),R(L,N)
C
C       DO 5 I=1,L
C       DO 5 J=1,N
C       R(I,J)=0.D0
C       DO 5 K=1,M
C    5  R(I,J)=R(I,J)+A(I,K)*B(K,J)
C
        DO 10 I=1,L
         DO 10 J=1,N
          IND1=L*(J-1)+I
           R(IND1)=0.D0
           DO 10 K=1,M
            IND2=L*(K-1)+I
            IND3=M*(J-1)+K
            R(IND1)=R(IND1)+A(IND2)*B(IND3)
   10   CONTINUE
C
        RETURN
        END
