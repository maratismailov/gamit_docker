Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE ATPA(AT,IPNTAT,IROWAT,NRAT,NCAT,P,W,WORK)

      implicit none

      INTEGER  IPNTAT(*),IROWAT(*) 
      integer ndim,nrat,ncat,is,it,ir,iv,iu,ix,iy,iz,ii,ia,jelf
     .      , j,i,k,l,m

      REAL*8   AT(*),P(*),W(*),WORK(*),x
C
      NDIM=NRAT*(NRAT+1)/2
      DO 100 I=1,NDIM
 100  W(I)=0.0D0
      DO 10 I=1,NRAT
      IS=IROWAT(I+1)+1
      IT=IROWAT(I+2)
      IF(IS.GT.IT) GO TO 10
      DO 11 II=1,NCAT
   11 WORK(II)=0.D0
      DO 15 J=IS,IT
      IR=IPNTAT(J)
      IV=(IR*(IR-1))/2
      X=AT(J)
      DO 15 K=1,NCAT
      IF(IR-K) 16,17,17 
c     this statement fails with HP fort77 (v ?) with FFLAGS = +O3 +U77 -K
c  16 IU=IR+(K*(K-1))/2  
c     so substitute the following:
c   16 k1 = k*(k-1)
c     iu = ir + k1/2  
c     go one step further and substitute a subroutine call (has checking) rwk 980121
   16 iu = ir + jelf(k)
      GO TO 18
   17 IU=K+IV
   18 WORK(K)=WORK(K)+X*P(IU)
   15 CONTINUE
      DO 30 L=I,NRAT
      IY=IROWAT(L+1)+1
      IZ=IROWAT(L+2)
      IF(IY.GT.IZ) GO TO 30
      IX=(L*(L-1))/2+I
      DO 50 M=IY,IZ
      IA=IPNTAT(M)
      W(IX)=W(IX)+WORK(IA)*AT(M)
   50 CONTINUE
   30 CONTINUE
C
   10 CONTINUE
      RETURN
      END
