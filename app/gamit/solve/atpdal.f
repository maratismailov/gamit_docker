Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE ATPDAL(AT,IPNTAT,IROWAT,NRAT,NCAT,P,W,WORK)

      implicit none

      INTEGER  IPNTAT(*),IROWAT(*)
      integer nrat,ncat,ia,ii,ir,is,it,ix,iy,iz,i,j,l,m
      REAL*8   AT(*),W(*),WORK(*), p
C
      DO 5 I=1,NRAT*(NRAT+1)/2
    5 W(I)=0.D0
C
      DO 10 I=1,NRAT
      IS=IROWAT(I+1)+1
      IT=IROWAT(I+2)
      IF(IS.GT.IT) GO TO 10
      DO 11 II=1,NCAT
   11 WORK(II)=0.D0
      DO 15 J=IS,IT
      IR=IPNTAT(J)
      WORK(IR)=WORK(IR)+AT(J)*P
   15 CONTINUE
      DO 30 L=I,NRAT
      IY=IROWAT(L+1)+1
      IZ=IROWAT(L+2)
      IF(IY.GT.IZ) GO TO 30
      IX=(L*(L-1))/2+I
      DO 31 M=IY,IZ
      IA=IPNTAT(M)
      W(IX)=W(IX)+WORK(IA)*AT(M)
   31 CONTINUE
   30 CONTINUE
   10 CONTINUE
C
      RETURN
      END
