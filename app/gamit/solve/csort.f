Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE CSORT(DISTC,IBSTAT,NSTAT)
C
C     SORT DISTANCES FORM CENTROID BY SELECTION AND EXCHANGE (FLORES, P.32)

      implicit none

      real*8 DISTC(*),c1,c2,temp
      integer*4 IBSTAT(*),inext,nstat,ilow,itemp,i
C
C    FIND MINIMUM DISTANCE FROM CENTROID
      DO 50 INEXT=1,NSTAT-1
        ILOW=INEXT
        DO 100 I=INEXT,NSTAT-1
                C1=DISTC(ILOW)
                C2=DISTC(I+1)
                IF(C1.LE.C2) GO TO 100
                ILOW=I+1
  100   CONTINUE
                ITEMP=IBSTAT(INEXT)
                IBSTAT(INEXT)=IBSTAT(ILOW)
                IBSTAT(ILOW)=ITEMP
                TEMP=DISTC(INEXT)
                DISTC(INEXT)=DISTC(ILOW)
                DISTC(ILOW)=TEMP
C
   50 CONTINUE
C
C
      RETURN
      END
