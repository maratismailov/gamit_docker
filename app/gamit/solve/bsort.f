Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE BSORT(NCHI,SORT,ISORT,XLOW,IBSLOT,NFIX,NN,KEY)
C
C     sort arrays from the smallest to the bigest
C     key : control variable.

      implicit none
              
      integer*4 nfix,nchi,isort,nn,ibslot,key,iswtch,i,k

      REAL*8 XLOW,X1,sort,temp 

      dimension nfix(*),ISORT(NCHI),NN(NCHI),IBSLOT(NCHI)
     .        , XLOW(NCHI),sort(nchi) 

C
                IF(NCHI.LE.1) GO TO 1417
                DO 1416 K=1,32000
                  ISWTCH=0
                  DO 1412  I=2,NCHI
               GO TO (10,20,30,40,50,60), KEY
 10                IF(SORT(I).GE.SORT(I-1)) GO TO 1412
                     GOTO 80
 20                IF(ISORT(I).GE.ISORT(I-1)) GO TO 1412
                     GOTO 80
 30                IF(XLOW(I).GE.XLOW(I-1)) GO TO 1412
                     GOTO 80
 40                IF(IBSLOT(I).GE.IBSLOT(I-1)) GO TO 1412
                     GOTO 80
 50                IF(NFIX(I).GE.NFIX(I-1)) GO TO 1412
                     GOTO 80
 60                IF(NN(I).GE.NN(I-1)) GO TO 1412
 80                     TEMP=SORT(I)
                        SORT(I)=SORT(I-1)
                        SORT(I-1)=TEMP
C
                        X1=XLOW(I)
                        XLOW(I)=XLOW(I-1)
                        XLOW(I-1)=X1
C
                        ISWTCH=ISORT(I)
                        ISORT(I)=ISORT(I-1)
                        ISORT(I-1)=ISWTCH
C
                        ISWTCH=NN(I)
                        NN(I)=NN(I-1)
                        NN(I-1)=ISWTCH
C
                        ISWTCH=NFIX(I)
                        NFIX(I)=NFIX(I-1)
                        NFIX(I-1)=ISWTCH
C
                        ISWTCH=IBSLOT(I)
                        IBSLOT(I)=IBSLOT(I-1)
                        IBSLOT(I-1)=ISWTCH
1412                    CONTINUE
                IF(ISWTCH.EQ.0) GO TO 1417
1416            CONTINUE
C
1417            CONTINUE
C
      RETURN
      END
