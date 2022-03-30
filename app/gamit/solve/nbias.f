      Subroutine NBIAS( lane,biskut,mopt )

c      lane    : 1= NL (NBIAS called by BISOPT for LC_HELP or L1L2_IND
c                2= WL (NBIAS called by GETWL for LC_HELP,
c                    by BISOPT for L1L2_IND and L1&L2)

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h' 
      include 'parameters.h'

C MAX NUMBER OF BIAS PARAMS = MAXSIT * MAXSAT * 2 (L1 AND L2)
C MODIFIED 4/26/85 TO SORT BIASES IN ORDER OF ERROR IF
C THERE ARE FEWER THAN THE MAXIMUM NUMBER OF BIASES
C MODIFIED 1/24/86 TO AUTOMATE BIAS FIXING
C PASS IN BISKUT SO QUESTION IS ONLY ASKED ONCE DND 870909
C Modified to fix BISOPT.FTN. -DND- 871217
C---- Change a lot to fit new solution-update algorithm.
C                                -DND- 880409
C     Be careful about the type of variables.
c Add explicit declaration of variables (new are lowercase)- rwk 930223


      REAL*8 CHIREL
C 10 LOWEST CHI-SQUARES
      REAL*8 DSAVE(10,20),DTEMP(10,20),XLOW(20),sort(20)
      REAL*8 MINCHI(10),MINTMP(10),BISKUT
      real*8 fnc

      integer*4 ipuse,isort,id,ii,jj,lh,kk,lk,lm,nn,mlive0
     .        , lane,nchi,isig,ifix,mopt,iautbs,ibslot,lmj
     .        , nchio,ifast,ncsav,mlive,i,j
      integer jelf
      INTEGER IDIAG
      DIMENSION ISORT(20)
      DIMENSION ID(20),NN(20)
      DIMENSION IBSLOT(20)
      logical pseu
C
*     we are looping
      integer iter    ! Iteration counter

      IFAST = 0
      pseu = .false.
      isig = 0
      if (mopt.eq.2) pseu = .true.
      ipuse = 0
C
C     FOR AUTOMATIC BIAS-FIXING (ALWAYS EXIT THIS WAY TO AUTOMATE
C      THE LAST SOLUTION)
      IAUTBS = 1
C
 999  CONTINUE
C
c====================== bias rounding =======================
c
* MOD TAH 031202: goto's come back to the 123 line, set iteration 
*     counter before hand
      iter = 0
 123  MLIVE0 = NLIVE
c
C---- Option 7: decision function with pseudo-range priority
c
      if (noptin.eq.7.and.pseu) then
         call nbiasp
         ipuse = ipuse +1
         pseu = .false.
         mlive = nlive
         if (mlive.lt.mlive0) goto 120
      endif
           
c      print *,'NBIAS calling NBIASR lane nlive,nfix ',lane,nlive
c      write(*,'(10i5)') (nfix(i),i=1,5)
      CALL NBIASR(lane,mopt )  
c      print *,'  after NBIASR nlive ',nlive
* MOD TAH 031202: Incremnent iteration counter and report number
*     of biases fixed
      iter = iter + 1
        
      MLIVE = NLIVE

* NOTE: TAH 031202: Option to use pseudorange seems to be only
*     called if no bias parameters are fixed.
      IF (MLIVE.EQ.MLIVE0.AND.MOPT.EQ.2) then
         if (.not.pseu) GOTO 284
         call nbiasp
         if( logprt ) write(6,5050) iter, 'NBIAS', mlive-nlive
         write(10,5050) iter, 'NBIAS', mlive-nlive
5050     format('NBIAS: Iteration ',i4,' Biases fixed in ',a,1x,i4)
   
         mlive = nlive
         ipuse = ipuse +1
* MOD TAH 031202: Commented to two lines below. Want to skip
*        multiple bias fix option associated with line 284 on.
C         if (mlive.eq.mlive0) goto 284
C         pseu = .false.
      endif

C----  update solution with fixed integer biases
c
 120  IFIX = MLIVE0 - MLIVE
* MOD TAH 031202: Only execute code below is ifix is > 0
*     (Added because goto 284 commented out).
      if( ifix.gt.0 ) then 
         CALL BNEW(MLIVE0,IFIX,FREE,BDEV,NFIX)
         DO 510 I = 1,MLIVE  
* MOD TAH 980310: Replaced calc by jelf call
            IDIAG = jelf(I+1)            !  = I*(I+1)/2
            SIGMA(I) = DSIGN(DSQRT(DABS(A(IDIAG))),A(IDIAG))
 510     CONTINUE

C SET ISIGMA ARRAY - WITH SLOTS OF LIVE PARAMETERS 

c      print *,'NBIAS 1 ntpart isigma ',ntpart
c      write(*,'(10i7)') (isigma(i),i=1,ntpart)
 
         ISIG = 0
         DO 98 I = 1,NTPART
            IF(FREE(I).EQ.0) GO TO 98
            ISIG = ISIG+1
            ISIGMA(ISIG) = I
   98    CONTINUE  
c      print *,'NBIAS 2 ntpart isigma ',ntpart
c      write(*,'(10i7)') (isigma(i),i=1,ntpart)
      endif
* MOD TAH 980310: Removed the ipuse test in goto so that
*     algorithm will be iterated while ever new biases are
*     being fixed.
C     IF (MLIVE.LT.MLIVE0.and.ipuse.le.1) GOTO 123
      if( mlive.lt.mlive0 ) goto 123
C
c====================== bias searching =======================
c
 284  CONTINUE
* MOD TAH 980310: Skip this whole code.  I think nchi is the
*     number of biases fixed by this routine (so in existing
*     code, if it is zero we goto 1988.  To skip code, go
*     directly to 1988.
*     Report number of biases fixed.
      if( logprt ) write(6,5055) iter-1
      write(10,5055) iter-1
 5055 format('NBIAS: ',i4,' Biases fixed to integers')
      goto 1988
      CALL NBIAS1(XLOW,NN,IBSLOT,NCHI,
     1        NCSAV,FNC,IAUTBS,ISIG,MOPT)
      IF (NCHI.EQ.0) GOTO 1988

      DO 10 I = 1,10
      MINCHI(I) = 3.0D+10
   10 CONTINUE
      CALL ZERO1I(1,20,ID)
C
C REMOVE CURRENT BIAS SUBSET FROM PARAMETER LIST
      DO 26 LMJ = 1,NCHI
      LM = IBSLOT(LMJ)
      BDEV(LMJ) = XLOW(LMJ)-ADJUST(LM)
   26 CONTINUE
C      WRITE (6,'(I9,F9.3,I9)') (I,BDEV(I),NFIX(I),I = 1,NCHI)
C
      MLIVE0 = NLIVE
      MLIVE = MLIVE0-NCHI
C
C LEAST SQUARES SOLUTION WITHOUT CURRENT BIAS SUBSET
C---- Calculate inverse(Q22)
      IFAST = 0

      CALL BNEW1(MLIVE0,MLIVE,BDEV,NFIX,CHI2,IFAST)
C
C---------------------------------------------------------------------
C
C HEAD OF BIAS SEARCH LOOP
C
6543  CONTINUE
C
C.....NOW CALCULATE CHI SQUARED FOR DIFFERENT BIAS VALUES
C
      IFAST = 1
      DO 29 I = 1,NCHI
         LM = IBSLOT(I)
         LH = HALF(LM-LPART)
         BDEV(I) = XLOW(I)+ID(I)/DBLE(LH)-ADJUST(LM)
29    CONTINUE
      CALL BNEW1(MLIVE0,MLIVE,BDEV,NFIX,CHI2,IFAST)
      CHIREL = CHI2
C
C.....AUXILIARY ARRAYS FOR TRANSFER
      DO 40 II = 1,10
         MINTMP(II) = MINCHI(II)
         DO 30 JJ = 1,NCHI
            DTEMP(II,JJ) = DSAVE(II,JJ)
   30    CONTINUE
   40 CONTINUE
C.....CHECK IF CHI SQUARED IS OF THE 10 SMALLEST; IF YES SAVE
C.....APPROPRIATE INTEGER BIAS VALUES
      DO 50 II = 1,10
         IF (CHIREL.GE.MINCHI(II)) GO TO 50
         MINCHI(II) = CHIREL
         DO 49 JJ = 1,NCHI
            LM = IBSLOT(JJ)
            DSAVE(II,JJ) = BDEV(JJ)+ADJUST(LM)
49       CONTINUE
         GO TO 55
   50 CONTINUE
      GO TO 100
   55 CONTINUE
      IF (II.EQ.10) GO TO 100
C.....SHIFT VALUES
      II = II+1
      DO 70 JJ = II,10
         MINCHI(JJ) = MINTMP(JJ-1)
         DO 60 KK = 1,NCHI
            DSAVE(JJ,KK) = DTEMP(JJ-1,KK)
   60    CONTINUE
   70 CONTINUE
  100 CONTINUE
C
C.....................................................
C
      ID(1) = ID(1)+1
      IF(NCHI.LE.1) GO TO 201
      DO 200 I = 2,NCHI
         J = I-1
         IF(ID(J).LT.NN(J)) GO TO 200
         ID(J) = 0
         ID(I) = ID(I)+1
200   CONTINUE
C
201   IF(ID(NCHI).LT.NN(NCHI)) GO TO 6543
C    END OF LOOP
C-----------------------------------------------------------------

C
      WRITE (6,1000) (MINCHI(JJ),JJ = 1,10)
      WRITE (10,1000) (MINCHI(JJ),JJ=1,10)
 1000 FORMAT (/,' Lowest chi-squares: ',10(1X,F7.0))

C
C COMPUTE CONTRASTS
      DO 7001 JJ = 2,10
       MINCHI(JJ)=((MINCHI(JJ)/MINCHI(1))-1.)*DSQRT(DBLE(NOBS-MLIVE))
7001  CONTINUE
      MINCHI(1) = 0.D0
C
      WRITE (6,1002) (MINCHI(JJ),JJ=1,10)
      WRITE (10,1002) (MINCHI(JJ),JJ=1,10)
1002  FORMAT (/,' (chi2/chi20 - 1.)*sqrt(df): ',10(1X,F7.2))

      WRITE(6,7020)
      WRITE (10,7020)
7020  FORMAT(/,' Parameter biases: ')
C
C---- CHANGE FORMAT TO OUTPUT 10 COMBINATIONS. -DND- 871216
      DO 90 KK = 1,NCHI
       LM = IBSLOT(KK)
       XLOW(KK) = DSAVE(1,KK)
       WRITE (6,1001) LM,(DSAVE(JJ,KK),JJ = 1,10)
       WRITE (10,1001) LM,(DSAVE(JJ,KK),JJ=1,10)
 1001 FORMAT (I5,10(1X,F7.1))
C
   90 CONTINUE
c      WRITE(10,1009)
c 1009 FORMAT(/)
C
C------------
C FOR AUTOMATIC BIAS SEARCH
      NCHIO=NCHI
      IF((IAUTBS.EQ.1).AND.(MINCHI(2).LT. BISKUT))
     1  CALL BISET(BISKUT,NCHI,NFIX,MINCHI,DSAVE)
C------------
C NO BIASES FIXED - TURN OFF AUTO
      IF(NCHI.EQ.0) IAUTBS=0
      IF (NCHI.EQ.0) GOTO 1399
C---- All biases have enough contrast. Original order keeps same.
      IF (NCHIO.EQ.NCHI) GOTO 210
C---- Kick out biases with smaller contrast
      CALL BSORT(NCHIO,SORT,ISORT,XLOW,IBSLOT,NFIX,NN,5)

C---- Fix biases with enough contrast
 210  NLIVE = MLIVE0-NCHI
      DO 220 J = 1,NCHI
         LK = IBSLOT(J)
         FREE(LK) = 0
         BDEV(J) = XLOW(J)-ADJUST(LK)
         ADJUST(LK) = XLOW(J)
 220  continue

C---  Remove fixed bias index from IDXB.
      DO 112 J = 1,MBIAS
         LK = IDXB(J)
         IF (LK.LE.0) GOTO 112
         IF (FREE(LK).EQ.0) IDXB(J) = -LK
 112  CONTINUE

C FINAL SOLUTION WITH FIXED INTEGER BIASES
      IFIX = MLIVE0 - NLIVE
      CALL BNEW(MLIVE0,IFIX,FREE,BDEV,NFIX)
      DO 610 I = 1,MLIVE
         IDIAG = I*(I+1)/2
         SIGMA(I) = DSIGN(DSQRT(DABS(A(IDIAG))),A(IDIAG))
 610  CONTINUE

1399  CONTINUE
C------------
C FOR AUTOMATIC BIAS SEARCH
      IF((IAUTBS.EQ.1).AND.(NCHI.LT.NCSAV)) GO TO 999
1988  continue
C
C  COUNT REMAINING LIVE PARAMETERS
      NLIVE = 0
      DO 443 I = 1,NTPART
         IF(FREE(I).EQ.1) NLIVE = NLIVE+1
  443 CONTINUE
C
      RETURN
      END
