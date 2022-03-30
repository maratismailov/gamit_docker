C
      Subroutine NBIAS1( XLOW,NN,
     1   IBSLOT,NCHI,NCSAV,FNC,IAUTBS,ISIG,MOPT)
C
C     CALLED BY NBIAS (BROKEN UP FOR OVERLAYING)

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer*4 ntemp,isort,ih,nn,nchi,isig,ia1,ia2,l1b,ld1
     .        , ld2,ncut,issn,mopt,iautbs,ibslot
     .        , llbad,ncsav,k1,k2,i,j,k

      REAL*8 BL1,BL2,BL,cr,cj
      REAL*8 SRC,SRCRGE,sfctr,sfct
      REAL*8 XLOW(20),sort(20)
      real*8 fnc

      DIMENSION ISORT(20)
      DIMENSION NN(20)
      dimension ibslot(20)
      logical away
c
      SFCTR = 3.0d0
      IF(IAUTBS.EQ.1) GO TO 999
      call report_stat('FATAL','SOLVE','nbias1',' '
     .            , 'Interaction expected--something wrong',0)
      SFCT = 3.0d0
      WRITE (6,6000) SFCT
 6000 FORMAT (1X,'ENTER MULTIPLE OF BIAS UNCERTAINTY FOR SEARCH ',
     1 F3.1,'>: ')
      READ (5,6500) SFCTR
 6500 FORMAT (F3.0)
      IF (SFCTR.EQ.0.) SFCTR = SFCT
C
  999 CONTINUE
C
      NCHI = 0
      ncut = 5
      call bfwork(mopt,issn,k,k2,ld1,ld2,llbad,ia1,ia2,
     .     l1b,i,k1,away,1)
         call bfwork(mopt,issn,k,k2,ld1,ld2,llbad,ia1,ia2,
     .        l1b,i,k1,away,2)
         DO 2 I = IA1,IA2
            call bfwork(mopt,issn,k,k2,ld1,ld2,llbad,ia1,ia2,
     .           l1b,i,k1,away,3)
            if (away) goto 2
c---- We don't need keep big arrays to store all informations.  Only
c---- 11 with smallest uncertainties are necessary.
            IF (NCHI.LT.11) NCHI = NCHI+1
            IF(NCHI.GE.MBIAS) GO TO 2
            IBSLOT(NCHI) = K1
            NFIX(NCHI) = K+K2
C SEARCH RANGE
            IH = HALF(I)
            CR = 1.0/dble(IH)
            SRC = 3.D0*SFCTR*SCLERR*SIGMA(K+K2)*DBLE(IH)
            SRCRGE = DINT(SRC+.5D0*DSIGN(1.0D0,SRC))/DBLE(IH)
            CJ = SRCRGE+CR
C NEAREST INTEGER VALUE
            BL = ADJUST(K1)*DBLE(IH)
            BL1 = DINT(BL+.5D0*DSIGN(1.D0,BL))/DBLE(IH)
            XLOW(NCHI) = BL1-CJ
            CJ = 2.0*CJ+CR
            NN(NCHI) = NINT(CJ)
C
            BL2 = DABS(ADJUST(K1)-BL1)
            SORT(NCHI) = SRCRGE+BL2
            ISORT(NCHI) = NCHI
c---- reorder according to the uncertainties
            IF (NCHI.GE.11)
     .      CALL BSORT(NCHI,SORT,ISORT,XLOW,IBSLOT,NFIX,NN,1)

  2      CONTINUE
         IF (IBAND.EQ.2) ISSN = ISSN+L1BIAS
c
CD        WRITE (6,'(3X,I8,3X,F12.4)') (ISIGMA(I),SIGMA(I),I = 1,ISIG)
C
C EXIT IF NO BIASES REMAIN
       IF(NCHI.EQ.0) RETURN

          NCSAV = NCHI
          NTEMP = NCHI
CD               WRITE(6,*) NCHI,IAUTBS
                IF(IAUTBS.EQ.0) GO TO 1275
C SEARCH 5 BIASES AT A TIME
                NTEMP = MIN0(NCHI,ncut)
                GO TO 1278
1275            CONTINUE
          IF(NCHI.LE.ncut) GO TO 1278
c** This interactive code should be removed ---rwk 961123
          WRITE(6,1300) NCHI
1300      FORMAT(I5,' BIASES WILL BE SEARCHED. ENTER C.R. OR NUMBER ')
          READ(5,1302) NTEMP
1302      FORMAT(I5)
          IF(NTEMP.LE.0) NTEMP = NCHI
C-----------------------------------------------
1278      CONTINUE
C
C   SORT IN ORDER OF ERROR
      CALL BSORT(NCHI,SORT,ISORT,XLOW,IBSLOT,NFIX,NN,1)
C

C   ALL SORTED. NOW WE CAN REDUCE THE NUMBER WE SEARCH
        NCHI = NTEMP
cd        WRITE(6,333) NCHI,(IBSLOT(I),I = 1,NCHI)   
cd  333   FORMAT(2X,'NCHI = ',I3,'  INDEX = ',6I5)

C
1419    CONTINUE
C
        FNC = 1.
        DO 1303 I = 1,NCHI
           FNC = FNC*NN(I)
1303    CONTINUE
C
C
cd       WRITE(6,3) FNC
cd 3     FORMAT(F15.0,' CHI2 EVALUATIONS. PROCEED?')
       IF(IAUTBS.EQ.0) GO TO 1279
                IF (FNC.LE.5000) GO TO 1305
                NCHI = NCHI-1
                GO TO 1419
1279            CONTINUE
       call report_stat('FATAL','SOLVE','nbias1',' '
     .                 ,'Interaction expected--something wrong',0)
c**        CALL CHECK(IFLAG)
c          interactive check Y (=1) / N (=2)  
c**        IF(IFLAG.NE.1) GO TO 999
1305    CONTINUE
C
C REMOVE FIXED BIASES FROM SIGMA ARRAY OF LIVE PARAMETERS
C         
       if( isig.gt.0 ) then
        DO 10 I = 1,NCHI
           DO 15 J = 1,ISIG
              IF(IBSLOT(I).EQ.ISIGMA(J)) ISIGMA(J) = -ISIGMA(J)
   15      CONTINUE
   10   CONTINUE           
       endif
C
C   SORT IN ORDER OF index.
      CALL BSORT(NCHI,SORT,ISORT,XLOW,IBSLOT,NFIX,NN,4)
C
       RETURN
       END
