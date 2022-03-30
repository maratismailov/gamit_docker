      Subroutine USELSS( phi,iphi,igood )
C
C   CHECK IF THERE ARE ANY DOUBLE DIFFERENCES
C    THIS PART OF CODE FOR THIS SCENARIO:
C
C        X   -   X
C        X   X   -
C        -   X   -
C    NO DOUBLE DIFFERENCE IN THE ABOVE CASE.
C    SCAN COLUMNS OF ONEWAY PHASE MATRIX TO FIND DOUBLE MATCHES
C    ELIMINATE USELESS OBSERVATIONS


      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      real*8 PHI(MAXSIT,MAXSAT),sm

      integer IPHI(MAXSIT,MAXSAT),IPCHK(MAXSIT,MAXSAT)
     .      , igood,igood0,i,j,k,l

      SM=1.0D-8
C
C INITIALIZE ARRAY THAT DENOTES USEFUL OBSERVATIONS
      DO 123 J=1,nsat
         DO 123 I=1,nsite
         IPCHK(I,J)=0
  123 CONTINUE
C
      IGOOD=0
      DO 10 I=1,nsat-1
        DO 11 J=I+1,nsat
          DO 12 K=1,nsite-1
          IF (DABS(PHI(K,I)).LT.SM.OR.DABS(PHI(K,J)).LT.SM) GOTO 12
            IGOOD0=IGOOD
            DO 13 L=K+1,nsite
            IF (DABS(PHI(L,I)).LT.SM.OR.DABS(PHI(L,J)).LT.SM) GOTO 13
C A MATCH
            IPCHK(L,I)=1
            IPCHK(L,J)=1
            IGOOD=IGOOD+1
   13      CONTINUE
           IF (IGOOD.GT.IGOOD0) IPCHK(K,I)=1
           IF (IGOOD.GT.IGOOD0) IPCHK(K,J)=1
   12     CONTINUE
   11   CONTINUE
   10 CONTINUE
C
C REMOVE USELESS OBSERVATIONS
      DO 133 J=1,nsat
         DO 133 I=1,nsite
         IF (IPCHK(I,J).EQ.1) GOTO 133
         PHI(I,J)=0.D0
         IPHI(I,J)=0
  133 CONTINUE
C
CD     WRITE(6,80)((PHI(ISTAT,ISAT),ISAT=1,nsat),ISTAT=1,nsite)
cd   80 FORMAT(1X,6F12.3)
C
      RETURN
      END
