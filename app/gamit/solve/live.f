      Subroutine LIVE(FREE,nsite,ISTAT,ICOORD,ILIVE) 
C
C DETERMINE POSITION OF PARAMETER IN LIVE LIST
C     FREE    : FREE = 0, DEAD PARAMETER
C     FREE = 1, LIVE PARAMETER
C     NURT    : TOTAL NUMBER OF STATIONS
C     ISTAT   : STATION NUMBER
C     ICOORD  : COORDINATE NUMBER (ALWAYS 3)
C     ILIVE   : LIVE PARAMETER INDEX

      implicit none

      include '../includes/dimpar.h'

      integer free(maxprm),ilive,nsite,istat,icoord,mlive
     .      , indx,indx1,i,j

      dimension ilive(icoord)


      CALL ZERO1I(1,ICOORD,ILIVE)
      MLIVE = 0
      DO 10 I = 1,nsite
         IF (I.GT.ISTAT) GOTO 100
         INDX = 3*(I-1)
         DO 20 J = 1,ICOORD
            INDX1 = INDX+J
            IF(FREE(INDX1).EQ.0) GO TO 20
            MLIVE = MLIVE+1
            IF(I.EQ.ISTAT) ILIVE(J) = MLIVE
   20    CONTINUE
   10 CONTINUE
C
  100 CONTINUE
      RETURN
      END
