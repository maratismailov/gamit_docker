      SUBROUTINE HEIGHT( COORD,COORD2,A,FINV,SHFT )
C
C      Transform from geocentric to geodetic
C
      implicit none
C
      integer*4 i,j

      REAL*8 COORD(3),COORD2(3),TEMP(3),SHFT(3),DUMMY(9)
     .     , finv,a,del,delta

      DATA DELTA/0.00000001D0/
C
      DO 1 I=1,3
1     COORD2(I)=COORD(I)
C
      DO 10 J=1,150
      CALL GDETIC( COORD2,TEMP,DUMMY,A,FINV,SHFT )
      DO 5 I=1,3
      DEL=COORD(I)-TEMP(I)
      COORD2(I)=COORD2(I)+DEL
5     CONTINUE
      IF (DABS(DEL).LT.DELTA) GO TO 20
C      JCOUNT=J
   10 CONTINUE
C
20    CONTINUE
C      WRITE(6,30) JCOUNT
C30    FORMAT(I10,' ITERATIONS')
C
      RETURN
      END
