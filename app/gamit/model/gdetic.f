      SUBROUTINE GDETIC(FIN,OUT,PART,A,FINV,SHFT )
C
C
      implicit none
C
C     Given geodetic coordinates in FIN return geocentric
C     coordinates in OUT and the Jacobian in PART
C     Formulae are closed form:
C
C       TAN(LAT')=((1-F)**2 + 2*F*H/A) * TAN(LAT)
C
C      S.A. Gourevitch 6/25/81

      REAL*8 FIN(3),OUT(3),PART(3,3),SHFT(3)
     .     , ZERO,ONE,TWO,PIOTWO,ONEFSQ,H,LAT,HILAT
     .     , fden,finv,fnum,temp,a,f,fac,del

      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
C
C HIGH LATITUDE CUTOFF
      DATA HILAT/5.D-4/
C
C PI/2
      PIOTWO = TWO*DATAN(ONE)
C
      LAT=FIN(1)
      H=FIN(3)
      F= ONE/FINV
      ONEFSQ = (ONE-F)**2
C
C LATITUDE : TAKE INTO ACCOUNT LAT=90
      DEL=DABS(LAT)-PIOTWO
      IF(DABS(DEL).GT.HILAT) THEN
         OUT(1)=DATAN((ONEFSQ+TWO*F*H/A) * DTAN(LAT))
      ELSE
         OUT(1)=DSIGN(PIOTWO + DEL/(ONEFSQ+TWO*F*H/A),LAT)
      ENDIF
C
C LONGITUDE
      OUT(2)=FIN(2)
C
C RADIUS
      FNUM=DCOS(LAT)**2 + ONEFSQ**2 * DSIN(LAT)**2
      FDEN=DCOS(LAT)**2 + ONEFSQ    * DSIN(LAT)**2
      FAC=DSQRT(FNUM/FDEN)
      OUT(3)= A*FAC + H*(ONE-F*F*DSIN(TWO*LAT)/TWO)
C
      CALL SHFTOR(SHFT,OUT)
C
C PARTIALS: D GEOCENTRIC(I)/D GEODETIC(J)
      TEMP = DCOS(LAT)**2 + ((ONEFSQ+TWO*F*H/A) * DSIN(LAT))**2
C
      PART(1,1) = (ONEFSQ+TWO*F*H/A) / TEMP
      PART(1,2) = ZERO
      PART(1,3) = (F*DSIN(TWO*LAT)/A) / TEMP
C
      PART(2,1) = ZERO
      PART(2,2) = ONE
      PART(2,3) = ZERO
C
      PART(3,1) = A*ONEFSQ*(ONEFSQ-ONE)*DSIN(TWO*LAT)/(TWO*FAC*FDEN**2)
      PART(3,2) = ZERO
      PART(3,3) = ONE-F*F*DSIN(TWO*LAT)/TWO
C
      RETURN
      END
