      Subroutine legndr_p(z,zz,nzone,ntess,leg,gleg) 
        
c     R. King March 2017: Shortened version of legendr.f with derivatives omitted.
c     Eventually replace this with a single routine with a flag for derivatives,
c     and also incorporate a flag for normalized or not normalized functions.

Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
cc      SUBROUTINE LEGNDR(Z,ZZ,NZONE,NTESS,LEG,LEG1,GLEG,GLEG1)
C
C        F AMUCHASTEGUI - SEPTEMBER 1969 - SUBROUTINE LEGNDR
C        R ABBOT OCTOBER 1984, SMALL MODIFICATIONS FOR GPS
C
C        EVALUATION OF LEGENDRE POLYNOMIALS AND LEGENDRE
C        FUNCTIONS USING RECURSION FORMULAS
C
C        LEGNDR AND LEGND2 EVALUATE THE LEGENDRE POLYNOMIALS, ITS
C        DERIVATIVES AND ALSO THE ASSOCIATED LEGENDRE FUNCTIONS
C        THE LEGENDRE POLYNOMIALS ARE DEFINED AS:
C        PN(Z) = (1 / (N|*2**N) ) * (NTH DERIVATIVE OF (Z*Z-1)**N)
C        AND THE ASSOCIATED LEGENDRE FUNCTIONS AS
C        P(N,H) = ((DSQRT(1-Z*Z))**H) * PH(N) ,
C        WHERE PH(N) IS THE  HTH DERIVATIVE OF  PN(Z)
C        WE CAN SEE THAT FOR H=0  P(N,0)=PN(Z)
C
C        THE FORMULAE ARE FROM  "TABLES OF INTEGRALS, SERIES
C        AND PRODUCTS" - GRADSTEYN & RYZHIK,PAGS 1004-27
CNOTE    THE FACTOR (-1)**H  FOR THE LEGENDRE FUNCTION USED IN THIS
C        REFERENCE HAS BEEN SET = 1 TO USE THE CONVENTION OF
C        THE SMITHSONIAN ASTR. OBSERVATORY JOURNAL
C        THE FORMULAE FOR  P'(Z)  AND  P''(Z)  ARE FROM "DIFFERENTIAL
C        EQUATIONS WITH APPLICATIONS", RITGER & ROSE, PAGE 223
C
cc      REAL*8 LEG,LEG1,GLEG,GLEG1,Z,ZZ
      real*8 leg(*),gleg(*),z,zz
      REAL*8 A1,A2,DIK,SIGN
      integer nzone, ntess,nsize,i,k,l,m
c      DIMENSION LEG(*),LEG1(*),GLEG(*),GLEG1(*)
C
      NSIZE = NZONE-1
      IF( NTESS.GT.NZONE ) NSIZE = NTESS-1
C
C        EVALUATION OF  P(Z)
C        THE FIRST POLYNOMIAL IN THE ARRAY  LEG(N)  IS THE
C        SECOND ORDER POLYNOMIAL, SINCE THE FIRST ORDER IS NOT
C        USED IN P.E.P.
C
      LEG(1) = 1.5D0*Z**2 - 0.5D0
      LEG(2) = 2.5D0*Z**3 - 1.5D0*Z
cc      IF( NSIZE.LE.2 ) GO TO 20    
      if( nsize.le.2 ) go to 100
      DO 10 I=3,NSIZE
      LEG(I) = ( DBLE(2*I+1)*Z*LEG(I-1) - DBLE(I)*LEG(I-2) ) /
     1           DBLE(I+1)
  10  CONTINUE
C
C        EVALUATION OF  P'(Z)
c 20  LEG1(1) = 3.0D0*Z
c     LEG1(2) = 7.5D0*Z**2 - 1.5D0
c     IF( NSIZE.LE.2 ) GO TO 100
c     DO 50 I=3,NSIZE
c     LEG1(I) = LEG1(I-2) + DBLE(2*I+1)*LEG(I-1)
c 50  CONTINUE
C
C        EVALUATION OF P(N,H)
C        THE ORDER OF THE LEGENDRE FUNCTIONS IN THE ARRAY GLEG,
C        AS WELL AS THE ORDER OF THE PARTIALS IN GLEG1, GLEG2 IS
C        P(2,1), P(2,2), P(3,1), P(3,2), P(3,3), ....... , P(N,N)
C
 100  IF( NTESS.LE.1 ) GO TO 999
      GLEG(1) = 3.0D0*Z*ZZ
      GLEG(2) = 3.0D0*ZZ**2
      IF( NTESS.LE.2 ) GO TO 150
      L = 2
      DO 140 I=3,NTESS
      SIGN = 1.0D0
      M = ((I-1)*(I-2))/2
      DO 135 K=1,I
      SIGN = -SIGN
      L = L+1
      DIK = DBLE(I+K-1)
      IF( K.EQ.1 ) GO TO 120
      IF( K.NE.I ) GO TO 110
      GLEG(L) = SIGN*DIK*ZZ*GLEG(M-1)
      GO TO 130
 110  GLEG(L) = (Z*GLEG(M) + DIK*ZZ*GLEG(M-1)) * SIGN
      GO TO 130
 120  GLEG(L) = -Z*GLEG(M) - DIK*ZZ* LEG(I-2)
 130  M = M+1
 135  GLEG(L) = GLEG(L)*SIGN
 140  CONTINUE
 150  continue                            

C        EVALUATION OF  P'(N,H)
cc 150  GLEG1(1) = -3.0D0*(Z**2)/ZZ + 3.0D0*ZZ
cc      GLEG1(2) = -6.0D0*Z
cc     IF( NTESS.LE.2 ) GO TO 999
cc      L = 2
cc      DO 190 I=3,NTESS
cc      SIGN = 1.0D0
cc      DO 180 K=1,I
cc      SIGN = -SIGN
cc      L = L+1
cc      A1 = DBLE(K)
cc      A2 = DBLE((I-K+1)*(I+K))
cc      IF( K.EQ.1 ) GO TO 160
cc      GLEG1(L) = (A1*Z*GLEG(L) - A2*ZZ*GLEG(L-1)) * SIGN
cc      GO TO 170
cc 160  GLEG1(L) = -A1*Z*GLEG(L) + A2*ZZ* LEG(I-1)
cc 170  GLEG1(L) = GLEG1(L)*SIGN/(ZZ**2)
cc 180  CONTINUE
cc 190  CONTINUE
C
 999  RETURN
      END
