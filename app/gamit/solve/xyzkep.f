Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE XYZKEP(X,COND)
C
C     CONVERTS STATE VECTOR TO KEPLERIAN ELEMENTS
C     MODIFIED FROM PEP ROUTINE CHNCNC
C
C GOOSE       = INPUT VALUE OF SQUARE ROOT OF GRAVITATIONAL CONSTANT
C TIMES MASS OF CENTRAL BODY. UNITS ARE L**3/2*(T**-1)
C X(1-6)      = INPUT CARTESIAN COORDINATES
C     COND(1) = A    = SEMI-MAJOR AXIS
C     COND(2) = E    = ECCENTRICITY
C     COND(3) = INC  = INCLINATION
C     COND(4) = ASC  = RIGHT ASCENSION OF ASCENDING NODE
C     COND(5) = PER  = ARGUMENT OF PERIGEE
C     COND(6) = ANOM = INITIAL MEAN ANOMALY
C
      implicit none  

      real*8 X(6),COND(6),RCV(3),Y(2),DY(2),inc,mu,twopi,convd
     .     , goose,g2,g,p,p2rt,r,r2,v,v2,cinc,sinc,casc,sasc
     .     , asc,a,e2,absa,a2rt,quan3,quan4,quan5,e,secc,cecc,ecc
     .     , anom,sper,cper,per

      DATA MU/398603.46D0/
      TWOPI=8.D0*DATAN(1.D0)
      CONVD=TWOPI/360.D0
C
      GOOSE=DSQRT(MU)
      CALL CROSS(X(1),X(4),RCV)
      CALL DOT(RCV,RCV,G2)
C      G2 = DOT(RCV(1),RCV(1),3)
      G  = DSQRT(G2)
      P  = G2/MU
      P2RT = G/GOOSE
      QUAN5 = G
      CALL DOT(X(1),X(1),R2)
C      R2 = DOT(X(1),X(1),3)
      R = DSQRT(R2)
      CALL DOT(X(4),X(4),V2)
C      V2 = DOT(X(4),X(4),3)
      V = DSQRT(V2)
C
C	 CALCULATE INCLINATION AND ASCENDING NODE
      CINC = RCV(3)/G
      SINC =1.0D0-CINC**2
      IF(SINC.GT.0.0D0) GO TO 105
      CINC =DSIGN(1.0D0,CINC)
  105 INC = DACOS(CINC)
      CASC = -RCV(2)
      SASC = RCV(1)
      IF((CASC.NE.0.0D0).AND.(SASC.NE.0.0D0)) GO TO 111
      ASC  = 0.0D0
      GO TO 121
  111 ASC = DATAN2(SASC,CASC)
      IF( ASC.LT.0.0D0 ) ASC=ASC+TWOPI
  121 CONTINUE
C
C     SEMI-MAJOR AXIS AND ECCENTRICITY
      A = 2.0D0/R-V2/MU
      E2 = 1.0D0-P*A
      A = 1.0D0/A
      ABSA = DABS(A)
      A2RT = DSQRT(ABSA)
C
C     ELLIPTIC MOTION
      IF( E2.GT.1.0D-16   ) GO TO 341
C  CIRCULAR MOTION
      E2 = 0.0D0
      Y(1) = R
      Y(2) = 0.0D0
      DY(1) = 0.0D0
      DY(2) = V
  341 QUAN3 = P2RT*A2RT
      QUAN4 = A2RT*GOOSE
      E = DSQRT(E2)
      CALL DOT(X(1),X(4),SECC)
      SECC=SECC/(QUAN4*E)
C      SECC = DOT(X(1),X(4),3)/(QUAN4*E)
      CECC = (1.0D0-R/ABSA)/E
      Y(1) = A*(CECC-E)
      Y(2)   =QUAN3*SECC
      DY(1)= -QUAN4*SECC/R
      DY(2)  =QUAN5*CECC/R
      ECC = DATAN2(SECC,CECC)
      ANOM   = (ECC-E*SECC)
      IF(ANOM.LT.0.0D0) ANOM=ANOM+TWOPI
C
C     ARGUMENT OF PERIGEE
      IF(SINC.GT.0D0) GOTO 521
      SPER = X(1)*DY(1)-X(4)*Y(1)
      CPER = X(1)*DY(2)-X(4)*Y(2)
      GO TO 525
  521 SPER = X(3)*DY(2)-X(6)*Y(2)
      CPER = X(6)*Y(1)-X(3)*DY(1)
  525 PER = DATAN2(SPER,CPER)
      IF( PER.LT.0.0D0 ) PER=PER+TWOPI
C
C
C    MOVE ORBITAL ELEMENTS INTO OUTPUT VECTOR
      COND(1) =A
      COND(2) =E
      COND(3) =INC/CONVD
      COND(4) =ASC/CONVD
      COND(5) =PER/CONVD
      COND(6) =ANOM/CONVD
C
      RETURN
      END
