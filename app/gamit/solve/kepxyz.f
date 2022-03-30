Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
       SUBROUTINE KEPXYZ(T,COND,NV,Y,DY)
C
C      CONVERTS KEPLERIAN ELEMENTS TO STATE VECTOR (CARTESIAN)
C      COMPUTES JACOBIAN FOR COVARIANCE PROPAGATION
C      MODIFICATION OF PEP SUBROUTINE JNITL
C
C COND - ARRAY OF INITIAL CONDITIONS (KEPLERIAN ELEMENTS)
C USAGE OF ARRAY 'SETP' --- STORAGE FOR SAVED QUANTITIES FROM
C SETTING UP.  IT IS EFFECTIVELY EQUIVALENCED TO THE FOLLOWING:
C      B(3,2),TSV,A,E,ANOM0,MOTPI,MU2,SECC,CECC,QUAN2,QUAN3,QUAN4,QUAN5,
C      1      7   8 9  10    11    12  13   14   15    16    17    18
C
C      QUAN12(3),QUAN13(3),QUAN11,MOTION,MU,SASC,CASC,SPCI,CPCI
C            19    22       25      26   27  28   29   30   31
C
C      SINC,CINC,SPER,CPER
C       32   33   34   35
C
C NOTE: ALL VARIABLES UP THRU QUAN5 (SETP(18)) ARE USED IN ALL
C ELLIPTIC CALCULATIONS, BUT THOSE FROM THERE ON ARE NEEDED ONLY
C IF PARTIALS WITH RESPECT TO THE INITIAL CONDITIONS ARE NEEDED.
C
C NV - OPERATION INDICATOR.
C      1 DO POSITION AND VELOCITY
C     >1 DO PARTIALS AS WELL (UP TO DY(1-6,NV-1) )
C      0 DO JUST POSITION
C Y - ARRAY FOR OUTPUT COORDINATES (CARTESIAN)
C RY - OUTPUT RADIAL DISTANCE
C DY - OUTPUT ARRAY OF PARTIALS
C
C INPUT ELLIPTIC ELEMENTS
C  [units: km, -, dec. deg, dec. deg, dec. deg, dec. deg]
C     EQUIVALENCE (COND(1),AP),(COND(2),EP),(COND(3),INCP),
C    1 (COND(4),ASCP),(COND(5),PERP),(COND(6),ANOMP)

      implicit none
  
      character*256 message
              
      integer nv,int,iyb,i,j

      real*8 motion,mu,cond(6),setp(35),y(6),ybar(6),dy(6,6),ry,cf
     .     , t,goose,twopi,convd,quan1,th,cinc,sinc,casc,sasc,sper
     .     , cper,spci,cpci,anom,anoms,qq1,qq2,qq4,qq5,qq6,qq7,qq8


C  [units: km**3/sec**2]
      DATA MU/398603.46D0/
C  [units: km**(3/2)/sec]
      GOOSE=DSQRT(MU)
      TWOPI=8.D0*DATAN(1.D0)
C  Conversion degrees to radians
      CONVD=TWOPI/360.D0
C
C  SETP(7) NOT USED
      SETP(7)=0.D0
      IF(T.GT.0.D0) GO TO 700
C  [units: km]
      SETP(8)=COND(1)
C  [units: -]
      SETP(9)=COND(2)
C  [units: dec. deg/360]
      SETP(10)=COND(6)/360.D0
C  [units: -]
      QUAN1=1.D0-SETP(9)**2
c     test if quan1 is negative before taking square root!
      if (quan1 .ge. 0.d0) then
C  [units: -]
         SETP(15)=DSQRT(QUAN1)
      else
c        write error message to Q-file
         write(message,50) quan1
         write(15,'(a)') message
  50     format(1x,'Negative sqrt of quan1--check for cycle slips')
         call report_stat('FATAL','SOLVE','kepxyz',' ',message,0)
      endif
c  [units: km]
      SETP(16)=SETP(8)*SETP(15)
c  [units: km**(1/2)]
      SETP(17)=DSQRT(SETP(8))
c  [units: km/sec]
      MOTION=GOOSE/SETP(8)/SETP(17)
c  [units: km/sec/2pi
      SETP(11)=MOTION/TWOPI
c  [units: km**3/sec**2]
      SETP(12)=GOOSE**2
c  [units: km**2/sec]
      SETP(17)=SETP(17)*GOOSE
c  [units: km**2/sec]
      SETP(18)=SETP(17)*SETP(15)
c  [units: km**2/sec]
      SETP(17)=-SETP(17)
c  [units: radians]
      TH=COND(3)*CONVD
      CINC=DCOS(TH)
      SINC=DSIN(TH)
c  [units: radians]
      TH=COND(4)*CONVD
      CASC=DCOS(TH)
      SASC=DSIN(TH)
c  [units: radians]
      TH=COND(5)*CONVD
      SPER=DSIN(TH)
      CPER=DCOS(TH)
      SPCI=SPER*CINC
      CPCI=CPER*CINC
c  [units: -]
      SETP(1)= CASC*CPER-SASC*SPCI
c  [units: -]
      SETP(4)=-CASC*SPER-SASC*CPCI
c  [units: -]
      SETP(2)= SASC*CPER+CASC*SPCI
c  [units: -]
      SETP(5)= CASC*CPCI-SASC*SPER
c  [units: -]
      SETP(3)= SPER*SINC
c  [units: -]
      SETP(6)= CPER*SINC
C DECIDE IF PARTIALS SETUP NEEDED
C      IF(KIND.EQ.0) RETURN
C FILL REST OF SETP ARRAY FOR PARTIALS
c  [units: -]
      SETP(25)=SETP(9)/QUAN1
c  [units: km/sec]
      SETP(26)=MOTION
c  [units: km**(3/2)/sec]
      SETP(27)=GOOSE
c  [units: -]
      SETP(28)=SASC
c  [units: -]
      SETP(29)=CASC
c  [units: -]
      SETP(30)=SPCI
c  [units: -]
      SETP(31)=CPCI
      DO 100 I=1,3
c  [units: km]
      SETP(I+18)=SETP(I)*SETP(8)
c  [units: km/sec]
  100 SETP(I+21)=SETP(I+3)*SETP(25)
      DY(3,4)=0.D0
      DY(6,4)=0.D0
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  700 CONTINUE
C COMPUTE ELLIPTIC QUANTITIES FOR TIME T FROM INITIAL EPOCH.
C CORRESPONDING SETP ARRAY MUST ALREADY HAVE BEEN SET UP.
C
C GET MEAN ANOMALY
      ANOM=SETP(11)*T+SETP(10)
      INT=ANOM
      IF(ANOM.LT.0.D0) INT=INT-1
      ANOMS=TWOPI*(ANOM-DBLE(INT))
C SOLVE KEPLER'S EQN FOR ECCENTRIC ANOMALY
      QQ1=ANOMS+SETP(9)*(DSIN(ANOMS)+SETP(9)*DSIN(2.D0*ANOMS)*0.5D0)
CD     WRITE(26,*) ANOMS,SETP(9),QQ1
  150 QQ2=(ANOMS-QQ1+SETP(9)*DSIN(QQ1))/(1.D0-SETP(9)*DCOS(QQ1))
      QQ1=QQ1+QQ2
      IF(DABS(QQ2).GT.5D-15) GOTO 150
C SAVE RESULTS OF SOLUTION AS SINE,COSINE
      SETP(13)=DSIN(QQ1)
      SETP(14)=DCOS(QQ1)
C GET RADIAL DISTANCE
      RY=SETP(8)*(1.D0-SETP(9)*SETP(14))
C COMPUTE VECTOR IN ORBIT PLANE COORDINATES
      CONTINUE
      YBAR(1)=SETP(8)*(SETP(14)-SETP(9))
      YBAR(2)=SETP(16)*SETP(13)
      YBAR(3)=SETP(17)*SETP(13)/RY
      YBAR(4)=SETP(18)*SETP(14)/RY
      IF(NV.GT.0) GOTO 300
C COMPUTE ONLY POSITION
      IYB=1
C ROTATE VECTOR TO STANDARD FRAME
      CALL MMPLY(SETP,YBAR(IYB),Y,3,2,1)
C  250 CALL PRODCT(SETP,YBAR(IYB),Y,3,2,1)
      RETURN
C COMPUTE POSITION AND VELOCITY
C GET POSITION AND VELOCITY
  300 CALL MMPLY(SETP,YBAR(1),Y(1),3,2,1)
      CALL MMPLY(SETP,YBAR(3),Y(4),3,2,1)
C  300 CALL PRODCT(SETP,YBAR(1),Y(1),3,2,1)
C      CALL PRODCT(SETP,YBAR(3),Y(4),3,2,1)
      IF(NV.EQ.1) GOTO 400
      CF=-SETP(12)/RY**3
      QQ1=T*1.5D0
      QQ2=CF*QQ1
      QQ4=SETP(13)/SETP(26)
      QQ5=SETP(8)*SETP(14)/RY
      QQ6=CF*QQ4
      QQ7=-YBAR(4)*SETP(25)
      DO 330 I=1,3
      J=I+3
C PARTIALS W.R.T. A
C UNITS (UNITLESS -- YB)
      DY(I,1)=(Y(I)-QQ1*Y(J))/SETP(8)
C UNITS (1/SEC -- YB)
      DY(J,1) = (-0.5D0*Y(J)-QQ2*Y(I))/SETP(8)
C PARTIALS W.R.T. E (UNITLESS -- YB)
      DY(I,2)=QQ4*Y(J)-SETP(I+18)-SETP(I+21)*YBAR(2)
      DY(J,2)=QQ5*Y(J)+QQ6*Y(I)+QQ7*SETP(J)
  330 CONTINUE
      IF(NV.LE.3) GOTO 400
C PARTIALS W.R.T. INC (RADIANS)
      DY(1,3)= SETP(28)*Y(3)
      DY(2,3)=-SETP(29)*Y(3)
      DY(3,3)= SETP(30)*YBAR(1)+SETP(31)*YBAR(2)
      DY(4,3)= SETP(28)*Y(6)
      DY(5,3)=-SETP(29)*Y(6)
      DY(6,3)= SETP(30)*YBAR(3)+SETP(31)*YBAR(4)
      IF(NV.LE.4) GOTO 400
C PARTIALS W.R.T. ASC (RADIANS)
      DY(1,4)=-Y(2)
      DY(2,4)= Y(1)
      DY(4,4)=-Y(5)
      DY(5,4)= Y(4)
      IF(NV.LE.5) GOTO 400
      QQ8=CF/SETP(26)
      DO 340 I=1,3
      J=I+3
C PARTIALS W.R.T. PER (RADIANS MEASURED FROM NODE)
      DY(I,5)=SETP(J)*YBAR(1)-SETP(I)*YBAR(2)
      DY(J,5)=SETP(J)*YBAR(3)-SETP(I)*YBAR(4)
C PARTIALS W.R.T. ANOM0 (RADIANS FROM PERIAPSE)
      DY(I,6)=Y(J)/SETP(26)
      DY(J,6)=Y(I)*QQ8
  340 CONTINUE
  400 RETURN
      END
