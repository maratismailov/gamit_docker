C***** SAMPLE PROGRAM FOR SURFACE DEFORMATION                           00050001
      IMPLICIT REAL*8  (A-H,O-Z)                                        00060001
      DIMENSION  D1(3),D2(3),D3(3)                                      00070001
      DATA  D1 / 1.D0, 0.D0, 0.D0 /                                     00080001
      DATA  D2 / 0.D0, 1.D0, 0.D0 /                                     00090001
      DATA  D3 / 0.D0, 0.D0, 1.D0 /                                     00100001
      DATA  F0, F1, EPS / 0.D0, 1.D0, 1.D-3 /                           00110001
      RAD=6.283185307179586D0/3.6D2                                     00120001
      ALP=0.5D0                                                         00130001
C*****                                                                  00140001
    1 READ(5,500,END=90)  KTYPE,X,Y,DEP,DIP,AL,AW                       00150001
  500 FORMAT(I2,6F5.0)                                                  00160001
      IF(KTYPE.EQ.1)  WRITE(6,601)  X,Y,DEP,DIP                         00170001
      IF(KTYPE.EQ.2)  WRITE(6,602)  X,Y,DEP,DIP,AL,AW                   00180001
  601 FORMAT(//' *** POINT ***   X,Y,DEP,DIP =',4F8.2/)                 00190001
  602 FORMAT(//' *** FINITE ***   X,Y,DEP,DIP,AL,AW =',6F8.2/)          00200001
C*****                                                                  00210001
      SD=DSIN( DIP*RAD )                                                00220001
      CD=DCOS( DIP*RAD )                                                00230001
      IF(DABS(CD).GT.EPS)  GO TO 10                                     00240001
      CD=F0                                                             00250001
      IF(SD.GT.F0)  SD= F1                                              00260001
      IF(SD.LT.F0)  SD=-F1                                              00270001
C*****                                                                  00280001
C***** IN CASE OF POINT SOURCE                                          00290001
C*****                                                                  00300001
   10 IF(KTYPE.EQ.2)  GO TO 20                                          00310001
      DO 1111  I=1,3                                                    00320001
      CALL SPOINT(ALP,X,Y,DEP,SD,CD,D1(I),D2(I),D3(I),                  00330001
     *               U1,U2,U3,U11,U12,U21,U22,U31,U32)                  00340001
 1111 WRITE(6,610)   U1,U2,U3,U11,U12,U21,U22,U31,U32                   00350001
      GO TO 1                                                           00360001
C*****                                                                  00370001
C***** IN CASE OF FINITE FAULT                                          00380001
C*****                                                                  00390001
 20   DO 2222  I=1,3                                                    00400001
      CALL SRECTF(ALP,X,Y,DEP,AL,AW,SD,CD,D1(I),D2(I),D3(I),            00410001
     *               U1,U2,U3,U11,U12,U21,U22,U31,U32)                  00420001

 2222 WRITE(6,610)   U1,U2,U3,U11,U12,U21,U22,U31,U32                   00430001
      GO TO 1                                                           00440001
C*****                                                                  00450001
  610 FORMAT(  1H ,1P9E14.3)                                            00460001
   90 STOP                                                              00470001
      END                                                               00480001
      SUBROUTINE  SPOINT(ALP,X,Y,D,SD,CD,DISL1,DISL2,DISL3,             00490001
     *                   U1,U2,U3,U11,U12,U21,U22,U31,U32)              00500001
C*****                                                                  00510001
C*****    SURFACE DISPLACEMENT,STRAIN,TILT DUE TO BURIED POINT SOURCE   00520001
C*****    IN A SEMIINFINITE MEDIUM     CODED BY  Y.OKADA ... JAN 1985   00530001
C*****                                                                  00540001
C***** INPUT                                                            00550001
C*****   ALP   : MEDIUM CONSTANT  MYU/(LAMDA+MYU)                       00560001
C*****   X,Y   : COORDINATE OF STATION                                  00570001
C*****   D     : SOURCE DEPTH                                           00580001
C*****   SD,CD : SIN,COS OF DIP-ANGLE                                   00590001
C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)00600001
C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION      00610001
C*****                                                                  00620001
C***** OUTPUT                                                           00630001
C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL / AREA )   00640001
C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /          00650001
C*****   U31,U32         : TILT                 UNIT OF X,Y,D /AREA )   00660001
C*****                                                                  00670001
      IMPLICIT REAL*8 (A-H,O-Z)                                         00680001
      DATA  F0,F1,F2,F3,F4,F5,F8,F9                                     00690001
     *      /0.D0, 1.D0, 2.D0, 3.D0, 4.D0, 5.D0, 8.D0, 9.D0/            00700001
      PI2=6.283185307179586D0                                           00710001
C*****                                                                  00720001
      P =Y*CD + D*SD                                                    00730001
      Q =Y*SD - D*CD                                                    00740001
      S =P*SD + Q*CD                                                    00750001
      X2=X*X                                                            00760001
      Y2=Y*Y                                                            00770001
      XY=X*Y                                                            00780001
      D2=D*D                                                            00790001
      R2=X2 + Y2 + D2                                                   00800001
      R =SQRT(R2)                                                       00810001
      R3=R *R2                                                          00820001
      R5=R3*R2                                                          00830001
      QR=F3*Q/R5                                                        00840001
      XR =F5*X2/R2                                                      00850001
      YR =F5*Y2/R2                                                      00860001
      XYR=F5*XY/R2                                                      00870001
      DR =F5*D /R2                                                      00880001
      RD =R + D                                                         00890001
      R12=F1/(R*RD*RD)                                                  00900001
      R32=R12*(F2*R + D)/ R2                                            00910001
      R33=R12*(F3*R + D)/(R2*RD)                                        00920001
      R53=R12*(F8*R2 + F9*R*D + F3*D2)/(R2*R2*RD)                       00930001
      R54=R12*(F5*R2 + F4*R*D +    D2)/R3*R12                           00940001
C*****                                                                  00950001
      A1= ALP*Y*(R12-X2*R33)                                            00960001
      A2= ALP*X*(R12-Y2*R33)                                            00970001
      A3= ALP*X/R3 - A2                                                 00980001
      A4=-ALP*XY*R32                                                    00990001
      A5= ALP*( F1/(R*RD) - X2*R32 )                                    01000001
      B1= ALP*(-F3*XY*R33      + F3*X2*XY*R54)                          01010001
      B2= ALP*( F1/R3 - F3*R12 + F3*X2*Y2*R54)                          01020001
      B3= ALP*( F1/R3 - F3*X2/R5) - B2                                  01030001
      B4=-ALP*F3*XY/R5 - B1                                             01040001
      C1=-ALP*Y*(R32 - X2*R53)                                          01050001
      C2=-ALP*X*(R32 - Y2*R53)                                          01060001
      C3=-ALP*F3*X*D/R5 - C2                                            01070001
C*****                                                                  01080001
      U1 =F0                                                            01090001
      U2 =F0                                                            01100001
      U3 =F0                                                            01110001
      U11=F0                                                            01120001
      U12=F0                                                            01130001
      U21=F0                                                            01140001
      U22=F0                                                            01150001
      U31=F0                                                            01160001
      U32=F0                                                            01170001
C**************************************                                 01180001
C*****                            *****                                 01190001
C*****  STRIKE-SLIP CONTRIBUTION  *****                                 01200001
C*****                            *****                                 01210001
C**************************************                                 01220001
      IF(DISL1.EQ.F0)  GO TO 200                                        01230001
      UN=DISL1/PI2                                                      01240001
      QRX=QR*X                                                          01250001
      FX=F3*X/R5*SD                                                     01260001
      U1 =U1 - UN*( QRX*X + A1*SD )                                     01270001
      U2 =U2 - UN*( QRX*Y + A2*SD )                                     01280001
      U3 =U3 - UN*( QRX*D + A4*SD )                                     01290001
      U11=U11- UN*( QRX* (F2-XR)        + B1*SD )                       01300001
      U12=U12- UN*(-QRX*XYR      + FX*X + B2*SD )                       01310001
      U21=U21- UN*( QR*Y*(F1-XR)        + B2*SD )                       01320001
      U22=U22- UN*( QRX *(F1-YR) + FX*Y + B4*SD )                       01330001
      U31=U31- UN*( QR*D*(F1-XR)        + C1*SD )                       01340001
      U32=U32- UN*(-QRX*DR*Y     + FX*D + C2*SD )                       01350001
C**************************************                                 01360001
C*****                            *****                                 01370001
C*****    DIP-SLIP CONTRIBUTION   *****                                 01380001
C*****                            *****                                 01390001
C**************************************                                 01400001
  200 IF(DISL2.EQ.F0)  GO TO 300                                        01410001
      UN=DISL2/PI2                                                      01420001
      SDCD=SD*CD                                                        01430001
      QRP=QR*P                                                          01440001
      FS=F3*S/R5                                                        01450001
      U1 =U1 - UN*( QRP*X - A3*SDCD )                                   01460001
      U2 =U2 - UN*( QRP*Y - A1*SDCD )                                   01470001
      U3 =U3 - UN*( QRP*D - A5*SDCD )                                   01480001
      U11=U11- UN*( QRP*(F1-XR)        - B3*SDCD )                      01490001
      U12=U12- UN*(-QRP*XYR     + FS*X - B1*SDCD )                      01500001
      U21=U21- UN*(-QRP*XYR            - B1*SDCD )                      01510001
      U22=U22- UN*( QRP*(F1-YR) + FS*Y - B2*SDCD )                      01520001
      U31=U31- UN*(-QRP*DR*X           - C3*SDCD )                      01530001
      U32=U32- UN*(-QRP*DR*Y    + FS*D - C1*SDCD )                      01540001
C****************************************                               01550001
C*****                              *****                               01560001
C*****  TENSILE-FAULT CONTRIBUTION  *****                               01570001
C*****                              *****                               01580001
C****************************************                               01590001
  300 IF(DISL3.EQ.F0)  GO TO 900                                        01600001
      UN=DISL3/PI2                                                      01610001
      SDSD=SD*SD                                                        01620001
      QRQ=QR*Q                                                          01630001
      FQ=F2*QR*SD                                                       01640001
      U1 =U1 + UN*( QRQ*X - A3*SDSD )                                   01650001
      U2 =U2 + UN*( QRQ*Y - A1*SDSD )                                   01660001
      U3 =U3 + UN*( QRQ*D - A5*SDSD )                                   01670001
      U11=U11+ UN*( QRQ*(F1-XR)        - B3*SDSD )                      01680001
      U12=U12+ UN*(-QRQ*XYR     + FQ*X - B1*SDSD )                      01690001
      U21=U21+ UN*(-QRQ*XYR            - B1*SDSD )                      01700001
      U22=U22+ UN*( QRQ*(F1-YR) + FQ*Y - B2*SDSD )                      01710001
      U31=U31+ UN*(-QRQ*DR*X           - C3*SDSD )                      01720001
      U32=U32+ UN*(-QRQ*DR*Y    + FQ*D - C1*SDSD )                      01730001
C*****                                                                  01740001
  900 RETURN                                                            01750001
      END                                                               01760001
      SUBROUTINE  SRECTF(ALP,X,Y,DEP,AL,AW,SD,CD,DISL1,DISL2,DISL3,     01770001
     *                   U1,U2,U3,U11,U12,U21,U22,U31,U32)              01780001
C*****                                                                  01790001
C*****   SURFACE DISPLACEMENTS,STRAINS AND TILTS DUE TO RECTANGULAR     01800001
C*****   FAULT IN A HALF-SPACE       CODED BY  Y.OKADA ... JAN 1985     01810001
C*****                                                                  01820001
C***** INPUT                                                            01830001
C*****   ALP   : MEDIUM CONSTANT  MYU/(LAMDA+MYU)                       01840001
C*****   X,Y   : COORDINATE OF STATION                                  01850001
C*****   DEP   : SOURCE DEPTH                                           01860001
C*****   AL,AW : LENGTH AND WIDTH OF FAULT                              01870001
C*****   SD,CD : SIN,COS OF DIP-ANGLE                                   01880001
C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)01890001
C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION      01900001
C*****                                                                  01910001
C***** OUTPUT                                                           01920001
C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL     )      01930001
C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /          01940001
C*****   U31,U32         : TILT                 UNIT OF X,Y,,,AW )      01950001
C*****                                                                  01960001
C***** SUBROUTINE USED...SRECTG                                         01970001
C*****                                                                  01980001
      IMPLICIT REAL*8 (A-H,O-Z)                                         01990001
      DIMENSION  U(9),DU(9)                                             02000001
      DATA  F0, F1 / 0.D0, 1.D0 /                                       02010001
C*****                                                                  02020001
      P = Y*CD + DEP*SD                                                 02030001
      Q = Y*SD - DEP*CD                                                 02040001
C*****                                                                  02050001
      DO 1111  I=1,9                                                    02060001
 1111 U(I)=F0                                                           02070001
C*****                                                                  02080001
      DO 5555  K=1,2                                                    02090001
       IF(K.EQ.1)  ET=P                                                 02100001
       IF(K.EQ.2)  ET=P-AW                                              02110001
       DO 4444  J=1,2                                                   02120001
        IF(J.EQ.1)  XI=X                                                02130001
        IF(J.EQ.2)  XI=X-AL                                             02140001
        JK=J+K                                                          02150001
        IF(JK.NE.3)  SIGN= F1                                           02160001
        IF(JK.EQ.3)  SIGN=-F1                                           02170001
        CALL SRECTG(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3,                02180001
     *           DU(1),DU(2),DU(3),DU(4),DU(5),DU(6),DU(7),DU(8),DU(9)) 02190001
        DO 3333  I=1,9                                                  02200001
         U(I)=U(I)+SIGN*DU(I)                                           02210001
 3333   CONTINUE                                                        02220001
 4444  CONTINUE                                                         02230001
 5555 CONTINUE                                                          02240001
      U1 =U(1)                                                          02250001
      U2 =U(2)                                                          02260001
      U3 =U(3)                                                          02270001
      U11=U(4)                                                          02280001
      U12=U(5)                                                          02290001
      U21=U(6)                                                          02300001
      U22=U(7)                                                          02310001
      U31=U(8)                                                          02320001
      U32=U(9)                                                          02330001
      RETURN                                                            02340001
      END                                                               02350001
      SUBROUTINE  SRECTG(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3,           02360001
     *                   U1,U2,U3,U11,U12,U21,U22,U31,U32)              02370001
C*****                                                                  02380001
C*****   INDEFINITE INTEGRAL OF SURFACE DISPLACEMENTS, STRAINS AND TILTS02390001
C*****   DUE TO FINITE FAULT IN A SEMIINFINITE MEDIUM                   02400001
C*****                                    CODED BY  Y.OKADA ... JAN 198502410001
C***** INPUT                                                            02420001
C*****   ALP     : MEDIUM CONSTANT  MYU/(LAMDA+MYU)                     02430001
C*****   XI,ET,Q : FAULT COORDINATE                                     02440001
C*****   SD,CD   : SIN,COS OF DIP-ANGLE                                 02450001
C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)02460001
C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION      02470001
C*****                                                                  02480001
C***** OUTPUT                                                           02490001
C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL    )       02500001
C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /          02510001
C*****   U31,U32         : TILT                 UNIT OF XI,ET,Q )       02520001
C*****                                                                  02530001
      IMPLICIT REAL*8 (A-H,O-Z)                                         02540001
      DATA  F0,F1,F2/ 0.D0, 1.D0, 2.D0 /                                02550001
      PI2=6.283185307179586D0                                           02560001
C*****                                                                  02570001
      XI2=XI*XI                                                         02580001
      ET2=ET*ET                                                         02590001
      Q2=Q*Q                                                            02600001
      R2=XI2+ET2+Q2                                                     02610001
      R =DSQRT(R2)                                                      02620001
      R3=R*R2                                                           02630001
      D =ET*SD-Q*CD                                                     02640001
      Y =ET*CD+Q*SD                                                     02650001
      RET=R+ET                                                          02660001
      IF(RET.LT.F0)  RET=F0                                             02670001
      RD =R+D                                                           02680001
      RRD=F1/(R*RD)                                                     02690001
C*****                                                                  02700001
      IF( Q .NE.F0)  TT = DATAN( XI*ET/(Q*R) )                          02710001
      IF( Q .EQ.F0)  TT = F0                                            02720001
      IF(RET.NE.F0)  RE = F1/RET                                        02730001
      IF(RET.EQ.F0)  RE = F0                                            02740001
      IF(RET.NE.F0)  DLE= DLOG(RET)                                     02750001
      IF(RET.EQ.F0)  DLE=-DLOG(R-ET)                                    02760001
      RRX=F1/(R*(R+XI))                                                 02770001
      RRE=RE/R                                                          02780001
      AXI=(F2*R+XI)*RRX*RRX/R                                           02790001
      AET=(F2*R+ET)*RRE*RRE/R                                           02800001
      IF(CD.EQ.F0)  GO TO 20                                            02810001
C*****                                                                  02820001
C***** INCLINED FAULT                                                   02830001
C*****                                                                  02840001
      TD=SD/CD                                                          02850001
      X =DSQRT(XI2+Q2)                                                  02860001
      IF(XI.EQ.F0)  A5=F0                                               02870001
      IF(XI.NE.F0)                                                      02880001
     *A5= ALP*F2/CD*DATAN( (ET*(X+Q*CD)+X*(R+X)*SD) / (XI*(R+X)*CD) )   02890001
      A4= ALP/CD*( DLOG(RD) - SD*DLE )                                  02900001
      A3= ALP*(Y/RD/CD - DLE) + TD*A4                                   02910001
      A1=-ALP/CD*XI/RD        - TD*A5                                   02920001
      C1= ALP/CD*XI*(RRD - SD*RRE)                                      02930001
      C3= ALP/CD*(Q*RRE - Y*RRD)                                        02940001
      B1= ALP/CD*(XI2*RRD - F1)/RD - TD*C3                              02950001
      B2= ALP/CD*XI*Y*RRD/RD       - TD*C1                              02960001
      GO TO 30                                                          02970001
C*****                                                                  02980001
C***** VERTICAL FAULT                                                   02990001
C*****                                                                  03000001
   20 RD2=RD*RD                                                         03010001
      A1=-ALP/F2*XI*Q/RD2                                               03020001
      A3= ALP/F2*( ET/RD + Y*Q/RD2 - DLE )                              03030001
      A4=-ALP*Q/RD                                                      03040001
      A5=-ALP*XI*SD/RD                                                  03050001
      B1= ALP/F2*  Q  /RD2*(F2*XI2*RRD - F1)                            03060001
      B2= ALP/F2*XI*SD/RD2*(F2*Q2 *RRD - F1)                            03070001
      C1= ALP*XI*Q*RRD/RD                                               03080001
      C3= ALP*SD/RD*(XI2*RRD - F1)                                      03090001
C*****                                                                  03100001
   30 A2=-ALP*DLE - A3                                                  03110001
      B3=-ALP*XI*RRE - B2                                               03120001
      B4=-ALP*( CD/R + Q*SD*RRE ) - B1                                  03130001
      C2= ALP*(-SD/R + Q*CD*RRE ) - C3                                  03140001
C*****                                                                  03150001
      U1 =F0                                                            03160001
      U2 =F0                                                            03170001
      U3 =F0                                                            03180001
      U11=F0                                                            03190001
      U12=F0                                                            03200001
      U21=F0                                                            03210001
      U22=F0                                                            03220001
      U31=F0                                                            03230001
      U32=F0                                                            03240001
C**************************************                                 03250001
C*****                            *****                                 03260001
C*****  STRIKE-SLIP CONTRIBUTION  *****                                 03270001
C*****                            *****                                 03280001
C**************************************                                 03290001
      IF(DISL1.EQ.F0)  GO TO 200                                        03300001
      UN=DISL1/PI2                                                      03310001
      REQ=RRE*Q                                                         03320001
      U1 =U1 - UN*( REQ*XI +   TT    + A1*SD )                          03330001
      U2 =U2 - UN*( REQ*Y  + Q*CD*RE + A2*SD )                          03340001
      U3 =U3 - UN*( REQ*D  + Q*SD*RE + A4*SD )                          03350001
      U11=U11+ UN*( XI2*Q*AET - B1*SD )                                 03360001
      U12=U12+ UN*( XI2*XI*( D/(ET2+Q2)/R3 - AET*SD ) - B2*SD )         03370001
      U21=U21+ UN*( XI*Q/R3*CD + (XI*Q2*AET - B2)*SD )                  03380001
      U22=U22+ UN*( Y *Q/R3*CD + (Q*SD*(Q2*AET-F2*RRE)                  03390001
     *                            -(XI2+ET2)/R3*CD - B4)*SD )           03400001
      U31=U31+ UN*(-XI*Q2*AET*CD + (XI*Q/R3 - C1)*SD )                  03410001
      U32=U32+ UN*( D*Q/R3*CD + (XI2*Q*AET*CD - SD/R + Y*Q/R3 - C2)*SD )03420001
C**************************************                                 03430001
C*****                            *****                                 03440001
C*****    DIP-SLIP CONTRIBUTION   *****                                 03450001
C*****                            *****                                 03460001
C**************************************                                 03470001
  200 IF(DISL2.EQ.F0)  GO TO 300                                        03480001
      UN=DISL2/PI2                                                      03490001
      SDCD=SD*CD                                                        03500001
      U1 =U1 - UN*( Q/R             - A3*SDCD )                         03510001
      U2 =U2 - UN*( Y*Q*RRX + CD*TT - A1*SDCD )                         03520001
      U3 =U3 - UN*( D*Q*RRX + SD*TT - A5*SDCD )                         03530001
      U11=U11+ UN*( XI*Q/R3            + B3*SDCD )                      03540001
      U12=U12+ UN*( Y *Q/R3 - SD/R     + B1*SDCD )                      03550001
      U21=U21+ UN*( Y *Q/R3 + Q*CD*RRE + B1*SDCD )                      03560001
      U22=U22+ UN*( Y*Y*Q*AXI - (F2*Y*RRX + XI*CD*RRE)*SD + B2*SDCD )   03570001
      U31=U31+ UN*( D *Q/R3 + Q*SD*RRE + C3*SDCD )                      03580001
      U32=U32+ UN*( Y*D*Q*AXI - (F2*D*RRX + XI*SD*RRE)*SD + C1*SDCD )   03590001
C****************************************                               03600001
C*****                              *****                               03610001
C*****  TENSILE-FAULT CONTRIBUTION  *****                               03620001
C*****                              *****                               03630001
C****************************************                               03640001
  300 IF(DISL3.EQ.F0)  GO TO 900                                        03650001
      UN=DISL3/PI2                                                      03660001
      SDSD=SD*SD                                                        03670001
      U1 =U1 + UN*( Q2*RRE                       - A3*SDSD )            03680001
      U2 =U2 + UN*(-D*Q*RRX - SD*(XI*Q*RRE - TT) - A1*SDSD )            03690001
      U3 =U3 + UN*( Y*Q*RRX + CD*(XI*Q*RRE - TT) - A5*SDSD )            03700001
      U11=U11- UN*( XI*Q2*AET             + B3*SDSD )                   03710001
      U12=U12- UN*(-D*Q/R3 - XI2*Q*AET*SD + B1*SDSD )                   03720001
      U21=U21- UN*( Q2*(CD/R3 + Q*AET*SD) + B1*SDSD )                   03730001
      U22=U22- UN*((Y*CD-D*SD)*Q2*AXI - F2*Q*SD*CD*RRX                  03740001
     *                      - (XI*Q2*AET - B2)*SDSD )                   03750001
      U31=U31- UN*( Q2*(SD/R3 - Q*AET*CD) + C3*SDSD )                   03760001
      U32=U32- UN*((Y*SD+D*CD)*Q2*AXI + XI*Q2*AET*SD*CD                 03770001
     *                       - (F2*Q*RRX - C1)*SDSD )                   03780001
C*****                                                                  03790001
  900 RETURN                                                            03800001
      END                                                               03810001
