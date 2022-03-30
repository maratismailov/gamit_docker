      SUBROUTINE SMOOFT(Y,N,PTS)

c     Smooths an array Y of length N with a window whose full width
c     is of order PTS neigboring points, a user supplied value.
c     Y is returned modified.

c     routine from numerical recipes
c     converted to double precision

CKF      IMPLICIT REAL*8 (A-H,O-Z)
Ckf      PARAMETER(MMAX=4096)
ckf      DIMENSION Y(MMAX)
      real*8 y(*),pts,y1,yn,rn1,fac,const
      integer mo2,j,k,m,n,nmin

      M=2
      NMIN=N+2.*PTS
1     IF(M.LT.NMIN)THEN
        M=2*M
      GO TO 1
      ENDIF
      CONST=(PTS/M)**2
      Y1=Y(1)
      YN=Y(N)
      RN1=1./(N-1.)
      DO 11 J=1,N
        Y(J)=Y(J)-RN1*(Y1*(N-J)+YN*(J-1))
11    CONTINUE
c     zero padding here
      IF(N+1.LE.M)THEN
        DO 12 J=N+1,M
          Y(J)=0.
12      CONTINUE
      ENDIF
      MO2=M/2
      CALL REALFT(Y,MO2,1)
      Y(1)=Y(1)/MO2
      FAC=1.
      DO 13 J=1,MO2-1
        K=2*J+1
        IF(FAC.NE.0.)THEN
CKF       FAC=AMAX1(0.,(1.-CONST*J**2)/MO2)
          FAC=DMAX1(0.d0,(1.d0-CONST*J**2)/MO2)
          Y(K)=FAC*Y(K)
          Y(K+1)=FAC*Y(K+1)
        ELSE
          Y(K)=0.
          Y(K+1)=0.
        ENDIF
13    CONTINUE
CKF   FAC=AMAX1(0.,(1.-0.25*PTS**2)/MO2)
      FAC=DMAX1(0.0d0,(1.0d0-0.25*PTS**2)/MO2)
      Y(2)=FAC*Y(2)
      CALL REALFT(Y,MO2,-1)
      DO 14 J=1,N
        Y(J)=RN1*(Y1*(N-J)+YN*(J-1))+Y(J)
14    CONTINUE
      RETURN
      END


      SUBROUTINE REALFT(DATA,N,ISIGN)

c     routine from numerical recipes
c     converted to double precision

ckf      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,wis,wrs,c2,c1,h1i,h2i,h1r,h2r
      real*8 DATA(*)
      integer isign,i,i1,i2,i3,n2p3,n,i4

      THETA=6.28318530717959D0/2.0D0/DBLE(N)
      WR=1.0D0
      WI=0.0D0
      C1=0.5
      IF (ISIGN.EQ.1) THEN
        C2=-0.5
        CALL FOUR1(DATA,N,+1)
        DATA(2*N+1)=DATA(1)
        DATA(2*N+2)=DATA(2)
      ELSE
        C2=0.5
        THETA=-THETA
        DATA(2*N+1)=DATA(2)
        DATA(2*N+2)=0.0
        DATA(2)=0.0
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      N2P3=2*N+3
      DO 11 I=1,N/2+1
        I1=2*I-1
        I2=I1+1
        I3=N2P3-I2
        I4=I3+1
        WRS=SNGL(WR)
        WIS=SNGL(WI)
        H1R=C1*(DATA(I1)+DATA(I3))
        H1I=C1*(DATA(I2)-DATA(I4))
        H2R=-C2*(DATA(I2)+DATA(I4))
        H2I=C2*(DATA(I1)-DATA(I3))
        DATA(I1)=H1R+WRS*H2R-WIS*H2I
        DATA(I2)=H1I+WRS*H2I+WIS*H2R
        DATA(I3)=H1R-WRS*H2R+WIS*H2I
        DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
11    CONTINUE
      IF (ISIGN.EQ.1) THEN
        DATA(2)=DATA(2*N+1)
      ELSE
        CALL FOUR1(DATA,N,-1)
      ENDIF
      RETURN
      END

      SUBROUTINE FOUR1(DATA,NN,ISIGN)

c     routine from numerical recipes
c     converted to double precision

ckf      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,tempi,tempr
      real*8 DATA(*)
      integer isign,istep,i,j,m,n,mmax,nn
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END
