       SUBROUTINE ESTREF( NUMSIT,SITNAM,X1,X2,X2N,OMC,ERROR,DERIV)
C
C     Determine transformation parameters for a set of station
C     coordinates given in two different reference frames.  The
C     estimated parameters give the transformation from system 2
C     to system 1:  X1 = T + (1+SCALE)*R*X2 where T is a vector
C     of translation parameters, scale is a scalar, and R is the
C     linearized rotation matrix for rotations about the 1,2, & 3 axes.
C     The units of X and T are meters;the angles R1, R2, R2 are in arcs
C
C     R. King  3 February 1987
C     M. Murray 13 Jan 88 : non-linear inversion added
C
C
      implicit none
      include '../includes/tform.h'

      integer*4 ndf,iprint,numsit,mdim,iter,iset,i,j,l,m,mm,n,k
      real*8 x1,x2,x2n,omc,error,deriv,t,r,resid,adj,sigma
     .     , omc0,x10,x20,cmat,convds,eps,scale,chisqr,sigfac
      character*16 sitnam
C
      DIMENSION X1(3,NSDIM),X2(3,NSDIM),X2N(3,NSDIM)
     1        , T(3),R(3),SITNAM(NSDIM),RESID(3),SIGFAC(7)
     2        , OMC(NDIM),ERROR(NDIM),DERIV(NDIM,7),ADJ(7),SIGMA(7)
     3        , OMC0(NDIM),X10(3,NSDIM),X20(3,NSDIM),CMAT(7,7)
      DATA SIGFAC/1.D2,1.D2,1.D2,1.D3,1.D3,1.D3,1.D9/
      data iprint/6/
C
C     Radians to arcsec
      CONVDS= DATAN(1.D0)*4.D0/180.D0/3600.D0
C     Set convergence criterion = 1 mm (assumes ERROR= 1 m)
      EPS= .001
      MDIM= 7
      M= 7
      ITER= 0
C
      WRITE(ISCRN,10)
   10 FORMAT(/,1X,'Enter the transformation parameter set to be'
     1           ,' adjusted:',//,
     2 '   1 = scale only',/,
     3 '   2 = orientation only',/,
     4 '   3 = translation and orientation     ',/,
     5 '   4 = translation, orientation, and scale',/,
     6 '   5 = scale and orientation only')
      READ(ITERM,*) ISET
      MM=0
      IF( ISET.EQ.1 ) MM=1
      IF( ISET.EQ.2 ) MM=3
      IF( ISET.EQ.3 ) MM=6
      IF( ISET.EQ.4 ) MM=7
      IF( ISET.EQ.5 ) MM=4
C
C
C       Initial values of parameters
C
      DO 15 I=1,3
      T(I)= 0.
   15 R(I)= 0.
      SCALE= 0.
C     Move X2 coords into X2N array to get initial o-c correct
      DO 18 J=1,NUMSIT
      DO 18 I=1,3
      X10(I,J) = X1(I,J)
      X20(I,J) = X2(I,J)
   18 OMC0((J-1)*3+I) = X2(I,J) - X1(I,J)
C
C       Form observations and compute derivatives
C
   20 K= 0
C
C        Compute the new System 2 coordinates
C
      DO 70 J=1,NUMSIT
      CALL TF2REF( X2(1,J),X2N(1,J),SCALE,T,R )
   70 CONTINUE
C

      DO 50 J=1,NUMSIT
      DO 50 I=1,3
      K= K + 1
      OMC(K)= X2N(I,J) - X2(I,J) + OMC0(K)
      ERROR(K)= 1.D0
      DO 25 L=1,7
   25 DERIV(K,L)= 0.D0

C       Parameters 1-3 are translation - T(I)
      IF( ISET.EQ.3 .OR. ISET.EQ.4 ) THEN
       IF( I.EQ.1 ) DERIV(K,1)= 1.D0
       IF( I.EQ.2 ) DERIV(K,2)= 1.D0
       IF( I.EQ.3 ) DERIV(K,3)= 1.D0
      ENDIF

C       Parameters 4-6 are orientation - R(I)
      IF( ISET.NE.1 ) THEN
       IF( I.EQ.1 ) DERIV(K,5)= -X2(3,J) * (1.D0 + SCALE) * CONVDS
       IF( I.EQ.1 ) DERIV(K,6)=  X2(2,J) * (1.D0 + SCALE) * CONVDS
       IF( I.EQ.2 ) DERIV(K,4)=  X2(3,J) * (1.D0 + SCALE) * CONVDS
       IF( I.EQ.2 ) DERIV(K,6)= -X2(1,J) * (1.D0 + SCALE) * CONVDS
       IF( I.EQ.3 ) DERIV(K,4)= -X2(2,J) * (1.D0 + SCALE) * CONVDS
       IF( I.EQ.3 ) DERIV(K,5)=  X2(1,J) * (1.D0 + SCALE) * CONVDS
      ENDIF
C
C       Parameter 7 is scale
      IF( ISET.EQ.1 .OR. ISET.EQ.4 .OR. ISET.EQ.5 ) THEN
       IF( I.EQ.1 )
     1  DERIV(K,7)= X2(1,J) + (   R(3)*X2(2,J) - R(2)*X2(3,J) ) * CONVDS
       IF( I.EQ.2 )
     1  DERIV(K,7)= X2(2,J) + ( - R(3)*X2(1,J) + R(1)*X2(3,J) ) * CONVDS
       IF( I.EQ.3 )
     1  DERIV(K,7)= X2(3,J) + (   R(2)*X2(1,J) - R(1)*X2(2,J) ) * CONVDS
      ENDIF
C
C  Assume rotations and scale are small, observation partials
C  become rather simple :  B = [ I | -I ]
C
   50 CONTINUE
      N=K
C
C       Perform least squares estimation
C
      CALL ESTIMT( MDIM,N,M,OMC,ERROR,DERIV,ADJ,SIGMA,CMAT,ITER)
C
C
C       Update the parameters
C
      DO 60 I=1,3
      T(I) = T(I) - ADJ(I)
   60 R(I) = R(I) - ADJ(I+3)
      SCALE= SCALE - ADJ(7)
C
C
C       See if the solution has converged
C
      DO 80 I=1,M
   80 IF( DABS(ADJ(I)/SIGMA(I) ).GT.EPS ) GOTO 90
      GOTO 100
   90 ITER= ITER + 1
      IF( ITER.GE.10 ) GOTO 990
C
C  Update observations (only on System 2, due to the B-matrix
C  simplifications alluded to above).
C
      K = 0
      DO J = 1, NUMSIT
         DO I = 1, 3
            K = K + 1
            X2N(I,J) = 0.D0
            DO L = 1, 7
               X2N(I,J) = X2N(I,J) - DERIV(K,L)*ADJ(L)
            ENDDO
            X2(I,J) = X2(I,J) - (X2N(I,J)+OMC(K))/2.D0
         ENDDO
      ENDDO
C
      GOTO 20
C
C----------------------------------------------------------------------
C
  100 CONTINUE
C
C        Compute the new System 2 coordinates
C
      DO 71 J=1,NUMSIT
   71 CALL TF2REF( X20(1,J),X2N(1,J),SCALE,T,R )
C
C       Compute and print the postfit residuals
C
      WRITE(IPRNT,110)
  110 FORMAT(/,1X,'Postfit residuals:',/)
      K= 0
      CHISQR= 0.D0
      DO 140 J=1,NUMSIT
      DO 120 I=1,3
      K= K + 1
      RESID(I)= X10(I,J) - X2N(I,J)
      CHISQR= CHISQR + (RESID(I)/ERROR(K))**2
  120 CONTINUE
      WRITE(IPRNT,130) SITNAM(J),(RESID(I),I=1,3)
  130 FORMAT(/,1X,A16,3X,3F9.4)
  140 CONTINUE
      NDF= N - MM
      CHISQR= DSQRT(CHISQR/NDF)
      WRITE(IPRNT,150) CHISQR, NDF
  150 FORMAT(/,1X,'Sqrt( chi**2/df ) = ',E12.5,'   df=',I3)
C
C
C        Print the adjusted parameters and their uncertainties
C
C     Scale the uncertainties by sqrt(chi**2) pdf
      DO 160 I=1,M
  160 SIGMA(I)= SIGMA(I)*CHISQR
C
      DO I = 1, 3
         T(I) = T(I)*SIGFAC(I)
         R(I) = R(I)*SIGFAC(I+3)
         SIGMA(I) = SIGMA(I)*SIGFAC(I)
         SIGMA(I+3) = SIGMA(I+3)*SIGFAC(I+3)
      ENDDO
      SCALE = SCALE*SIGFAC(7)
      SIGMA(7) = SIGMA(7)*SIGFAC(7)
C
C        Print out a Covariance/Correlation matrix
C        Make upper triangle(including diagonal) covariance
C        Make lower triangle correlation
C
      DO I = 1, M
      DO J = 1, M
         IF (J.GE.I) THEN
            CMAT(I,J) = CMAT(I,J)*SIGMA(I)*SIGMA(J)
         ENDIF
      ENDDO
      ENDDO
C
      WRITE(IPRNT,170) T(1),SIGMA(1),T(2),SIGMA(2),T(3),SIGMA(3)
     1              , R(1),SIGMA(4),R(2),SIGMA(5),R(3),SIGMA(6)
     2              , SCALE,SIGMA(7)
  170 FORMAT(/,1X,'Parameter estimates (std dev):',//
     1        ,1X,'   T =',3(F9.4,3X,'(',F7.4,')',5X),//
     2        ,1X,'   R =',3(F9.5,1X,'(',F9.5,')',2X),//
     3        ,1X,'SCALE=',1X,E11.4,1X,'(',E11.4,')',// )
C
      DO I = 1, M
         WRITE(IPRINT,95) (CMAT(I,J),J = 1, M)
      ENDDO
   95 FORMAT(4X,7F9.3)
C
C
C         Print the new System 2 coordinates
C
      WRITE(IPRNT,180)
  180 FORMAT(/,1X,'Set 2 coordinates transformed to System 1:')
      WRITE(IPRNT,190) (SITNAM(J),(X2N(I,J),I=1,3),J=1,NUMSIT)
  190 FORMAT(/,1X,A16,3F13.3)
C
      GOTO 999
C
  990 WRITE(ISCRN,991)
  991 FORMAT(1X,'  More than 10 iterations in ESTREF, stop')
C
  999 RETURN
      END
