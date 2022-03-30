      SUBROUTINE ESTIMT( MDIM,N,M,OMC,ERROR,DERIV,ADJ,SIGMA,CMAT,ITER)

C  Routine to perform least squares estimation for program TFORM2
C  R.W. King  3 February 87 (from ESTIM 86.9.23)

      implicit none

      include '../includes/tform.h'

      integer*4 iuse,ierinv,mdim,iprint,i,j,k,n,m,iter  

      real*8 omc,error,deriv,adj,sigma,cmat,rhs,pvrow,pvcol,amat
     .     , ai,aj,pvrwb

      DIMENSION OMC(NDIM),ERROR(NDIM),DERIV(NDIM,7),PVROW(7),
     1          PVCOL(7),IUSE(7)
     1        , RHS(7),CMAT(7,7),SIGMA(7),AMAT(28),ADJ(7)
C

      data iprint/6/
C
C  Check diminstions
      IF( MDIM.NE.7 .OR. M.GT.7 ) GOTO 992
C
C
C  Loop over observations, calculating partials and forming normal eqns
C
      DO 25 I=1,M
      RHS(I)= 0.D0
      DO 25 J=1,M
   25 CMAT(I,J)= 0.D0
C
      DO  50 K=1,N
C
C
C -- DEBUG PRINTOUT
C     WRITE(IPRINT,21) K,OMC(K),(DERIV(K,L),L=1,M)
C  21 FORMAT(1X,'OBS #, OMC, PARTIALS=',I2,4D15.6)
C
C     Form RHS and design matrix
C
      DO 40 I=1,M
C
      AI= DERIV(K,I)
      RHS(I)= RHS(I) + ERROR(K)**(-2)*AI*OMC(K)
C
      DO 30 J=1,M
      AJ= DERIV(K,J)
      CMAT(I,J)= CMAT(I,J) + ERROR(K)**(-2)*AI*AJ
   30 CONTINUE
C
   40 CONTINUE
C
   50 CONTINUE
C     End loop on observations
C
C-- Debug printout
C     WRITE(IPRINT,51) RHS,CMAT
C  51 FORMAT(/,1X,'RHS=',3D20.9,//,1X,'CMAT='/,3(1X,3D20.10,/))
C
C     If a diagonal element is zero (not-adjusted), replace by 1.
C
      DO 55 I=1,M
   55 IF( CMAT(I,I).EQ.0.D0 ) CMAT(I,I)= 1.D0
C
C  Solve normal equations
C
C     Store design matrix as upper triangular part for inversion routine
C
      DO 60 I=1,M
      DO 60 J=1,I
   60 AMAT(J+I*(I-1)/2) = CMAT(I,J)
C
C     Invert the A-matrix
      CALL SYMINV( AMAT,RHS,M,1,PVROW,PVCOL,IUSE,PVRWB,IERINV)
C
C     Restore the square form of the design matrix
C
      DO 70 I=1,M
      DO 70 J=1,I
      CMAT(I,J)= AMAT(J+I*(I-1)/2)
   70 IF( I.NE.J ) CMAT(J,I)= CMAT(I,J)
C
C-- Debug printout
C     WRITE(IPRINT,51) RHS,CMAT
C
C
C  Compute adjustments, standard deviations, and oorrelation matrix
      DO 80 I=1,M
      ADJ(I)= RHS(I)
   80 SIGMA(I)= DSQRT( CMAT(I,I) )
      DO 90 J=1,M
      DO 90 I=1,M
   90 CMAT(I,J)= CMAT(I,J)/SIGMA(I)/SIGMA(J)
      IF( ITER.LE.0) WRITE(IPRINT,95) CMAT
   95 FORMAT(/,1X,'Correlation Matrix: ',//,7(4X,7F8.4,/))
C
C
      RETURN
C
C
C  Abnormal terminations
C
  992 WRITE(ITERM,993) MDIM,M
  993 FORMAT(/,1X,'MDIM=',I3,' gt 7  OR M=',I3,'gt 7, Stop in ESTIMT')
      STOP
C
      END

