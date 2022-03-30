      SUBROUTINE POLYFT(NOBS,X,T,T0,W,IORDER,XM,RES,SIG,KEY)
C
C     Polynomial fit for a time series
C     ----     DND  890524
C
C  INPUT:
C     NOBS  :  Number of data points
C     X     :  observation data
C     T     :  Time index of observation data
C     T0    :  center time tag for polynomial fitting
C     W     :  Weight matrix (the simplest case --- diagonal)
C     IORDER:  order of polynomials (minimum=0,maximum=29)
C              number of coefficients = IORDER + 1
C     KEY   :  control variable
C              0. coefficients only
C              1. smoothed curve and coefficients
C              2. smoothed curve, residuals and coefficients
C
C  OUTPUT:
C     XM    :  filtered data
C     RES   :  residuals
C     SIG   :  standard deviation for this fitting (KEY > 0)
C     W     :  coefficients
C
      implicit none
c      IMPLICIT REAL*8(A-H,O-Z)
c      IMPLICIT INTEGER*4(I-N)

      integer*4 index,index0,index1,it,icol,ierr,nobs,icoeff,irow,iorder
     .        , inxold,key,i1,i2,j1,i,j,jelf
C
      real*8    aa,bb,cc,dt,rr,xm,t,w,x,sig,res,t0,small

      INTEGER MAXCOF
      PARAMETER (MAXCOF=30)
      DIMENSION X(*),T(*),XM(*),RES(*),W(*)
      DIMENSION AA(MAXCOF*(MAXCOF+1)/2),BB(MAXCOF),CC(MAXCOF)
      SMALL=1.0D-10
C
C     KEEP THIS ROUTINE FROM BLOWING UP
      IF (IORDER + 1 .GT. MAXCOF) THEN
         call report_stat('STATUS','CVIEW','polyft',iorder,
     .   'Error, iorder is out of bounds. IORDER=',0)
      ENDIF
      IF (KEY.LT.0.OR.KEY.GT.3) THEN
         call report_stat('STATUS','CVIEW','polyft',key,
     .   'Error, key is out of bounds. KEY=',0)
      ENDIF
C
C---- Time offset moveout
C
      DO 5 IT=1,NOBS
         T(IT)=T(IT)-T0
 5    CONTINUE
C
C---- Form the normal matrix
C
      ICOEFF=IORDER+1
      INXOLD=0
      DO 10 IROW=1,ICOEFF
c        replace by subroutine call to avoid HP compiler bug  tah/rwk 980121
c         INDEX0=IROW*(IROW-1)/2
         index0 = jelf(irow)
         I2=(IROW-2)*2
         DO 20 ICOL=1,IROW
C----       Index
            INDEX=INDEX0+ICOL
            AA(INDEX)=0.0D0
C----       Special case
            IF (IROW.EQ.1) THEN
               DO 140 IT=1,NOBS
               AA(INDEX)=AA(INDEX)+W(IT)
 140           CONTINUE
               GOTO 20
            ENDIF
C----       Avoid redundant calculations
            I1=IROW+ICOL-2
            IF (IROW.GT.2.AND.I1.LE.I2) THEN
               INDEX1=INXOLD+ICOL+1
               AA(INDEX)=AA(INDEX1)
               GOTO 20
            ENDIF
C----       Epoch loop
            DO 30 IT=1,NOBS
               DT=T(IT)
               AA(INDEX)=AA(INDEX)+DT**I1*W(IT)
 30         CONTINUE
 20      CONTINUE
         INXOLD=INDEX0
C----    Right hand term
         I1=IROW-1
         BB(IROW)=0.0D0
         DO 40 IT=1,NOBS
            IF (IROW.EQ.1) THEN
               BB(IROW)=BB(IROW)+X(IT)*W(IT)
               GOTO 40
            ENDIF
            DT=T(IT)
            BB(IROW)=BB(IROW)+X(IT)*DT**I1*W(IT)
 40      CONTINUE
 10   CONTINUE
C
C---- Get the inverse of the normal matrix
C
      CALL GJ(AA,CC,ICOEFF,IERR)

c     if not positive definite, please tell us.
      IF (IERR.EQ.-1) then
         key = -key
         return
      endif
C
C---- Get the coefficients
C
      DO 50 I=1,ICOEFF
         CC(I)=0.0D0
         DO 60 J=1,ICOEFF
            IF (I.GE.J) THEN
c              replace by subroutine call to avoid HP compiler bug
c               INDEX=I*(I-1)/2+J    
                index = jelf(i) + j
            ELSE
c               INDEX=J*(J-1)/2+I 
                index = jelf(j) + i
            ENDIF
            CC(I)=CC(I)+AA(INDEX)*BB(J)
 60      CONTINUE
 50   CONTINUE
C
C---- Different options
C
         DO 110 I=1,ICOEFF
              W(I)=CC(I)
 110     CONTINUE
         IF (KEY.EQ.0) GOTO 130
C
C---- Calculate the polynomials fitted time series
C
      SIG=0.0D0
      DO 80 IT=1,NOBS
         XM(IT)=0.0D0
         DT=T(IT)
         IF (DABS(DT).LT.SMALL) THEN
            XM(IT)=CC(1)
            GOTO 78
         ENDIF
         DO 70 J=1,ICOEFF
            J1=J-1
            XM(IT)=XM(IT)+CC(J)*DT**J1
 70      CONTINUE
C78      RR=XM(IT)-X(IT)
c        residual = observed - calculated
 78      RR=X(IT)-XM(IT)
         SIG=SIG+RR**2
         IF (KEY.LE.1) GOTO 80
         RES(IT)=RR
 80   CONTINUE
      SIG=DSQRT(SIG/DBLE(NOBS-1))
C
C---- Change back the time offset
C
 130  DO 90 IT=1,NOBS
         T(IT)=T(IT)+T0
 90   CONTINUE
C
C
      END
