c---------------------------------------------------------------
      subroutine chi2_dis(n,alfa,up_bound,lo_bound)
c
c     Given number of freedom (n), and confidence level (alfa),
c     calculate the upper bound and lower bound of the chi square 
c     based on the probability distribution function
c                                   Dong 04/05/93
      integer i,j,i1,i2,itmax,jtmax,n,ft
      real alfa,up_bound,lo_bound,step,t1,t2,q,step2
      real gammq
      
      itmax = 1000
      jtmax = 1000
      step = 1.0
      t1 = 0.001-step
      ft = 0.5*n
      i1 = 0
      i2 = 0
c     upper bound
      do 20 i = 1,itmax
        t1 = t1+step
        q = gammq(ft,0.5*t1)
        if (q.gt.1.0-alfa) i1 = i
        if (q.le.1.0-alfa) i2 = i
        if (i2.gt.0) goto 10
 20   continue
 10   t1 = step*(i1-1)-0.099
      step = 0.1
      i1 = 0
      i2 = 0
      do 40 i = 1,itmax
        t1 = t1+step
        q = gammq(ft,0.5*t1)
        if (q.gt.1.0-alfa) i1 = i
        if (q.le.1.0-alfa) i2 = i
        if (i2.gt.0) then
          t2 = t1-step
          step2 = 0.01
          do 30 j = 1,jtmax
             t2 = t2+step2 
             if (q.le.1.0-alfa) then
               up_bound = t2
               goto 50
             endif
 30       continue
        endif
 40   continue
      up_bound = 99999.9
 50   step = 1.0
      if (n.le.6) step = 0.1
      if (n.le.3) step = 0.01
      t1 = 0.00001-step
      i1 = 0
      i2 = 0
c     lower bound
      do 60 i = 1,itmax
        t1 = t1+step
        q = gammq(ft,0.5*t1)
        if (q.gt.alfa) i1 = i
        if (q.le.alfa) i2 = i
        if (i2.gt.0) goto 70
 60   continue
 70   t1 = step*(i1-1)
      if (i1.gt.1) lo_bound = t1
      if (i1.le.1) t1 = 0.1d-5
      step = 0.1*step
      t1 = t1-step
      i1 = 0
      i2 = 0
      do 90 i = 1,itmax
        t1 = t1+step
        q = gammq(ft,0.5*t1)
        if (q.gt.alfa) i1 = i
        if (q.le.alfa) i2 = i
        if (i2.gt.0) then
          step2 = 0.1*step
          t2 = t1-step
          do 80 j = 1,jtmax
             t2 = t2+step2 
             if (q.le.alfa) then
               lo_bound = t2
               goto 100
             endif
 80       continue
        endif
 90   continue
      lo_bound = 0.0001  
 100  return
      end
c---------------------------------------------------------------
      FUNCTION GAMMP(A,X)
      REAL GAMSER,GAMMP,GLN,GAMMCF,A,X
      IF(X.LT.0..OR.A.LE.0.)PAUSE
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.-GAMMCF
      ENDIF
      RETURN
      END
c---------------------------------------------------------------
      FUNCTION GAMMQ(A,X)
      REAL GAMMQ,GAMSER,GLN,GAMMCF,A,X
      IF(X.LT.0..OR.A.LE.0.)PAUSE
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMQ=GAMMCF
      ENDIF
      RETURN
      END
c---------------------------------------------------------------
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      INTEGER ITMAX,N
      REAL GAMSER,GLN,GAMMLN,EPS,SUM,DEL,A,X,AP
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)PAUSE
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
c---------------------------------------------------------------
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      INTEGER ITMAX,N
      REAL X,A,GAMMCF,GLN,EPS,GOLD,GAMMLN,AN,ANA,A0,A1,B0,B1
      REAL FAC,ANF,G
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END
c---------------------------------------------------------------
      FUNCTION GAMMLN(XX)
      INTEGER J
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER,XX,GAMMLN
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
