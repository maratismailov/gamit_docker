      program fdist
c
c     Given two degrees of freedom, create a table of f-distributions
c
c
      print *,'Program to create a table of f-distributions'
c
      print*,'Input degrees of freedom'
      read(5,*) df1,df2
10    print*,'Input do loop specifications'
      read(5,*) begin,end,step

      do f = begin,end,step
         PROB = BETAI(0.5*DF2,0.5*DF1,DF2/(DF2+DF1*F))
c     *    +(1.-BETAI(0.5*DF1,0.5*DF2,DF1/(DF1+DF2/F)))
         write(6,'(2f10.4)') f,prob
      enddo
      goto 10
      end
c----------------------------------------------------------------
      FUNCTION BETA(Z,W)
      BETA=EXP(GAMMLN(Z)+GAMMLN(W)-GAMMLN(Z+W))
      RETURN
      END
c----------------------------------------------------------------
      FUNCTION BETACF(A,B,X)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      AM=1.
      BM=1.
      AZ=1.
      QAB=A+B
      QAP=A+1.
      QAM=A-1.
      BZ=1.-QAB*X/QAP
      DO 11 M=1,ITMAX
        EM=M
        TEM=EM+EM
        D=EM*(B-M)*X/((QAM+TEM)*(A+TEM))
        AP=AZ+D*AM
        BP=BZ+D*BM
        D=-(A+EM)*(QAB+EM)*X/((A+TEM)*(QAP+TEM))
        APP=AP+D*AZ
        BPP=BP+D*BZ
        AOLD=AZ
        AM=AP/BPP
        BM=BP/BPP
        AZ=APP/BPP
        BZ=1.
        IF(ABS(AZ-AOLD).LT.EPS*ABS(AZ))GO TO 1
11    CONTINUE
      PAUSE 'A or B too big, or ITMAX too small'
1     BETACF=AZ
      RETURN
      END
c----------------------------------------------------------------
      FUNCTION BETAI(A,B,X)
      IF(X.LT.0..OR.X.GT.1.)PAUSE 'bad argument X in BETAI'
      IF(X.EQ.0..OR.X.EQ.1.)THEN
        BT=0.
      ELSE
        BT=EXP(GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)
     *      +A*ALOG(X)+B*ALOG(1.-X))
      ENDIF
      IF(X.LT.(A+1.)/(A+B+2.))THEN
        BETAI=BT*BETACF(A,B,X)/A
        RETURN
      ELSE
        BETAI=1.-BT*BETACF(B,A,1.-X)/B
        RETURN
      ENDIF
      END
c----------------------------------------------------------------
      FUNCTION GAMMLN(XX)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
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
c----------------------------------------------------------------
      SUBROUTINE FTEST(DATA1,N1,DATA2,N2,F,PROB)
      DIMENSION DATA1(N1),DATA2(N2)
      CALL AVEVAR(DATA1,N1,AVE1,VAR1)
      CALL AVEVAR(DATA2,N2,AVE2,VAR2)
      IF(VAR1.GT.VAR2)THEN
        F=VAR1/VAR2
        DF1=N1-1
        DF2=N2-1
      ELSE
        F=VAR2/VAR1
        DF1=N2-1
        DF2=N1-1
      ENDIF
      PROB = BETAI(0.5*DF2,0.5*DF1,DF2/(DF2+DF1*F))
     *    +(1.-BETAI(0.5*DF1,0.5*DF2,DF1/(DF1+DF2/F)))
      RETURN
      END
c----------------------------------------------------------------
      SUBROUTINE AVEVAR(DATA,N,AVE,VAR)
      DIMENSION DATA(N)
      AVE=0.0
      VAR=0.0
      DO 11 J=1,N
        AVE=AVE+DATA(J)
11    CONTINUE
      AVE=AVE/N
      DO 12 J=1,N
        S=DATA(J)-AVE
        VAR=VAR+S*S
12    CONTINUE
      VAR=VAR/(N-1)
      RETURN
      END
