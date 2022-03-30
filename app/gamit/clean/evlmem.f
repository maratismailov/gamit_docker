      real*8 FUNCTION EVLMEM(FDT,COF,M,PM)
c     routine from numerical recipes
c     converted to double precision

c     given COF, M, PM as returnde by MEMCOF, this function returns
c     the power spectrum estimate p(f) as a function of FDT = f*delta

ckf   DIMENSION COF(M)
      real*8 COF(*)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      real*8 sumr,sumi,fdt,pm
      integer m,i

      THETA=6.28318530717959D0*FDT
      WPR=DCOS(THETA)
      WPI=DSIN(THETA)
      WR=1.D0
      WI=0.D0
      SUMR=1.
      SUMI=0.
      DO 11 I=1,M
        WTEMP=WR
        WR=WR*WPR-WI*WPI
        WI=WI*WPR+WTEMP*WPI
ckf     SUMR=SUMR-COF(I)*SNGL(WR)
ckf     SUMI=SUMI-COF(I)*SNGL(WI)
        SUMR=SUMR-COF(I)*WR
        SUMI=SUMI-COF(I)*WI
11    CONTINUE
      EVLMEM=PM/(SUMR**2+SUMI**2)
      RETURN
      END

