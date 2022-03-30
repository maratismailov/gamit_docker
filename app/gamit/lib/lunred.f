Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE LUNRED(ISPEED,FJD,RVEC)
C
C THIS VERSION OF LUNRED INTERPOLATES FROM TABULAR LUNAR
C EPHEMERIS. VELOCITIES ARE COMPUTED FROM POSITIONS.
C SERGEI GOUREVITCH, THE ORGIINAL EPHEMERIS INTERPOLATOR
C RICK ABBOT NOVEMBER 1984, MODIFICATIONS FOR LUNAR EPHEMERIS
C
C ISPEED = 0 DO ALL
C ISPEED = 1 POSTIONS ONLY
C   2 DO THE REST
C   3 VELOCITIES ONLY
C
      implicit none      

      include '../includes/dimpar.h' 
      include '../includes/units.h'
      include '../includes/arc.h'


      real*8 rvec(6),trun,fjd,p,q,p2,q2,sum1,sum2 
      integer ispeed,ir,il,j,k,kup,kdwn,insart,numper,ijjk,ijko,i

      logical debug/.false./

      TRUN = FJD - FJDBMN + 1.D0/SDELTM
C
C        AT THIS POINT TRUN IS TIME SINCE FIRST EPOCH
C
C        SDELTM = 1/INTERVAL
      P=TRUN*SDELTM
      JNOWM=P
C
  990 CONTINUE
      IF(iendfm.NE.0) GO TO 1000
      IF((JNOWM+5.LE.JLASTM).AND.(JLASTM.GE.10)) GO TO 1000
C
      JLASTM=JLASTM+1
      JI0M=MOD(JI0M,10)+1
      JILM=MOD(JI0M+3,10)+1
C
      READ(ilun,5000,END=995)
     $          (YYM(JILM,IR),IR=1,nintrsm)
 5000 FORMAT (6X,3F11.3)
C
      IF((JNOWM+8.LE.JLASTM).AND.(JLASTM.GE.12)) GO TO 990
C
      IY1M=IY2M
      IY2M=MOD(IY2M,2)+1
      if(debug) then
        do j=1,10
          write(*,'(a,i3,3f16.2)') 'yym j ',j,(yym(j,k),k=1,3)
        enddo
      endif
      DO 520 IL=1,nintrsm
      DO 500 J=1,5
      YTRPM(J,IY2M,IL)=0.0D0
      YTRPM(J,IY2M,IL+nintrsm)=0.0D0
      DO 400 K=5,2,-1
      KUP=MOD(JI0M+K+8,10)+1
      KDWN=MOD(JI0M-K+10,10)+1
      YTRPM(J,IY2M,IL)=YTRPM(J,IY2M,IL)
     $              +EVCF(K,J)*(YYM(KUP,IL)+YYM(KDWN,IL))
  400 CONTINUE
      YTRPM(J,IY2M,IL)=
     $               YTRPM(J,IY2M,IL)+EVCF(1,J)*YYM(JI0M,IL)
      YTRPM(J,IY2M,IL+nintrsm)=
     $               YTRPM(J,IY2M,IL)*FACT(J)*SDELTM
  500 CONTINUE
  520 CONTINUE  

      if(debug) then
        do j=1,2
          do i=1,5
            write(*,'(a,i2,3f16.2)') 'j i ytrpm',j,(ytrpm(i,j,k),k=1,3)
          enddo
        enddo
      endif


      GO TO 990
C
  995 iendfm=1
      JNOWM=JNOWM-1
C
 1000 CONTINUE
C
C        NOW THE INTERPOLATION VECTOR Y IS SET UP
C
      P=P-JLASTM+5
      Q=1.D0-P
      P2=P*P
      Q2=Q*Q

      if(debug) write(*,'(a,i3,4f10.5)') 
     .     'lunred jlasts p q p2 q2 ',jlasts,p,q,p2,q2

      NUMPER=nintrsm+3
C
C        ISPEED = 0	DO ALL
C        ISPEED = 1	POSITIONS ONLY
C                 2	DO THE REST
C                 3	VELOCITIES ONLY
C
      INSART=1
      IF(ISPEED.EQ.1) NUMPER=3
      IF(ISPEED.EQ.2) INSART=4
      IF(ISPEED.EQ.3) INSART=nintrsm+1
C
      DO 2500 IJJK=INSART,NUMPER
C
      SUM1=0.0D0
      SUM2=0.0D0
      IJKO=0
C
      IF(IJJK.LE.3) IJKO=IJJK
      IF((IJJK.GT.3).AND.(IJJK.LE.nintrsm)) IJKO=IJJK+3
      IF(IJJK.GT.nintrsm) IJKO=IJJK-nintrsm+3
      IF(IJJK.GT.nintrsm+3) IJKO=IJJK  
      if (ijko .eq. 0) then
        call report_stat('FATAL','ARC','lunred',' '
     .  ,'Error reading the luntab file. IJKO = 0',0)
      endif
C
      DO 2000 J=5,2,-1
      SUM1=P2*(SUM1+YTRPM(J,IY2M,IJJK))
      SUM2=Q2*(SUM2+YTRPM(J,IY1M,IJJK))
 2000 CONTINUE
C
      IF(IJJK.LE.nintrsm)
     $ RVEC(IJKO)=P*(YTRPM(1,IY2M,IJJK)+SUM1)
     $ + Q*(YTRPM(1,IY1M,IJJK)+SUM2)
      IF(IJJK.GT.nintrsm)
     $ RVEC(IJKO)=(YTRPM(1,IY2M,IJJK)+SUM1)
     $ - (YTRPM(1,IY1M,IJJK)+SUM2)
C
 2500 CONTINUE                            
      if(debug) then
        print *,'rvec ',(rvec(j),j=1,3)
cc        stop ' LUNRED stop ' 
      endif
C

C
      RETURN
      END
