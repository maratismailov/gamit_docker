Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE SOLRED(ISPEED,FJD,RVEC)
C
C THIS VERSION OF LUNRED INTERPOLATES FROM TABULAR SOLAR
C EPHEMERIS. VELOCITIES ARE COMPUTED FROM POSITIONS.
C SERGEI GOUREVITCH, THE ORGIINAL EPHEMERIS INTERPOLATOR
C RICK ABBOT NOVEMBER 1984, MODIFICATIONS FOR LUNAR EPHEMERIS
C
C ISPEED = 0	DO ALL
C ISPEED = 1	POSTIONS ONLY
C   2 DO THE REST
C   3 VELOCITIES ONLY
C
      implicit none  

      include '../includes/dimpar.h' 
      include '../includes/units.h'
      include '../includes/arc.h'
  
      real*8 rvec(6),trun,fjd,p,q,p2,q2,sum1,sum2
      integer ir,il,i,j,k,kup,insart,numper,ijjk,ijko,ispeed,kdwn

      logical debug/.false./

      TRUN = FJD - FJDBSN + 1.D0/SDELTS
C
C        AT THIS POINT TRUN IS TIME SINCE FIRST EPOCH
C
C        SDELTS = 1/INTERVAL
      P=TRUN*SDELTS
      JNOWS=P
      if( debug ) write(*,'(a,3f14.6,f4.0,f12.8,i4)') 
     .   'SOLRED fjd trun fjdbsn sdelts p jnows '
     .   ,        fjd,trun,fjdbsn,sdelts,p,jnows 

C
c     Burc points out that iendf must be initialized--kurt 901102
c     IENDF = 0
c     No, I think this is meant to be IENDFS, initialized in EPHDRD--rwk 950728
  990 CONTINUE
      IF(iendfs.NE.0) GO TO 1000
      IF((JNOWS+5.LE.JLASTS).AND.(JLASTS.GE.10)) GO TO 1000
C
      JLASTS=JLASTS+1
      JI0S=MOD(JI0S,10)+1
      JILS=MOD(JI0S+3,10)+1
C
      READ(isun,5000,END=995)
     $          (YYS(JILS,IR),IR=1,nintrss)
 5000 FORMAT (6X,3F11.0)
      if(debug) write(*,'(a,4i4,6f16.3)')
     .  ' jlasts ji0s jils nintrss READ yys '
     .  , jlasts,ji0s,jils,nintrss,(yys(jils,ir),ir=1,nintrss)

      IF((JNOWS+8.LE.JLASTS).AND.(JLASTS.GE.12)) GO TO 990
C
      IY1S=IY2S
      IY2S=MOD(IY2S,2)+1
                         
      if(debug) then
        do j=1,10
          write(*,'(a,i3,3f16.2)') 'yys j ',j,(yys(j,k),k=1,3)
        enddo
      endif
      DO 520 IL=1,nintrss
      DO 500 J=1,5
      YTRPS(J,IY2S,IL)=0.0D0
      YTRPS(J,IY2S,IL+nintrss)=0.0D0
      DO 400 K=5,2,-1
      KUP=MOD(JI0S+K+8,10)+1
      KDWN=MOD(JI0S-K+10,10)+1
cd      if(debug)  write(*,'(a,6i4,2f16.2)') 
cd     .      'il j k ji0s kdwn kup yys(dwn) yyys(up) '
cd     .     ,il,j,k,ji0s,kdwn,kup,yys(kdwn,il),yys(kup,il)
      YTRPS(J,IY2S,IL)=YTRPS(J,IY2S,IL)
     $              +EVCF(K,J)*(YYS(KUP,IL)+YYS(KDWN,IL))
cd      if(debug) write(*,'(a,7i3,2d16.2)')  
cd     .   'il j k ji0s kup kdwn iy2 evcf ytrps '
cd     .   ,il,j,k,ji0s,kup,kdwn,iy2s,evcf(k,j),ytrps(j,iy2s,il)
      if( debug.and.il.eq.1 ) 
     .  write(*,'(a,5i3,4f16.2)')
     .   'For il=1 iy2s j k kdwn kup yys-dwn yys-up evcf ytrps '
     .    ,iy2s,j,k,kup,kdwn,yys(kdwn,il),yys(kup,il),evcf(k,j)
     .    ,ytrps(j,iy2s,il)
      
  400 CONTINUE
      YTRPS(J,IY2S,IL)=
     $               YTRPS(J,IY2S,IL)+EVCF(1,J)*YYS(JI0S,IL)
      if( debug.and.il.eq.1 )   
     .   write(*,'(a,i3,3f16.2)')
     .    '   ji0s yys-0 evcf1j ytrps '
     .      , ji0s,yys(ji0s,il),evcf(1,j),ytrps(j,iy2s,il)
cd       if(debug) write(*,'(a,4i4,2f16.2)') 
cd     .     'il j iy2s ji0s evcf yys '
cd     .     ,il,j,iy2s,ji0s,evcf(1,j),ytrps(j,iy2s,il)
      YTRPS(J,IY2S,IL+nintrss)=
     $               YTRPS(J,IY2S,IL)*FACT(J)*SDELTS
cd      if(debug) write(*,'(a,3i4,2f16.2)') 
cd     .    'j iy2s il fact sdelts '
cd      ,    ,j,iy2s,il,fact(j),sdelts

  500 CONTINUE
  520 CONTINUE      

      if(debug) then
        do j=1,2
          do i=1,5
            write(*,'(a,i2,3f16.2)') 'j i ytrps',j,(ytrps(i,j,k),k=1,3)
          enddo
        enddo
      endif
         

      GO TO 990
C
  995 iendfs=1
      JNOWS=JNOWS-1
C
 1000 CONTINUE
C
C        NOW THE INTERPOLATION VECTOR Y IS SET UP
C
      P=P-JLASTS+5
      Q=1.D0-P
      P2=P*P
      Q2=Q*Q
      if(debug) write(*,'(a,i3,4f10.5)') 
     .     'jlasts p q p2 q2 ',jlasts,p,q,p2,q2
C
      NUMPER=nintrss+3
C
C        ISPEED = 0	DO ALL
C        ISPEED = 1	POSITIONS ONLY
C                 2	DO THE REST
C                 3	VELOCITIES ONLY
C
      INSART=1
      IF(ISPEED.EQ.1) NUMPER=3
      IF(ISPEED.EQ.2) INSART=4
      IF(ISPEED.EQ.3) INSART=nintrss+1
C
      DO 2500 IJJK=INSART,NUMPER
C
      SUM1=0.0D0
      SUM2=0.0D0
      IJKO=0
C
      IF(IJJK.LE.3) IJKO=IJJK
      IF((IJJK.GT.3).AND.(IJJK.LE.nintrss)) IJKO=IJJK+3
      IF(IJJK.GT.nintrss) IJKO=IJJK-nintrss+3
      IF(IJJK.GT.nintrss+3) IJKO=IJJK
      if (ijko .eq. 0) then
        call report_stat('FATAL','ARC','solred',' '
     .  ,'Error reading the soltab file. IJKO = 0',0)
      endif
C
      DO 2000 J=5,2,-1
      SUM1=P2*(SUM1+YTRPS(J,IY2S,IJJK))
      SUM2=Q2*(SUM2+YTRPS(J,IY1S,IJJK))         
      if(debug) write(*,'(a,4i3,2f18.2)') 'j iy1s iy2s ijjk sum1 sum2 '
     .      ,j,iy1s,iy2s,ijjk,sum1,sum2
 2000 CONTINUE
C
      IF(IJJK.LE.nintrss)
     $ RVEC(IJKO)=P*(YTRPS(1,IY2S,IJJK)+SUM1)
     $ + Q*(YTRPS(1,IY1S,IJJK)+SUM2)
      if(debug.and.ijko.eq.1) then
        write(*,'(a,3i3,2f16.2)') 
     .    'ijjk ijko iy2s p ytrps(1,iy2s) '
     .    ,ijjk,ijko,iy2s,ytrps(1,iy2s,ijjk)     
        write(*,'(a,3i3,2f16.2)') 
     .    'ijjk ijko iy1s q ytrps(1,iy1s) '
     .    ,ijjk,ijko,iy1s,ytrps(1,iy1s,ijjk)
      endif
      IF(IJJK.GT.nintrss)
     $ RVEC(IJKO)=(YTRPS(1,IY2S,IJJK)+SUM1)
     $ - (YTRPS(1,IY1S,IJJK)+SUM2)
C
 2500 CONTINUE
      if(debug) then 
        print *,'rvec ',(rvec(j),j=1,3)
cc        stop ' SOLRED stop ' 
      endif
C
      RETURN
      END
