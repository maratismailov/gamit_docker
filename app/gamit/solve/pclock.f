Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      SUBROUTINE PCLOCK(KL,ISTAT,ISAT,DELTA,TMPART,ICFILE)
C
C   COMPUTE CLOCK RATE AND ACCELERATION PARTIALS

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer kl,istat,isat,icfile,ii1

      real*8 delta,tmpart,day0,dayor,tmt0

                 CALL DAYNUM(IT0(1),IT0(2),IT0(3),DAY0)
                 II1 = ISTAT + ICFILE
                 CALL DAYNUM(ITOR(1,II1),ITOR(2,II1),
     1                       ITOR(3,II1),DAYOR)
                 DAY0=DAY0+T00(1)/24.D0+T00(2)/(24.D0*60.D0)
     1                    +T00(3)/(24.D0*60.D0*60.D0)
                 DAYOR=DAYOR+TOR(1,II1)/24.D0
     1                      +TOR(2,II1)/(24.D0*60.D0)
     2                      +TOR(3,II1)/(24.D0*60.D0*60.D0)
                 TMT0=DELTA+(DAY0-DAYOR)*86400.D0
CD                 WRITE(6,99) DAY0,DAYOR,DELTA,TMT0
CD                 WRITE(16,99) DAY0,DAYOR,DELTA,TMT0
CD   99            FORMAT(/,' TMT0 FOR RATE PARTIAL (SECS):',4F10.2)
C RATE
                 IF(KL.EQ.21)
     1           TMPART=TPART(5,ISTAT,ISAT)*TMT0
CD     WRITE(16,777) ISTAT,ISAT,TMPART
cd  777 FORMAT(1X,2I5,1X,D15.8)
C ACCELERATION
                 IF(KL.EQ.22)
     1           TMPART=0.5D0*TPART(5,ISTAT,ISAT)*TMT0**2/86400.D0
C
      RETURN
      END
