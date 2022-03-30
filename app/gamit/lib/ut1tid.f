Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE UT1TID( UT1,L,F,D,ASCM )
C
C       Compute UT1 from UT1R:  i.e, add the short-period
C       (9-, 14-, 30-day) terms due to tidal effects
C       which have been removed for smoothing
C       Ref:  Yoder et al., J. Geophys. Res. 86, 881-891, 1981
C         
      implicit none

C
      REAL*8 UT1PAR(2)
      REAL*8 L,KCF,KCM
      real*8 ascm,d,ut1,f,argmnt
C
      KCF= 0.94D0
      KCM= 0.94D0
      UT1PAR(1)= ( .0437D0*DSIN(ARGMNT(L + 2.D0*F + ASCM))
     1              + .1056D0*DSIN(ARGMNT(L + 2.D0*F + 2.D0*ASCM))
     2              + .0210D0*DSIN(ARGMNT(-L+2.D0*F+2.D0*D+2.D0*ASCM))
     3              + .0318D0*DSIN(ARGMNT(2.D0*F))
     4              + .3413D0*DSIN(ARGMNT(2.D0*F + ASCM))
     5              + .8252D0*DSIN(ARGMNT(2.D0*F + 2.D0*ASCM))
     6              + .0781D0*DSIN(ARGMNT(2.D0*D))
     7              + .0106D0*DSIN(ARGMNT(ASCM))*DCOS(ARGMNT(2.D0*D))
     8              + .0360D0*(1.D0 - .12D0*DCOS(ARGMNT(ASCM)))
     9                             *DSIN(ARGMNT(2.D0*L)) )*1.D-3
      UT1= UT1 - KCF*UT1PAR(1)
cd     WRITE (6,55555) UT1,UT1PAR(1)
cd55555 FORMAT (1X,'In UT1TID',2(D22.15,1X))
      UT1PAR(2)= (- .0188D0*DSIN(ARGMNT(-L + 2.D0*F + ASCM))
     1              - .0463D0*DSIN(ARGMNT(-L + 2.D0*F + 2.D0*ASCM))
     2              + .8788D0*(1.D0 -.131D0*DCOS(ARGMNT(ASCM)))
     3                              *DSIN(ARGMNT(L))
     4              + .1940D0*(1.D0 - .137D0*DCOS(ARGMNT(ASCM)))
     5                              *DSIN(ARGMNT(-L + 2.D0*D)) )*1.D-3
      UT1= UT1 - KCM*UT1PAR(2)
c     WRITE (6,55555) UT1,UT1PAR(2)
C
      RETURN
      END
