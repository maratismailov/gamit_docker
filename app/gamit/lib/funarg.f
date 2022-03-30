Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE FUNARG( TJD,L,F,D,ASCM )
C
C       M.E.ASH  SEPT 1967   SUBROUTINE FUNARG
C       FUNDAMENTAL ARGUMENTS OF BROWN LUNAR THEORY ARE EVALUATED AT
C       JULIAN EPHEMERIS DATE TJD

      implicit none

      REAL*8  L,ascm,funcof,d,f,tjd,t1
C
C.....TIMES FROM FUNDAMENTAL EPOCHS
      T1    =TJD-2415020.0D0
C
C.....MEAN LONGITUDE OF MOON ASCENDING NODE (BIG OMEGA)
      ASCM  =FUNCOF(T1,0.71995354167D0,-  1,-4.7094228332D-5,
     1                                         +0.432630D-14,+1.266D-22)
C.....LONGM-PERM (LITTLE L)
      L     =FUNCOF(T1,0.82251280093D0,+362,+9.1645684716D-5,
     1                                         +1.913865D-14,+8.203D-22)
C.....LONGS-PERS (LITTLE L PRIME)
C      LP    =FUNCOF(T1,0.99576620370D0,+ 27,+3.7778519279D-5,
C     1                                         -0.031233D-14,-1.900D-22)
C.....LONGM-ASCM (BIG F)
      F     =FUNCOF(T1,0.03125246914D0,+367,+4.8195691688D-5,
     1                                         -0.668609D-14,-0.190D-22)
C.....LONGM-LONGS (BIG D)
      D     =FUNCOF(T1,0.97427079475D0,+338,+6.3192198393D-5,
     1                                         -0.299023D-14,+1.077D-22)
C
      RETURN
      END
