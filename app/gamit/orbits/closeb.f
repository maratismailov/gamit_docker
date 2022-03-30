Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE CLOSEB
     1   ( IUTIN,IUTOUT,INUT,IUT1,IPOLE )
C
C      Close all rotation input and output files opened by OPENB for
c      subroutine TROT or its calling programs BCTOT, NGSTOT, TTONGS
C        Last mod for orbits directory   R King   18 Nov 91
C
      implicit none

      integer*4 iutin,iutout,inut,iut1,ipole
C
c
      close (iutin)
      close (iutout)
c     inut may not be open but I think this gives no error 
      close (inut)
      close (iut1)
      close (ipole)

      return
      end
