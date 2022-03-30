Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine DTWOPI
C
C     INITIALIZE VALUES ASSOCIATED WITH CIRCLES
C     RICK ABBOT - NOVEMBER 1984
C
      implicit none  

      include '../includes/dimpar.h'   
      include '../includes/arc.h'


c      COMMON/STWOPI/TWOPI,CDR,CASR
c      real*8 twopi,cdr,casr
c      DATA TWOPI/6.283185307179586D0/
c      DATA CDR  /1.74532925199433D-02/
c      DATA CASR /4.84813681109536D-06/

      twopi = 6.283185307179586D0
      cdr  =  1.74532925199433D-02
      casr =  4.84813681109536D-06

      return

      END
