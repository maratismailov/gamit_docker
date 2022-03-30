Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      FUNCTION ARGMNT(ANGLE)
C
C       M.E.ASH  SEPT 1967   FUNCTION SUBROUTINE ARGMNT
C       GET ANGLE IN REVOLUTIONS BETWEEN - TWOPI AND TWOPI RADIANS
C       ANGLE =  INPUT ANGLE IN REVOLUTIONS
C
      implicit none

      real*8 twopi, argmnt,angle
C
      TWOPI=8.D0*DATAN(1.D0)
C
      ARGMNT=  (DMOD(ANGLE,1.D0))*TWOPI
C
      RETURN
      END
