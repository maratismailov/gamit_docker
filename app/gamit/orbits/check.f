Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
c  mod pch mar 90 - allow lowercase responses.
C
      SUBROUTINE CHECK(IFLAG)
C.....CHECK FOR USER PROMPT OF (Y)ES OR (N)O

      implicit none

      CHARACTER*1 YESNO
      integer*4 iflag
C
   50 WRITE (6,1000)
 1000 FORMAT (30X,'<Y>ES OR <N>O: ',$)
      READ (5,2000) YESNO
 2000 FORMAT (A1)
      IF (YESNO.EQ.'Y '.OR.YESNO.EQ.'N ') GO TO 100
      IF (YESNO.EQ.'y '.OR.YESNO.EQ.'n ') GO TO 100
      WRITE (6,3000)
 3000 FORMAT (1X,'INAPPROPRIATE RESPONSE')
      GO TO 50
  100 CONTINUE
      IF (YESNO.EQ.'Y ' .OR. yesno .EQ. 'y ' ) IFLAG=1
      IF (YESNO.EQ.'N ' .OR. yesno .EQ. 'n ' ) IFLAG=2
      RETURN
      END
