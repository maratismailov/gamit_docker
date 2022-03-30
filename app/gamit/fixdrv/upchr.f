Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE  UPCHR( CHAR )
C
C     Update the session code (up to 10 sessions)
C     Match 1-10 to a-j (10 is 0)
C
      CHARACTER*1  CHAR
      CHARACTER*1  HEX(10), ALPHA(10)
      integer*4    i
C
C     HEX 1-9, 0
      DATA  HEX  /'1','2','3','4','5','6','7','8','9','0'/
      DATA  ALPHA/'a','b','c','d','e','f','g','h','i','j'/
C
      DO  10  I = 1, 10
         IF( CHAR .NE. HEX(I) )  GO TO 10
         CHAR = ALPHA(I)
         GO TO 11
   10 CONTINUE
   11 CONTINUE
C
      RETURN
      END
