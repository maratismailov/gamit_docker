Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE NEWFIL(FILEO,FILEN,NEWCHR)
C
C     CONVERT ONE TYPE OF FILE TO ANOTHER
C      CHANGE FIRST CHARACTER OF NAME

      implicit none

      integer icol

      CHARACTER*16 FILEO,FILEN
      CHARACTER*1 NEWCHR
C
      ICOL=INDEX(FILEO,'.')
      FILEN=FILEO
      FILEN(ICOL-6:ICOL-6)=NEWCHR
C
      RETURN
      END
