      SUBROUTINE UPNAM2(FILEO,FILEN)
C
C  UPDATE FILE NAME (Q-,M- FILES)
C    LOCATE CHARACTER BEFORE PERIOD (THE YEAR OR ITERATION CODE)
C     UPDATE NUMBER WITH CORRESPONDING LETTER

c     fix this routine so that input and output file name
c     can be the same  -kf 870914

C
      CHARACTER*16 FILEO,filen,hold
      CHARACTER*1 ICHAR

      integer*4 icol
C
      hold=FILEO
      ICOL = INDEX(FILEO,'.')
C CHARACTER BEFORE PERIOD
      ICHAR = FILEO(ICOL-1:ICOL-1)
C UPDATE CHARACTER   
c**   rwk 021022: change this since upchr no longer exists
c**      CALL UPCHR(ICHAR)         
      call newchr(ichar)
      hold(ICOL-1:ICOL-1) = ICHAR
C
      filen = hold

      RETURN
      END
