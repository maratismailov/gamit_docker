      INTEGER FUNCTION LUNIT()
C     FIND AN UNOPENED FORTRAN UNIT NUMBER BETWEEN 7 AND 100
c     from code by Geoff Abers
      LOGICAL LOPEN
      DO 1 LUNIT=7,99
         INQUIRE(UNIT=LUNIT,OPENED=LOPEN)
         IF (.NOT. LOPEN) RETURN
    1 CONTINUE
      LUNIT=0
      call suicid ('LUNIT: no unopenned unit between 7 and 99')
      RETURN
      END

