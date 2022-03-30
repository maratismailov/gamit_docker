      INTEGER FUNCTION LUNIT()
C     FIND AN UNOPENED FORTRAN UNIT NUMBER BETWEEN 7 AND 100
c     from code by Geoff Abers 

      implicit none
                   
      character*80 prog_name
      integer*4 len,rcpar 
      LOGICAL LOPEN

c     get calling program for report_stat
      len = rcpar(0,prog_name)

      DO 1 LUNIT=7,99
         INQUIRE(UNIT=LUNIT,OPENED=LOPEN)
         IF (.NOT. LOPEN) RETURN
    1 CONTINUE
      LUNIT=0   
      call report_stat('FATAL',prog_name,'comlib/lunit',' '
     .   ,'No unopened unit between 7 and 99',0)
      RETURN
      END

