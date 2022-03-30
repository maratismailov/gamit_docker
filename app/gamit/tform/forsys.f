      SUBROUTINE FORSYS(ITYP,IDMS,ILOC,IFILE)
C
C  Determine systems and formats for input and output
C  of data.
C
      implicit none

      include '../includes/tform.h'

      LOGICAL IEXIST
      integer*4 ityp,idms,iloc,ifile
      CHARACTER*16 FILNAM
C
C  System type
C
    9 WRITE(ISCRN,10)
   10 FORMAT(1X,'1=Cartesian  2=Spherical  3=Geodetic  4=Local  5=Cylind
     .rical: ',$)
      READ(ITERM,*) ITYP
      IF( ITYP.LT.1 .OR. ITYP.GT.5) THEN
         WRITE(ISCRN,1000)
         GOTO 9
      ENDIF

C  Format : note that deg/min/sec now available only with L-file format

      IF(ITYP.EQ.2 .or. ITYP.EQ.3) then
  21     WRITE(ISCRN,22)
  22     FORMAT(1X,'1=Decimal  2=L-file (Deg Min Sec): ',$)
         READ(ITERM,*) IDMS
         IF( IDMS.LT.1 .OR. IDMS.GT.2) THEN
            WRITE(ISCRN,1000)
            GOTO 21
         ENDIF 
      elseif (ityp.eq.5 ) then
  23     WRITE(ISCRN,24)
  24     FORMAT(1X,'1=Decimal : ',$)
         READ(ITERM,*) IDMS
         IF( IDMS.NE.1 ) THEN
            WRITE(ISCRN,1000)
            GOTO 23
         ENDIF 
 
      ENDIF

C  Local format

      IF(ITYP.EQ.4) THEN
   25  WRITE(ISCRN,26)
C   26    FORMAT(/,1X,'Format of local coordinates:',/,3x,
   26   FORMAT(/,1X,
     1   '1=North,East,Up (meters)',/,1x,
     2   '2=Delta Lat,Long (decimal degrees), Up (meters)',/,1x,
     4   'Format: ',$)
         READ(ITERM,*) ILOC
         IF( ILOC.LT.1 .OR. ILOC.GT.2) THEN
            WRITE(ISCRN,1000)
            GOTO 25
        ENDIF
      ENDIF

c  Enter file name for input or output

c**      IF (NUMSIT.NE.1) THEN
   39    IF (IFILE.EQ.1 .OR. IFILE.EQ.2) THEN
            WRITE(ISCRN,40)
         ELSEIF (IFILE.EQ.3) THEN
            WRITE(ISCRN,41)
         ELSEIF (IFILE.EQ.4) THEN
            WRITE(ISCRN,42)
         ENDIF
   40    FORMAT(/,1X,'Input Coordinate File ',/,1X,
     .   '(return for terminal input): ',$)
   41    FORMAT(/,1X,'Output Coordinate File ',/,1X,
     .   '(return for terminal input): ',$)
   42    FORMAT(/,1X,'Local Coordinate File ',/,1X,
     .   '(return for terminal input): ',$)
         READ(ITERM,'(A16)') FILNAM
         CALL LJUST(16,FILNAM)
         IF(FILNAM(1:1).NE.' ') THEN
            IF(IFILE.NE.3) THEN
               INQUIRE(FILE=FILNAM,EXIST=IEXIST)
               IF(.NOT.IEXIST) THEN
                  WRITE(ISCRN,45)
   45             FORMAT(/,1X,'File does not exist'/)
                  GOTO 39
               ENDIF
            ENDIF
            OPEN(UNIT=IFILE,FILE=FILNAM,status='unknown')
         ELSE
            IFILE = 0
         ENDIF
c**      ELSE
c**         IFILE = 0
c**      ENDIF
C
 1000 FORMAT(/,1X,'Not a valid option')
      RETURN
      END
