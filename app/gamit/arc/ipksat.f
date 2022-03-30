Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE IPKSAT

C     USER SELECTS THE SATELLITES TO ANALYZE, UP TO MAXSAT SATELLITES POSSIBLE
C     RICK ABBOT - NOVEMBER 1984

      implicit none

      include '../includes/dimpar.h'

      INTEGER*4 IPRN,NPKSAT,ISUNIT,IOPEP,I,ioerr

      CHARACTER* 3 UPPERC
      CHARACTER*16 SATNAM(MAXSAT),SNAME,STRING
      COMMON/PIKAUX/NPKSAT,ISUNIT,IOPEP
      COMMON/SATLNM/SATNAM,SNAME,IPRN
C
      DO 10 I=1,MAXSAT
   10 SATNAM(I)='                '
C
      NPKSAT=0
  100 NPKSAT=NPKSAT+1
      READ (5,150,iostat=ioerr) STRING
      if ( ioerr .ne. 0 ) then
        call report_stat('STATUS','ARC','ipksat','string',
     .  'Error reading satellites from the ARC input batch file',ioerr)
      endif
  150 FORMAT (A16)
      IF (UPPERC(STRING(1:3)).EQ.'END') GO TO 200
      SATNAM(NPKSAT)= STRING
      GO TO 100
  200 NPKSAT=NPKSAT-1
C
      RETURN
      END
