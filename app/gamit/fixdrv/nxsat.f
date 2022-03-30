      INTEGER*4  FUNCTION  NXSAT( XFILE, XORC, LXFIL, LCFIL )
C
C        S. Shimada  01/18/1990  at IGPP
C
C     function to read number of satellites of X-(C-)file
C     arguments
C        XFILE  : X-file name                               input
C        XORC   : First character of X-,C-file              input
C        LXFIL  : unit number of X-file                     input
C        LCFIL  : unit number of C-file                     input

      implicit none

      include '../includes/dimpar.h'
C

      CHARACTER* 1  XORC
      character* 3  rcvrsw
      CHARACTER*16  XFILE, CFILE
      character*20  rctype,rcvnum,anttyp,antnum
      CHARACTER*32  SITNAM
      CHARACTER*80  TEXT(MAXTXT)
      CHARACTER*120 LINE
C
      INTEGER*4     ISCHAN(MAXSAT),NTEXT, NDAT, LAMBDA(MAXSAT,MAXDAT)
     .            , NEXTRA, DATTYP(MAXDAT), lxfil, inter, nslip
     .            , id, im, iy, ihr, imin, nepoch, lfntop, nchan, lcfil
     .            , i0, mtime, ircint,isessn, ioerr

      INTEGER*2     ISLIP(MAXCSB),ISLPST(MAXCSB)

      REAL*4  SWVER

      REAL   * 8    EXTRA(MAXEXT), offarp(3),OFFSL1(3), OFFSL2(3)
     .            , sec
C
c initialise
      NXSAT = 0
c
C     X-file
C
      IF( XORC .EQ. 'x' )  THEN
C
C        find first epoch
C
         OPEN( LXFIL, FILE=XFILE, STATUS='OLD',err=90 )
   10    CONTINUE
         READ( LXFIL, '(A120)' )  LINE
         IF( LINE(1:3) .NE. 'END' )  GO TO 10
   20    CONTINUE
         READ( LXFIL, '(A120)' )  LINE
         IF( INDEX(LINE,'SATELLITES') .EQ. 0 )  GO TO 20
         READ( LINE(1:3), '(I3)' )  NCHAN
         CLOSE( LXFIL )
C
C     C-file
C
      ELSE
         I0 = LFNTOP( XFILE )
         CFILE = XFILE
         CFILE(I0:I0) = 'c'
c         I1 = NBLEN( XFILE )
c         CFILE = 'c'//XFILE(I0+1:I1)
c         CFILE(I1-I0-3:I1-I0-3) = PHSDT(2:2)
c         CALL  LOWERS( CFILE(I1-I0-3:I1-I0-3) )
         CALL  COPENS( CFILE, 'OLD', LCFIL, IOERR )
         OPEN( 8, STATUS='SCRATCH' )
         CALL  CHDRED( LCFIL, 8, NEPOCH, INTER, MTIME,
     .                 IY, IM, ID, IHR, IMIN, SEC,
     .                 NCHAN, ISCHAN, NDAT, DATTYP, LAMBDA,
     .                 offarp,OFFSL1, OFFSL2, SITNAM, RCVRSW, SWVER,
     .                 rctype,rcvnum,anttyp,antnum,
     .                 ircint,isessn,NTEXT, TEXT, NEXTRA, EXTRA,
     .                 NSLIP, ISLIP, ISLPST)
         CLOSE( 8 )
         CLOSE( LCFIL )
      ENDIF
C
      NXSAT = NCHAN
C
      RETURN

  90  call report_stat('FATAL','FIXDRV','nxsat','xfile'
     .             ,'Error opening X-file to get # SVs',0)

      END



