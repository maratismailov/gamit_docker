      SUBROUTINE QUERYB( NSITE, STBIAS)
C
C PURPOSE:
C     QUERY FOR BIAS TYPE: EXPLICIT OR IMPLICIT
C
      implicit none
C
      include '../includes/dimpar.h'

      INTEGER NSITE
      INTEGER I
      INTEGER MAXLEN
      CHARACTER*1  EXPLCT
      CHARACTER*1  STBIAS(maxsit)
      character*2  buf2
      CHARACTER*3  ANET
      CHARACTER*19 FORM1
      CHARACTER*15 FORM2
      CHARACTER*256 message
      character*6  form3
C
      PARAMETER (MAXLEN=5+maxsit*3)
      CHARACTER*(MAXLEN) JSCR,ISCR
      DATA EXPLCT/'E'/
C
      IF(NSITE.GT.maxsit) then
        write(message,'(3(a,i2))')
     .       'Number of sites in the solution: (',nsite
     .       , '), exceeds maxsit: (',maxsit,')'
        call report_stat('FATAL','CFMRG','queryb',' ',message,0)
      endif
C
      call uppers(explct)
      DO 10 I = 1, MAXLEN
         JSCR(I:I) = ' '
         ISCR(I:I) = ' '
   10 CONTINUE

      WRITE(ANET,'(I3)') maxsit
      FORM1 = '(1X,A,' //ANET// '(1X,A2),/)'
      FORM2 = '(5X,'   //ANET// '(2X,A1))'

      form3 = '(50a1)'
      buf2 = '  '
      write(buf2,'(i2)') nsite
      write(form3(2:3),'(a2)') buf2
      read(5,form3) (stbias(i),i=1,nsite)

      RETURN
      END
