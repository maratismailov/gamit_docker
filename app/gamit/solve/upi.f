C*** This subroutine is incomplete
      SUBROUTINE UPI

C     Updates an I-file
C     Modified from UPL -- Yehuda Bock
      
      implicit none

      include '../includes/dimpar.h'    
      include 'solve.h'

      CHARACTER*4 SNAME(maxsit),STRBK,lowerc
      CHARACTER*16 TMPNAM
      CHARACTER*144 BUF144

      LOGICAL FCHECK

      integer ioerr,istat,i,j


C     LOOP THRU ALL STATIONS AND UPDATE CLOCK PARAMETERS
C***** fill in here

C     LOOP THRU I-FILE
      TMPNAM=iinf
      call lowers(tmpnam)
      if (fcheck(tmpnam)) then
         OPEN (UNIT=16,FILE=tmpnam,STATUS='OLD',IOSTAT=IOERR)
         IF (IOERR.NE.0) THEN
           call report_stat('WARNING','SOLVE','upi',tmpnam
     .          , 'Cannot open old I-file--no upodate',ioerr)
            RETURN
         ENDIF
      else
         call report_stat('WARNING','SOLVE','update',tmpnam
     .          , 'Cannot find I-file--no upodate',0)
         RETURN
      endif

      OPEN (UNIT=17,STATUS='SCRATCH',ERR=601)

C     SEARCH FOR ALL OCCURRENCES OF STATION
      DO 220 J=1,32000
         READ(16,'(A140)',END=209) BUF144
         STRBK = BUF144(1:4)
C
C        SORT TO FIND STATION
         DO 500 ISTAT=1,nsite
c            write(6,*) 'strbk,sname',strbk,sname(istat)
            IF(lowerc(STRBK).NE.lowerc(SNAME(ISTAT))) GO TO 500
C           WE FOUND IT
            WRITE(6,400) SNAME(ISTAT),STRBK
 400        FORMAT(1X,A4,1X,A4)
            GO TO 501

 500     CONTINUE

C        NO MATCH
         WRITE(17,'(A144)') BUF144
         go to 220

 501     CONTINUE
C
C        A MATCH - UPDATE COORDINATES
C******  fill in here
            WRITE(17,201) BUF144
201         FORMAT(A144)


220   CONTINUE
C     End global search on I-file

209   CLOSE(UNIT=16)
      REWIND 17
      TMPNAM = ioutf
      CALL LOWERS(TMPNAM)
      OPEN (UNIT=18,FILE=TMPNAM,STATUS='UNKNOWN')
C     COPY SCRATCH FILE TO NEW I-FILE
      DO 600 I=1,32000
         READ(17,201,END=601) BUF144
         WRITE(18,201) BUF144
  600 CONTINUE

  601 CLOSE(UNIT=18)
      CLOSE(UNIT=17)

      RETURN
      END
