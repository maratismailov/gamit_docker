      SUBROUTINE INVER2(A,B,IJOB,N,rcond,IER)

      include '../includes/dimpar.h'

C     Modified to use LAPACK routines by Pascal Gegout and Simon McClusky March 2006.

      REAL*8 A(*),B(*),WORK(3*MAXBIS),RCOND,ANORM
      INTEGER IJOB,N,IER,IWORK(MAXBIS)
c** uncomment for debug
c      CHARACTER*256 MESSAGE

      CHARACTER*1 UPLO
      INTEGER NRHS,LDB,INFO 

      logical debug/.false./
      integer i

C     WRITE(MESSAGE,'(A5,I5)') 'IJOB=',IJOB
C     CALL REPORT_STAT('STATUS','SOLVE','INVER2',' ',MESSAGE,0)
              
      if(debug) then 
        print *,'INVER2 ijob n maxbis ',ijob,n,maxbis
c        stop 
        print *,   'A ',(a(i),i=1,n) 
      endif

      UPLO='U'
      NRHS=1
      LDB=N
      INFO=0
      ANORM=1
      RCOND = 0.D0

      CALL DPPTRF(UPLO,N,A,INFO)

      IF (INFO.EQ.0 .AND. IJOB.EQ.4) THEN
         CALL DPPCON(UPLO,N,A,ANORM,RCOND,WORK,IWORK,INFO)
      ENDIF

C     WRITE(MESSAGE,'(A12,I5)') 'DPPTRF INFO=',INFO
C     CALL REPORT_STAT('STATUS','SOLVE','INVER2',' ',MESSAGE,0)

      IF (INFO.EQ.0) THEN
        IF (IJOB.EQ.2.OR.IJOB.EQ.3) THEN
        CALL DPPTRS(UPLO,N,NRHS,A,B,LDB,INFO)
C     WRITE(MESSAGE,'(A12,I5)') 'DPPTRS INFO=',INFO
C     CALL REPORT_STAT('STATUS','SOLVE','INVER2',' ',MESSAGE,0)
        END IF
      END IF

      IF (INFO.EQ.0) THEN
        IF (IJOB.EQ.1.OR.IJOB.EQ.3) THEN
        CALL DPPTRI(UPLO,N,A,INFO)
C     WRITE(MESSAGE,'(A12,I5)') 'DPPTRI INFO=',INFO
C     CALL REPORT_STAT('STATUS','SOLVE','INVER2',' ',MESSAGE,0)
        END IF
      END IF
      if(debug) print *,' info ',info       

      IER=INFO

c Set ier=130 if matrix is singular
      IF (IJOB.EQ.4.AND.INFO.NE.0) IER=130

      END
