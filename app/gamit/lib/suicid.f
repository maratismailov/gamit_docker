      SUBROUTINE SUICID(STRING)
      character*(*) string
c      character*1 letter,lowerc
c
c     UNIX version
C
C     AS THE NAME IMPLIES!! COME TO THIS SUBROUTINE ON AN ABORT
C     S.A. GOUREVITCH   6/81
c
c     To get a good trace back, it is necessary to dump core.
C
      write (6,*) string
c      write (6,'(a,1x,$)') 'SUICID: Do you want to dump core?'
c      read (5,*) letter
c
c      if (lowerc(letter) .ne. 'n') then
c         write (6,*) 'SUICID: Please wait while I dump core...'
c         call abort
c      else
         stop 'SUICID'
c      endif

      END
