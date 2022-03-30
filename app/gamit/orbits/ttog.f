      Program TTOG

c     R. King  3 April 1987
c     Modified to increase non-grav. force coeffs. to 3 -- rwk 880808
c     Incorporated into orbits directory, with mods -- rwk 911227
c     Reduced to a call to GMAKE - rwk 940719

      implicit none

      include '../includes/dimpar.h'

      integer*4 iscrn,itsat(maxsat)
c*      integer*4 iutin

      character*6  prog
      character*16 tname,gname
      character*120 version
      character*256 message

      data iscrn/6/
c*      data iutin/10/
C
c          Print the version and machine

      call oversn(version)
      write(iscrn,'(a)')' '
      write(message,5)version
    5 format('Started TTOG ',a120)
      call report_stat('STATUS','TTOG','ttog',' ',message,0)

C          Open the T-File

      write(iscrn,'(/,1x,a)')'Enter the input T-file name: '
      read(5,'(a16)') tname

c          Read the T-file and write the G-file

      prog = 'TTOG  '
      call gmake( prog,tname,gname,itsat )

      call report_stat('STATUS','TTOG','ttog',' ',
     .'Normal end in TTOG',0)

      stop
      end

