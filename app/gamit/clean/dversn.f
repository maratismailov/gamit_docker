      subroutine DVERSN(VERS)
      character*10 getmac,mach
      CHARACTER*80 VERS
      integer nblen

      mach = getmac(1)
      write (vers,5) mach(1:nblen(mach))
   5  format ('DBLCLN v. 8.18 of 92/2/23 11:15:00 (',a,')')

      WRITE(6,10) vers
   10 FORMAT(a80,//)

      WRITE(6,20)
   20 FORMAT (
     .   1X,'This program is for cleaning double difference',
     .   1X,'cycle slips                                    ',//,
     .   1X,'New features include:                                 ',/,
c    version 1.1 kurt June 91
     .   1X,'Read M-file name from command line.                   ',/,
     .   1X,'Overwrite old C-files on output.                      ',/,
     .   /)

c   Update comments merged with CVERSN, and call to DVERSN removed
c   from DBLCLN.   King 4/22/92


      return
      end





