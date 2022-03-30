      subroutine RVERSN(VERS)
      character*10 getmac,mach
      CHARACTER*80 VERS
      integer nblen

      mach = getmac(1)
      write (vers,5) mach(1:nblen(mach))
   5  format ('RMBIAS v. 8.13 of 91/12/19 12:05:00 (',a,')')

      WRITE(6,10) vers
   10 FORMAT(a80,//)

      WRITE(6,20)
   20 FORMAT (
     .   1X,'This program is for removing biases in 1-way phase',
     .   1X,'except after big gaps (input as in SINCLN)            ',//,
c     .   1X,'New features include:                                 ',/,
c
c version 8.13 Fixed multi-session read from M-file problem in GETFIL
c                     Yehuda Bock 12/19/91


     .   /)
      return
      end





