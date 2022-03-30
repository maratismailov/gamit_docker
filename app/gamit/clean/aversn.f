      subroutine AVERSN(VERS)
      character*10 getmac,mach
      CHARACTER*80 VERS
      integer nblen

      mach = getmac(1)
      write (vers,5) mach(1:nblen(mach))
   5  format ('SINCLN v. 8.18 of 92/02/23 11:15:00 (',a,')')

      WRITE(6,10) vers
   10 FORMAT(a80,//)

      WRITE(6,20)
   20 FORMAT (
     .   1X,'This program is for automatic cycle-slip cleaning',
     .   1X,'                                                 ',//,
c     .   1X,'New features include:                                 ',/,
c
c   Update comments merged with CVERSN, and call to AVERSN removed
c   from SINCLN.   King 4/22/92.
c
c

     .   /)
      return
      end





