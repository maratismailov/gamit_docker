      subroutine filnmb( bcfile,xfile,tfile,tfilef )

c     Input file names for BCTOT
c     Last mod for compatibility with OPENB - R King  18 Nov 91

      implicit none

      integer*4 iper
      character*16 bcfile,tfile,xfile,tfilef

C     Enter Broadcast, X-,  and T-file names

c     E-file
      write(6,10)
   10 format(' Enter broadcast ephemeris (E-) file name >: ')
      read(5,11) bcfile
   11 format(a16)


c     X-file (optional)
      write(6,40)
   40 format(' Enter X-file name (if none enter CR) >: ')
      read(5,60) xfile
   60 format(a16)


c     output T-file
      write(6,20)
   20 format(' Enter output T-file name >: ')
      read(5,11) tfile


C     Define earth-fixed T-file name
c     it is opened later, in OPENB
      tfilef=tfile
      iper=index(tfile,'.')
      tfilef(iper-1:iper-1)='e'

      return
      end
