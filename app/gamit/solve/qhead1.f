c
      subroutine qhead1(mode)
c
c     mode = 1: print effective observation number
c     mode = 2: print parameter number information

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer mode,nded

      if (mode.eq.1) then
         if( logprt ) write(6,70) nones,nd1obs,nd2obs
         write(10,70) nones,nd1obs,nd2obs
      endif
c
      if (mode.eq.2) then
         nded = ntpart-nlive
         if( logprt ) write(6,80) ntpart,nlive,nded
      endif
c
 70   format(/,' Number of good oneway phases:',i7,/
     1  ' Number of single differences:',i7,/
     2  ' Number of double differences:',i7,/)
c
 80   format(/,' Total number of parameters :',i5,/,
     1         '  Number of live parameters :',i5,/,
     2         '  Number of dead parameters :',i5,/)

      return
      end


