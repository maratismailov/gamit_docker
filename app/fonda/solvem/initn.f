      subroutine initn
c
c     initialize working arrays
c
      implicit real*8(a-h,o-z)
      include 'solvem.fti'
      integer ipa,i,ino
c
c     normal matrix
      ipa = nsit*6
      ino = ipa*(ipa+1)/2
      do 10 i = 1,ino
         anorm(i) = 0.0d0
 10   continue
c
c     chi square
      chi2 = 0.0d0
c
c     right hand term array
      do 20 i = 1,ipa
         bnorm(i) = 0.0d0
 20   continue
c
c     index array
      do 30 i = 1,nsit
         iexc(i) = 0
 30   continue
c
      return
      end
