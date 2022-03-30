      subroutine reodrd(nlen,ist)
c
c     Reorder D-optimal by length of baselines
c
c     ist   : shift index
c     nlen  : baseline number
c                              
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer nlen,ist,i1,i2,jj,i,j

        do 30 j = 1,32000
           jj = 0
           do 20 i = 2,nlen
              if(ipntdt(i).ge.ipntdt(i-1)) go to 20
              i1 = ipntdt(i)
              ipntdt(i) = ipntdt(i-1)
              ipntdt(i-1) = i1
c
              i2 = (i-1)*ist
              do 10 i1 = i2+1,i2+ist
                 jj = ipntd(i1)
                 ipntd(i1) = ipntd(i1-ist)
                 ipntd(i1-ist) = jj
 10           continue
 20        continue
           if(jj.eq.0) go to 40
 30     continue

 40   continue

      return
      end
