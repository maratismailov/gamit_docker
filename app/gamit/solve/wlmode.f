
      Subroutine WLMODE(free,mode)
c
c     free or fix WL bias parameters
c     mode = 1: fix WL biases
c     mode = 2: free WL biases
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'

      integer free(maxprm),mode,ifirst,i1,i2,i3,isn,i

      ifirst = idxb(1)
      if (ifirst.lt.0) ifirst = -ifirst
      i1 = ifirst
      i2 = i1-1
      i3 = 0
c
         isn = l1bias
         i3 = i3+isn
         i1 = i1+isn
         i2 = i2+isn*2
         do 20 i = i1,i2
            i3 = i3+1
            if (mode.eq.1) then
               free(i) = 0
               idxb(i3) = -i
            endif
            if (mode.eq.2) then
               free(i) = 1
               idxb(i3) = i
            endif
 20      continue
         i1 = i1+isn
c
c     recalculate the live parameters
c
      i1 = 0
      do 30 i = 1,ntpart
         if (free(i).eq.0) goto 30
         i1 = i1+1
 30   continue
      nlive = i1
c      print *,'End of WL mode new nlive ',nlive
c
      return
      end


