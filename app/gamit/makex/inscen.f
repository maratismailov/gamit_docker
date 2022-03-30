      integer*4 function inscen(nprn,isatcount,satarray)
c     return the number of the channel in the array

      include '../includes/dimpar.h'
      include '../includes/makex.h'

C     sat PRN number to check
      integer*4        nprn
C     number of sats in scenario
      integer*4        isatcount
C     list of scenario
      integer*4        satarray(maxsat)
      integer*4     i

      inscen = 0

      if (isatcount .le. maxsat) then
         do 20 i=1,isatcount
            if (nprn .eq. satarray(i) .and. nprn .gt. 0) then
               inscen = i
               return
            else
               inscen = -1
            endif
 20      continue
      else
         inscen = -1
      endif

      return
      end
