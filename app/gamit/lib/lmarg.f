      logical function lmarg (iflag)

c     return .true. if data is marginally OK
c     these flags are not used in solutions, but
c     should be plotted with a different symbol

      include '../includes/errflg.h'

      integer iflag,itemp

      if (iflag .gt. 100) then
         itemp = iflag - 200
      else
         itemp = iflag
      endif

      if (itemp .eq. igloel .or.
     .    itemp .eq. igunwt .or.
     .    itemp .eq. igoutl .or.
     .    itemp .eq. iglamp .or.
     .    itemp .eq. ig2few) then
         lmarg = .true.
      else
         lmarg = .false.
      endif

      return
      end

