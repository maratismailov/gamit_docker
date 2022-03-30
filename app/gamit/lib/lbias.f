      logical function lbias (iflag)

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

      if (itemp .eq. igbias) then
         lbias = .true.
      else
         lbias = .false.
      endif

      return
      end

