      logical function lgood (iflag)

c     return .true. if data is OK for use in a solution

      include '../includes/errflg.h'

      integer iflag,itemp

      if (iflag .gt. 100) then
         itemp = iflag - 200
      else
         itemp = iflag
      endif

      if (itemp .eq. iggood .or.
     .    itemp .eq. igisok .or.
     .    itemp .eq. igrewt .or.
     .    itemp .eq. igbias) then
         lgood = .true.
      else
         lgood = .false.
      endif

      return
      end

