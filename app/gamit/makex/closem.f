      subroutine closem
c     close the files open for current day

      implicit none

      include '../includes/makex.h'

      if (qclock) close(uclock)
      if (qxfile) close(uxfile)
      if (qsvclk) close(usvclk)
      if (qsp3) close(usp3)
      if (qnav  ) close(unav)
      if (qrinex) close(urinex)  
      if (qanthi) close(uanthi)

      return
      end

