      subroutine bound (ichan1,ichan2,isite1,isite2,jj0,jj1)

c     return the index of the first and last epoch for the N-tuple

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer*4 ichan1,ichan2,isite1,isite2,jj0,jj1

c     first epoch
      jj0 = 1
c     last epoch
      jj1  = maxepc

c     check if the series has any data
      if (ichan1 .ne. 0 .and. isite1 .ne. 0) then
         jj0 = max(int(kk0(ichan1,isite1)),jj0)
         jj1 = min(int(kk1(ichan1,isite1)),jj1)
         if (ichan2 .ne. 0) then
            jj0 = max(int(kk0(ichan2,isite1)),jj0)
            jj1 = min(int(kk1(ichan2,isite1)),jj1)
         endif
         if (isite2 .ne. 0) then
            jj0 = max(int(kk0(ichan1,isite2)),jj0)
            jj1 = min(int(kk1(ichan1,isite2)),jj1)
         endif
         if (ichan2 .ne. 0 .and. isite2 .ne. 0) then
            jj0 = max(int(kk0(ichan2,isite2)),jj0)
            jj1 = min(int(kk1(ichan2,isite2)),jj1)
         endif
      endif

      return
      end
