c************************************************************************
      subroutine findslip_scn (isite1,isite2,ichan1,ichan2,idd,
     .                      nsat,nobs,ii0,ii1,rms1,rms2,rms3)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer         isite1,isite2,ichan1,ichan2,idd
      integer         nsat
      integer         ii0,ii1,islip
      integer         nseq,nslip,nobs
      logical         slip
      real*8          yadd(2)
      real*8          ytot
      real*8          rms1,rms2,rms3

c     maximum number of automatic patching attempts and slip
      integer   maxslip
      parameter (maxslip  = 25)

      character*4 bias_count
      common /biascount/bias_count(maxepc)

c     **************************************************************
      nslip = 0
      islip = ii0
      slip = .true.
      ytot = 0.d0

      do 100 while (slip .and. nslip. lt. maxslip)
           nslip = nslip + 1

           call search_slip(
     .                   isite1,isite2,ichan1,ichan2,nobs,
     .                   ii0,ii1,islip,nseq,yadd,slip)

           if (slip)
     .          write(idd,1000)islip,rms2,ichan2,ichan1,
     .                    isite2,isite1,yadd(1),yadd(2),
     .                                 bias_count(islip)

  100 continue

 1000 format (i5,4x,f7.2,3x,4i5,4x,f6.2,4x,f6.2,7x,a4)

      return
      end
c
c************************************************************************

