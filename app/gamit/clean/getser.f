      logical function getser
     .         (ichan1,ichan2,isite1,isite2,l12b,imode,domarg,
     .          wv,kw,nobs,ngood,iig0,iig1,iim0,iim1)

c     retrieve a time series from the big data arrays
c
c     input:
c       ichan1,ichan2,isite1,isite2
c       l12b
c          -2 for azimuth to sat (radians CW from North)
c          -1 for elevation angle (radians above horizon)
c           1 for L1
c           4 for L2
c           5 for P1
c           6 for P2
c           disable these two to save time and space: kurt 900307
c          12 for L1 multipath phase delay in cycles
c          13 for L2 multipath phase delay in cycles
c
c          16 for station clock offset in seconds
c       imode = 0 for prefit, 1 for postfit
c       domarg = .true. to show marginal points
c     output
c       wv = observations
c       kw = error flags
c       nobs = number of observations
c       ngood = number of good observations
c       iig0,iig1 first and last good observations
c       iim0,iim1 first and last marg observations
c       getser = .true. for no error

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      real*8 y

      integer ngood,nobs
     .,       ichan1,ichan2,isite1,isite2
     .,       ierrcm
     .,       i,itemp
     .,       imode,idat,icount,jcount

      integer iig0,iig1,iim0,iim1,l12b
      logical lgood,lmarg,domarg

      real*8 wv(maxepc)
      integer kw(maxepc)
*******modified  ********START:  for bias_count  m.burc oral 2/26/1992
      character*4 bias_count, char_bias
      common /biascount/bias_count(maxepc)
*******modified  ********    END :  for bias_count  m.burc oral 2/26/1992

      ngood = 0
      getser = .false.
      iig0 = nobs
      iig1 = 1
      iim0 = nobs
      iim1 = 1
      icount = 0

c     what kind of data?
      if (l12b .eq. -2 .or.
     .    l12b .eq. -1 .or.
     .    l12b .eq.  1 .or.
     .    l12b .eq. 12 .or.
     .    l12b .eq. 16) then
         idat = 1
      else if (l12b. eq. 4 .or. l12b .eq. 13) then
         idat = 2
      else if (l12b. eq. 5) then
         idat = 3
      else if (l12b. eq. 6) then
         idat = 4
      else
         print *,'GETSER requested unknown data type = ',l12b
         idat = 1
      endif

c     determine if this observable can be formed
      icount = 0
      jcount = 0
      if (ichan1 .ne. 0 .and. isite1 .ne. 0) then
         jcount = jcount + 1
         if (lambds(isite1,ichan1,idat) .ne. 0) then
            icount = icount + 1
         endif
      endif
      if (ichan1 .ne. 0 .and. isite2 .ne. 0) then
         jcount = jcount + 1
         if (lambds(isite2,ichan1,idat) .ne. 0) then
            icount = icount + 1
         endif
      endif
      if (ichan2 .ne. 0 .and. isite1 .ne. 0) then
         jcount = jcount + 1
         if (lambds(isite1,ichan2,idat) .ne. 0) then
            icount = icount + 1
         endif
      endif
      if (ichan2 .ne. 0 .and. isite2 .ne. 0) then
         jcount = jcount + 1
         if (lambds(isite2,ichan2,idat) .ne. 0) then
            icount = icount + 1
         endif
      endif

      if (icount .eq. jcount) then
         getser = .true.
      else
         getser = .false.
         ngood = 0
         return
      endif

      do 100 i=1,nobs
*********   modified  ******START:  for bias_count  m.burc oral 2/26/1992
         char_bias(1:4) = "0000"
         if (ichan1*isite1 .ne. 0) then
         if( ierr(i,ichan1,isite1) .eq. igbias ) char_bias(1:1) = "1"
         endif
         if (ichan1*isite2 .ne. 0) then
         if( ierr(i,ichan1,isite2) .eq. igbias ) char_bias(2:2) = "1"
         endif
         if (ichan2*isite1 .ne. 0) then
         if( ierr(i,ichan2,isite1) .eq. igbias ) char_bias(3:3) = "1"
         endif
         if (ichan2*isite2 .ne. 0) then
         if( ierr(i,ichan2,isite2) .eq. igbias ) char_bias(4:4) = "1"
         endif
         bias_count(i) = char_bias
*********   modified  ******START:  for bias_count  m.burc oral 2/26/1992
c        figure out the error codes
c        assume that (isite1,ichan1) always has something in it.
c        A bad assumption?
         if (ichan1*ichan2*isite1*isite2 .ne. 0) then
            kw(i) = ierr(i,ichan1,isite1)
            itemp = ierr(i,ichan2,isite1)
            kw(i) = ierrcm( kw(i),itemp )
            itemp = ierr(i,ichan1,isite2)
            kw(i) = ierrcm( kw(i),itemp )
            itemp = ierr(i,ichan2,isite2)
            kw(i) = ierrcm( kw(i),itemp )
         else if (ichan1*ichan2*isite1 .ne. 0) then
            kw(i) = ierr(i,ichan1,isite1)
            itemp = ierr(i,ichan2,isite1)
            kw(i) = ierrcm( kw(i),itemp )
         else if (ichan1*ichan2*isite2 .ne. 0) then
            kw(i) = ierr(i,ichan1,isite2)
            itemp = ierr(i,ichan2,isite2)
            kw(i) = ierrcm( kw(i),itemp )
         else if (ichan1*isite1*isite2 .ne. 0) then
            kw(i) = ierr(i,ichan1,isite1)
            itemp = ierr(i,ichan1,isite2)
            kw(i) = ierrcm( kw(i),itemp )
         else if (ichan2*isite1*isite2 .ne. 0) then
            kw(i) = ierr(i,ichan2,isite1)
            itemp = ierr(i,ichan2,isite2)
            kw(i) = ierrcm( kw(i),itemp )
         else if (ichan1*isite1 .ne. 0) then
            kw(i) = ierr(i,ichan1,isite1)
         else if (ichan1*isite2 .ne. 0) then
            kw(i) = ierr(i,ichan1,isite2)
         else if (ichan2*isite1 .ne. 0) then
            kw(i) = ierr(i,ichan2,isite1)
         else if (ichan2*isite2 .ne. 0) then
            kw(i) = ierr(i,ichan2,isite2)
         else
            kw(i) = ignone
         endif

c        count the number of points we want to use.
         if (lgood(kw(i)) .or. (domarg.and.lmarg(kw(i)))) then
            ngood = ngood + 1
         endif

c        Find first and last good points.
         if (lgood(kw(i))) then
            iig0 = min(iig0,i)
            iig1 = max(iig1,i)
         endif

c        Find first and last marginal points, whether or not
c        we wish to plot them.
         if (lmarg(kw(i)) .or. lgood(kw(i))) then
            iim0 = min(iim0,i)
            iim1 = max(iim1,i)
         endif

c        Think about this.
c        Do we want to always use the same time tag for
c        observables differenced between stations?
         tx(i) = tag(i,isite1)

c        continue processing if point is good or even marginal
         if (lgood(kw(i)) .or. lmarg(kw(i))) then

c           form difference (double, single or none)

c           initialize
            y = 0.0d0

            if (ichan1 .ne. 0 .and. isite1 .ne. 0) then
               if (icount .lt. 4) icount = icount + 1
               if (l12b .eq. -2) y = y + aza(i,ichan1,isite1)
               if (l12b .eq. -1) y = y + ela(i,ichan1,isite1)
               if (l12b .eq.  1) y = y + yl1(i,ichan1,isite1)
               if (l12b .eq.  4) y = y + yl2(i,ichan1,isite1)
               if (l12b .eq.  5) y = y + pr1(i,ichan1,isite1)
               if (l12b .eq.  6) y = y + pr2(i,ichan1,isite1)
ckf            if (l12b .eq. 12) y = y + rd1(i,ichan1,isite1)
ckf            if (l12b .eq. 13) y = y + rd2(i,ichan1,isite1)
               if (l12b .eq. 12) y = 0.0d0
               if (l12b .eq. 13) y = 0.0d0
               if (l12b .eq. 16) y = y + clk(i,ichan1,isite1)
CD              print *,'y11= ',y
            endif
            if (ichan1 .ne. 0 .and. isite2 .ne. 0) then
               if (icount .lt. 4) icount = icount + 1
               if (l12b .eq. -2) y = y - aza(i,ichan1,isite2)
               if (l12b .eq. -1) y = y - ela(i,ichan1,isite2)
               if (l12b .eq.  1) y = y - yl1(i,ichan1,isite2)
               if (l12b .eq.  4) y = y - yl2(i,ichan1,isite2)
               if (l12b .eq.  5) y = y - pr1(i,ichan1,isite2)
               if (l12b .eq.  6) y = y - pr2(i,ichan1,isite2)
ckf            if (l12b .eq. 12) y = y - rd1(i,ichan1,isite2)
ckf            if (l12b .eq. 13) y = y - rd2(i,ichan1,isite2)
               if (l12b .eq. 12) y = 0.0d0
               if (l12b .eq. 13) y = 0.0d0
               if (l12b .eq. 16) y = y - clk(i,ichan1,isite2)
CD              print *,'y12= ',y
            endif
            if (ichan2 .ne. 0 .and. isite1 .ne. 0) then
               if (icount .lt. 4) icount = icount + 1
               if (l12b .eq. -2) y = y - aza(i,ichan2,isite1)
               if (l12b .eq. -1) y = y - ela(i,ichan2,isite1)
               if (l12b .eq.  1) y = y - yl1(i,ichan2,isite1)
               if (l12b .eq.  4) y = y - yl2(i,ichan2,isite1)
               if (l12b .eq.  5) y = y - pr1(i,ichan2,isite1)
               if (l12b .eq.  6) y = y - pr2(i,ichan2,isite1)
ckf            if (l12b .eq. 12) y = y - rd1(i,ichan2,isite1)
ckf            if (l12b .eq. 13) y = y - rd2(i,ichan2,isite1)
               if (l12b .eq. 12) y = 0.0d0
               if (l12b .eq. 13) y = 0.0d0
               if (l12b .eq. 16) y = y - clk(i,ichan2,isite1)
CD              print *,'y21= ',y
            endif
            if (ichan2 .ne. 0 .and. isite2 .ne. 0) then
               if (icount .lt. 4) icount = icount + 1
               if (l12b .eq. -2) y = y + aza(i,ichan2,isite2)
               if (l12b .eq. -1) y = y + ela(i,ichan2,isite2)
               if (l12b .eq.  1) y = y + yl1(i,ichan2,isite2)
               if (l12b .eq.  4) y = y + yl2(i,ichan2,isite2)
               if (l12b .eq.  5) y = y + pr1(i,ichan2,isite2)
               if (l12b .eq.  6) y = y + pr2(i,ichan2,isite2)
ckf            if (l12b .eq. 12) y = y + rd1(i,ichan2,isite2)
ckf            if (l12b .eq. 13) y = y + rd2(i,ichan2,isite2)
               if (l12b .eq. 12) y = 0.0d0
               if (l12b .eq. 13) y = 0.0d0
               if (l12b .eq. 16) y = y + clk(i,ichan2,isite2)
CD              print *,'y22= ',y
            endif

            wv(i)=y

CD           print *,'GETSER ',i,kw(i),wv(i)
         endif
  100 continue



CD     print *,'GETSER iig0,iig1',iig0,iig1
CD     print *,'GETSER ngood',ngood
      return

      end




