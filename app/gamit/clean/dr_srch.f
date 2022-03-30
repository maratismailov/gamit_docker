c************************************************************************
      subroutine search_slip (
     .                   isite1,isite2,ichan1,ichan2,nobs,
     .                   ii0,ii1,islip,icyc,yadd,slip)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer         isite1,isite2,ichan1,ichan2,nobs
      integer         ii0,ii1,islip,lastslip,icyc
      integer         nlong,ifind(20)
      logical         lgood,backp,frontp,getgn
      logical         slip,slipwl
      real*8          dwl,yadd(2)
      real            epsilon(5),gnoise(5)

c     lower bound of signal to noise ratio
      parameter (nlong = 5)

      lastslip = islip
      islip = islip - 1
      slip = .false.
      getgn = .true.

      do 100 while((.not. slip) .and. (islip .lt. ii1))
           islip = islip + 1

           if (islip .le. 0) islip = 1
           if (lgood(kww(islip)) .and. (.not. slip)) then
                call find5(
     .            lastslip,islip,ii0,ii1,ifind,nlong,backp,frontp)

                call get_eps(
     .                epsilon,islip,ifind,backp,frontp)
                call dr_check_slip(ii0,ii1,
     .                islip,epsilon,ifind,backp,frontp,slip)

                if ((.not. slip) .and. (islip .gt. ii0))
     .                 call check_global(
     .                    ii0,ii1,islip,slip,getgn,gnoise,ichan2)

c               check in WL for slips
                if (dabs(pc2(islip)) .gt. 0.1) then
                     call check_wl(ii0,ii1,islip,slipwl,dwl)
                     if(slipwl) slip = .true.
                else
                     dwl = 0.d0
                endif
           endif

           if(slip .and. (islip .eq. lastslip)) slip = .false.

c          confirm that PATCH will improve the situation
           if(slip) call confirm_slip(
     .               isite1,isite2,ichan1,ichan2,nobs,
     .               ii0,ii1,islip,yadd,slip,slipwl,ifind,dwl)

  100 continue

      return
      end
c************************************************************************
      subroutine find5(
     .        lastslip,ipoint,ii0,ii1,ifind,nlong,backp,frontp)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         ipoint,ii0,ii1,lastslip
      integer         i,j,nlong,ifind(20)
      integer         jpoint,jtot,jback,jfront
      logical         lgood,backp,frontp

      backp  = .false.
      frontp = .false.
      jpoint = 0 

c     see if there are enough backward point to start a backward search
      if ((ipoint - lastslip) .gt. nlong) then

          i = ipoint
          j = 0

          do 100 while ((j .lt. nlong) .and. (i .gt. lastslip))
              i = i - 1
              if(lgood(kww(i))) then
                  j = j + 1
                  ifind(j) = i
              endif
  100     continue

c         signal the end of the sequance
          if (j .ge. nlong) then
              backp = .true.
              ifind(j + 1) = 0
          endif
      endif

c     see if there are enough forward point to start a forward search
      if ((.not. backp) .and. ((ii1 - ipoint) .gt. nlong)) then

          i = ipoint
          j = 0

          do 200 while ((j .lt. nlong) .and. (i .lt. ii1))
              i = i + 1
              if(lgood(kww(i))) then
                  j = j + 1
                  ifind(j) = i
              endif
  200     continue

c         signal the end of the sequance
          if (j .ge. nlong) then
              frontp = .true.
              ifind(j + 1) = 0
          endif

      endif

c     see if there are enough backward point to start a forward search
c     ignoring the previous slip
      if (((.not. backp) .and. (.not. frontp)) .and.
     .             ((ipoint - ii0) .gt. nlong)) then

          i = ipoint
          j = 0

          do 250 while ((j .lt. nlong) .and. (i .gt. ii0))
              i = i - 1
              if(lgood(kww(i))) then
                  j = j + 1
                  ifind(j) = i
              endif
  250     continue

c         signal the end of the sequance
          if (j .ge. nlong) then
              backp = .true.
              ifind(j + 1) = 0
          endif
      endif

c     see if there are enough point for a middle search
      if ((.not. backp) .and. (.not. frontp)) then
          j = 0

          do 300 i= ii0,ii1
              if(lgood(kww(i))) then
                  j = j + 1
                  ifind(j) = i
                  if (i .eq. ipoint) jpoint = j
              endif

  300     continue
          jtot = j

          jback  = jpoint - 1
          jfront = jtot - jpoint

          if (jback .ge. jfront) then
               ifind(jpoint) = 0
               backp = .true.
          else
               do 400 j=1,jfront
                  ifind(j) = ifind(jpoint + j)
  400          continue
               ifind(jfront + 1) = 0
               frontp = .true.
          endif
      endif


      return
      end
c************************************************************************
      subroutine get_eps(
     .           epsilon,islip,ifind,backp,frontp)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         i,j,k,islip,ifind(20)
      logical         backp,frontp
      real            epsilon(5),signal
      real*8          ymin(5),ymax(5),yk(5)

c     lower bound of signal to noise ratio

      signal = 5.0d0

c     initialize min/max arrays
      do 100 i=1,4
           ymax(i) = -1.0d25
           ymin(i) =  1.0d25
  100 continue

      j = 0
      k = 1

c     find min/max
      do 300 while (k .ne. 0)

          j = j + 1
          k = ifind(j)

          if (k .ne. 0) then

              yk(1) = wl1(k)
              yk(2) = wl2(k)
              yk(3) = wl1(k)
     .              - faclr(jsat1)*(wl2(k)-gear(jsat1)*wl1(k))
              yk(4) = wl2(k) - gear(jsat1)*wl1(k)

               do 200 i=1,4
                  ymax(i) = dmax1(ymax(i),yk(i))
                  ymin(i) = dmin1(ymin(i),yk(i))
  200         continue

          endif

  300 continue

      do 400 i=1,4
           epsilon(i) = (ymax(i) - ymin(i)) * signal
  400 continue

      return
      end
c************************************************************************
      subroutine dr_check_slip(ii0,ii1,
     .           islip,epsilon,ifind,backp,frontp,slip)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         ii0,ii1,islip,ifind(20)
      integer         j,j1,j2,k
      logical         backp,frontp,slip
      real            epsilon(5)
      real*8          y1(5),y2(5)
      real*8          yobs(5),yexp(5),dy,slop


c     get observed values
      yobs(1) = wl1(islip)
      yobs(2) = wl2(islip)
      yobs(3) = wl1(islip)
     .      - faclr(jsat1)*(wl2(islip)-gear(jsat1)*wl1(islip))
      yobs(4) = wl2(islip) - gear(jsat1)*wl1(islip)

c     get location of first and last points
      j1 = ifind(1)

      j = 0
      k = 1
      do while (k .ne. 0)
         j = j + 1
         k = ifind(j)
      end do
      j2 = ifind(j-1)

      if(j1 .eq. j2) then
         j1 = ii0
         j2 = ii1
      endif

c     get values of first and last points
      y1(1) = wl1(j1)
      y1(2) = wl2(j1)
      y1(3) = wl1(j1) - faclr(jsat1)*(wl2(j1)-gear(jsat1)*wl1(j1))
      y1(4) = wl2(j1) - gear(jsat1)*wl1(j1)

      y2(1) = wl1(j2)
      y2(2) = wl2(j2)
      y2(3) = wl1(j2) - faclr(jsat1)*(wl2(j2)-gear(jsat1)*wl1(j2))
      y2(4) = wl2(j2) - gear(jsat1)*wl1(j2)

      j = 0
      do 100 while ((.not.slip) .and. (j .lt. 4))
           j= j + 1
           slop = (y2(j) - y1(j))/(j2 - j1)

           yexp(j) = y2(j) + slop*(islip-j2)

c          check if the difference between observed and expected
c          are above the noise level
           dy = abs(yobs(j)-yexp(j))

           if (dy .gt. epsilon(j)) then
                slip = .true.
                if (frontp) islip = ifind(1)
           endif

  100 continue


      return
      end
c************************************************************************
      subroutine confirm_slip(
     .           isite1,isite2,ichan1,ichan2,nobs,
     .           ii0,ii1,islip,yadd,slip,slipwl0,ifind,dwl0)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         isite1,isite2,ichan1,ichan2
      integer         ii0,ii1,nobs,islip,ifind(20)
      integer         j,j1,j2,k
      logical         slip,slipwl0,slipwl1,unwt
      real*8          yn,yadd(2),oldrms,newrms
      real*8          ynew,yold,yexp,slop
      real*8          ylc1,ylc2,dlc
      real*8          dwl0,dwl1


c     buffer the series for undo
      call vcopy (kww,wl1,kbb,bl1,nobs)
      call vcopy (kww,wl2,kbb,bl2,nobs)

      unwt = .false.
c      if(ichan2 .eq. 0) unwt = .true.
c     determines the slip amount by patching
      call dr_patch (
     .         isite1,isite2,ichan1,ichan2,unwt,
     .         yadd,isite1,ichan1,nobs,ii0,ii1,islip)

c     define what is a significant patch
      if(ichan2 .eq. 0) then
           yn = 2.6
           yn = 0.1
      else
           yn = 0.1
      endif

c     CHECK 1
c     check if patch change significantly the slip
      if ((dabs(yadd(1)) .lt. yn) .and.
     .    (dabs(yadd(2)) .lt. yn)) slip = .false.

c     CHECK 2
c     checks if patch contribute to LC
      if (slip) then
           dlc = yadd(1) - faclr(jsat1)*(yadd(2)-gear(jsat1)*yadd(1))
           if (dabs(dlc) .lt. 0.1) slip = .false.
      endif

c     CHECK 3
c     check if LC is improved
      if (slip) then
c         get location of first and last points
          j1 = ifind(1)

          j = 0
          k = 1
          do while (k .ne. 0)
             j = j + 1
             k = ifind(j)
          end do
          j2 = ifind(j-1)

          if(j1 .eq. j2) then
             j1 = ii0
             j2 = ii1
          endif

c         get values of first and last points
          ylc1 = wl1(j1) - faclr(jsat1)*(wl2(j1)-gear(jsat1)*wl1(j1))
          ylc2 = wl1(j2) - faclr(jsat1)*(wl2(j2)-gear(jsat1)*wl1(j2))

          slop = (ylc2 - ylc1)/(j2 - j1)
          yexp = ylc2 + slop*(islip-j2)

c         get old and new observed values
          ynew = wl1(islip) - faclr(jsat1)*(wl2(islip)
     .               - gear(jsat1)*wl1(islip))
          yold = bl1(islip) - faclr(jsat1)*(bl2(islip)
     .               - gear(jsat1)*bl1(islip))

c         checks the difference between new/old observed and expected
c         if slip increases because of patching skip that slip
          if(dabs(yexp-ynew) .ge. dabs(yexp-yold)) slip = .false.
      endif

c     CHECK 4
c     checks if slip is reduecd at WL observable
      if((slip) .and. (dabs(pc2(islip)).gt.0.1)) then
            if(dwl0 .le. 5) then
                call check_wl(ii0,ii1,islip,slipwl1,dwl1)
                if (dwl1 .ge. dwl0) slip = .false.
                if (slipwl1 .and. (.not. slipwl0)) slip = .false.
            endif
      endif

c     CHECK 5
c     checks if LC-RMS increases
      if((slip) .and. (isite2 .ne. 0)) then
             call cal_rms
     .            (ichan1,ichan2,isite1,isite2,
     .             newrms,kww,wl1,wl2,ii0,ii1)
             call cal_rms
     .            (ichan1,ichan2,isite1,isite2,
     .             oldrms,kbb,bl1,bl2,ii0,ii1)
             if (newrms .gt. oldrms) slip = .false.
       endif



c     return to the original series
      call vcopy (kbb,bl1,kww,wl1,nobs)
      call vcopy (kbb,bl2,kww,wl2,nobs)

      return
      end
c************************************************************************
      subroutine check_global (
     .                 ii0,ii1,islip,slip,getgn,gnoise,ichan2)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer         ii0,ii1,islip,ichan2
      integer         i,i1,npoints
      real            gnoise(5),signal,dy(5)
      real*8          ylc0,ylc1,ylg0,ylg1
      logical         lgood,slip,getgn,good

c     one-way value
      if(ichan2 .eq. 0) then
           signal = 8.0d0
      else
           signal = 4.0d0
      endif


c     calculate global noise
      if(getgn) then

           getgn = .false.
           npoints = 0
           gnoise(1) = 0.d0
           gnoise(2) = 0.d0
           gnoise(3) = 0.d0
           gnoise(4) = 0.d0

           i1 = ii0
           do 200 i=ii0+1,ii1

               if (lgood(kww(i))) then

                    if((i-i1) .lt. biggap) then

                       ylc0 = wl1(i) - faclr(jsat1)*(wl2(i)
     .                               - gear(jsat1)*wl1(i))
                       ylc1 = wl1(i1) - faclr(jsat1)*(wl2(i1)
     .                                - gear(jsat1)*wl1(i1))
                       ylg0 = wl2(i) - gear(jsat1)*wl1(i)
                       ylg1 = wl2(i1) - gear(jsat1)*wl1(i1)

                       gnoise(1) = gnoise(1) + dabs(wl1(i)-wl1(i1))
                       gnoise(2) = gnoise(2) + dabs(wl2(i)-wl2(i1))
                       gnoise(3) = gnoise(3) + dabs(ylc0-ylc1)
                       gnoise(4) = gnoise(4) + dabs(ylg0-ylg1)

                       npoints = npoints + 1
                    endif

                    i1 = i
               endif
  200      continue

           if (npoints .gt. 0) then
                gnoise(1) = gnoise(1)*signal/float(npoints)
                gnoise(2) = gnoise(2)*signal/float(npoints)
                gnoise(3) = gnoise(3)*signal/float(npoints)
                gnoise(4) = gnoise(4)*signal/float(npoints)
           endif

      endif


c     check if dy is greater than global noise
      i  = islip
      i1 = islip
      good = .false.
      do 300 while (.not. good)
          i1 = i1 -1
          if (lgood(kww(i1))) good = .true.
  300 continue

      if((i-i1) .lt. biggap) then
          ylc0 = wl1(i) - faclr(jsat1)*(wl2(i)-gear(jsat1)*wl1(i))
          ylc1 = wl1(i1) - faclr(jsat1)*(wl2(i1)-gear(jsat1)*wl1(i1))
          ylg0 = wl2(i) - gear(jsat1)*wl1(i)
          ylg1 = wl2(i1) - gear(jsat1)*wl1(i1)

          dy(1) =  dabs(wl1(i)-wl1(i1))
          dy(2) =  dabs(wl2(i)-wl2(i1))
          dy(3) =  dabs(ylc0-ylc1)
          dy(4) =  dabs(ylg0-ylg1)

          if ((dy(1) .gt. gnoise(1)) .or.
     .        (dy(2) .gt. gnoise(2)) .or.
     .        (dy(3) .gt. gnoise(3)) .or.
     .        (dy(4) .gt. gnoise(4))) slip = .true.
      endif

      return
      end
c************************************************************************
      subroutine find_gap (islip,igap)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer         islip,igap,ipoint
      logical         lgood,still


      still = .true.
      ipoint = islip
      igap = 0

      do 100 while (still .and. (ipoint .gt. 1))
            ipoint  = ipoint - 1

            if(.not.lgood(kww(ipoint))) then
                 igap = igap + 1
            else
                 still = .false.
            endif
  100 continue

      return
      end
c************************************************************************

