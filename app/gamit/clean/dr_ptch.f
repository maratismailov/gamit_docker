c************************************************************************
      subroutine dr_patch (
     .                   isite1,isite2,ichan1,ichan2,unwt,
     .                   yadd,jsite,jchan,nobs,ii0,ii1,islip)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer         nobs,iibl,iibr,ii0,ii1,islip
      integer         isite1,isite2,ichan1,ichan2
      integer         jsite,jchan,idd
      real*8          yadd(2)
      logical         lpatch,doit,ldumb,unwt

      integer          nmlfit
      integer          indfl(15), indfr(15), errl,errr
c     parameter ( mxfitl = 15 )

c     for compatibility with other callers of PATCH
      integer*2 iwin(4)
      real*8 confid	


C     There are 2 flavors of patch, depending on observable:
c     L1,L2:  simple linear interpolation
c     LC,LG:  Use Herring's PATCH routine to fix L1 and L2 using LC and LG
c
c     Assume user wants PATCH to apply to interval from iibl to iibr inclusive
c     and that this has been determined from the 2 brackets

      lpatch = .true.
      iibl = islip
      iibr = nobs
      yadd(1) = 0.d0
      yadd(2) = 0.d0
      idd = ichan1*ichan2*isite1*isite2
      idd = 0


c     check if there are some trouble point near margin - if so delete them
      if (ichan2 .eq. 0) then

           lpatch = .false.
           nmlfit = 3
           call scanlgood(iibl,kww, nobs, nmlfit, -1, indfl, errl )
           call scanlgood(iibl,kww, nobs, nmlfit, +1, indfr, errr )

c          unweight at the left end
           if( errl.ne.0 .and. errr.eq.0 ) then

               call unwt_left (
     .                   isite1,isite2,ichan1,ichan2,
     .                   yadd,jsite,jchan,nobs,ii0,ii1,islip)

c          unweight at the right end
           else if( errl.eq.0 .and. errr.ne.0 ) then

               call unwt_right (
     .                   isite1,isite2,ichan1,ichan2,
     .                   yadd,jsite,jchan,nobs,ii0,ii1,islip)

           else if (unwt) then
c              check if there are trouble point near the patching point
c              and unweight if necessary.
               call check_unwt (islip,isite1,ichan1,lpatch)

           else
               lpatch = .true.
           endif
      else
           lpatch = .true.
      endif


      if (lpatch) then

c         use if possible 5 point at each side
          nmlfit = 5

          call scanlgood(iibl,kww, nobs, nmlfit, -1, indfl, errl )
          call scanlgood(iibl,kww, nobs, nmlfit, +1, indfr, errr )

c         New call sequence for patch: kurt 910620
          doit = .false.
          if( errl.eq.0 .and. errr.eq.0 ) then
             call patch (iibl,iibr,ii0,ii1,iwin,jsite,jchan,
     .               idd,nobs,ldumb,nmlfit,yadd,confid,doit)
          else
             nmlfit = 3
             call patch (iibl,iibr,ii0,ii1,iwin,jsite,jchan,
     .               idd,nobs,ldumb,nmlfit,yadd,confid,doit)
          endif
      endif


 1000 format (4x,'------------- Unweighted data point ',2i6)
 2000 format (4x,'------------- Unflagged   data point ',2i6)
 3000 format (4x,'------------- Flagged     data point ',2i6)

      return
      end
c************************************************************************
      subroutine unwt_left (
     .                   isite1,isite2,ichan1,ichan2,
     .                   yadd,jsite,jchan,nobs,ii0,ii1,islip)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer         nobs,ii0,ii1,islip
      integer         isite1,isite2,ichan1,ichan2
      integer         jsite,jchan,i,iichan1,iisite1
      real*8          yadd(2)
      logical         lgood,lbias,putbias,found
               
c     initialization to satisfy some compilers
      iisite1 = 0
      iichan1 = 0

c     unweight at the left end

      if (ichan2 .ne. 0) then

           if (ii0 .gt. 5*inext) then
                 i = 0
           else         
                 i = 5
           endif
           found = .false.
c          look for the shorter sequance
           do 100 while ((.not. found) .and. (i .le. 4))

               i = i + 1
               iichan1 = ichan1
               iisite1 = isite1
               if ((i .eq. 2) .or. (i .eq. 4)) iichan1 = ichan2
               if ((i .eq. 3) .or. (i .eq. 4)) iisite1 = isite2

          if  ((.not. lgood(ierr(ii0-1*inext,iichan1,iisite1))) .and.
     .         (.not. lgood(ierr(ii0-2*inext,iichan1,iisite1))) .and.
     .         (.not. lgood(ierr(ii0-3*inext,iichan1,iisite1))) .and.
     .         (.not. lgood(ierr(ii0-4*inext,iichan1,iisite1))) .and.
     .         (.not. lgood(ierr(ii0-5*inext,iichan1,iisite1))))
     .                                                 found = .true.

               if ((i .eq. 2) .and. (isite2 .eq. 2)) i = 5

  100      continue

       else
           iichan1 = ichan1
           iisite1 = isite1
           found = .true.
       endif

       if (found) then
           putbias = .false.
           do 400 i=ii0,islip-inext,inext
               if(lbias(ierr(i,iichan1,iisite1)))putbias = .true.
               if(lgood(kww(i))) then
                     ierr(i,iichan1,iisite1) = igunwt
                     kww(i) = igunwt
                     kbb(i) = igunwt
                     write(6,1000)iisite1,iichan1,i
               endif
  400      continue

           ii0 = islip

           if (putbias) then
                ierr(ii0,iichan1,iisite1) = igbias
                kww(ii0) = igbias
                kbb(ii0) = igbias
               write(6,2000)iisite1,iichan1,ii0
           endif
       endif


 1000 format (4x,'------------- Unweighted data point ',3i6)
 2000 format (4x,'------------- Flagged     data point ',3i6)

      return
      end
c************************************************************************
      subroutine unwt_right (
     .                   isite1,isite2,ichan1,ichan2,
     .                   yadd,jsite,jchan,nobs,ii0,ii1,islip)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer         nobs,ii0,ii1,islip
      integer         isite1,isite2,ichan1,ichan2
      integer         jsite,jchan,i,iichan1,iisite1
      real*8          yadd(2)
      logical         good,lgood,lbias,putbias,found

             
c     initialization to satisfy some compilers
      iisite1 = 0
      iichan1 = 0

c     unweight at the left end

      if (ichan2 .ne. 0) then

           if (ii1 .lt. nobs-6*inext) then
                 i = 0
           else
                 i = 5
           endif
           found = .false.
c          look for the shorter sequance
           do 100 while ((.not. found) .and. (i .le. 4))

               i = i + 1
               iichan1 = ichan1
               iisite1 = isite1
               if ((i .eq. 2) .or. (i .eq. 4)) iichan1 = ichan2
               if ((i .eq. 3) .or. (i .eq. 4)) iisite1 = isite2

          if  ((.not. lgood(ierr(ii1+1*inext,iichan1,iisite1))) .and.
     .         (.not. lgood(ierr(ii1+2*inext,iichan1,iisite1))) .and.
     .         (.not. lgood(ierr(ii1+3*inext,iichan1,iisite1))) .and.
     .         (.not. lgood(ierr(ii1+4*inext,iichan1,iisite1))) .and.
     .         (.not. lgood(ierr(ii1+5*inext,iichan1,iisite1))))
     .                                            found = .true.

               if ((i .eq. 2) .and. (isite2 .eq. 2)) i = 5

  100      continue

       else
           iichan1 = ichan1
           iisite1 = isite1
           found = .true.
       endif

       if (found) then
           putbias = .false.
           do 400 i=islip,ii1,inext
               if(lbias(ierr(i,iichan1,iisite1)))putbias = .true.
               if(lgood(kww(i))) then
                     ierr(i,iichan1,iisite1) = igunwt
                     kww(i) = igunwt
                     kbb(i) = igunwt
                     write(6,1000)iisite1,iichan1,i
               endif
  400      continue

           good = .false.
           i = islip

           do 500 while ((.not.good) .and. (i .gt. ii0))
             i = i - 1
             if(lgood(kww(i))) good = .true.
  500      continue
           ii1 = i

           if (putbias) then
               do 600 while ((.not.good) .and. (i .lt. nobs))
                    i = i + 1
                    if(lgood(ierr(i,iichan1,iisite1))) then
                        good = .true.
                        if(.not. lbias(ierr(i,iichan1,iisite1)))
     .                                                     then
                            ierr(i,iichan1,iisite1) = igbias
                            kww(i) = igbias
                            kbb(i) = igbias
                            write(6,2000)iisite1,iichan1,i
                        endif
                    endif
  600           continue
           endif
       endif


 1000 format (4x,'------------- Unweighted data point ',3i6)
 2000 format (4x,'------------- Flagged     data point ',3i6)

      return
      end
c************************************************************************
      subroutine dr_save (
     .                 isite1,isite2,ichan1,ichan2,
     .                 islip,jsite,jchan,nobs,yadd)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

c     FUNCTIONS:
      real*8 rnhalf
      real*8 rnqrtr

      integer         isite1,isite2,ichan1,ichan2
      integer         islip,jsite,jchan,nobs,i
      integer*2       iabs2

      real*8 yadd(2)

c     Change yval for all points, not just good ones
c     Allow down to quarter cycle slips for lambda = +/- 2.
c     Allow up to half cycle slips for lambda = +/- 2.

      if (iabs2(lambds(jsite,jchan,1)) .eq. 1) then
            yadd(1) = rnhalf(yadd(1))
      else if (iabs2(lambds(jsite,jchan,1)) .eq. 2) then
            yadd(1) = rnqrtr(yadd(1))
      else
            yadd(1) = 0.0d0
      endif

      if (iabs2(lambds(jsite,jchan,2)) .eq. 1) then
            yadd(2) = rnhalf(yadd(2))
      else if (iabs2(lambds(jsite,jchan,2)) .eq. 2) then
            yadd(2) = rnqrtr(yadd(2))
      else
            yadd(2) = 0.0d0
      endif

      if((jsite .eq. isite1 .and. jchan .eq. ichan2) .or.
     .   (jsite .eq. isite2 .and. jchan .eq. ichan1)) then
             yadd(1) = -yadd(1)
             yadd(2) = -yadd(2)
      endif

      do 3200 i = islip, nobs

          yl1(i,jchan,jsite) = yl1(i,jchan,jsite) + yadd(1)
          yl2(i,jchan,jsite) = yl2(i,jchan,jsite) + yadd(2)

 3200 continue

c     make sure there is bias at patching point
c      ierr(islip,jchan,jsite) = igbias


      return
      end
c************************************************************************
      subroutine check_unwt (islip,isite1,ichan1,lpatch)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer         islip,isite1,ichan1,i,i1,i0
      logical         lgood,lbias,lpatch,unwt


c     check to the left of the slip
      if(lbias(kww(islip-inext))) then
           ierr(islip-inext,ichan1,isite1) = igunwt
           kww(islip-inext) = igunwt
           kbb(islip-inext) = igunwt
           write(6,1000)ichan1,islip-inext

      else if(lgood(kww(islip-inext))) then

           if(.not. lgood(kww(islip-2*inext))) then
               ierr(islip-inext,ichan1,isite1) = igunwt
               kww(islip-inext) = igunwt
               kbb(islip-inext) = igunwt
               write(6,1000)ichan1,islip-inext

           else if(lbias(kww(islip-2*inext))) then
               ierr(islip-inext,ichan1,isite1) = igunwt
               kww(islip-inext) = igunwt
               kbb(islip-inext) = igunwt
               write(6,1000)ichan1,islip-inext

               ierr(islip-2*inext,ichan1,isite1) = igunwt
               kww(islip-2*inext) = igunwt
               kbb(islip-2*inext) = igunwt
               write(6,1000)ichan1,islip-2*inext
          endif
      endif

c     check to the right of the slip
      unwt = .false.
      i0 = islip
      i1 = i0 + 5*inext
      i1 = min(i1,maxepc)
      do 100 while ((.not. unwt) .and. (i0 .lt. i1))
           i0 = i0 + inext
           if(.not. lgood(kww(i0))) unwt = .true.
  100 continue

      if(unwt) then
           do 200 i=islip,i0-inext,inext

               ierr(i,ichan1,isite1) = igunwt
               kww(i) = igunwt
               kbb(i) = igunwt
               write(6,1000)ichan1,i
               lpatch = .false.
  200       continue

      else
            lpatch = .true.
      endif


 1000 format (4x,'------------- Unweighted data point ',2i6)
 2000 format (4x,'------------- Unflagged   data point ',2i6)

      return
      end
c************************************************************************
      subroutine cal_rms
     .         (ichan1,ichan2,isite1,isite2,
     .          yssd,kw,wv1,wv2,ii0,ii1)

c     retrieve a time series from the big data arrays
c
c     input:
c       ichan1,ichan2,isite1,isite2
c       imode = 0 for prefit, 1 for postfit
c
c     output
c       wv = observations
c       kw = error flags
c       nobs = number of observations
c       ngood = number of good observations
c       ii0,ii1 first and last good observations

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer ngood,ichan1,ichan2,isite1,isite2,i

      integer ii0,ii1
      logical lgood

      real*8 wv1(maxepc),wv2(maxepc)
      real*8 yssd,ylc,yavg
      integer kw(maxepc)

      if (ii1 .gt. ii0) then

c         sum all LC data
          yavg = 0.d0
          ngood = 0
          do 200 i = ii0,ii1
             if (lgood(kw(i))) then
               yavg = yavg+wv1(i)
     .            - faclr(jsat1)*(wv2(i)-gear(jsat1)*wv1(i))
               ngood = ngood+1
             endif
 200      continue

c         get average Lc
          if (ngood .gt. 1) then
             yavg = yavg/dble(ngood)
          endif

          yssd = 0.d0
          ngood = 0
          do 300 i = ii0,ii1
              if (lgood(kw(i))) then
                 ylc = wv1(i) - faclr(jsat1)*(wv2(i)-gear(jsat1)*wv1(i))
                 yssd = yssd + (ylc-yavg)**2
                 ngood = ngood + 1
              endif
 300     continue

         if (ngood .gt. 1) then
            yssd = dsqrt(yssd/dble(ngood-1))
         else
            yssd = 0.0d0
         endif
      else
         yssd = 0.0d0
      endif

      return

      end

c**********************************************************************************
