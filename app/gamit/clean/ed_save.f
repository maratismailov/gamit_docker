      subroutine ed_save (
     .                 nsave,nobs,imenu,mwarn,jsite,jchan,
     .                 replot,newbox,redraw,lfound,
     .                 isite1,isite2,jsite1,jsite2,
     .                 ichan1,ichan2,jchan1,jchan2,
     .                 isvmode)

c     save the working array into the large arrays

c     combined with ed_save1 by including isvmode in call, kurt June 91

c     isvmode  = 0 do error flags and buffer site and chan pointers
c     isvmode  = 1 do not do error flags, or buffer pointers

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

c     FUNCTIONS:
      real*8 rnhalf
      real*8 rnqrtr

      logical        redraw,replot,newbox,lfound

      integer*2       mwarn(4),iabs2
      integer         imenu,i

      integer         ichan1,ichan2,jchan1,jchan2
      integer         isite1,isite2,jsite1,jsite2
      integer         nobs,nsave,jsite,jchan,isvmode

      real*8 yadd1,yadd2


      if (jsite .le. 0 .or. jchan .le. 0) then
         call gmsg (mwarn,'ED_SAVE error!')
         print *,'ED_SAVE: zero! ',jsite,jchan
         return
      endif

      if (isvmode .eq. 0) then
         do 3000 i=1,nobs
c           User can only apply certain error flags.
c           Do not clobber other values set elsewhere.
c           Flags assigned by the conjunction of observations
c           in double or single differences will have flags equal
c           to 100 + these values, and thus not be propogated.
            if (kww(i) .eq. igunwt .or.
     .         kww(i) .eq. igrewt .or.
     .         kww(i) .eq. igbias .or.
     .         kww(i) .eq. igchop .or.
     .         kww(i) .eq. iggood) then
              ierr(i,jchan,jsite) = kww(i)
            endif
 3000     continue
       endif


       do 3200 i=1,nobs
c        Change yval for all points, not just good ones
c        Allow down to quarter cycle slips for lambda = +/- 2.
c        Allow up to half cycle slips for lambda = +/- 2.
         if (iabs2(lambds(jsite,jchan,1)) .eq. 1) then
            yadd1=rnhalf(wl1(i)-cl1(i))
         else if (iabs2(lambds(jsite,jchan,1)) .eq. 2) then
            yadd1=rnqrtr(wl1(i)-cl1(i))
         else
            yadd1 = 0.0d0
         endif

         if (iabs2(lambds(jsite,jchan,2)) .eq. 1) then
            yadd2=rnhalf(wl2(i)-cl2(i))
         else if (iabs2(lambds(jsite,jchan,2)) .eq. 2) then
            yadd2=rnqrtr(wl2(i)-cl2(i))
         else
            yadd2 = 0.0d0
         endif

c        flip sign if necessary
         if (isvmode .eq. 1) then
            if(jsite .eq. isite1 .and. jchan .eq. ichan2) then
               yadd1 = -yadd1
               yadd2 = -yadd2
            else if(jsite .eq. isite2 .and. jchan .eq. ichan1) then
               yadd1 = -yadd1
               yadd2 = -yadd2
            endif
         endif

         yl1(i,jchan,jsite) = yl1(i,jchan,jsite) + yadd1
         yl2(i,jchan,jsite) = yl2(i,jchan,jsite) + yadd2
 3200 continue

      call gmsg (mwarn,'SAVEd your changes.')
      replot = .false.
      newbox = .true.
      redraw = .true.
      imenu = -1

      if (isvmode .eq. 0) then
c        if the site and chan has changed due to a FIND, revert.
         if (lfound) then
            ichan1 = jchan1
            ichan2 = jchan2
            isite1 = jsite1
            isite2 = jsite2
            lfound = .false.
         else
            jchan1 = ichan1
            jchan2 = ichan2
            jsite1 = isite1
            jsite2 = isite2
         endif
      endif

      return
      end
