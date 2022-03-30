c************************************************************************
      subroutine ed_bias (
     .        key,nbrack,nobs,iibr,iibl,ii0,ii1,
     .        isite1,isite2,jsite,ichan1,ichan2,jchan,
     .        ntemp,itemp,mwarn,amsg,replot,domarg)

c     **** BIAS ***
c     Handle biases

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         nbrack,nobs,iibl,iibr,i,ii0,ii1
      character*1     key
      integer         isite1,isite2,ichan1,ichan2
      integer         jsite,jchan,ntemp,itemp
      logical         replot,domarg,lbias,lmarg,lgood
      integer*2       mwarn(4)
      character*96    amsg

      if (nbrack .ge. 0 .and. nbrack .le. 2) then
c        first, buffer the old values into b vectors
         call vcopy (kww,wl1,kbb,bl1,nobs)
         call vcopy (kww,wl2,kbb,bl2,nobs)

         if (key .eq. MOUSE1) then
c           remove biases
            if (nbrack .eq. 0) then
c              no brackets: remove all biases
               iibl = ii0
               iibr = ii1
            else if (nbrack .eq. 1) then
c              1 bracket: remove bias on bracketed point
               iibr = iibl
            endif

            ntemp = 0
            do 3065 i = iibl,iibr
c              Remove all the bias flags in span.
c              But only if point is good.
               if (lbias(ierr(i,ichan1,isite1)) .and.
     .            lgood(kww(i))) then
                  kww(i)=iggood
                  ntemp = ntemp + 1
               else
                  if (nbrack .eq. 1) then
                     call gmsg (mwarn,'No bias here!')
                  endif
               endif
 3065       continue
            if (ntemp .gt. 0) then
               write (amsg,3067) ntemp
 3067          format ('Removed ',i5,' bias flags.')
               call gmsg (mwarn,amsg)
            endif
         else if (key .eq. MOUSE2 .and. nbrack .eq. 1) then
c           find the bias
            ntemp = 0
            if (lbias(ierr(iibl,ichan1,isite1))) then
               jchan = ichan1
               jsite = isite1
               ntemp = ntemp + 1
            endif
            if (isite2 .gt. 0) then
               if (lbias(ierr(iibl,ichan1,isite2))) then
                  jchan = ichan1
                  jsite = isite2
                  ntemp = ntemp + 1
               endif
            endif
            if (ichan2 .gt. 0) then
               if (lbias(ierr(iibl,ichan2,isite1))) then
                  jchan = ichan2
                  jsite = isite1
                  ntemp = ntemp + 1
               endif
            endif
            if (ichan2 .gt. 0 .and. isite2 .gt. 0) then
               if (lbias(ierr(iibl,ichan2,isite2))) then
                  jchan = ichan2
                  jsite = isite2
                  ntemp = ntemp + 1
               endif
            endif
            if (ntemp .eq. 0) then
               call gmsg (mwarn,'No bias flag here.')
            else if (ntemp .eq. 1) then
               write (amsg,3066) jchan,jsite
 3066          format ('Bias in chan ',i2,' site ',i2)
               call gmsg (mwarn,amsg)
            else
               call gmsg (mwarn,'No bias flag here.')
            endif
         else if (key .eq. MOUSE3 .and. nbrack .ge. 1) then
c           Add bias flags to whole span
            if (nbrack .eq. 1) iibr = iibl
            do 3068 i = iibl,iibr
c              Add a bias flag.
               itemp = ierr(i,ichan1,isite1)
               if (lgood(itemp) .or.
     .            (domarg.and.lmarg(itemp))) then
                  kww(i)=igbias
               else
                  if (nbrack .eq. 1) then
                     call gmsg (mwarn,'Cannot add bias here!')
                  endif
               endif
3068        continue
         endif
      else
         call gmsg (mwarn,'Invalid bracket for BIAS.')
      endif
      replot = .true.
      nbrack = 0

      return
      end
