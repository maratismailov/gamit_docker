      subroutine ed_elim (
     .                 nelim,nbrack,nobs,imenu,mwarn,l12e,
     .                 iibl,iibr,ii0,ii1,iibeg,iiend,
     .                 replot,newbox,redraw,lfound,
     .                 isite1,isite2,jsite1,jsite2,
     .                 ichan1,ichan2,jchan1,jchan2)

c     **** ELIM ***
c     Eliminate all biases in all combinations in the span.
c     By its nature, this cannot be udone.


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      logical        redraw,replot,newbox,lfound,lbias

      integer*2       mwarn(4)
      integer         imenu,i

      integer         ichan1,ichan2,jchan1,jchan2
      integer         isite1,isite2,jsite1,jsite2
      integer         nobs,nelim,l12e,ii0,ii1
      integer         nbrack,iibl,iibr,iibeg,iiend



c     Update observations and data flags
      nelim = nelim + 1
      if (nelim .eq. 1) then

c         first, buffer the old values into b vectors
          if (l12e .eq. 1) then
                call vcopy (kww,wl1,kbb,bl1,nobs)
          else if (l12e .eq. 4) then
                call vcopy (kww,wl2,kbb,bl2,nobs)
          else
                call vcopy (kww,wl1,kbb,bl1,nobs)
                call vcopy (kww,wl2,kbb,bl2,nobs)
          endif

           if (nbrack .eq. 2) then
                iibeg = iibl
                iiend = iibr
           else if (nbrack .eq. 1) then
                iibeg = iibl
                iiend = iibl
           else
                iibeg = ii0
                iiend = ii1
           endif

           call gmsg (mwarn,'Confirm ELIM or UNDO')
           do 3201 i=iibeg,iiend
             if (lbias(kww(i))) then
                 kww(i) = iggood
              endif
 3201     continue
          redraw = .true.
          replot = .true.
      else if (nelim .ge. 2) then
          do 3202 i=iibeg,iiend
c
             if (lbias(kbb(i))) then

                if(ichan1.ne.0 .and. isite1.ne.0 .and.
     .             lbias(ierr(i,ichan1,isite1)))
     .                   ierr(i,ichan1,isite1) = iggood
                if(ichan2.ne.0 .and. isite1.ne.0 .and.
     .             lbias(ierr(i,ichan2,isite1)))
     .                   ierr(i,ichan2,isite1) = iggood
                if(ichan1.ne.0 .and. isite2.ne.0 .and.
     .             lbias(ierr(i,ichan1,isite2)))
     .                   ierr(i,ichan1,isite2) = iggood
                if(ichan2.ne.0 .and. isite2.ne.0 .and.
     .             lbias(ierr(i,ichan2,isite2)))
     .                   ierr(i,ichan2,isite2) = iggood
             endif

 3202      continue

           call gmsg (mwarn,'SAVEd your changes.')
           replot = .false.
           newbox = .true.
           redraw = .true.
           imenu = -1
           nelim = 0
c          if the site and chan has changed due to a FIND, revert.
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



