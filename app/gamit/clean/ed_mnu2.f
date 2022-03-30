c************************************************************************
      subroutine ed_abort (redraw,newbox,replot,lfound,
     .            imenu,mwarn,isite1,isite2,
     .          ichan1,ichan2,jsite1,jsite2,jchan1,jchan2)

c     abort the edits on current set of modifications.

      logical    redraw,newbox,replot,lfound
      integer  imenu,isite1,isite2,ichan1,ichan2
      integer  jsite1,jsite2,jchan1,jchan2
      integer*2         mwarn(4)

      call gmsg (mwarn,'Aborting without saving.')
      imenu = -1
      redraw = .true.
      newbox = .true.
      replot = .false.
c     if the chan or site indexes have changed due to a FIND, revert.
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

      return
      end
c************************************************************************
      subroutine ed_span (nbrack,lspan,allpts,redraw,replot,
     .                     nl,ii0,ii1,iibl,nobs,iibr)
      logical     lspan,allpts,redraw,replot
      integer     nbrack,nl,ii0,ii1,iibl,nobs,iibr


      if (nbrack .eq. 1) then
c        Center span on point
         nl = max(ii1-ii0,10)
         ii0 = max(iibl-nl/2,1)
         ii1 = min(iibl+nl/2,nobs)
         lspan = .false.
         allpts = .false.
      else if (nbrack .eq. 2) then
c        Show range between brackets
         ii0 = iibl
         ii1 = iibr
         lspan = .false.
         allpts = .false.
      else
c        Toggle SPAN
         lspan = (.not. lspan)
         if (lspan) allpts = .false.
      endif
      redraw = .true.
      replot = .true.

      return
      end
c************************************************************************
      subroutine ed_all (allpts,lspan,replot,redraw)
      logical      allpts,lspan,redraw,replot

      allpts = (.not. allpts)
      if (allpts) lspan = .false.
      replot = .true.
      redraw = .true.

      return
      end
c************************************************************************
      subroutine ed_marg (domarg,redraw,replot)
         logical      domarg,redraw,replot

c              toggle plotting of marginal points
               domarg = (.not. domarg)
               redraw = .true.
               replot = .true.


      return
      end
c************************************************************************
      subroutine ed_hide (hide_b,redraw,replot)
         logical      hide_b,redraw,replot

c              toggle plotting of marginal points
               hide_b = (.not. hide_b)
               redraw = .true.
               replot = .true.


      return
      end
c************************************************************************
      subroutine ed_lt (nl,ii1,ii0,replot)
      logical      replot
      integer      nl,ii0,ii1

      nl = max(ii1 - ii0,10)
      ii1 = ii0 + nint(2.0 * nl)
      replot = .true.

      return
      end
c************************************************************************
      subroutine ed_gt (nl,ii1,ii0,replot)
      logical      replot
      integer      nl,ii0,ii1

      nl = max(ii1 - ii0, 10)
      ii1 = ii0 + nint(0.5 * nl)
      replot = .true.

      return
      end
c************************************************************************
      subroutine ed_mvleft (key,MOUSE2,nl,inter,ii0,ii1,replot)
      logical      replot
      character*1  MOUSE2,key
      integer      nl,inter,ii0,ii1

c     move back 4 minutes.
      if (key .eq. MOUSE2) then
         nl = 240/inter
      else
         nl = nint(0.8 * max(ii1 - ii0, 10))
      endif
      ii0 = ii0 - nl
      ii1 = ii1 - nl
      replot = .true.

      return
      end
c************************************************************************
      subroutine ed_mvright (key,MOUSE2,nl,inter,ii0,ii1,replot)
      logical      replot
      character*1  MOUSE2,key
      integer      nl,inter,ii0,ii1

c     move back 4 minutes.
      if (key .eq. MOUSE2) then
         nl = 240/inter
      else
         nl = nint(0.8 * max(ii1 - ii0, 10))
      endif
      ii0 = ii0 + nint(0.8 * nl)
      ii1 = ii1 + nint(0.8 * nl)
      replot = .true.

      return
      end
c************************************************************************

