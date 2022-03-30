      subroutine ed_data (
     .             nbrack,nobs,igraph,mplot,mrect,mtaxi,mwarn,
     .             csize,lbx,rbx,ix1,ix2,iy1,iy2,ix0,iwin,nwins,
     .             pos,iibl,iibr,i0,ii0,ii1,l12p,l12e,l12v,dt,
     .             xs,ttt,yold,plt,amsg,igerr,replot,domarg)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         nobs,nbrack,igraph,j
      integer         iibl,iibr,i0,ii0,ii1
      integer         l12p,l12e,l12v(maxwin),iwin,nwins
      integer*2       lbx,rbx,ix1,iy1,ix2,iy2,ix0,pos(2)
      real*8          dt,xs,ttt,yold,plt(maxepc)
      logical         replot,lgood,ltemp,lmarg,domarg
      integer*2       mplot(maxwin,4),mrect(4),mtaxi(4)
      integer*2       mwarn(4),csize(2)
      integer*4       igerr
      character*96    amsg


c     figure out which plotting window
      do j=1,nwins
         if  (pos(1) .ge. mplot(j,1) .and.
     .        pos(1) .le. mplot(j,1)+mplot(j,3) .and.
     .        pos(2) .ge. mplot(j,2) .and.
     .        pos(2) .le. mplot(j,2)+mplot(j,4)) then
            l12p = l12v(j)
            iwin = j
         endif
      enddo

      nbrack = nbrack + 1
      l12e = l12p

c     squeeze the brackets inward
      if (nbrack.gt.0 .and. mod(nbrack,2) .eq. 1) then
c        first bracket click, assumed on left
         lbx = pos(1)
         if (igraph .eq. 1 .or. igraph .eq. 2) then
            iibl = nint((lbx-ix0)/(dt*xs)) + ii0
            iibl = min(nobs,iibl)
         else
            ttt = (tx(maxepc)-tx(1))/(maxepc-1)
            iibl = nint((lbx-ix0)/(ttt*xs)) + 1
            iibl = min(maxepc,iibl)
         endif

         if (igraph .eq. 1) then
c           look to the right for a good point
 3405       continue
            if (iibl .lt. ii1 .and. iibl .gt. 0) then
               ltemp=(lgood(kww(iibl)) .or.
     .              (domarg.and.lmarg(kww(iibl))))
               if (.not. ltemp) then
                  iibl = iibl + 1
                  goto 3405
               endif
            endif
         endif
         iibl = max(1,iibl)
         i0 = iibl
      else if (nbrack.gt.0 .and. mod(nbrack,2) .eq. 0) then
c        second bracket click; assumed on right
         rbx = pos(1)
         if (igraph .eq. 1 .or. igraph .eq. 2) then
            iibr = nint((rbx-ix0)/(dt*xs)) + ii0
            iibr = min(nobs,iibr)
         else
            ttt = (tx(maxepc)-tx(1))/(maxepc-1)
            iibr = nint((rbx-ix0)/(ttt*xs)) + 1
            iibr = min(maxepc,iibr)
         endif
         if (igraph .eq. 1) then
c           look to the left for a good point
 3410       continue
            if (iibr .gt. 1) then
               ltemp=(lgood(kww(iibr)) .or.
     .              (domarg.and.lmarg(kww(iibr))))
               if (.not. ltemp) then
                  iibr = iibr - 1
                  goto 3410
               endif
            endif
         endif
         iibr = max(1,iibr)
         i0 = iibr
      endif

      if (nbrack .gt. 2) then
         if (igraph .eq. 1 .or. igraph .eq. 2) then
            call gmsg (mwarn,'Invalid bracket.')
            nbrack = 0
            replot = .true.
         endif
      endif

      if (nbrack .eq. 2) then
c        make sure smaller index is on left
         if (rbx.lt.lbx .and. igraph .ne. 3) then
            nbrack = 0
            replot = .true.
            call gmsg(mwarn,'Bracket left to right.')
         endif
      endif

      if (nbrack .ge. 1) then
c        draw the bracket
         if (igraph .eq. 3) then
            ttt = tx(i0)-tx(1)
         else
            ttt = tx(i0)-tx(ii0)
         endif
         ix1 = ix0 + int(xs*ttt)
         ix2 = ix0 + int(xs*ttt)
c        go over all windows
         iy1 = mplot(1,2)
         iy2 = mplot(nwins,2) + mplot(nwins,4)
         call gline (ix1,iy1,ix2,iy2,igerr)

c        label the base of the bracket
         if (igraph .eq. 1) then
c           for time series, use the epoch number
            write (amsg,'(i4)') i0
         else if (igraph .eq. 3) then
c           for power spectrum, use the frequency
            write (amsg,'(1pe8.1,''Hz'')') tx(i0)
         else
c           don't plot anything
            write (amsg,'(1x)')
         endif
         mrect(1) = ix1 - int(2.1*csize(1))
         mrect(2) = mtaxi(2)
         mrect(3) = int(4.2*csize(1))
         mrect(4) = int(1.1*csize(2))
         call gmsg (mrect,amsg)

c        print value.  This only works if there is only one window.
c        Otherwise, the value in plt is not nescessarily current.
c        It would be nice if this worked better in the future.
         if (igraph .eq. 1 .and. nwins .eq. 1) then
            if (mod(nbrack,2) .eq. 1) then
              write (amsg,'(''Y = '',1pg20.2)') plt(i0)
            else
              write (amsg,'(''Ydiff = '',1pg20.2)') plt(i0)-yold
              yold = plt(i0)
            endif
         else if (igraph .eq. 3) then
            write(amsg,'(''Period = '',f10.2,''s'')')
     .      1.d0/tx(i0)
         else
            write (amsg,'(1x)')
         endif
         call gmsg (mwarn,amsg)
      endif

      return
      end

