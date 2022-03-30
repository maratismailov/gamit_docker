c************************************************************************
      subroutine ed_line (
     .             mplot,mrect,mtaxi,mwarn,csize,
     .             ix1,ix2,iy1,iy2,ix0,nwins,
     .             iibl,iibr,i0,ii0,ii1,dt,
     .             xs,ttt,amsg,igerr,slip)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         iibl,iibr,i0,ii0,ii1,nwins
      integer*2       ix1,iy1,ix2,iy2,ix0
      real*8          dt,xs,ttt
      integer*2       mplot(maxwin,4),mrect(4),mtaxi(4)
      integer*2       mwarn(4),csize(2)
      integer*4       igerr
      character*96    amsg
      logical         slip


      if (slip) then
         i0 = iibl
         ttt = tx(i0)-tx(ii0)
         ix1 = ix0 + int(xs*ttt)
         ix2 = ix0 + int(xs*ttt)
c        go over all windows
         iy1 = mplot(1,2)
         iy2 = mplot(nwins,2) + mplot(nwins,4)
         call gline (ix1,iy1,ix2,iy2,igerr)

c        label the base of the bracket
         write (amsg,'(i4)') i0

         mrect(1) = ix1 - int(2.1*csize(1))
         mrect(2) = mtaxi(2)
         mrect(3) = int(4.2*csize(1))
         mrect(4) = int(1.1*csize(2))
         call gmsg (mrect,amsg)
      endif

      return
      end
c************************************************************************

