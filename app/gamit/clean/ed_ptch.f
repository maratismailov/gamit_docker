c************************************************************************
      subroutine ed_ptch (
     .        nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .        ngood,yend,iiend,ybeg,iibeg,
     .        yadd,yslp,itemp,jsite,jchan,idd,
     .        key,mwarn,mhelp,amsg,replot,domarg,nmlfit)

c     **** PATCH *** invoked directly

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         nbrack,nobs,iibl,iibr,l12e,ii0,ii1
      character*1     key
      integer         jsite,jchan,idd,nmlfit,ngood
      integer         itemp,iibeg,iiend
      logical         replot,domarg,doit
      real*8          ybeg,yend,yslp,yadd
      integer*2       mwarn(4),mhelp(4)
      character*96    amsg
      real*8          yadd2(2),confid

C     There are 2 flavors of patch, depending on observable:
c     L1,L2:  simple linear interpolation
c     LC,LG:  Use Herring's PATCH routine to fix L1 and L2 using LC and LG

c     Assume user wants PATCH to apply to interval from iibl to iibr inclusive
c     and that this has been determined from the 2 brackets

      if (nbrack .eq. 1 .or. nbrack .eq. 2) then
c      will occur on all points to the right of the bracket
            if (nbrack .eq. 1) iibr = nobs

c           first, buffer the old values into b vectors
            call vcopy (kww,wl1,kbb,bl1,nobs)
            call vcopy (kww,wl2,kbb,bl2,nobs)

            if (l12e .eq. 1) then
                 call ptch_L1(
     .                nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .                ngood,yend,iiend,ybeg,iibeg,
     .                yadd,yslp,itemp,jsite,jchan,key,
     .                mwarn,mhelp,amsg,replot,domarg,nmlfit)
                  replot = .true.

            else if (l12e .eq. 4) then
                 call ptch_L2(
     .                nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .                ngood,yend,iiend,ybeg,iibeg,
     .                yadd,yslp,itemp,jsite,jchan,key,
     .                mwarn,mhelp,amsg,replot,domarg,nmlfit)
                  replot = .true.

            else if (l12e .eq. 3) then
c                 LC window: use LC and LG to clean up L1 and L2

                  nmlfit = 5
                  doit = .true.
                  if (jsite .ne. 0 .and. jchan .ne. 0) then
                     call PATCH (iibl,iibr,ii0,ii1,mhelp,jsite,jchan,
     .                    idd,nobs,replot,nmlfit,yadd2,confid,doit)
                     write (amsg,100) jsite,jchan,yadd2, confid
 100                 format(2(i2,1x),2(f5.1,1x),f4.0,'% confid')
                  else
                     write (amsg,'(''ED_PTCH: '',2(i2,1x))') jchan,jsite
                  endif
                  call gmsg(mhelp, amsg)
                  replot = .true.
            else if (l12e .eq. 2) then
c                 LG window: use LG and WL to clean up L1 and L2

                  nmlfit = 5
                  doit = .true.
                  if (jsite .ne. 0 .and. jchan .ne. 0) then
                     call PATCH_LGWL (iibl,iibr,ii0,ii1,mhelp,jsite,
     .                    jchan,nobs,replot,nmlfit,yadd2,confid,doit)
                  else
                     write (amsg,'(''ED_PTCH: '',2(i2,1x))') jchan,jsite
                     call gmsg(mhelp, amsg)
                  endif
                  replot = .true.
            else if ((l12e .eq. 7 .or. l12e .eq. 8) .and.
     .                        (dabs(pc2(iibl)) .gt. 0.1)) then
c                 If user clicks in WL or W* window, then
c                 use WL and LG to clean up L1 and L2

                   call ptch_WL(
     .                nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .                yadd2,jsite,jchan,mwarn,mhelp,amsg)
                  replot = .true.
            else
                  call gmsg(mwarn,'May not PATCH this observable.')
            endif

      else
            call gmsg (mwarn,'Need at least 1 bracket for PATCH.')
            replot = .true.
      endif
      nbrack = 0

      return
      end
c************************************************************************
      subroutine patch_save (
     .              isite1,isite2,ichan1,ichan2,jsite,jchan,
     .              nbrack,nobs,iibr,iibl,l12e,ii0,ii1,ngood,
     .              yend,iiend,ybeg,iibeg,ytot,yadd,yslp,itemp,
     .              key,mwarn,mhelp,amsg,replot,domarg,nmlfit,nsave,
     .              imenu,newbox,redraw,lfound,imode,
     .              jsite1,jsite2,jchan1,jchan2)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         iibl,iibr,nobs,nbrack,l12e,ii0,ii1,ngood
      integer         iiend,iibeg,itemp,nmlfit,nsave,imenu,imode
      integer         isite1,isite2,ichan1,ichan2,jsite,jchan
      integer         jsite1,jsite2,jchan1,jchan2
      logical         lfound,redraw,newbox,replot,domarg
      real*8          ytot,yadd,yend,ybeg,yslp
      integer*2       mwarn(4),mhelp(4)
      character*96    amsg
      character*1     key
      integer         isvmode

      call ed_ptch1 (
     .       nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .       ngood,yend,iiend,ybeg,iibeg,
     .       ytot,yadd,yslp,itemp,jsite,jchan,
     .       key,mwarn,mhelp,amsg,replot,domarg,nmlfit)

      isvmode = 1
      call ed_save(
     .        nsave,nobs,imenu,mwarn,jsite,jchan,
     .        replot,newbox,redraw,lfound,
     .        isite1,isite2,jsite1,jsite2,
     .        ichan1,ichan2,jchan1,jchan2,isvmode)

      write (amsg,'(''Patched chan '',i2," site ",i2)') jchan,jsite

      call gmsg (mwarn,amsg)

      return
      end
c
c************************************************************************
      subroutine ed_ptch1(
     .        nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .        ngood,yend,iiend,ybeg,iibeg,
     .        ytot,yadd,yslp,itemp,jsite,jchan,
     .        key,mwarn,mhelp,amsg,replot,domarg,nmlfit)

c     PATCH as invoked by SLIP function

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         nbrack,nobs,iibl,iibr,l12e,ii0,ii1
      character*1     key
      integer         jsite,jchan,idd,nmlfit,ngood
      integer         itemp,iibeg,iiend
      logical         replot,domarg,doit
      real*8          ybeg,yend,yslp,yadd,ytot
      integer*2       mwarn(4),mhelp(4)
      character*96    amsg
      real*8          yadd2(2),confid

c     Assume user wants PATCH to apply to interval from iibl to iibr inclusive
c     and that this has been determined from the 2 brackets

      iibr = nobs

      nmlfit = 5
      doit = .true.
      call PATCH (iibl,iibr,ii0,ii1,mwarn,jsite,jchan,
     .            idd,nobs,replot,nmlfit,yadd2,confid,doit)

      ytot  = ytot + dabs(wl1(iibl) - bl1(iibl))
     .             + dabs(wl2(iibl) - bl2(iibl))
      replot = .true.

      return
      end
c
c************************************************************************
      subroutine ptch_L1 (
     .        nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .        ngood,yend,iiend,ybeg,iibeg,
     .        yadd,yslp,itemp,jsite,jchan,key,
     .        mwarn,mhelp,amsg,replot,domarg,nmlfit)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         nbrack,nobs,iibl,iibr,l12e,ii0,ii1
      character*1     key
      integer         jsite,jchan,nmlfit,ngood
      integer         i,itemp,iibeg,iiend
      logical         replot,domarg,lgood
      real*8          ybeg,yend,yslp,yadd
      integer*2       mwarn(4),mhelp(4),iabs2
      character*96    amsg
      real*8          rnhalf


c     if user clicks in window for L1 or L2, PATCH is
c     a dumb linear interpolation.

c     look to the left for the last good points
      ngood = 0
      do 3400 i = (iibl-1), 1, -1
           if (lgood(kww(i))) then
              if (ngood .eq. 0) then
                  yend = wl1(i)
                  iiend = i
                  ngood = 1
               else if (ngood .eq. 1) then
                   ybeg = wl1(i)
                   iibeg = i
                   ngood = 2
                endif
            endif
 3400  continue

c      determine move amount
       if (ngood .eq. 0) then
             yadd = 0.0d0
       else if (ngood .eq. 1) then
             yadd = yend - wl1(iibl)
       else if (ngood .eq. 2) then
             yslp = (yend - ybeg)/(iiend - iibeg)
             yadd = yend + (iibl-iiend)*yslp - wl1(iibl)
       endif

c      For lamba = +/- 2 cycle slip may be be a half-integer.
c      Otherwise, it should be a whole integer.

       if (iabs2(lambds(jsite,jchan,1)) .eq. 2) then
c           half-cycle -slip
            yadd = rnhalf(yadd)
       else
c            cycle slip must be an integer
             yadd = dnint(yadd)
       endif

c      add yadd to all points
       do i=iibl,iibr
             wl1(i) = wl1(i) + yadd
       enddo
       write (amsg,'(''Patched L1: '',2(i2,1x),f8.1)') jsite,jchan,yadd
       call gmsg(mhelp, amsg)


      return
      end
c************************************************************************
      subroutine ptch_L2 (
     .        nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .        ngood,yend,iiend,ybeg,iibeg,
     .        yadd,yslp,itemp,jsite,jchan,key,
     .        mwarn,mhelp,amsg,replot,domarg,nmlfit)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         nbrack,nobs,iibl,iibr,l12e,ii0,ii1
      character*1     key
      integer         jsite,jchan,nmlfit,ngood
      integer         i,itemp,iibeg,iiend
      logical         replot,domarg,lgood
      real*8          ybeg,yend,yslp,yadd
      integer*2       mwarn(4),mhelp(4),iabs2
      character*96    amsg
      real*8          rnhalf


c     if user clicks in window for  L2, PATCH is
c     a dumb linear interpolation.

c     look to the left for the last good points
      ngood = 0
      do 3400 i = (iibl-1), 1, -1
           if (lgood(kww(i))) then
              if (ngood .eq. 0) then
                  yend = wl2(i)
                  iiend = i
                  ngood = 1
               else if (ngood .eq. 1) then
                   ybeg = wl2(i)
                   iibeg = i
                   ngood = 2
                endif
            endif
 3400  continue

c      determine move amount
       if (ngood .eq. 0) then
             yadd = 0.0d0
       else if (ngood .eq. 1) then
             yadd = yend - wl2(iibl)
       else if (ngood .eq. 2) then
             yslp = (yend - ybeg)/(iiend - iibeg)
             yadd = yend + (iibl-iiend)*yslp - wl2(iibl)
       endif

c      For lamba = +/- 2 cycle slip may be be a half-integer.
c      Otherwise, it should be a whole integer.

       if (iabs2(lambds(jsite,jchan,2)) .eq. 2) then
c           half-cycle -slip
            yadd = rnhalf(yadd)
       else
c            cycle slip must be an integer
             yadd = dnint(yadd)
       endif

c      add yadd to all points
       do i=iibl,iibr
             wl2(i) = wl2(i) + yadd
       enddo
       write (amsg,'(''Patched L2: '',2(i2,1x),f8.1)') jsite,jchan,yadd
       call gmsg(mhelp, amsg)


      return
      end
c************************************************************************
      subroutine check_patch (
     .                        nobs,ii0,ii1,domarg,plt,yflag0,yflag1,
     .                        ngood,lcando,dosave)



      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'


      integer      nobs,ngood,l12p,ii0,ii1,nmarg,istep,ijump
      real*8       ymin,ymax,yavg,yssd
      real*8       yflag1,yflag0,yjump,ystep
      logical      domarg,dosave
      real*8       plt(maxepc)
      integer*4    kpl(maxepc)
      logical      lcando(-3:ntypes)



      l12p = 3
      call combo (nobs,ii0,ii1,l12p,domarg,kpl,plt,
     .          ymin,ymax,yavg,yssd,yjump,ystep,
     .          ngood,nmarg,ijump,istep,lcando)

c     **** check if new LC-RMS is smaller then the old one
      yflag0 = yflag1
      yflag1 = yssd

      if (yflag1 .lt. yflag0) then
             dosave = .true.
      else
             dosave = .true.
      endif

      return
      end
c************************************************************************
      subroutine ptch_WL (
     .        nbrack,nobs,iibr,iibl,l12e,ii0,ii1,
     .        yadd2,jsite,jchan,mwarn,mhelp,amsg)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         nbrack,nobs,iibl,iibr,l12e,ii0,ii1,i
      integer         jsite,jchan
      logical         slipwl
      real*8          yadd2(2),dwl,dy1,dy2
      integer*2       mwarn(4),mhelp(4),iabs2
      character*96    amsg
      real*8          rnhalf


c     if user clicks in window for WL, PATCH according
c     to WL and dLG = 0;
c            dWL = L2 - L1
c            dLG = L2 - gL1 = 0


        call check_wl (ii0,ii1,iibl,slipwl,dwl)
        dy1 = dwl/(gear(jsat1)-1.0d0)
        dy2 = dwl+dy1

c      For lamba = +/- 2 cycle slip may be be a half-integer.
c      Otherwise, it should be a whole integer.

       if (iabs2(lambds(jsite,jchan,2)) .eq. 2) then
c           half-cycle -slip
            dy1 = rnhalf(dy1)
            dy2 = rnhalf(dy2)
       else
c            cycle slip must be an integer
             dy1 = dnint(dy1)
             dy2 = dnint(dy2)
       endif

       yadd2(1) = dy1
       yadd2(2) = dy2

c      add yadd to all points
       do i=iibl,iibr
             wl1(i) = wl1(i) + dy1
             wl2(i) = wl2(i) + dy2
       enddo

       write (amsg,100) jsite,jchan,dy1,dy2
 100   format('WL: ',2(i2,1x),2(f8.1,1x),f4.0)
       call gmsg(mhelp, amsg)

      return
      end
c************************************************************************

