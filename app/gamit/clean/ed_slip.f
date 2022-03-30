c************************************************************************
c           **** SLIP ***
      subroutine ed_slip (
     .               replot,nextslip,lslip,npatch,nslip,
     .                                  simple,key,key0)



      logical       nextslip,lslip,replot,simple
      integer       npatch,nslip
      character*1   key,key0

      nextslip = .true.
      lslip  = .true.
      simple = .true.
      npatch = 0
      nslip  = 0
      key0   = key

      return
      end
c************************************************************************
      subroutine autoslip (
     .         lfound,seekob,newser,replot,redraw,newbox,domarg,lslip,
     .         nslip,nbrack,nobs,ntemp,npatch,nwins,iiend,iibeg,itemp,
     .         jsite1,jsite2,jchan1,jchan2,isite1,isite2,ichan1,ichan2,
     .         jsite,jchan,key,pos,l12e,l12v,l12p,ngood,nsave,nmlfit,
     .         ymin2,ymax2,yadd,yend,ybeg,yslp,ymult,yold,slip,nextslip,
     .         ii0,ii1,iibl,iibr,iibl0,epsilon,istart,ilast,imenu,
     .         mplot,mrect,mtaxi,menus,mclip,csize,mwarn,mhelp,
     .         amsg,igerr,
     .         amenus,lbx,rbx,ix1,ix2,iy1,iy2,ix0,i0,dt,xs,ttt,ytot,
     .         iprob,nprob,pbias,simple,nsat,ncfls,imode)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      logical       lfound,replot,slip,nextslip,lslip
      logical       seekob,newser,redraw,newbox,domarg,simple,mbias
      integer       npatch,nslip,nwins,istart,ilast,iibl0
      integer       jsite,jchan,jchan1,jchan2,itemp,nmlfit
      integer       nbrack,nobs,ntemp,iiend,iibeg,ngood,nsave
      integer       isite1,isite2,ichan1,ichan2,jsite1,jsite2,imenu
      integer       iibl,iibr,i0,ii0,ii1,l12v(maxwin),l12e,l12p
      integer       iprob,nprob,pbias(maxepc),nsat,ncfls,imode
      integer*2     lbx,rbx,ix1,iy1,ix2,iy2,ix0,pos(2)
      real*8        ytot,yadd,yend,ybeg,yslp,ymult,yold,dt,xs
      real*8        ymin2(maxwin),ymax2(maxwin),epsilon(4),ynoise(4),ttt
      character*1   key
      integer*2     mplot(maxwin,4),mrect(4),mtaxi(4),menus(nmenus,4)
      integer*2     mwarn(4),csize(2),mclip(4),mhelp(4)
      integer*4     igerr
      character*16  amenus(nmenus)
      character*96  amsg

c     maximum number of automatic patching attempts and slip
      integer   maxpatch,maxslip
      parameter (maxslip  = 15)
      parameter (maxpatch = 2)

c     **************************************************************
c     *** Cycle-slip finder
c     **************************************************************
      nslip = nslip + 1
      if (nslip .eq. 1)
     .     call initial(
     .         lslip,npatch,nslip,nbrack,mwarn,amsg,iibl0,ytot,
     .         nobs,iibeg,iiend,istart,ilast,ii0,ii1,iibl,iibr,
     .         iprob,nprob)

      if (lslip) then
           call findslip_ed(
     .          nbrack,nobs,npatch,nslip,iibl0,ynoise,
     .          istart,ilast,iibr,iibl,l12e,ii0,ii1,
     .          slip,ngood,yend,iiend,ybeg,iibeg,epsilon,
     .          ymin2,ymax2,yadd,yslp,itemp,isite1,ichan1,
     .          key,mwarn,amsg,replot,domarg,nmlfit)  

           if (slip) then
                call ed_line (
     .               mplot,mrect,mtaxi,mwarn,csize,
     .               ix1,ix2,iy1,iy2,ix0,nwins,
     .               iibl,iibr,i0,ii0,ii1,dt,
     .               xs,ttt,amsg,igerr,slip)
                if (key .eq. MOUSE2 .or. key .eq. MOUSE3) then
                     call ed_find1(
     .                   nbrack,ntemp,iibl,itemp,jsite,jchan,
     .                   redraw,replot,newser,lfound,mwarn,
     .                   npatch,amsg,isite1,isite2,ichan1,ichan2,
     .                   nprob,pbias,mbias)

                     if (lfound) then
                        call patch_save (
     .                     isite1,isite2,ichan1,ichan2,jsite,jchan,
     .                     nbrack,nobs,iibr,iibl,l12e,ii0,ii1,ngood,
     .                     yend,iiend,ybeg,iibeg,ytot,yadd,yslp,itemp,
     .                     key,mwarn,mhelp,amsg,replot,domarg,nmlfit,
     .                     nsave,imenu,newbox,redraw,lfound,imode,
     .                     jsite1,jsite2,jchan1,jchan2)  
                     else if (mbias) then
                         call gmsg (mwarn,'More than 1 bias')
                         call multibias (
     .                       lfound,seekob,newser,replot,redraw,
     .                       newbox,domarg,lslip,nslip,nbrack,nobs,
     .                       ntemp,npatch,nwins,iiend,iibeg,itemp,
     .                       jsite1,jsite2,jchan1,jchan2,
     .                       isite1,isite2,ichan1,ichan2,
     .                       jsite,jchan,key,pos,l12e,l12v,l12p,
     .                       ngood,nsave,nmlfit,ymin2,ymax2,yadd,
     .                       yend,ybeg,yslp,ymult,yold,slip,nextslip,
     .                       ii0,ii1,iibl,iibr,iibl0,epsilon,istart,
     .                       ilast,imenu,mplot,mrect,mtaxi,menus,mclip,
     .                       csize,mwarn,amsg,igerr,amenus,lbx,rbx,
     .                       ix1,ix2,iy1,iy2,ix0,i0,dt,xs,ttt,ytot,
     .                       iprob,nprob,pbias,simple,nsat,ncfls,imode)
                         call patch_save (
     .                      isite1,isite2,ichan1,ichan2,jsite,jchan,
     .                      nbrack,nobs,iibr,iibl,l12e,ii0,ii1,ngood,
     .                      yend,iiend,ybeg,iibeg,ytot,yadd,yslp,itemp,
     .                      key,mwarn,mhelp,amsg,replot,domarg,nmlfit,
     .                      nsave,
     .                      imenu,newbox,redraw,lfound,imode,
     .                      jsite1,jsite2,jchan1,jchan2)
                     else
                            call gmsg (mwarn,'No bias/gap')
                     endif
              write (amsg,'("Patch:  site",i4,"  chan",i4)')jsite,jchan
                     call gmsg (mwarn,amsg)
                     call ed_plot(newser,nobs)
                endif
           else
                write (amsg,'(''ytot #  '',1pe9.2," epsilon ",1pe9.2)')
     .                      ytot,epsilon(2)
                call gmsg (mwarn,amsg)
                if (ytot .gt. 0.01d0) then
                     npatch = npatch + 1
                     nslip = 0
                else
                     npatch = maxpatch
                endif
           endif
      else
          npatch = maxpatch
      endif

      nbrack = 0
      lslip  = .false.
      nextslip = .false.

      if (nslip .eq. maxslip) then
                write (amsg,'(''nslip = maxslip  '')')
                call gmsg (mwarn,amsg)
           if (ytot .gt. 0.01d0) then
                npatch = npatch + 1
                nslip = 0
                iibl = iibeg
            else
                npatch = maxpatch
            endif
       endif
       if (key .eq. MOUSE2 .or. key .eq. MOUSE3) then
            if (npatch .lt. maxpatch) nextslip = .true.
       endif
       if (key .eq. MOUSE3 .and. npatch .ge. maxpatch) then
             call ed_seek (newser,seekob,nobs)
             nextslip = .true.
             npatch = 0
             nslip  = 0
        endif


      return
      end
c

c************************************************************************
      subroutine initial (
     .             lslip,npatch,nslip,nbrack,mwarn,amsg,iibl0,ytot,
     .             nobs,iibeg,iiend,istart,ilast,ii0,ii1,iibl,iibr,
     .             iprob,nprob)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'


      logical      lslip,good,lgood
      integer      nobs,npatch,nslip,nbrack,iibl0,nprob,iprob,i
      integer      iibeg,iiend,istart,ilast,ii0,ii1,iibl,iibr
      real*8       ytot
      integer*2    mwarn(4)
      character*96    amsg


      ytot = 0.d0
      nprob = 0
      iprob = 0
      if (npatch .eq. 0)  iibl0 = 1

      ilast  = ii1
      if (nbrack .eq. 2) ilast = iibr

c     search for first good observation
      good = .false.
      i = ii0 - 1
      do 200 while ((.not.good) .and. (i .le. ilast))
          i = i + 1
          if(lgood(kww(i))) good = .true.
  200 continue

      if (good) then
          iibeg = i

c         search for last good observation
          good = .false.
          i = ii1 + 1
          do 300 while ((.not.good) .and. (i .gt. iibeg+2))
              i = i - 1
              if(lgood(kww(i))) good = .true.
  300     continue
          if (good) then
               iiend = i
               ilast = min(ilast,iiend)
          else
               write (amsg,'(''Only one good data point'')')
               call gmsg (mwarn,amsg)
               lslip = .false.
          endif
      else
          write (amsg,'(''No data in this span'')')
          call gmsg (mwarn,amsg)
          lslip = .false.
      endif

      return
      end
c************************************************************************
      subroutine findslip_ed (
     .        nbrack,nobs,npatch,nslip,iibl0,ynoise,
     .        istart,ilast,iibr,iibl,l12e,ii0,ii1,
     .        slip,ngood,yend,iiend,ybeg,iibeg,epsilon,
     .        ymin2,ymax2,yadd,yslp,itemp,isite1,ichan1,
     .        key,mwarn,amsg,replot,domarg,nmlfit)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         nbrack,nobs,iibl,iibl0,iibr,l12e,ii0,ii1
      character*1     key
      integer         isite1,ichan1,nmlfit,ngood,npatch,nslip
      integer         j,istart,ilast,itemp,iibeg,iiend,iyp
      logical         replot,domarg
      logical         slip
      real*8          ymin2(maxwin),ymax2(maxwin)
      real*8          ybeg,yend,yslp,yadd
      real*8          signal,epsilon(4),ynoise(4)
      integer*2       mwarn(4)
      character*96    amsg

c     lower bound of signal to noise ratio
      parameter (signal = 5.0d0)

      if (slip) then
c          slip = .false.
          istart = iibl + 1
      else
          istart = max(iibeg,iibl)
          istart = iibeg
          call noise_lev (ynoise,iiend,iibeg,ymin2,ymax2,iyp)
          do 100 j=1,4
              epsilon(j) =  ynoise(j) * signal
  100     continue
      endif

      call ed_check_slip(
     .       iibl,iibl0,istart,ilast,iiend,iibeg,epsilon,
     .       ymin2,ymax2,mwarn,amsg,slip,npatch,nslip)


      return
      end
c************************************************************************
      subroutine ed_check_slip (
     .        iibl,iibl0,istart,ilast,iiend,iibeg,epsilon,
     .        ymin2,ymax2,mwarn,amsg,slip,npatch,nslip)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         npatch,nslip,iibl,iibl0,iobs
      integer         i,j,istart,ilast,iibeg,iiend,iprev
      logical         lgood,slip,good
      real*8          ymin2(maxwin),ymax2(maxwin)
      real*8          yzero(4),ybeg(4),yend(4),yexpb(4),yexpe(4)
      real*8          yobs(4),yprev(4),slopb(4),slope(4)
      real*8          dy,dyb,dye,epsilon(4)
      integer*2       mwarn(4)
      character*96    amsg



      yzero(1)  =  ymin2(1) + (ymax2(1) - ymin2(1))/2.d0
      yzero(2)  =  ymin2(2) + (ymax2(2) - ymin2(2))/2.d0
      yzero(3)  =  ymin2(3) + (ymax2(3) - ymin2(3))/2.d0
      yzero(4)  =  ymin2(4) + (ymax2(4) - ymin2(4))/2.d0

      ybeg(1) = wl1(iibeg) - yzero(1)
      yend(1) = wl1(iiend) - yzero(1)
      ybeg(2) = wl2(iibeg) - yzero(2)
      yend(2) = wl2(iiend) - yzero(2)
      ybeg(3) = wl1(iibeg) 
     .          - faclr(jsat1)*(wl2(iibeg)-gear(jsat1)*wl1(iibeg)) 
     .          - yzero(3)
      yend(3) = wl1(iiend) 
     .          - faclr(jsat1)*(wl2(iiend)-gear(jsat1)*wl1(iiend)) 
     .          - yzero(3)
      ybeg(4) = wl2(iibeg) - gear(jsat1)*wl1(iibeg) - yzero(4)
      yend(4) = wl2(iiend) - gear(jsat1)*wl1(iiend) - yzero(4)


c     get an initial slop if possible from previous points
      good = .false.
      i = istart
      do 100 while ((.not.good) .and. (i .gt. iibeg))
           i = i - 1
           if(lgood(kww(i))) good = .true.
  100 continue

      if (good) then
           iprev = i
c     get the initial slop from forward points
      else

           i = istart - 1
           do 200 while ((.not.good) .and. (i .lt. ilast))
                i = i + 1
                if(lgood(kww(i))) good = .true.
  200      continue
           iprev = i
      endif

      yprev(1)  = wl1(iprev) - yzero(1)
      yprev(2)  = wl2(iprev) - yzero(2)
      yprev(3)  = wl1(iprev) 
     .            - faclr(jsat1)*(wl2(iprev)-gear(jsat1)*wl1(iprev)) 
     .            - yzero(3)
      yprev(4)  = wl2(iprev) - gear(jsat1)*wl1(iprev) - yzero(4)

      do 300 j=1,4
          if (iprev .ne. iibeg) then
             slopb(j) = (yprev(j) - ybeg(j))/(iprev - iibeg)
          else
             slopb(j) = 1.0d8
          endif
          if (iprev .ne. iiend) then
              slope(j) = (yprev(j) - yend(j))/(iprev - iiend)
          else
              slope(j) = 1.0d8
          endif
  300 continue

c     start looking for out of range observation
      if (istart .gt. iprev) then
           i = istart - 1
      else
           i = iprev
      endif
      slip = .false.

      do 500 while ((.not.slip) .and. (i .lt. ilast))
          i = i + 1

          if(lgood(kww(i))) then
               iobs = i
               yobs(1)  = wl1(iobs) - yzero(1)
               yobs(2)  = wl2(iobs) - yzero(2)
               yobs(3)  = wl1(iobs) 
     .                 -  faclr(jsat1)*(wl2(iobs)-gear(jsat1)*wl1(iobs))
     .                 - yzero(3)
               yobs(4)  = wl2(iobs) - gear(jsat1)*wl1(iobs) - yzero(4)

               j = 0
               do 400 while ((.not.slip) .and. (j .lt. 4))
                   j= j + 1
                   if (iobs .ne. iibeg) then
                       yexpb(j) = ybeg(j) + slopb(j)*(iobs - iibeg)
                       slopb(j) = (yobs(j) - ybeg(j))/(iobs - iibeg)
                   else
                       yexpb(j) = 1.0d25
                   endif
                   if (iobs .ne. iiend) then
                       yexpe(j) = yend(j) - slope(j)*(iiend -iobs)
                       slope(j) = (yend(j) - yprev(j))/(iiend - iprev)
                   else
                       yexpe(j) = 1.0d25
                   endif

c                  get the smallest difference between observed and expected
                   dyb = abs(yobs(j)-yexpb(j))
                   dye = abs(yobs(j)-yexpe(j))
                   dy = min(dyb,dye)
                   if (dy .gt. epsilon(j)) then
                        slip = .true.
                   endif
  400          continue
            endif

  500   continue

        if(slip) then
              iibl = i
        else
              iibl = iibeg
        endif

      return
      end
c************************************************************************
      subroutine noise_lev (ynoise,iiend,iibeg,ymin2,ymax2,iyp)

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         iibeg,iiend,i,j,iyp,iobs
      logical         lgood
      real*8          yzero(4),ybeg(4),yend(4),yexpb(4),yexpe(4)
      real*8          ymin2(maxwin),ymax2(maxwin)
      real*8          yobs(4),slopb(4),slope(4)
      real*8          dy,sumdy(4),dyb,dye,ynoise(4)


c     get an initial estimate for epsilon
      yzero(1)  =  ymin2(1) + (ymax2(1) - ymin2(1))/2.d0
      yzero(2)  =  ymin2(2) + (ymax2(2) - ymin2(2))/2.d0
      yzero(3)  =  ymin2(3) + (ymax2(3) - ymin2(3))/2.d0
      yzero(4)  =  ymin2(4) + (ymax2(4) - ymin2(4))/2.d0

      ybeg(1) = wl1(iibeg) - yzero(1)
      yend(1) = wl1(iiend) - yzero(1)
      ybeg(2) = wl2(iibeg) - yzero(2)
      yend(2) = wl2(iiend) - yzero(2)
      ybeg(3) = wl1(iibeg) 
     .          - faclr(jsat1) *(wl2(iibeg)-gear(jsat1)*wl1(iibeg)) 
     .          - yzero(3)
      yend(3) = wl1(iiend) 
     .          - faclr(jsat1)*(wl2(iiend)-gear(jsat1)*wl1(iiend)) 
     .          - yzero(3)
      ybeg(4) = wl2(iibeg) - gear(jsat1)*wl1(iibeg) - yzero(4)
      yend(4) = wl2(iiend) - gear(jsat1)*wl1(iiend) - yzero(4)

      do 100 j=1,4
         slope(j) = (yend(j) - ybeg(j))/(iiend - iibeg)
         slopb(j) = 2.0d10 * slope(j)
         sumdy(j) = 0.0d0
  100 continue

c     start looking for noise level
      i = iibeg
      iyp = 0
      do 300 while (i .lt. iiend)

           i = i + 1
           if(lgood(kww(i))) then
                iobs  = i
                iyp = iyp + 1
                yobs(1)  = wl1(iobs) - yzero(1)
                yobs(2)  = wl2(iobs) - yzero(2)
                yobs(3)  = wl1(iobs) 
     .                  - faclr(jsat1)*(wl2(iobs)-gear(jsat1)*wl1(iobs))
     .                  - yzero(3)
                yobs(4)  = wl2(iobs) - gear(jsat1)*wl1(iobs) - yzero(4)

                do 200 j=1,4
                    yexpb(j) = ybeg(j) + slopb(j)*(iobs - iibeg)
                    yexpe(j) = yend(j) - slope(j)*(iiend -iobs)

                    slopb(j) = (yobs(j) - ybeg(j))/(iobs - iibeg)
                    if (iobs .ne. iiend) then
                        slope(j) = (yend(j) - yobs(j))/(iiend - iobs)
                    else
                        yexpe(j) = 1.0d25
                    endif

c                   get the smallest difference between observed and expected
                    dyb = abs(yobs(j)-yexpb(j))
                    dye = abs(yobs(j)-yexpe(j))
                    dy = min(dyb,dye)
                    sumdy(j) = sumdy(j) + dy
  200           continue

           endif

  300 continue

        do 400 j=1,4
            if(iyp .gt. 1) then
                 ynoise(j) = sumdy(j)/iyp
            else
                 ynoise(j) = 1.0d20
            endif
  400 continue

      return
      end
c************************************************************************

