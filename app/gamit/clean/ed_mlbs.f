c************************************************************************
      subroutine multibias (
     .         lfound,seekob,newser,replot,redraw,newbox,domarg,lslip,
     .         nslip,nbrack,nobs,ntemp,npatch,nwins,iiend,iibeg,itemp,
     .         jsite1,jsite2,jchan1,jchan2,isite1,isite2,ichan1,ichan2,
     .         jsite,jchan,key,pos,l12e,l12v,l12p,ngood,nsave,nmlfit,
     .         ymin2,ymax2,yadd,yend,ybeg,yslp,ymult,yold,slip,nextslip,
     .         ii0,ii1,iibl,iibr,iibl0,epsilon,istart,ilast,imenu,
     .         mplot,mrect,mtaxi,menus,mclip,csize,mwarn,amsg,igerr,
     .         amenus,lbx,rbx,ix1,ix2,iy1,iy2,ix0,i0,dt,xs,ttt,ytot,
     .         iprob,nprob,pbias,simple,nsat,ncfls,imode)


      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      logical       lfound,replot,slip,nextslip,lslip
      logical       seekob,newser,redraw,newbox,domarg,simple
      logical       s1,s2,c1,c2
      integer       npatch,nslip,nwins,istart,ilast,iibl0
      integer       jsite,jchan,jchan1,jchan2,itemp,nmlfit
      integer       imode
      integer       nbrack,nobs,ntemp,iiend,iibeg,ngood,nsave
      integer       isite1,isite2,ichan1,ichan2,jsite1,jsite2,imenu
      integer       iibl,iibr,i0,ii0,ii1,l12v(maxwin),l12e,l12p
      integer       iprob,nprob,pbias(maxepc),nsat,ncfls
      integer*2     lbx,rbx,ix1,iy1,ix2,iy2,ix0,pos(2)
      real*8        ytot,yadd,yend,ybeg,yslp,ymult,yold,dt,xs
      real*8        ymin2(maxwin),ymax2(maxwin),epsilon,ttt
      character*1   key
      integer*2     mplot(maxwin,4),mrect(4),mtaxi(4),menus(nmenus,4)
      integer*2     mwarn(4),csize(2),mclip(4)
      integer*4     igerr
      character*16  amenus(nmenus)
      character*96  amsg
                                              

c     dummy statement for avoid warning for unused calling argument
      if( itemp.eq.0 ) then
        continue
      endif


c     **************************************************************
c     *** Multi-bias cycle-slip finder
c     **************************************************************


       call ed_find2(
     .       iibl,ntemp,s1,s2,c1,c2,mwarn,
     .       isite1,isite2,ichan1,ichan2)

       if (ntemp .eq. 2) then
           if (s1 .and. .not. s2) then
                    jsite = isite1
                    jchan = 0
                    call formsd (
     .                   isite1,isite2,ichan1,ichan2,jsite,jchan,
     .                   nobs,iibl,imode,domarg)
           else if (s2 .and. .not. s1) then
                    jsite = isite2
                    jchan = 0
                    call formsd (
     .                   isite1,isite2,ichan1,ichan2,jsite,jchan,
     .                   nobs,iibl,imode,domarg)
           else if (c1 .and. .not. c2) then
                    jsite = 0
                    jchan = ichan1
                    call formsd (
     .                   isite1,isite2,ichan1,ichan2,jsite,jchan,
     .                   nobs,iibl,imode,domarg)
           else if (c2 .and. .not. c1) then
                    jsite = 0
                    jchan = ichan2
                    call formsd (
     .                   isite1,isite2,ichan1,ichan2,jsite,jchan,
     .                   nobs,iibl,imode,domarg)
           else
                    write (amsg,'(''I still do not know how to do it '',
     .                                      i4)')iibl

           endif

           write (amsg,'(''Patched site '',i2," chans",i3)')
     .                                                    jsite,jchan
           call gmsg (mwarn,amsg)
       else if (ntemp .gt. 2) then
             write (amsg,'(''at '',2i4," bias - too complicated")')
     .              iibl,ntemp
             call gmsg (mwarn,amsg)
       else
             write (amsg,'(''at '',2i4," bias ?????????????")')
     .              iibl,ntemp
             call gmsg (mwarn,amsg)
       endif
       call gmsg (mwarn,amsg)

      return
      end
c

c************************************************************************
      subroutine newdd (
     .                   isite1,isite2,ichan1,ichan2,jsite,jchan,
     .                   nbrack,nobs,iibr,iibl,l12e,ii0,ii1,ngood,
     .                   yend,iiend,ybeg,iibeg,ytot,yadd,yslp,itemp,
     .                   key,mwarn,amsg,replot,domarg,nmlfit,nsave,
     .                   imenu,newbox,redraw,lfound,imode,
     .                   jsite1,jsite2,jchan1,jchan2)



      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         iibl,iibr,nobs,nbrack,l12e,ii0,ii1,ngood
      integer         iiend,iibeg,itemp,nmlfit,nsave,imenu,imode
      integer         isite1,isite2,ichan1,ichan2,jsite,jchan,it
      integer         jsite1,jsite2,jchan1,jchan2,iig0,iig1,iim0,iim1
      logical         lfound,redraw,newbox,replot,domarg,dummy,getser
      real*8          ytot,yadd,yend,ybeg,yslp
      integer*2       mwarn(4)
      character*96    amsg
      character*1     key



c     dummy statement for avoid warning for unused calling argument
      if( itemp.eq.0 ) then
        continue
      endif


c     get L1, which is always needed
      it = 1
      dummy = getser (ichan1,ichan2,isite1,isite2,it,imode,
     .         domarg,cl1,kcc,nobs,ngood,iig0,iig1,iim0,iim1)
      call vcopy (kcc,cl1,kww,wl1,nobs)
c     get L2 if it is there
      it = 4
      dummy = getser (ichan1,ichan2,isite1,isite2,it,imode,
     .         domarg,cl2,kcc,nobs,ngood,iig0,iig1,iim0,iim1)
      call vcopy (kcc,cl2,kww,wl2,nobs)


      return
      end
c************************************************************************
      subroutine formsd (
     .                   isite1,isite2,ichan1,ichan2,jsite,jchan,
     .                   nobs,iibl,imode,domarg)



      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'
      include '../includes/macdep.h'

      integer         iibl,iprev,i,nobs,imode
      integer         is11,is12,is21,is22,ic11,ic12,ic21,ic22
      integer         ngood,iig0,iig1,iim0,iim1
      integer         isite1,isite2,ichan1,ichan2,jsite,jchan,it
      logical         lgood,good,dummy,getser,domarg
      real*8          dl1(maxepc),dl2(maxepc),d1l1,d1l2,d2l1,d2l2
      integer*2       mwarn(4)
c     Error flags must be integer*4! (Kurt June 91)
ckf   integer*2       kdd(maxepc)
      integer*4       kdd(maxepc)
      character*96    amsg

      iprev = 0

      if(jchan .eq. 0) then
         is11 = isite1
         is21 = isite2
         ic11 = ichan1
         ic21 = 0
         is12 = isite1
         is22 = isite2
         ic12 = 0
         ic22 = ichan2
      else
         is11 = isite1
         is21 = 0
         ic11 = ichan1
         ic21 = ichan2
         is12 = 0
         is22 = isite2
         ic12 = ichan1
         ic22 = ichan2
      endif

c     get L1, which is always needed
      it = 1
      dummy = getser (ic11,ic21,is11,is21,it,imode,
     .         domarg,dl1,kdd,nobs,ngood,iig0,iig1,iim0,iim1)
c     get L2 if it is there
      it = 4
      dummy = getser (ic11,ic21,is11,is21,it,imode,
     .         domarg,dl2,kdd,nobs,ngood,iig0,iig1,iim0,iim1)

c     get location of previous good point
      good = .false.
      i = iibl
      do 100 while ((.not.good) .and. (i .gt. 1))
           i = i - 1
           if(lgood(kdd(i))) good = .true.
c           if(lgood(kww(i))) good = .true.
  100 continue

      if (good) then
           iprev = i
      else
           write (amsg,'(''Wrong single diff. - chan'',i3)') ichan1
           call gmsg (mwarn,amsg)
      endif

      d1l1 = abs(dl1(iibl) - dl1(iprev))
      d1l2 = abs(dl2(iibl) - dl2(iprev))


c     get L1, which is always needed
      it = 1
      dummy = getser (ic12,ic22,is12,is22,it,imode,
     .         domarg,dl1,kdd,nobs,ngood,iig0,iig1,iim0,iim1)
c     get L2 if it is there
      it = 4
      dummy = getser (ic12,ic22,is12,is22,it,imode,
     .         domarg,dl2,kdd,nobs,ngood,iig0,iig1,iim0,iim1)
c      call vcopy (kcc,cl2,kww,wl2,nobs)

c     get location of previous good point
      good = .false.
      i = iibl
      do 300 while ((.not.good) .and. (i .gt. 1))
           i = i - 1
           if(lgood(kdd(i))) good = .true.
  300 continue

      if (good) then
           iprev = i
      else
           write (amsg,'(''Wrong single diff. - chan'',i3)') ichan2
           call gmsg (mwarn,amsg)
      endif

      d2l1 = abs(dl1(iibl) - dl1(iprev))
      d2l2 = abs(dl2(iibl) - dl2(iprev))

      write (amsg,'(''dls'',4f7.2)') d1l1,d1l2,d2l1,d2l2
      call gmsg (mwarn,amsg)

      if(d1l1 .gt. d2l1 .and. d1l2 .gt. d2l2) then
          if(jchan .eq. 0) then
              jchan = ichan1
          else
              jsite = isite1
          endif
      else if(d1l1 .lt. d2l1 .and. d1l2 .lt. d2l2) then
          if(jchan .eq. 0) then
              jchan = ichan2
          else
              jsite = isite2
          endif
      else
c           write (amsg,'(''Ambigous '',i3)') iibl
c           call gmsg (mwarn,amsg)
      endif

      return
      end
c************************************************************************

