c**************************************************************
      subroutine scannerdd (nobs,nsat,ncfls,iv,ir,idd)
c
c     SCAN  driver
c
c     Shimon Wdowindi, March 1991
c
c     based on the CVIEW editor
c
      implicit none

c     standard GAMIT array dimensions
      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      logical        seekob,domarg,getser
c     can an observable be formed?
      logical        lcando(-3:ntypes)

      integer       isite1,isite2,ichan1,ichan2,ijump,istep
      integer       nobs,ii0,ii1
     .,             iig01,iig11,iim01,iim11
     .,             iig02,iig12,iim02,iim12
      integer       nsat,ncfls,imode,ngood,l12p,iv,it,nmarg
     .,             kpl(maxepc),ngood1,ngood2,ir,idd
     .,             isite1_last,n1

      real*8        rms1,rms2,rms3
      real*8        yssd,ymin,ymax,yjump,ystep,yavg
     .,             plt(maxepc),percent_done,total,togo

      character*256 message

c     ignore marginal points
      domarg = .false.

c     **** initializing parameters *****
      seekob = .true.

      isite1_last = 0
      isite1 = 1
      ichan1 = 1
      ichan2 = 1
      if (ncfls .gt. 1) then
        isite2 = 2
      else
        isite2 = 0
      endif

      write (iv,40)
 40   FORMAT(1X,'RMS      = LC RMS of full solution (cycles)',/,
     .       1X,'FLAGGED  = Largest flagged jump (cycles)',/,
     .       1X,'UNFLAGD  = Largest UNFLAGGED offset',/,
     .       1X,'EPOCH(F) = Epoch of largest flagged jump',/,
     .       1X,'EPOCH(U) = Epoch of largest UNFLAGGED jump')
      write (iv,50)
 50   FORMAT(1X,'CHAN1',1x,'CHAN2',1X,'SIT1',2X,'SIT2',3X,'NDATA',6X,
     *       'RMS',4X,'EPOCH(F)',5X,'FLAGGED',2X,'EPOCH(U)',5X,
     *       'UNFLAGD')

      write (ir,60)
 60   FORMAT(5X,'quick  RMS',3x,'full   RMS',3x,'total RMS',3x,
     .        'ndata',3x,'CHAN1',1x,'CHAN2',1X,'SITE1',2X,'SITE2',/)

c     **** loop over all combinations *****
      do 100 while (seekob)

c         *** look for next series with data
          call nextob (ichan1,ichan2,nsat,isite1,isite2,ncfls,seekob)

c         check how far we have got thru the scanning
          if ( isite1 .gt. isite1_last ) then
             total = ncfls*(ncfls-1)/2
             togo =  (ncfls-isite1_last)*(ncfls-(isite1_last+1))/2
             percent_done = 100.d0 * (1.d0 - (togo / total))
             write(message,70)percent_done,isite1
 70          format('Scandd approx ',f5.1,'% completed. '
     .       ,'Now scanning against site: ',i3)
             call report_stat('STATUS','SCANDD','scannerdd',' ',
     .       message,0)
             isite1_last = isite1
          endif

          if (seekob) then
c             get L1, which is always needed
              it = 1
              lcando(it) = getser
     .        (ichan1,ichan2,isite1,isite2,it,imode,
     .        domarg,cl1,kcc,nobs,ngood1,
     .        iig01,iig11,iim01,iim11)
c             load working vectors
              call vcopy (kcc,cl1,kww,wl1,nobs)

c             get L2 if it is there
              it = 4
              lcando(it) = getser
     .        (ichan1,ichan2,isite1,isite2,it,imode,
     .        domarg,cl2,kcc,nobs,ngood2,
     .        iig02,iig12,iim02,iim12)
              call vcopy (kcc,cl2,kww,wl2,nobs)

c             use LC observable
              l12p = 3

c             choose end points wisely

              ii0 = min(iig01,iig02)
              ii1 = max(iig11,iig12)

c             form the observable, get extrema, and count good points
              call combo (nobs,ii0,ii1,l12p,domarg,kpl,plt,
     .        ymin,ymax,yavg,yssd,yjump,ystep,
     .        ngood,nmarg,ijump,istep,lcando)

              if (ngood .gt. 2) then

                   call findslip_scn (
     .                   isite1,isite2,ichan1,ichan2,idd,
     .                   nsat,nobs,ii0,ii1,rms1,rms2,rms3)
              endif

c             pushs bias and updates statistics
              call push_bias (
     .             isite1,isite2,ichan1,ichan2,kpl,plt,
     .             nobs,ii0,ii1,ymin,ymax,yavg,yssd,
     .             yjump,ystep,ngood,nmarg,ijump,istep)

              if (ngood .ge. 1) then

                  call get_rms (ir,ichan2,ichan1,isite2,isite1,
     .                              ii0,ii1,n1,rms1,rms2,rms3)
                  write(ir,4000) rms3,rms2,rms1,n1
     .                         , ichan2,ichan1,isite2,isite1
 4000             format (1x,f10.2,1x,3x,f10.2,1x,2x,f10.2,4x,i5,2x,4i6)
                  write (iv,5000) ichan1,ichan2,isite1,isite2
     .                          , ngood,rms2,istep,ystep,ijump,yjump
 5000             format (2x,4(i2,4x),i4,f11.2,2(i8,sp,f14.2,ss))

              endif

          endif

  100 continue

      return
      end
c
c**************************************************************
      subroutine push_bias (
     .             isite1,isite2,ichan1,ichan2,kpl,plt,
     .             nobs,ii0,ii1,ymin,ymax,yavg,yssd,
     .             yjump,ystep,ngood,nmarg,ijump,istep)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

c     functions to sort out error codes
      logical lgood,lbias,pbias

      integer nobs,ngood,ii0,ii1,i,nmarg,istep,ijump,ichan1,ichan2
     .      , isite1,isite2
      real*8 ymin,ymax,yavg,yssd,ylast
     .      ,yjump,ystep

      character*4 bias_count
      common /biascount/bias_count(maxepc)

c     what we are plotting
      real*8  plt(maxepc)
      integer*4 kpl(maxepc)

C     Bias pushing
      pbias = .false.
      do 99 i = ii0,ii1
         if(lgood(kww(i)) .and. pbias) then
               kww(i) = igbias
               pbias = .false.
CCCC           write(*,*) 'PUSHED BIAS at ',i,
CCCC .                     isite1,isite2,ichan1,ichan2
         else if((.not.lgood(kww(i))) .and.
     .          (bias_count(i) .ne. "0000")) then
               pbias = .true.
CCCC           write(*,*) 'BIAS at ',i
         endif
   99 continue


c     loop through again to get sample standard deviation
c     Also find largest jump and largest step
      ijump = 0
      istep = 0
      yssd =  0.0d0
      yjump = 0.0d0
      ystep = 0.0d0
      ylast = plt(ii0)
      do 200 i = ii0,ii1
         if (lgood(kww(i))) then
            yssd = yssd + (plt(i)-yavg)**2
            if (.not. lbias(kww(i))) then
               if (dabs(plt(i)-ylast) .gt. dabs(yjump)) then
                  yjump = plt(i) - ylast
                  ijump = i
               endif
            else if (lbias(kww(i))) then
               if (dabs(plt(i)-ylast) .gt. dabs(ystep)) then
                  ystep = plt(i) - ylast
                  istep = i
               endif
            endif
            ylast = plt(i)
         endif
 200  continue

      if (ngood. gt. 1) then
          yssd = dsqrt(yssd/dble(ngood-1))
      else
          yssd = 0.0d0
      endif

      return
      end
c
c*************************************************************************
      subroutine get_rms (ir,ichan2,ichan1,isite2,isite1,
     .                              ii0,ii1,n1,rms1,rms2,rms3)

      implicit none

      include '../includes/dimpar.h'
      include '../includes/makex.h'
      include '../includes/cview.h'
      include '../includes/errflg.h'

      integer   isite1,isite2,ichan1,ichan2,ir
      integer   ii0,ii1,i,j
      integer   n1,n2,n3,igap
      integer   n1a,n1b,n2a,n2b,n3a,n3b
      real*8    rms1,rms2,rms3
      real*8    r1,r2,r3,ylc
      logical   lgood,lbias

      if (ii1 .gt. ii0) then

c         finds number of epoches between two data points
          inext = max(jnext(isite1),jnext(isite2))

c         sum all LC data
          rms1 = 0.0
          rms2 = 0.0
          rms3 = 0.0

          r1 = 0.0
          r2 = 0.0
          r3 = 0.0

          n1 = 0
          n2 = 0
          n3 = 0
          igap = 0

          n1a = ii0
          n2a = ii0
          n3a = ii0

          do 300 i = ii0,ii1,inext

c            sum each segment of QUICK solution data points
c             if (lbias(kww(i)) .or. (i .eq. ii1)) then
              if (i .gt. ii0) then
              if (lbias(kww(i))            .or.
     .           (.not.lgood(kww(inext))) .or.
     .           (i .eq. ii1))             then

                 if (n3 .gt. 1) then
                    r3 = r3/dble(n3)
                    n3b = i - inext
                    if((i .eq. ii1) .and.
     .                (.not. lbias(kww(i)))) n3b = i

                    do 100 j = n3a,n3b
                       if (lgood(kww(j))) then
                       ylc = wl1(j) - faclr(jsat1)
     .                     * (wl2(j)-gear(jsat1)*wl1(j))
                          rms3 = rms3 + (ylc-r3)**2
                       endif
 100                 continue
                 endif

                 r3 = 0.0
                 n3 = 0
                 n3a = i
             endif
             endif

c            sum each segment of FULL solution data points
              if (i .gt. ii0) then
              if (lbias(kww(i)) .or.
     .           (i .eq. ii1))  then

                 if (n2 .gt. 1) then
                    r2 = r2/dble(n2)
                    n2b = i - 1
                    if((i .eq. ii1) .and.
     .                (.not. lbias(kww(i)))) n2b = i

                    do 200 j = n2a,n2b
                       if (lgood(kww(j))) then
                       ylc = wl1(j) - faclr(jsat1) 
     ,                     * (wl2(j)-gear(jsat1)*wl1(j))
                          rms2 = rms2 + (ylc-r2)**2
                       endif
 200                 continue
                 endif

                 r2 = 0.0
                 n2 = 0
                 n2a = i
             endif
             endif

c            sum TOTAL data points
             if (lgood(kww(i))) then
                 ylc = wl1(i)-faclr(jsat1)*(wl2(i)-gear(jsat1)*wl1(i))

                 r1 = r1 + ylc
                 r2 = r2 + ylc
                 r3 = r3 + ylc

                 n1 = n1 + 1
                 n2 = n2 + 1
                 n3 = n3 + 1

                 igap = 0
             else

                 igap = igap + 1
             endif

 300      continue

          if (n1 .gt. 1) then
              r1 = r1/dble(n1)
          else
              r1 = 0.0
          endif

          n1b = ii1
          do 400 j = n1a,n1b
              if (lgood(kww(j))) then
                 ylc = wl1(j) - faclr(jsat1)*(wl2(j)-gear(jsat1)*wl1(j))
                 rms1 = rms1 + (ylc-r1)**2
               endif
 400     continue

         if (n1 .gt. 1) then
            rms1 = dsqrt(rms1/dble(n1-1))
            rms2 = dsqrt(rms2/dble(n1-1))
            rms3 = dsqrt(rms3/dble(n1-1))
         else
            rms1 = 0.0d0
            rms2 = 0.0d0
            rms3 = 0.0d0
         endif
      else
         rms1 = 0.0d0
         rms2 = 0.0d0
         rms3 = 0.0d0
      endif

      return

      end


