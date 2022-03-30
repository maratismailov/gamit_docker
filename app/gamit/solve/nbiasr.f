      Subroutine NBIASR( lane,mopt )

c      lane    : 1=NL (NBIAS called by BISOPT)   2= WL (NBIAS called by GETWL)

c     round biases to nearest integer:
c      1. decision function approach (see Dong & Bock, 1989)
c      2. chi-square contrast testing approach
c       (1) add 2-sigma uncertainty to real-valued bias
c           (uncertainty is 3 times scaled formal errors,
c            so actually 6xsigma)
c       (2) check if nearest integer bias is spanned within a half cycle
c       (3) if passes test, fix integer bias
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer lane

c       input criteria for wide-lane (wl) and narrow-lane (nl) bias-fixing
c       -- all real*8
c       common/bcrit/wldev,wlsig,wlcut,wlrat,wldmax,
c                    nldev,nlsig,nlcut,nlrat,nldmax,
c                   , prdev,prsig,prcut
c       decision-function parameters:
c           dev:  threshold deviation for taper function (= 0.4 in Dong and Bock;
c                 default now 0.15 in BDECI and in FIXDRV)
c           sig:  threshold sigma for taper function (=0.33 in Dong and Bock;
c                 default now 0.15 in BDECI and in FIXDRV)
c           cut:  decision function d(x,sigma) = FT/Q, inverse of
c                 1. - allowable rate for a type 1 error (alpha),
c        chi-square search parameters:
c           rat:  threshhold ratio for accepting search
c
      logical away

      integer ipseu,ih,llbad,ia1,ia2,ifix,l1b,icut
     .      , ld1,ld2,junk(15),issn,mopt
     .      , junk1(15),junk2(15),j1,j2,j3,k1,k2
     .      , i,j,k

      real*8 xintp,xdiff1,xdiff2,deci,dmax(15),prob
     .     , cr1,cutdec,xint,xlow,srcrge,cutdev,pright(15),cutsig
     .     , add,dev,sma,awl,est,awv,xhigh,a1,a2
* MOD TAH 031203: Added sma_save (This is the sigma of the fixed
*     bias)
      real*8 sma_save

c---- test fix bias one by one approach.  880923

c     initialization
      icut=1
      ifix=0
      prob=1.0d0
      ipseu  =  0             
c     the following necessary to avoid compiler warnings
      cutdec = 0.d0

c     set parameters for wide-lane or narrow-lane
      if (lane.eq.1 ) then
         cutdev = nldev
         cutsig = nlsig
         cutdec = nlcut
      else if (lane.eq.2 ) then
         cutdev = wldev
         cutsig = wlsig
         cutdec = wlcut
      else
         call report_stat('FATAL','SOLVE','nbiasr',' '
     .                      ,'Bad value of lane',0)
      endif
      if (cutdec.lt.1.0d0) cutdec = 1000.0d0
      call zero1d(1,icut+1,dmax)
      call bfwork(mopt,issn,k1,k2,ld1,ld2,llbad,ia1,ia2,
     .     l1b,i,k,away,1)
         call bfwork(mopt,issn,k1,k2,ld1,ld2,llbad,ia1,ia2,
     .        l1b,i,k,away,2)    
c         print *,' NBIASR ia1 ia2 ',ia1,ia2
         do 2 i=ia1,ia2                              
c            print *,' At 2 i = ',i
            if (i.gt.issn) ipseu = ipseu + 1
            call bfwork(mopt,issn,k1,k2,ld1,ld2,llbad,ia1,ia2,
     .           l1b,i,k,away,3)  
c            print *,' aft BFWORK ipseu away ',ipseu,away
            if (away) goto 2
c
            ih=half(i)   
            est=adjust(k)
            sma=sigma(k1+k2)
c **new code
            if( lane.eq.2 ) sma = sclerr*sigma(k1+k2)
c----    bias decision function approach.  -dnd- 880403
            if (noptin.eq.6.or.noptin.eq.7) then
c              old code:
c              Should no longer do this since we don't scale full solution anymore:
c              sma=sma*3.0d0
c              call bdeci(adn,ih,sclerr,sma,prob,deci,factor)
c              new code:
               call bdeci( est,sma,ih,cutdev,cutsig,prob,deci )

               if (deci.le.cutdec) goto 2
* MOD TAH 031203: Output for fixed bias? These are biases that could
*              fixed.  The values still need to be ranked and the
*              highest deci value is fixed.  Skip output.
C              write(10,505) k, est, sma, deci, ifix
C505           format('Bias Fixd? Number ',i4,' Est ',F9.2,' +/- ',
C    .                 F9.2,' Deci ',F12.1,' Ifix ',i4)
               goto 60
            endif
c for the time being, choose factor 6.0 -dnd- 871210
            add=est*dble(ih)
            xint=dint(add+.5d0*dsign(1.d0,add))/dble(ih)
            srcrge=6.d0*sclerr*sma
            xhigh=est+srcrge
            xlow=est-srcrge     
c            print *,'xint xhigh xlow ',xint,xhigh,xlow
            if(xint.gt.xhigh.or.xint.lt.xlow) go to 2
            xdiff1=xhigh-xint
            xdiff2=xint-xlow
            cr1=0.5d0/dble(ih) 
c            print *,'cr1 xdiff1 xdiff2 ',cr1,xdiff1,xdiff2
            if(xdiff1.gt.cr1.or.xdiff2.gt.cr1) go to 2
            deci=sma
 60         continue                   
c            print *,'icut dmax deci ',icut,dmax(icut),deci
            if (deci.le.dmax(icut)) goto 2
            if (ifix.lt.icut) ifix=ifix+1
            do 70 j=1,ifix   
c            print *,'j dmax deci ',j,dmax(j),deci
            if (deci.le.dmax(j)) goto 70  
c            print *,'NBIASR ifix ',ifix
            if (ifix.eq.1) goto 85
            do 80 j1=j,ifix
               j3=ifix-j1+j
               j2=j3+1
               dmax(j2)=dmax(j3)
               nfix(j2)=nfix(j3)
               junk(j2)=junk(j3)
               junk1(j2)=junk1(j3)
               junk2(j2)=junk2(j3)
               pright(j2)=pright(j3) 
c               print *,' 80 loop  j j1 j2 j3 junk junk1 '
c     .              ,j,j1,j2,j3,junk(j2),junk1(j2)
  80        continue
  85        dmax(j)=deci
            nfix(j)=k1+k2  
            junk(j)=k   
            junk1(j)=i
            junk2(j) = ipseu   
c            print *,' aft 85 k1 k2 i j junk junk1 '
c     .              ,k1,k2,i,j,junk(j),junk1(j)

c            print *,'NBIASR i j k ',i,j,k
            pright(j)=prob
* MOD TAH 031203: Code above is executed for the high decision value
*           found so far.  Save the sigma of the bias at this point
            sma_save = sma
            goto 2
  70        continue
   2     continue
         if (iband.eq.2) issn=issn+l1bias
C---- Reorder array by junk1 to fit BNEW (dmax is working array only)
      if (ifix.eq.1) goto 1417
      call bsort(ifix,pright,junk1,dmax,junk,nfix,junk2,2)
c
1417  continue
c     check fixed biases against pseudorange estimates
c     (NBIASP called later by NBIAS to use pseudoranges with biases
c      unresolved by the decision function here)
      do 90 i=1,ifix
         k=junk(i)
         k1=junk1(i)
         est=adjust(k)
         add=est*dble(half(k1))
c        nearest integer (or half-integer) for phase/ion estimate
         xint=dint(add+.5d0*dsign(1.d0,add))/dble(half(k1))
c        deviation of phase estimate from integer
         bdev(i)=xint-est
c        store phase-estimate integer in 'adjust'
         adjust(k)=xint
         free(k)=0
         idxb(k1)=-k
         nlive=nlive-1
c        mopt=2 means searching wide-lane (why test here?)
         if (mopt.eq.2) then
            k2 = junk2(i)
            awl = ddwl(k2)
            awv = ddwv(k2)
c           nearest integer from PR estimate
            xintp=dint(awl+.5d0*dsign(1.d0,awl))
c           difference between nearest integers for phase and PR estimates
            dev = dabs(xint-xintp)
c           if the pseudorange sigma is < prsig (default=0.1) and the
c           phase and PR integers are different, see if PR is better
            if (awv.lt.prsig .and. dev.gt.0.9d0) then
               write (6,120) est,awl
               write (10,120) est,awl
c  I am bypassing this since I found two cases where wrong widelane was
c  fixed according to PR when it was obvious that phase estimate was
c  correct -- YB 94/7/2
c  This code now restored but can be bypassed by inputting small values
c  of the PR deviation and sigma -- RWK 94/11/16
c              compare the two estimates
                a1 = bdev(i)
                a2 = xintp-awl
c              see if TWICE the PR deviation is less than the phase deviation
c              (and make sure the two are within 0.5 cycles to guard against
c               a bad PR estimate)
                if (dabs(a2*2.0d0).lt.dabs(a1).and.(dev.lt.0.5d0)) then
                   bdev(i) = xintp-est
                   write(6,130) xintp,adjust(k)
                   write(6,130) xintp,adjust(k)
                   adjust(k) = xintp
                endif
            endif
         endif
* MOD TAH 031202: Added additional output to better see what is happening.
*        It seems in the above code that ifix is never more than 1 for the
*        rounding code.
         if( logprt ) write (6,110) k,rlabel(k)
     .        , est,sma_save, adjust(k), dmax(1)
         write (10,110) k,rlabel(k),est,sma_save, adjust(k), dmax(1)
  90   continue
 110  FORMAT ('Fix No.',I5,1x,a20,' bias  from',F9.2,' +- ',F9.2,
     .        ' to',F8.1,' Decision Function ',E12.4)
 120  format (2x,'Warning: big difference between 2 WLs (phase ',
     .   f11.2, '   PR ',f11.2,')')
 130  format (2x,'Warning: PR WL (',f11.2,') overriding phase est ('
     .   ,f11.2,')')
c
       return
       end
