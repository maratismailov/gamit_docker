c
      Subroutine BISOPT( mopt )

c     Narrow_lane ambiguity resolution for L1_ONLY,L2_ONLY - mopt=1
c     Narrow_lane ambiguity resolution for L1L2_INDEPENDENT - mopt=3
c     Narrow-lane ambiguity resolution for LC_HELP or LC_RANGE solutions
c          called with mopt=1 from SOLVE
c     Wide-lane and narrow-lane for L1&L2 (with ionospheric constraint) solutions
c           BISOPT called twice from SOLVE, first with mopt=2 for widelane,
c           then with mopt=1 for narrowlane.  BISOPT is not called for LC_HELP or
c           LC_AUTCLN wide-lane resolution.

c     Approach options :

c        1. take nearest integers
c        2. take nearest integers (WL using pseudo-range) (no use)
c        3. take expective real values (using cumulative probability)
c        4. take expective real values sequentialy (JPL)
c        5. bias rounding + bias searching (LSQX default)
c        6. decision function (GAMIT default)
c        7. decision function with PR priority (Bob's favor)

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'                                
      include 'parameters.h'

c     Control for searching, set according to method:

c      mopt = 1  search L1 (or L2) biases only
c           = 2  search L2-L1 biases only
c           = 3  search both L2-L1 and L1 biases

      integer mopt,lane,k1,k2,k3,issn,isla
     .      , ia1,ia2,ifix,nlive0,idiag,i1,ii,i,k

      real*8 bint,s1,adn,smax,bepc,slarge,sgm,adb,biskut

      call report_stat('STATUS','SOLVE','bisopt',' '
     .          , 'Resolving narrow-lane ambiguities',0)

      
      goto (100,200,300,400,500,500,500,600), noptin
c
 100  continue
      do 110 i = 1,mbias
         i1 = idxb(i)
         if (i1.le.0) goto 110
         bint = dint(adjust(i1)+.5d0*dsign(1.d0,adjust(i1)))
         adjust(i1) = bint
c reduce the number of live perameters.
         nlive = nlive-1
         free(i1) = 0
         idxb(i) = -i1
 110  continue
      goto 2000
c
 200  continue
      goto 2000
c
 300  continue
      do 320 i = 1,mbias
         ii = idxb(i)
c----    dead bias parameters
         if (ii.le.0) goto 320
         ii = ii+1
         s1 = sigma(ii+msig-1)
         adn = adjust(ii)
         call bexpec(ii,adn,sclerr,s1,smax,bepc,0)
         free(ii) = 0
         nlive = nlive-1
         adjust(ii) = bepc
         idxb(i) = -ii
 320  continue
      goto 2000
c
 400  continue
 405  k1 = msig-1
      issn = 0
      k2 = 0
      k3 = 0 
      isla = 0
c---- calculate number of live bias parameters for searching
c---- find bias with maximum probability.
      slarge = -1.0d0
         k1 = k1+k2
         ia1 = issn*iband+1
         issn = issn+l1bias
         ia2 = issn*iband
         k2 = 0
         do 410 i = ia1,ia2
            k = idxb(i)
            if(k.le.0) go to 410
            k2 = k2+1
            if (mopt.eq.1.and.i.gt.issn) goto 410
            if (mopt.eq.2.and.i.le.issn) goto 410
            k3 = k3+1
            s1 = sigma(k2+k1)
            adn = adjust(k)
            call bexpec(k,adn,sclerr,s1,smax,bepc,1)
            if (smax.le.slarge) goto 410
            slarge = smax
            sgm = s1
            isla = i
            ii = k
            adb = adn
 410     continue
      if (k3.le.0) goto 2000
c---- fix this bias parameter
      call bexpec(ii,adb,sclerr,sgm,smax,bepc,0)
c      write (6,*) 'bias index ', isla,idxb(isla)
      free(ii) = 0
      nlive0 = nlive
      nlive = nlive-1
      bdev(1) = bepc-adjust(ii)
      adjust(ii) = bepc
c      The following cause and a non-understood error on the HP:   
c    Compiler error line 125 of bisopt.f: Undefined NSYM for ucode address :   (5726)
c      if( isla.le.0 ) call report_stat('FATAL','SOLVE',bisopt,' '
c     .                   ,'Bias index = 0; something wrong',0)
      idxb(isla) = -ii
c---- solve normal equation with the bias fixed
      ifix  =  nlive0 - nlive
      call bnew( nlive0,ifix )
      if(nlive.le.0) go to 2000
c---- update sigma array
      do 450 i = 1,nlive
         idiag = i*(i+1)/2
         sigma(i) = dsign(dsqrt(dabs(a(idiag))),a(idiag))
 450  continue
      goto 405
c=============================================================
 500  continue

      if( mopt.eq.1 .or. mopt.eq.3) then
        if( logprt ) write(6,540) nldev,nlsig,nlcut,nlrat,nldmax
        write(10,540) nldev,nlsig,nlcut,nlrat,nldmax
 540  format(/,' Narrow-lane bias-fixing criteria: ',
     1           ' deviation : ',f5.2,'  sigma : ',f5.2,
     2           '  decision func. : ',f8.1,' ratio : ',f6.1,/,
     3           36x,'maximum distance :',f8.1)
        if (nlrat.gt.1.0d-3.and.nlrat.lt.1.0d6)  biskut = nlrat
c       narrowlane
        lane = 1
      elseif (mopt.eq.2 ) then
        if( logprt ) write(6,550) 
     .      wldev,wlsig,wlcut,wlrat,wldmax,prdev,prsig,prcut
        write(10,550) wldev,wlsig,wlcut,wlrat,wldmax,prdev,prsig,prcut
 550  format(/,' Wide-lane bias-fixing criteria'/,
     .           ' Phase:       deviation  ',f5.2,'  sigma  ',f5.2,
     .           '  decision func.  ',f8.1,' ratio  ',f6.1,/,
     .           '  maximum distance  ',f8.1,/,
     .           ' Pseudorange: deviation  ',f5.2,'  sigma  ',f5.2,
     .           '  decision func.  ',f8.1)
        if (wlrat.gt.1.0d-3.and.wlrat.lt.1.0d6)  biskut = wlrat
c       widelane
        lane = 2
      else
        call report_stat('FATAL','SOLVE','bisopt',' '
     .                  ,'mopt must be 1, 2 or 3 ',0)
      endif
c      write(*,*) 'lane,biskut,mopt,noptin',lane,biskut,mopt,noptin
      call nbias( lane,biskut,mopt )
      goto 2000
 600  continue

c-----print the number of live and dead parameters

 2000 if( mopt.eq.1 .or. mopt.eq.3) then
        if( logprt ) write(6,2001)
        write(10,2001)
 2001   format(/,' Narrow-lane bias-fixing complete')  
      elseif( mopt.eq.2 ) then
        if( logprt ) write(6,2002)
        write(10,2002)
 2002   format(' Wide-lane bias-fixing complete') 
      endif
      call qhead1 (2)

      return
      end
