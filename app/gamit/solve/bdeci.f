       subroutine bdeci( est,sigma,ih,cutdev,cutsig
     .                , prob,deci )
c
c     Calculate the decision-function value as per Appendix A
c     of Dong and Bock [1989].
c        Da-nan Dong  880403
c     Add damping factor. DND  880602
c     Change definition of decision function based on hypothesis
c        test theory.      DND 880929
c     Remove arbitrary factor of 3., scale wide-lane but not narrow-lane
c        sigmas.    King 930322

c     Input:
c       est:     estimated (real) bias value
c       sigma:   estimated uncertainty of bias value, scaled by nrms for
c                  widelane, unscaled for narrow lane
c       ih :     control for receiver ambiguity
c                  = 1  unit=one cycle
c                  = 2  unit=half cycle
c       cutdev:  threshold deviation for taper function (= 0.4 in Dong and Bock;
c                  default still 0.4 here and in FIXDRV)
c       cutsig:  threshold sigma for taper function (=0.33 in Dong and Bock;
c                  default still 0.4 here and in FIXDRV)

c     Output:
c       prob : area of decision-function region, set = 1.0 for
c              each bias on input to BDECI, reduced by the probability
c              of an error in rounding each bias. (Since we now search
c              only one bias at a time, the cumulative probability is
c              not calculated by NBIASR or NBIASP.)
c       deci : decision function d(x,sigma) = FT/Q, inverse of
c              1. - allowable rate for a type 1 error (alpha),
c              compared with input wlcut or nlcut in NBIASR.
c              (But F is no longer used because we search one bias at
c              a time.)
c
      implicit none

      real*8 est,sigma,cutdev,cutsig,prob,dev,deci
     .     , term1,term2,c,d1,a1,bint,taper,add
     .     , cdev,csig,trun,s1,s2
      real*8 b1,b2,erfc,erfcb1,erfcb2
      integer*4 ih,j                  
      external erfc
         
c** Flag for debug printout
      integer*4 idebug
c      idebug = 0  : none 
c               1  : print to q- and log file decision values for fixed and not-fixed biases
c               2  : print to log the fixing criteria
c               3  : full print of logic 
       idebug = 1    
   


      data s2/1.414213562373095d0/

c
c     compute the deviation of the estimated value from an integer or half-integer
      add = est*dble(ih)
      bint = dint(add+.5d0*dsign(1.d0,add))/dble(ih)
      dev = dabs(bint-est)

c     set the cutoff deviation from the input or default value
      cdev = cutdev
c     default was 0.4 prior to 930315, now 0.15
      if (cutdev.lt.1.0d-3) cdev = 0.15d0
c     this was 0.6d0*cdev by mistake, prior to 930319
      if (ih.eq.2) cdev = 0.5d0*cdev
          
      if( idebug.ge.2 ) 
     .   print *,'BDECI cutdev dev cdev sig cutsig '
     .     ,cutdev,dev,cdev,sigma,cutsig
c     if the deviation is greater than the cutoff, set prob and deci and exit
      if (dev.ge.cdev) then
         prob =   1.0d0
         deci = 0.0d0
         goto 100
      endif

c**old code:
c     scale the estimated sigma by the nrms (sclerr) for both widelane and narrowlane
c     (these should be treated differently)
c     sigtemp = sclerr*sigma
c     s1 = 1.0d0/(sigtemp*s2)
c**new code
      s1 = 1.0d0/(sigma*s2)
c     numerical truncation tolerance
      trun = 1.d-9

c     compute the taper (T)
c     this term is (1 - dev/0.4) in Dong and Bock;
      term1 = 1.0d0-dev/cdev
c     this term is (1. - 3*scaled_sigma) in Dong and Bock; since cutsig is
c     0.4 now by default, term2 = 1.3 - 3*scaled_sigma
c     Dong and Bock:  term2 = (.333 - sigtemp)*3.
c     New:
      csig = cutsig
c     default changed from 0.4 to 0.15 930319
      if( cutsig.lt.1.0d-3) csig = 0.15d0
c**old code:
c     term2 = (csig-sigtemp)*3.d0
c**new code
      term2 = (csig-sigma)*3.d0
c**old code      if (bcigma.lt.1.0d-3) term1 = (0.4d0-c)*3.0d0
      if (term2.lt.0.0d0) term2 = 0.0d0
c     we now square the first term (linear in Dong and Bock) to
c     achieve a greater taper
      taper = term1**2 * term2

c     compute Q according to equation A-12 in Dong and Bock
      c = 0.0d0
      do 430 j = 1,50
         a1 = dble(j)
c         b1 = sngl((a1-dev)*s1)
c         b2 = sngl((a1+dev)*s1)
c         d1 = dble(erfc(b1)-erfc(b2))
         b1 = (a1-dev)*s1
         b2 = (a1+dev)*s1
c        limit the range of erf to  avoid underflows
         if (b1.lt.0d0 .or. b1.gt.15.d0 ) then
             erfcb1 = 0.d0
         else
             erfcb1 = erfc(b1)
         endif
         if (b2.lt.0d0 .or. b2.gt.15.d0 ) then
             erfcb2 = 0.d0
         else
             erfcb2 = erfc(b2)
         endif
c*         d1 = erfc(b1)-erfc(b2)
         d1 = erfcb1 - erfcb2
         c = c+d1
         if (d1.lt.trun) goto 440
 430  continue

c     return the decision function and reduced probability
 440  prob = 1.0d0-c
      if (c.lt.1.0d-9) c = 1.0d-9
      deci = taper/c  
      if( idebug.ge.3) print *,'csig term1 term2 taper prob deci '
     .       , csig,term1,term2,taper,prob,deci

 100  continue
      return
      end

c--------------------------------------------------------------
      function erf(x)
      real*8 erf,gammp,x,half
      half=0.5d0
      if(x.lt.0.)then
        erf = -gammp(half,x**2)
      else
        erf = gammp(half,x**2)
      endif
      return
      end
c--------------------------------------------------------------
      function erfc(x)
      real*8 erfc,x,gammp,gammq,half
      half=0.5d0
      if(x.lt.0.)then
        erfc = 1.+gammp(half,x**2)
      else
        erfc = gammq(half,x**2)
      endif
      return
      end
c--------------------------------------------------------------
      function gammp(a,x)
      real*8 gammp,a,x,gln,gammcf
c  rwk 070328: commented out because Solaris doesn't like 'pause'
c      if(x.lt.0..or.a.le.0.)pause
      if(x.lt.a+1.)then
        call gser(gammp,a,x,gln)
      else
        call gcf(gammcf,a,x,gln)
        gammp = 1.-gammcf
      endif
      return
      end
c--------------------------------------------------------------
      function gammq(a,x)
      real*8 gammq,a,x,gamser,gln   
c  rwk 070328: commented out because Solaris doesn't like 'pause'
c      if(x.lt.0..or.a.le.0.)pause
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq = 1.-gamser
      else
        call gcf(gammq,a,x,gln)
      endif
      return
      end
c--------------------------------------------------------------
      function gammln(xx)
      real*8 cof(6),stp,half,one,fpf,x,xx,tmp,ser,gammln
      integer*4 j
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     *    -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      x = xx-one
      tmp = x+fpf
      tmp = (x+half)*log(tmp)-tmp
      ser = one
      do 11 j = 1,6
        x = x+one
        ser = ser+cof(j)/x
11    continue
      gammln = tmp+log(stp*ser)
      return
      end
c--------------------------------------------------------------
      subroutine gcf(gammcf,a,x,gln)
      real*8 gammcf,gammln,a,anf,x,gln,gold,a0,a1,b0,b1,fac,an,ana,g,eps
      integer*4 n,itmax
      parameter (itmax = 100,eps = 3.e-7)
      gln = gammln(a)
      g = 0.
      gold = 0.
      a0 = 1.d0
      a1 = x
      b0 = 0.
      b1 = 1.d0
      fac = 1.d0
      do 11 n = 1,itmax
        an = float(n)
        ana = an-a
        a0 = (a1+a0*ana)*fac
        b0 = (b1+b0*ana)*fac
        anf = an*fac
        a1 = x*a0+anf*a1
        b1 = x*b0+anf*b1
        if(a1.ne.0.)then
          fac = 1.d0/a1
          g = b1*fac
          if(abs((g-gold)/g).lt.eps)go to 1
          gold = g
        endif
11    continue
ckf   pause 'a too large, itmax too small'
      call report_stat('FATAL','SOLVE','bdeci/gcf',' '
     .                ,'a too large, itmax too small',0)

1     gammcf = dexp(-x+a*dlog(x)-gln)*g
      return
      end
c--------------------------------------------------------------
      subroutine gser(gamser,a,x,gln)
      real*8 gamser,a,x,gln,gammln,ap,sum,del,eps
      integer*4 n,itmax
      parameter (itmax = 100,eps = 3.e-7)
      gln = gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)  call report_stat('FATAL','SOLVE','bdeci',' '
     .     ,'x < 0',0)
        gamser = 0.d0
        return
      endif
      ap = a
      sum = 1.d0/a
      del = sum
      do 11 n = 1,itmax
        ap = ap+1.d0
        del = del*x/ap
        sum = sum+del
        if(abs(del).lt.abs(sum)*eps)go to 1
11    continue
ckf   pause 'A too large, ITMAX too small'
      call report_stat('FATAL','SOLVE','bdeci/gser',' '
     .                ,'a too large, itmax too small',0)
1     gamser = sum*exp(-x+a*log(x)-gln)
      return
      end

