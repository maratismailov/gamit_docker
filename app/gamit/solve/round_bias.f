      Subroutine ROUND_BIAS( lane,mopt )   

c     Resolve ambiguities using the rounding/decision-function algorithm from
c     Dong and Bock [1989]. Called by RSOLVE_WL and RESOLVE_NL. 
c     R. King 9 January 2004

c     This code moved from, and supersedes the code in NBIASR, called in the 
c     original scheme by NBIAS.

c      lane    : 1=NL  2= WL   

c      mopt    : 1= LC_HELP or LC_AUTCLN for NL, L1_ONLY, L2_ONLY 
c                2= LC_HELP or LC_AUTCLN for WL, L1&L2 with ion constraint (this may not work)
c                3= L1L2 independent

      implicit none  
                      

c-- Global commons

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

c-- Calling arguments

      integer lane
c        = 1 NL   =2 WL

      integer mopt  

c-- Other commons

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

c     -- in parameters.h: parameter values, slots, and uncertainties
c      real*8 preval(maxprm),sigma_com(maxprm)
c      integer*4 idms(maxprm),islot1(maxprm),free_com(maxprm)
c      character*20 rlabel(maxprm)
c      common /params/preval,sigma_com,idms,islot1,free_com,rlabel  
*      MOD TAH 031218: Added the params common so that we output and check
*      label on the bias parameters.  Need to change some names because
*      some of these variables are passed to the routine.  (_com added
*      to name to denote common variable)


c-- Local variables

      real*8 dmax(15),pright(15),xintp
     .     , deci,prob
     .     , cutdec,xint,cutdev,cutsig
     .     , add,dev,sma,awl,est,awv,a1,a2
*        MOD TAH 031203: Added sma_save (This is the sigma of the fixed bias)
      real*8 sma_save

      integer ipseu,ih,ia1,ia2,ifix,icut,ld1,ld2,llbad,l1b
     .      , issn,junk(15),junk1(15),junk2(15),idum
     .      , j1,j2,j3,k1,k2,i,j,k
 
      logical away,ldum
           

c** Flag for debug printout
      integer*4 idebug
c      idebug = 0  : none 
c               1  : print to q- and log file decision values for fixed and not-fixed biases
c               2  : print to log the fixing criteria
c               3  : full print of logic 
       idebug = 1
                 

c-- Initialization

      icut=1
      ifix=0
      prob=1.0d0
      ipseu  =  0             
c     the following necessary to avoid compiler warnings
      cutdec = 0.d0   
      call zero1d(1,icut+1,dmax) 
      idum = 0      
      call bfwork(mopt,issn,k1,k2,ld1,ld2,llbad,ia1,ia2
     .           , l1b,idum,idum,ldum,1)   
      if(idebug.ge.3 ) then
        print *,'ROUND_BIAS initial call to BFWORK '
         print *,' mopt,issn k1 k2 ld1 ld2  '
     .       ,  mopt,issn,k1,k2,ld1,ld2
       endif


c-- Set parameters for wide-lane or narrow-lane

      if (lane.eq.1 ) then
         cutdev = nldev
         cutsig = nlsig
         cutdec = nlcut   
      else if (lane.eq.2 ) then
         cutdev = wldev
         cutsig = wlsig
         cutdec = wlcut
      else
         call report_stat('FATAL','SOLVE','round_bias',' '
     .                      ,'Bad value of lane',0)
      endif                
      if (cutdec.lt.1.0d0) cutdec = 1000.0d0    


c---- Get the bias indices                       

       call bfwork(mopt,issn,k1,k2,ld1,ld2,llbad,ia1,ia2,
     .        l1b,i,k,away,2)     
       if( idebug.ge.3 ) then
         print *,'ROUND_BIAS call BFWORK to get indices '
         print *,' mopt issn k1 k2 ld1 ld2 llbad ia1 ia2 l1b away '
     .       ,  mopt,issn,k1,k2,ld1,ld2,llbad,ia1,ia2,l1b,away     
        endif
           
c---  Find the highest value of the decision function among the currently free biases
                 
         do 10 i=ia1,ia2   
c            print *,' Start loop i = ',i
            if (i.gt.issn) ipseu = ipseu + 1
            call bfwork(mopt,issn,k1,k2,ld1,ld2,llbad,ia1,ia2,
     .           l1b,i,k,away,3) 
            if(idebug.ge.3) then
              print *,'ROUND_BIAS bias to BFWORK '
              print *,' i k k1 k2 ld1 ld2 llbad ia1 ia2 l1b away '
     .       ,  i,k,k1,k2,ld1,ld2,llbad,ia1,ia2,l1b,away
              print *,' aft BFWORK ipseu away ',ipseu,away
             endif
            if (away) goto 10
            ih=half(i)   
            est=adjust(k)
            sma = bscale(i)*sigma(k1+k2)       
            if( idebug.ge.2 ) then
               print *,'Calling BDECI i k k1k2 est sma cutdev cutsig '
     .                ,               i,k,k1+k2,est,sma,cutdev,cutsig
            endif
            call bdeci( est,sma,ih,cutdev,cutsig,prob,deci )    
            if( idebug.ge.2 ) then
             write (6,105) k,rlabel(k), est,sma,deci 
 105         format('Decision for bias No.',i5,2x,a20,2x,f9.2,' +- '
     .            ,f9.2,' Decision Function ',d12.4)
            endif   
            if (deci.le.cutdec) goto 10     
            if (deci.le.dmax(icut)) goto 10
            if (ifix.lt.icut) ifix=ifix+1
            do 15 j=1,ifix    
              if(idebug.ge.2) print *,' In fix loop j icut ',j,icut 
c             save the values only if deci is now the maximum value
              if (deci.le.dmax(j)) goto 15  
c             logic here clumsy; since ifix always = 1, this loop never executed
              if (ifix.ne.1) then
                do  j1=j,ifix
                  j3=ifix-j1+j
                  j2=j3+1
                  dmax(j2)=dmax(j3)
                  nfix(j2)=nfix(j3)
                  junk(j2)=junk(j3)
                  junk1(j2)=junk1(j3)
                  junk2(j2)=junk2(j3)
                  pright(j2)=pright(j3) 
                enddo 
              endif
              dmax(j)=deci
              nfix(j)=k1+k2  
              junk(j)=k   
              junk1(j)=i
              junk2(j) = ipseu   
              pright(j)=prob
*             MOD TAH 031203: Code above is executed for the highest decision 
*             value found so far.  Save the sigma of the bias at this point
              sma_save = sma
              goto 10
   15       continue 
c           ---end of loop over j = 1,ifix (only executed once so far as we can tell)
   10    continue
c        ---end of loop over biases in parameter list

         if (iband.eq.2) issn=issn+l1bias


c-- Reorder array by junk1 to fit BNEW (dmax is working array only)   

c     (never execulted since ifix always = 1)  
      if(idebug.ge.2 ) print *,'   End of deci loop ifix ',ifix
      if (ifix.gt.1)  call bsort( ifix,pright,junk1,dmax,junk
     .                          , nfix,junk2,2)


c-- Attempt to fix the bias (if WL, check fixed biases PR est)

c     (in original scheme NBIASP called later by NBIAS to use pseudoranges with 
c      biases unresolved by the decision function here; we may want to restore that option)  

c     (this loop is awkward in that it combines two different functions)       

c     (In the current code, ifix is always =1, and junk(i) is the global
c      parameter index for the bias with the largest decision function, as
c      a result of the logic of the last loop)
      
      do 30 i=1,ifix
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
c        if NL, set a flag for LCLOOS  to indicate that this was a resolved
c        bias and not a dependent one
         if( mopt.eq.1 ) then 
           nlres(k1) = k
c        mopt=2 means searching wide-lane (why test here?)
         elseif (mopt.eq.2) then
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
               if( logprt ) write (6,120) est,awl
               write (10,120) est,awl
c              I am bypassing this since I found two cases where wrong widelane 
c              was fixed according to PR when it was obvious that phase estimate 
c              was correct -- YB 94/7/2
c              This code now restored but can be bypassed by inputting small 
c              values of the PR deviation and sigma -- RWK 94/11/16
c              compare the two estimates
                a1 = bdev(i)
                a2 = xintp-awl
c              see if TWICE the PR deviation is less than the phase deviation
c              (and make sure the two are within 0.5 cycles to guard against
c               a bad PR estimate)
                if (dabs(a2*2.0d0).lt.dabs(a1).and.(dev.lt.0.5d0)) then
                   bdev(i) = xintp-est
                   if( logprt ) write(6,130) xintp,adjust(k)
                   if( logprt ) write(6,130) xintp,adjust(k)
                   adjust(k) = xintp
                endif
            endif
         endif    
*        MOD TAH 031202: Added additional output to better see what is happening.
*        It seems in the above code that ifix is never more than 1 for the rounding code.
         if( logprt.or. idebug.ge.2 ) 
     .      write (6,110) k,rlabel(k), est,sma_save, adjust(k), dmax(1) 
         write (10,110) k,rlabel(k),est,sma_save, adjust(k), dmax(1)
  30   continue      
c--- end of loop to fix bias (actually executed only once in current code)

c** debug print when bias not fixed
 110  FORMAT ('Fix No.',I5,1x,a20,' bias  from',F9.2,' +- ',F9.2,
     .        ' to',F8.1,' Decision Function ',E12.4)
 120  format (2x,'Warning: big difference between 2 WLs (phase ',
     .   f11.2, '   PR ',f11.2,')')
 130  format (2x,'Warning: PR WL (',f11.2,') overriding phase est ('
     .   ,f11.2,')')

       return
       end

