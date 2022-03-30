      Subroutine BFWORK(mopt,issn,k1,k2,ld1,ld2,llbad,
     .   ia1,ia2,l1b,i,k,away,job)
                    
c     Calculate indices for bias-fixing.  D. Dong March 1991
                       
c      job = 1: initialization
c      job = 2: get nessesary indexes
c      job = 3: check if the bias belongs to the fixing catalogue

c      mopt = 1 : L1 (or L2) biases only
c           = 3 : Both L2-L1 and L1 biases (not used?)

c      issn   : Number of L1 biases 

c      k1 : Number of live non-bias parameters
c      k2 : Number of live bias parameters    

c      llbad : Index (or counter) for biases corresponding to bad stations; this
c              code may not work and/or is bypassed.  See comments in LSQUAR. 
c      ld1, ld2: Indices associated with llbad and lbfre.

c      ia1 : Always = 1 in single-session mode
c      ia2 : Number of total biases (=1 for single-frequency, =2 for dual-freq)
c      l1b : Number of L1 biases (= common/bbi/ l1bias)
          
c      i  : Bias number (in common/bbii/ idxb array), dummy w/ ijob =1,2 

c      away = false : bias available for fixing
c           = true  : bias not available for fixing (negative in idxb array)      


      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      logical away 

      integer mopt,issn,k1,k2,ld1,ld2,llbad,ia1,ia2,l1b,job
     .      , lim1,ih,i,k

c      real*8 bl,bl1,bl2
       real*8 adjwl  

c      Function to deterimine is the WL for an associated NL has been resolved
       logical wl_fixed

      if (job.eq.1) then  
c         print *,'BFWORK msig ',msig
         k1=msig-1
         issn=0
         k2=0
         ld1=0
         ld2=0
         goto 100
      endif

      if (job.eq.2) then   
         l1b=l1bias
         llbad=lbfre(1)
         k1=k1+k2
         ia1=issn+1
         ia2=issn+l1b*iband
         issn=issn+l1b
         k2=0
         if (llbad.gt.0) then
           if (mopt.eq.1.or.mopt.eq.3) then
               ld1=lbfre(2)
               ld2=lbfre(3)
           endif
         endif
         goto 100
      endif

      if (job.eq.3) then   
c         print *,'BFWORK 3, i idxb issn llbad ld1 ld2 iband '
c     .                    , i,idxb(i),issn,llbad,ld1,ld2,iband
         away=.false.
         k=idxb(i)
         if(k.le.0) go to 2
         k2=k2+1   
         if (mopt.eq.1.and.i.gt.issn) goto 2
         if (mopt.eq.2.and.i.le.issn) goto 2
c----    force all continental scale biases free
         lim1=ia1+limitb-1         
c         print *,'  ia1 limitb lim1 '
c     .             ,ia1,limitb,lim1   
         if(i.le.issn.and.i.ge.lim1) goto 2
         if(i.gt.issn.and.i.ge.lim1+l1b) goto 2
c----    keep biases connected with bad stations always free
         if (llbad.gt.0) then
               if (k.ge.ld1.and.k.le.ld2) goto 2
         endif
c----    in lc after l1,l2 mode, keep l1 bias (which correspond
c           l2-l1 bias is unfixed) free.
         if (mopt.eq.1.and.iband.eq.2) then   
c** rwk 070412: replace this code with a function call, housed here but
c               also used by get_bias_scale.
c            i1=i+l1bias
c            ih=half(i1)
c            bl1=adjust(k+l1bias
c            bl=bl1*dble(ih)
c            bl2=dint(bl+.5d0*dsign(1.d0,bl))/dble(ih)
c            bl2=dabs(bl1-bl2)       
cc            print *,  'i1 ih bl1 bl2 ',i1,ih,bl1,bl2
cc **        CHANGE this to a smaller tolerance
c            if (bl2.ge.0.01) goto 2   
c             if (bl2.ge.0.0001) goto 2            
           ih = half(i + l1bias)
           adjwl = adjust(k+l1bias)   
           if( .not.wl_fixed(ih,adjwl) ) goto 2
         endif
         goto 100
      endif

 2    away=.true.
 100  continue

      return
      end   



