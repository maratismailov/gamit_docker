c
      Subroutine FXBIAS(insat,istat,itag,r2lc)
c
c     When a bias flag is encountered, solve for the biases implicitly 
c     and remove them from normal equations

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'    
      include 'parameters.h'

      integer*4 insat,istat,itag,ib,nsate,loop,lbias1,k,j
      integer*4 jk,i,jp,ifirst,lb1,iapp,n0,imax,imin,iaip
      integer*4 n1,ipp,n2,ill
      real*8 r2lc,alc21,facc,fac
      character*7 status
      logical adlc
       
       logical debug/.false./
C
C----------------------------------------------------
C METHOD:
C The simultaneous equations we need to solve, A|x> = |B>, written out:
C
C A11 X1 + A12 X2 + A13 X3 .....     =B1
C A21 X1 + A22 X2 + A23 X3 .....     =B2
C .........
C
C Just to make life easy, let X1 be the parameter we wish to implicitly
C  solve for and eliminate from the equations:
C
C X1 = (B1 - A12 X2 - A13 X3....)/A11
C
C If we plug this back in, the first line becomes B1 = B1 or all 0.
C  The second line:
C
C   (A22 -A21 A12/A11) X2  + (A23 - A21 A13/A11) X3 .   = B2-A21 B1/A11
C
C and further rows follow:
C
C Thus; if we wish to eliminate the kth parameter, all Aij, i,j .ne.k
C get replaced by:
C  Aij-Aik Akj/Akk
C
C  All Bj become Bj- Bk Ajk/Akk
C
C  and then all A's or B's with index k get zeroed.
C
C  What happens to chis**2 ? chi2 = r2sum + <x|A|x> - 2 <B|x>
C  After the parameters are adjusted it's r2sum-<B|x>
C
C  The difference in <B|x> before the parameter is eliminated and
C  after is B1**2/A11
C----------------------------------------------------
c                 
c     initialization for G77 compiler--logic not checked
      alc21 = 0.d0
      facc = 0.d0

c if nsate < 0 means do just that satellite
c    
      ib=1
      nsate=insat
      if(nsate.lt.0) then
         nsate=-nsate
         ib=nsate
      endif
       if( logprt.or.debug ) then
      if( itag.eq.10) write(*,889) iepoch,istat,ib,nsate
  889 format(' Update the bias parameter at epoch ',i5,' for site',
     *  i4,3x,'sat',i3,' to',i3)
       endif
c
      loop=1
      if (l2flag.ge.2) loop=2
      lbias1=nsite*nsat
      adlc = .false.
      if (l2flag.ge.1 ) adlc = .true.
c                                         
      if( debug ) print *,'FXBIAS loop adlc ib nsate '
     .                   ,loop,adlc,ib,nsate
      do 1500 k=ib,nsate
c
c     first find the bias params
      do 110 j=1,npartm(istat)
      jk=islot2(j,istat)
      i=islot1(jk)
c
c   skip if not bias parameter
c
c modify for multiple zenith delay parameters
cd      if(debug) print *,'FXBIAS chk slot jk i ',jk,i 
      if(i.le.2900 .or. i.gt.11900) go to 110
c
c  or bias flag
c 
      if(debug) print *,'FXBIAS itag ',itag  
      if(itag.eq.10) go to 112
c
110   continue
      if(debug) write(6,111) iepoch,istat,k,i
111   format(/,' Explicit bias at epoch ',i5,'for station',i3
     .   ,' sat',i3,i5)
      go to 1500
c------------------------------------------------------------
  112 continue
c
c solve the implicit bias parameters
c it's the first one we came to. so:
      jp=jk-1+k
c     write(*,*) 'jp',jp
      ifirst=jk-(istat-1)*nsat
c     write(*,*) 'ifirst',ifirst
      lb1=ifirst-1+lbias1
c     write(*,*) 'lb1,lbias1',lb1,lbias1
c           
      do 200 ill=1,loop
         if (ill.eq.2) then
            jp=jp+lbias1
            adlc = .false.
         endif
c
c we wish to eliminate the jp'th param from the normal equations
c
c for a(jp,jp):
         iapp=jp*(jp+1)/2
         n0=jp-ifirst+1   
cd          print *,'jp iapp a ',jp,iapp,a(iapp)

         if(dabs(a(iapp)).le.1.0d-6) then
            if (ill.eq.1) status='(  L1 )'
            if (ill.eq.2) status='(L2-L1)'
c           write (6,30) a(iapp),status,jp
c30         format ('The diagonal term is ',e8.1,'.',a7,'  row=',i4)
            go to 50
         endif
         do 300 i=1,ntpart
            imax=max0(i,jp)
            imin=min0(i,jp)
c for a(i,jp):
            iaip=imin+imax*(imax-1)/2
            fac=a(iaip)/a(iapp)
            if (adlc.and.i.le.lb1) then
               if (i.gt.lpart.and.i.lt.ifirst) goto 320
               if (i.le.lpart) then
                  n1=(n0-1)*lpart+i
                  alc21=clc(n1)
               endif
               if (i.ge.ifirst) then
                  n1=(imax-ifirst+1)*(imax-ifirst)/2+imin-ifirst+1
                  alc21=an22(n1)
               endif
               n1=n0*(n0+1)/2
c if the bias parameter happened at the beginning of effective
c double-difference combination, its diagonal term is zero.
               if (dabs(an22(n1)).lt.1.0d-8) then
                  facc=0.0d0
               else
                  facc=alc21/an22(n1)
               endif
               if(i.ne.jp) blc(i)=blc(i)-facc*blc(jp)
            endif
c
c eliminate:
 320        if(i.ne.jp) b(i)=b(i)-fac*b(jp)
c
            ipp=i*(i-1)/2
            do 310 j=1,i
            ipp=ipp+1
            if((j.eq.jp).or.(i.eq.jp)) go to 310
c a(i,j)=a(i,j)-fac*a(jp,j)
c       =a(i,j)- (a(i,jp)/a(jp,jp)) * a(jp,j)
            imax=(max0(j,jp))
            imin=(min0(j,jp))
c for a(jp,j):
            iaip=imin+imax*(imax-1)/2
            a(ipp)=a(ipp)-fac*a(iaip)
            if (i.gt.lpart.and.i.lt.ifirst) goto 310
            if (j.gt.lpart.and.j.lt.ifirst) goto 310
            if (adlc.and.i.le.lb1.and.dabs(facc).gt.1.0d-12) then
               if (i.le.lpart) then
                  n1=(n0-1)*lpart+j
                  alc(ipp)=alc(ipp)-facc*clc(n1)
               endif
               if (i.ge.ifirst) then
                  if (j.ge.ifirst) then
                     n1=i-ifirst+1
                     n2=n1*(n1-1)/2+j-ifirst+1
                     n1=(imax-ifirst+1)*(imax-ifirst)/2+imin-ifirst+1
                     an22(n2)=an22(n2)-facc*an22(n1)
                  endif
                  if (j.le.lpart) then
                     n1=(i-ifirst)*lpart+j
                     n2=(n0-1)*lpart+j
                     clc(n1)=clc(n1)-facc*clc(n2)
                  endif
               endif
             endif  

310         continue
300      continue
                    
c------------------------------------------------
c  now fix chisq             
        if(debug) print *,'FXBIAS jp old r2sum r2lc ',jp,r2sum,r2lc 
         r2sum=r2sum-b(jp)*b(jp)/a(iapp)
         if (adlc) then
            n1=n0*(n0+1)/2
            if (dabs(an22(n1)).lt.1.0d-8) then
               facc=0.0d0
            else
               facc=1.0d0/an22(n1)
            endif
            r2lc=r2lc-blc(jp)*blc(jp)*facc
         endif  
         if(debug) then   
           print *,'FXBIAS adlc r2sum r2lc ',adlc,r2sum,r2lc 
           print *,'FXBIAS jp n0 n1 an22 facc b blc r2sum r2lc '
     .               , jp,n0,n1,an22(n1),facc,b(jp),blc(jp),r2sum,r2lc  
         endif 
c
c  we left the jp'th rows and cols untouched.
c  now zero them
 50      do 309 j=1,ntpart
            imax=(max0(j,jp))
            imin=(min0(j,jp))
c for a(jp,j):
            iaip=imin+imax*(imax-1)/2
            a(iaip)=0.0d0
            if (adlc.and.j.le.lb1) then
               if (j.gt.lpart.and.j.lt.ifirst) goto 309
               if (j.le.lpart) then
                  n1=(n0-1)*lpart+j
                  clc(n1)=0.0d0
               else
                  n1=(imax-ifirst+1)*(imax-ifirst)/2+imin-ifirst+1
                  an22(n1)=0.0d0
               endif
            endif
309      continue
         b(jp)=0.0d0
         if (adlc) blc(jp)=0.0d0
c
 200  continue
c------------------------------------------------------------
1500  continue
c
      return
      end
