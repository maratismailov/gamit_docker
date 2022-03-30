c      Program to compute a cummulative chi distribution to compare with histograms
c      of velocity residuals.  Follows formula given by Press et al. in Numerical
c      Recipes, p. 537, with check against the table on p. 536.
c      R. King 12 February 2004

          
c The calling sequence is

c    chicurve [ndf] [n] 

c      where [ndf] is the degrees of freedom (usually 1 or 2)
c            [n]  is the number of points to calculate (default 100)

c      The range is hardwired to 0-100%, 0- to 3-sigma   
c      The output is to the screen

      implicit none

      integer iarg,iclarg,numpts,ndf,i

      real a,x,d,chi,p,gammq                 
     
      character*8 arg
  
c  Get the run-string arguements

       iarg = iclarg(1,arg) 
       if( iarg.eq.0 ) then
         print *,'Missing ndf in calling arguments'
         stop
       else
         read(arg,'(i1)') ndf  
         a = float(ndf)/2.
       endif

       iarg = iclarg(2,arg)
       if( iarg.eq.0 ) then
         numpts = 100
       else
          read(arg,'(i8.0)') numpts
       endif

c Compute the values

c    p = 1 - gammq(a,x)

c       where p = probability
c             a = df/2.
c             x = dchi2/2.
c             gammq is the incomplete gamma function (Num. Rec.p. 162)


c  Get the increment from the range and number of pts.

      d = 3.0/float(numpts)      
             
c  Compute the values for plotting

       chi = 0. 
       do i=1,numpts    
         chi = chi + d     
         x = chi**2/2.
         p = 1 - gammq(a,x) 
         p = p*100. 
         write(*,'(2f9.3)') chi,p
       enddo
       stop
       end
     
c---------------------------------------------------------------------------

      Function GAMMQ(a,x)  

c      From Numerical recipes, p. 162  the incomplete Gamma function 

      implicit none
      real a,x,gamser,gln,gammcf,gammq

      if( x.lt.0 .or .a.le.0 ) then
          print *,'Problem x a ',x,a
          stop
       endif
      if( x.lt.a+1. ) then 
         call gser(gamser,a,x,gln)
         gammq = 1. - gamser        
      else
         call gcf(gammcf,a,x,gln)  
         gammq = gammcf
      endif
      return
      end

      Subroutine GSER(gamser,a,x,gln)  

c      From Numerical Recipes, p. 162, the incomplete gamma function Q(a,x)
c      evaluated by its series representation.  Also returns ln (gamma(a) as GLN).

      implicit none
      integer itmax,n
      real eps,gln,a,x,gamser,ap,sum,del,gammln

      parameter(itmax=100, eps=3.e-7)
      gln = gammln(a) 
      if( x.le.0. ) then
        if( x.lt.0.) then
           print *,'Problem in GSER x < 0 ',x
           stop
         endif
         gamser = 0.
         return
      endif
      ap = a
      sum = 1./a
      del = sum  
      do 11 n=1,itmax
        ap = ap + 1.
        del = del*x/ap
        sum = sum + del
        if( abs(del).lt.abs(sum)*eps ) goto 1
   11 continue
      print *,'GSR A too large, ITMAX too small ',a,itmax
      stop
    1 gamser = sum*exp(-x+a*log(x)-gln)
      return
      end 

c-----------------------------------------------------------
   
      Subroutine GCF(gammcf,a,x,gln)          

c      From Numerical Recipes, p. 162, the incomplete gamma function Q(a,x) 
c      evaluated by it scontinued fraxction representation as GAMMCF.
c      Also returns gamma(a) as gln.

      implicit none
      integer itmax,n
      real eps,gammcf,gammln,a,x,gln,gold,a0,a1,b0,b1,fac,an,ana,anf,g
   
      parameter (itmax=100,eps=3.e-7) 
      gln = gammln(a) 
      gold = 0.
      a0 = 1.
      a1 = x
      b0 = 0.
      b1 = 1.
      fac = 1.
      do 11 n=1,itmax
        an = float(n)
        ana = an -a 
        a0 = (a1+a0*ana)*fac
        b0 = (b1+b0*ana)*fac
        anf = an*fac
        a1 =  x*a0 + anf*a1
        b1 = x*b0 + anf*b1
        if( a1.ne.0. ) then
          fac = 1./a1
          g = b1*fac
          if( abs((g-gold)/g).lt.eps) go to 1  
          gold = g
        endif
   11 continue
      print *,'A too large, ITMAX too small '
      stop
    1 gammcf = exp(-x+a*alog(x)-gln)*g
      return
      end 

c-------------------------------------------------------------------
            
      Function GAMMLN(xx)

c      From Numerical Recipes, returns the value ln (gamma(xx))).

      implicit none

      integer j

      real gammln,xx
      real*8 cof(6),stp,half,one,fpf,x,tmp,ser

      data cof/76.18009173d0,-86.50532033d0,24.01409822d0,
     .     -1.231739516d0,0.120858003d-2,-0.536382d-5/ 
      data stp/2.50662827465d0/ 
      data half,one,fpf/0.5d0,1.0d0,5.5d0/  

      x = xx-one
      tmp = x + fpf
      tmp = (x+half)*log(tmp) - tmp
      ser = one    
      do j=1,6
        x = x + one
        ser = ser + cof(j)/x
      enddo   
      gammln = tmp + dlog(stp*ser)
      return
      end


      
