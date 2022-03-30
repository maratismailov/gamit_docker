      subroutine nrmscl(a,b,scale,ncoef,iflag,iflag2) 
c
      real*8 a,b,scale(ncoef)
      integer iflag,iflag2,ncoef,i,idiag,ij,j
      dimension a((ncoef*(ncoef+1))/2),b(ncoef)
c
c.....
      go to (10,20), iflag
c
   10 continue
c      write(6,750)
  750 format(' normalization in progress')
c
c normalize normal equations
c.....compute scale factors for the normal equations
      do 100 i=1,ncoef
      idiag=i*(i+1)/2
      if (a(idiag).gt.0.0d0) go to 50
      print*,' ill-diagonal term at: ',i,idiag,a(idiag)
      scale(i)=1.d0
      go to 100
   50 scale(i)=1.d0/dsqrt(a(idiag))
  100 continue
c.....scale the coefficient matrix and right sides of the
c.....normal equations
   20 continue
c      write(6,*) iflag
      ij=0
      do 600 i=1,ncoef
      do 500 j=1,i
      ij=ij+1
      a(ij)=a(ij)*scale(i)*scale(j)
c      write (6,1000) scale(i),a(ij)
c      write (*,*) scale(i),'a',a(ij)
 1000 format (1x,2(e22.15,1x))
  500 continue
      if(iflag2.eq.1) b(i)=b(i)*scale(i)
c      write (6,1000) scale(i),b(i)
c      write (*,*) scale(i),'b',b(i)
  600 continue
c
      return
      end
