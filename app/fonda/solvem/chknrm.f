      subroutine chknrm(ncoef,a,scale)
c
      real*8 a,scale(ncoef),small,work
      dimension a((ncoef*(ncoef+1))/2)
      integer ncoef,i,j,idiag,ij
c
      small = 1.0d-12
c
c      write(6,750)
  750 format(' check the quality of covariance matrix')
c
c.....compute scale factors for the covariance matrix
      do 100 i=1,ncoef
         idiag=i*(i+1)/2
         if (a(idiag).lt.small) then
            scale(i)=1.d0
            print*,' ill diagonal at: ',i,a(idiag)
         else
            scale(i)=1.d0/dsqrt(a(idiag))
         endif
  100 continue
c.....scale the covariance matrix 
      ij=0
      do 600 i=1,ncoef
         do 500 j=1,i
            ij=ij+1
            work=a(ij)*scale(i)*scale(j)
            if (dabs(work).ge.1.0d0+small)
     .      write (6,'(a,2i5,f15.7)') ' error at: ',i,j,work
  500    continue
  600 continue
c
      return
      end
