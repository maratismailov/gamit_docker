      subroutine trans(l,m,a,at)
c	
c     transpose matrix
c
      implicit real*8(a-h,o-z)
      integer l,m,il,i1,i2,im
  
      dimension a(l*m),at(l*m)
  
      do 30 il = 1,l
         do 20 im = 1,m
            i1 = (il-1)*m+im
            i2 = (im-1)*l+il
            at(i2) = a(i1)
 20      continue
 30   continue

      return
      end
