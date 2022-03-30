      subroutine ata(l,n,at,atxa)
c	
c     matrix multiplication (A~)x(A)
c
      implicit real*8(a-h,o-z)
      integer l,n,i1,i2,ii,in
  
      dimension at(l*n),atxa(l*(l+1)/2)
  
      do 30 i1 = 1,l
         do 20 i2 = i1,l
            ii = i2*(i2-1)/2+i1
            temp = 0.0d0
            do 10 in = 1,n
            temp = temp+at((i1-1)*n+in)*at((i2-1)*n+in)
 10         continue
            atxa(ii) = temp
 20      continue
 30   continue

      return
      end
