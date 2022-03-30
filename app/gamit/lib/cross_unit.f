      Subroutine cross_unit(a,b,c)

c Compute the unit vector C orthogonal to the vectors A and B such
c that C = A x B  
c  --R. King 98/7/2

      implicit none 
         
      integer i
      real*8 a,b,c,rmag,dot
      dimension a(3),b(3),c(3)
                       
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)

      rmag = dsqrt( dot(c,c) )
      do i = 1, 3
       c(i) = c(i)/rmag
      enddo
      
      return
      end
