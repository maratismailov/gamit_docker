      program tmod

      implicit none

      integer*4 i,j,k, n
      
      do i = 1,32
         j = mod(i,15)
         n = i + 1
         k = mod(i+1,15)
         print *,' MOD ',i,j, n, k
      enddo
      end
