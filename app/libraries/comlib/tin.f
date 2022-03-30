      program int

      implicit none

*     Test if inkey works 

      integer*2 inkey, ikey
      integer*4 arg
      
      integer*4 i,j

      do i = 1, 100
         j = inkey(ikey)
         if( ikey.gt.0 ) then
         write(*,100) i,j, ikey, ikey, arg
  100    format(3i7, a2, o8)
         end if
      end do

      end
