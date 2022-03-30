c
      subroutine copy1i(istart,iend,ishift,irray1,irray2)
c
c copy 1-d integer array from irray1 to irray2


      implicit none
                      
      integer*4 iend,ishift,istart,irray1,irray2,i

      dimension irray1(iend),irray2(iend+ishift)
c
      do 20 i = istart,iend
         irray2(i+ishift) = irray1(i)
 20   continue
c
      return
      end
