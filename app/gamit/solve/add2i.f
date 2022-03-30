c
      subroutine add2i(irow,icolmn,irow1,iclmn1,irray1,irray2)
c
c add 2-d integer array irray1 to irray2
c            
      implicit none
            
      integer irow,icolmn,irow1,iclmn1,irray1,irray2,i,j
      dimension irray1(irow,icolmn),irray2(irow,icolmn)
c
      do 20 j = 1,iclmn1
         do 10 i = 1,irow1
            irray2(i,j) = irray1(i,j) + irray2(i,j)
 10      continue
 20   continue
c
      return
      end
