      subroutine indexi ( n, iarrin,indx )

c     indexes an integer array IARRIN of length N; ie outputs the array
c     INDX such that IARRIN(INDX(J)) is in ascending order for
c     J=1,2,...N.  The input quantities N and IARRIN are not changed
c     Routine from Press et al., Numerical Recipes
c     Added to GAMIT lib directory by R King  11 February 1997
            
      implicit none

      integer*4 n,iarrin(*),indx(*),i,j,l,iq,ir,indxt

      if( n.eq.1 ) then
        indx(1) = 1
        return
      endif

      do 11 j=1,n
         indx(j) = j
   11    continue
      l = n/2 + 1
      ir = n

   10 continue
        if( l.gt.1 ) then
           l= l-1
           indxt=indx(l)
           iq = iarrin(indxt)
        else
           indxt=indx(ir)
           iq= iarrin(indxt)
           indx(ir)=indx(1)
           ir=ir-1
           if(ir.eq.1) then
             indx(1)=indxt
             return
            endif
         endif
         i=l
         j=l+l

   20 if(j.le.ir) then
         if(j.lt.ir) then
             if(iarrin(indx(j)).lt.iarrin(indx(j+1))) j=j+1
         endif
         if(iq.lt.iarrin(indx(j))) then
             indx(i)=indx(j)
             i=j
             j=j+j
         else
             j=ir+1
         endif
       goto 20
       endif
       indx(i)=indxt
      goto 10

      end

