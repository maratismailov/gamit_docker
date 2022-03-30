      subroutine secsum (iwkn1,sow1,delta,iwkn2,sow2)
c
c     add delta to epoch 1 to make epoch 2
c     can handle going over weeks
c
      integer*4   iwkn1,iwkn2
      real*8      sow1,sow2,delta

      real*8      spw

      data        spw/604800.d0/


      sow2 = sow1 + delta
      iwkn2 = iwkn1

 100  continue
      if (sow2 .ge. spw) then
         sow2 = sow2 - spw
         iwkn2 = iwkn2 + 1
         goto 100
      else if (sow2 .lt. 0.d0) then
         sow2 = sow2 + spw
         iwkn2 = iwkn2 - 1
         goto 100
      endif

      return
      end
