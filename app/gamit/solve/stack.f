c
      subroutine stack(istack,isame,ipflg,ieflg)
c
c     determine if the data should be stacked
c
      logical istack
      integer isame,ipflg,ieflg

c
c     reorder normal matrix and copy it to unit 27
c
      if (isame.gt.0) then
         if (ipflg.eq.0.and.ieflg.eq.0) then
            isame = isame+1
         else
            istack = .false.
         endif
         goto 10
      endif

      if (ipflg.gt.0) isame = 1

 10   continue
c
      return
      end


