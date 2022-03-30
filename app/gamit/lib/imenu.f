      subroutine imenu (i,n)

c     handle a list of integer choices

c     index
      integer i

c     max allowed choice
      integer n

      integer ncount

      ncount = 0
 10   continue
      ncount = ncount+1
      if (ncount .lt. 5) then
         write (6,'(1x,a,1x)') 'Pick a number.'
c        the double read below is a kluge to avoid a pointer problem
c        created by a previous call to pickfn
         read (5,*,err=10) i

         if (i .le. 0) then
            write (6,'(1x,a,1x,$)') 'Too small.'
            goto 10
         else if (i .gt. n) then
            write (6,'(1x,a,1x,$)') 'Too large.'
            goto 10
         else
            return
         endif
      endif
      return
      end





