      Function IC4ARRAY( ctest,carray,n  )

C     Determine the position in the CARRAY of the token CTEST   
c     (assume case matched)

C        R. King 31 December 2003
         
      implicit none

      integer ic4array, n, i
      character*4 carray(*),ctest

      ic4array = 0
      do i=1,n
        if( ctest.eq.carray(i) ) then
           ic4array = i
        endif
      enddo
      return
      end

