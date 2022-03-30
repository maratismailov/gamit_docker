CTITLE 'DWVMV'
 
      subroutine dwvmv( from, incf, to, inct, iter)
c
c
*     Routine to copy one vector to another one
c
*          if, it       - Indexes in from and to vectors
*   i            - Loop counte
*   incf         - Increment for from
*   inct         - Increment to
*   iter         - Number of iterations
      integer*4 if, it, i, incf, inct, iter
*
c
C     real*8
c    .    from(incf,1)      ! the from vector
c    .,   to(inct,1)        ! the two vector
 
      real*8 from(1), to(1)
*
c
c
***** Copy
      if = 1
      it = 1
      do i = 1, iter
         to(it) = from(if)
c        to(1,i) = from(1,i)
         it = it + inct
         if = if + incf
      end do
c
c     write(*,*) if,it, iter, incf,inct, to(1), from(1)
c
***** Thats all
      return
      END
 
