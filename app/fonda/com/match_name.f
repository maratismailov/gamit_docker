c-----------------------------------------------------------------
      integer function match_name(n,m,name,code)
c
      integer n,m,i
      character*(*) name(n),code
c
      match_name = 0
      do i = 1,n
         if (code(1:m).eq.name(i)(1:m)) then
            match_name = i
            goto 50
         endif
      enddo
c
 50   continue
      return
      end
