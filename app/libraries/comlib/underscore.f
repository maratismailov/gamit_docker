CTITLE   ............................................................
 
      subroutine underscore(char_list,n)

      implicit none
c
c     Routine to replace blanks (but not trailing blanks) with
c     an underscore '_'
c
c
      integer*4 n
 
c
      character*(*) char_list(n)
 
c
      integer*4 TrimLen, i, length, index
 
c
c.... Loop over character array
      do i = 1, n
c
c....   Get length of string, not including trailing blanks
        length = TrimLen(char_list(i))
c
c....   If nonzero length, loop over characters
        if (length .gt. 0) then
c
c....     Loop over characters
          do index = 1, length
c
c....       Check for blanks
            if (char_list(i)(index:index) .eq. ' ')
     .          char_list(i)(index:index) =    '_'
c
          end do
c
        end if
c
      end do
c
      end
 
