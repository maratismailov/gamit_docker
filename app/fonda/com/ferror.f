      subroutine ferror(ios,lu)
c     print the translation for error code

c     kurt feigl for sun

      integer ios,lu
      character*127 string
                          
      call gerror(string)
      write (lu,*) string

      return
      end

