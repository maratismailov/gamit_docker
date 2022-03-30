      integer function get_path(name,pathname)

C     get the pathname of the specified "name"

      character*(*) name,pathname
      integer nblen

      call getenv(name,pathname)
      get_path = nblen(pathname)

      return
      end
