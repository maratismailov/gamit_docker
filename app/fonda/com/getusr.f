      subroutine getusr(uname)

C     GET THE USER'S  LOGIN NAME

      character*(*) uname

      call getlog(uname)

      return
      end
