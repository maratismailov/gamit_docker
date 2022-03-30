      subroutine getusr(uname)

C     GET THE USER'S  LOGIN NAME

      character*(*) uname
      
      

c Getlog will work only in foreground
      call getlog(uname)
c Try getting the environment variable instead
c SCM check for nulls and blanks in uname. 10/29/03
      if( uname(1:3).eq.'   '.or.ichar(uname(1:1)).eq.0) then
          call getenv('USER',uname)
       endif

      return
      end
