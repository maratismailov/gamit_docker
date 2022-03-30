CTITLE    ..................................................................
 
      subroutine command_error(buffer,ierr)
c
c     Routine to print out error meassage when command can not be
c     interpretted
c
c Variables
c --------
c buffer -- a character buffer to be written to the user's lu
c ierr   -- the error associated with the command error
c
c iterm  -- the users logical unit number (found by call to loglu)
c
      character*(*) buffer
 
c
      integer*4 iterm, ierr
 
c
c Functions used
c --------------
c trimlen -- HP utility
c loglu   -- HP utility
c
      integer*4 trimlen
 
c
      iterm = 6
c.... write out buffer
      if( ierr.eq.-1 ) then
         write(iterm,110) buffer(1:trimlen(buffer))
 110     format(" Unknown command: ",a)
         return
      end if
c
c     ierr = -2 is empty string and no error message is printed
c
      if( ierr.eq.-3 ) then
         write(iterm,130) buffer(1:trimlen(buffer))
 130     format(" String too long in get_cmd: ",a)
         return
      end if
c
      if( ierr.eq.-4 ) then
         write(iterm,140) buffer(1:trimlen(buffer))
 140     format(" Ambiguous command: ",a)
         return
      end if
c
      return
      end
 
