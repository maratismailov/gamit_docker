CTITLE    ..................................................................
 
      subroutine report(buffer)
c
c     Routine to report messages to the user lu
c
c Variables
c ---------
c buffer -- a character buffer to be written to the user's lu
c
      character*(*) buffer
 
c
c Functions used
c --------------
c trimlen -- HP utility
c loglu   -- HP utility
c
      integer*4 trimlen
 
c
c.... write out buffer
      if( trimlen(buffer).gt.0 ) then
         write(*,'(a)') buffer(1:trimlen(buffer))
      end if
c
      return
      end
 
