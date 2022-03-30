CTITLE    ..................................................................
 
      subroutine get_cmd(buffer, commands, num, iel)
c
c     Routine to extract the command from a user supplied buffer.
c     The iel returned is the sequential number of the command
c     in the 'commands' string
c
c Variables
c ---------
c buffer -- user supplied buffer with complete command plus variables
c     in it.
c commands -- a list of available commands given in PLOBD
c num      -- number of commands in 'commands'
c iel    -- the position of command from buffer in command list
c     Also used to return error messages if negative
c
c iblank  -- the position of the last character of the command
c new_buffer -- a copy of the buffer for working
c blank -- a blank string used to pad commands
* ierr  -- iostat error
* i,j,k -- Loop counters
c
      character*(*) buffer, commands(1)
 
c
      character*80 new_buffer
 
c
      character*8 blank
 
c
      integer*4 num, iel, iblank, ierr, i
 
c
c.... Initialize blank
      blank = '        '
c
c.... Remove the leading blanks from buffer
      call trim_lead(buffer,ierr)
c
c.... Make sure no error ocurred
*                             ! there was an error in trimming the string
      if( ierr.lt.0 ) then
c
         iel = ierr
c
*                             ! try to extract command
      else
c
c
c....    Now find the first blank
         iblank = index(buffer,' ')  - 1
c
c....    add blanks to string
         if( iblank.lt.8 ) then
            new_buffer = buffer(1:iblank)//blank(1:9-iblank)//
     .         buffer(iblank+1:)
            buffer  = new_buffer
         end if
c
c....    Now convert command part to upper case
         call casefold(buffer(1:8))
c
c....    loop over available commands
         iel = -1
         do i = 1, num
            if( index(buffer,commands(i)(1:iblank)).eq.1 ) then
c
c....          We found the command (see if any before have been found
*                                     ! ambiguous command
               if( iel.ne. -1 ) then
                  iel = -4
*                                     ! OK first occurence of command
               else
                  iel = i
               end if
            end if
         end do
c
      end if
c
      return
      end
 
