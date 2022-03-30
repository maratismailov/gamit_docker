CTITLE BAD_OPTION
 
      subroutine bad_option(option, name)

      implicit none
 
c     Writes a message to the CRT that we could not understand
c     the option passed to a "RW" subroutine.
 
 
*   option          - OPtion or block name used
 
 
      character*(*) name, option
 
***** Write the message
      write(*,100) option, name
  100 format(' *** Bad option ',a,' passed to routine ',A)
 
      end
 
