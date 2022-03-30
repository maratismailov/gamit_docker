
CTITLE rcpar4
 
      integer*4 function rcpar4( iel, arg )

      implicit none

 
*     Routine to emulate rcpar4 using the igetarg UNIX subroutine
*     Modified to use getarg
 
*         iel       - Element of runstring to get
*       igetarg     - UNIX routine to read runstring
*       len_arg     - Length of arg.
*       trimlen     - Get length of string
 
      integer*4 iel, len_arg, trimlen
 
*             arg   - Arg of runstring
 
      character*(*) arg
 
****  Get length of argument and runstring
 
      len_arg = LEN(arg)
      call getarg( iel, arg )
      rcpar4 = trimlen( arg )
 
***** Thats all
      return
      end
 
