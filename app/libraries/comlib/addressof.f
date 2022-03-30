CTITLE AddressOf
 
 
      integer*8 function addressOf( arg )

      implicit none

* MOD TAH 190511: Increased to I*8 for 64-bit machines
 
*     routine to use Sun function LOC to get the byte address
*     of the first byte in argumnet arg.  Arg is defined to be
*     any type.  Here we say it is interger*4
 
 
      integer*4 arg
 
      addressof = loc( arg )
 
***** Thats all
      return
      end
 
 
