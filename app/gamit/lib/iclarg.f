      integer*4 function iclarg ( iel, arg )

c     return iel'th argument from the command line in arg
c     length of the string is returned in iclarg

c     This is the Sun version

c     equivalent to rcpar of HP1000 library

      integer*4 iel, rcpar
      character*(*) arg

c     call getarg( iel, arg )
c     iclarg = nblen( arg )
* MOD TAH 010610: Use rcpar which has been modfied to
*     account for different starting argument numbers
      iclarg = rcpar( iel, arg)
      
      return
      end
