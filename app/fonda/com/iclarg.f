      integer*4 function iclarg ( iel, arg )

c     return iel'th argument from the command line in arg
c     length of the string is returned in iclarg

c     This is the Sun version

c     equivalent to rcpar of HP1000 library
 
      integer*4 iel, nblen    
      character*(*) arg
 
      call getarg( iel, arg )
      iclarg = nblen( arg )
 
      return
      end
