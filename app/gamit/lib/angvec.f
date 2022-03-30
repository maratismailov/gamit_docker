      Function angvec(a,b,c)
c
c     Subroutine to compute the angle (radians) between two vectors measured
c     from A to B, where right-hand positive is defined by A x B aligning
c     with C (where C might be, e.g,. the direction of orbital angular momentum).  
c     In the diagram below, if C is the (right-hand) positive Z axis, and A and B
c     are each inclined 60 degrees wrt the X-axis, then the angle  A->B 
c     is 240 degrees.  R. King 980706
c     
c
c
c                                          
c           B     Y
c            .    |
c             .   |
c              .  |
c               . |
c       --------------------X              
c               . |
c              .  |
c             .   |
c            .    |
c           A
c 

      real*8 a(3),b(3),c(3),au(3),bu(3),sinab,cosab,xvec(3)
     .     , dot,amag3,angvec
      
c  Convert the inputs to unit vectors
      call normalise(a,au)
      call normalise(b,bu)

c  Get the cosine of the angle (will be properly signed)
      cosab = dot(au,bu) 

c  Get the sine of the angle
      call cross(au,bu,xvec)
      sinab = dsign( amag3(xvec), dot(xvec,c) ) 
      
c  Get the angle
      angvec = datan2(sinab,cosab)

      return
      end


