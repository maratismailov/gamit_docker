      real*8 a(3),b(3),c(3),ang

      a(1) = 1.
      a(2) = 2.
      a(3) = 0.
      b(1) = -1.
      b(2) = 4.
      b(3) = 0.
      c(1) = 0.
      c(2) = 0.
      c(3) = 1.
      call angvec(a,b,c,ang)
      print *,'a,b,c,ang ',a,b,c,ang

      a(1) = 1.
      a(2) = 2.
      a(3) = 0.
      b(1) = -3.
      b(2) = .3
      b(3) = 0.
      c(1) = 0.
      c(2) = 0.
      c(3) = 1.
      call angvec(a,b,c,ang) 
      print *,'a,b,c,ang ',a,b,c,ang

      a(1) = 1.
      a(2) = 2.
      a(3) = 0.
      b(1) = -1.
      b(2) = -3.
      b(3) = 0.
      c(1) = 0.
      c(2) = 0.
      c(3) = 1.
      call angvec(a,b,c,ang) 
      print *,'a,b,c,ang ',a,b,c,ang

      a(1) = 1.
      a(2) = 2.
      a(3) = 0.
      b(1) = 1.
      b(2) = .3
      b(3) = 0.
      c(1) = 0.
      c(2) = 0.
      c(3) = 1.
      call angvec(a,b,c,ang)
      print *,'a,b,c,ang ',a,b,c,ang

      stop
      end
               
   
      Subroutine angvec(a,b,c,ang)
c
c     Subroutine to compute the angle (radians) between two vectors measured
c     from A to B, where right-hand positive is defined by A x B aligning
c     with C (e.g. the direction of orbital angular momentum).  In the dia-
c     gram below, if C is the (right-hand) positive Z axis, and A and B
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
     .     , dot,amag3,ang
      
c  Convert the inputs to unit vectors
      call normalise(a,au)
      call normalise(b,bu)

c  Get the cosine of the angle (will be properly signed)
      cosab = dot(au,bu) 

c  Get the sine of the angle
      call cross(au,bu,xvec)
      sinab = dsign( amag3(xvec), dot(xvec,c) ) 
      
c  Get the angle
      ang = datan2(sinab,cosab)
      print *,'au,bu,sinab cosab ang ',au,bu,sinab,cosab,ang

      return
      end

      SUBROUTINE CROSS(A,B,C)
C
C COMPUTE CROSS PRODUCT OF TWO VECTORS A AND B WITH RESULTS IN C
C S. A. GOUREVITCH	6/81
C
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 a,b,c
      DIMENSION A(3),B(3),C(3)
C
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
C
      RETURN
      END

    
      FUNCTION DOT(A,B)
C
C COMPUTE DOT PRODUCT OF VECTORS A AND B
C S. A. GOUREVITCH	6/81
C
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 a,b,dot
      DIMENSION A(3),B(3)
      DOT=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
      RETURN
      END

       subroutine normalise (vect,norm_vect)

c
c  Purpose: to normalise a 3 x 1 vector (vect) into norm_vect
c
c  IN: vect - vector to be normalised   R*8 (3)
c
c OUT: norm_vect - normalised vector    R*8 (3)
c
c SUBROUTINES CALLED: amag3
c
c CREATED: 3rd April, 1995                     LAST MODIFIED: 3rd April, 1995
c
c AUTHOR: P. Tregoning
c
c COPYRIGHT: DEPARTMENT OF EARTH AND PLANETRY SCIENCES
c            M.I.T. 1995
c
       implicit none
c
       real*8 vect(3),norm_vect(3),len,amag3
       integer i

       len = amag3(vect)
       do 10 i=1,3
         norm_vect(i) = vect(i)/len
10     continue

       return
       end


      function amag3(a)
c
c computes the magnitude of a vector using a dot product
c P Tregoning 3/95
c
      implicit none
      real*8 amag3,a(3)

      amag3 = dsqrt(a(1)**2 +a(2)**2 + a(3)**2)
      return
      end

