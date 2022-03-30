C	linalg.f
C
Clinear algebra functions and subroutines
C
C
Cfunctions
************************************************************************
      function lengthvect(vect)
      real vect(3), lengthvect

      lengthvect = sqrt((vect(1)**2) + (vect(2)**2) + 
     .                   (vect(3)**2))

      return
      end

************************************************************************

      function dotprod(vecta,vectb)

      real
     +      dotprod, vecta(3), vectb(3)

      dotprod = (vecta(1) * vectb(1)) + (vecta(2) * vectb(2)) + 
     .          (vecta(3) * vectb(3))

      return
      end

Csubroutines
************************************************************************

      subroutine crossprod(vecta, vectb, vectc)

      real vecta(3), vectb(3), vectc(3)

      vectc(1) = (vecta(2) * vectb(3)) -  (vecta(3) *  vectb(2))
      vectc(2) = (vecta(3) * vectb(1)) -  (vecta(1) *  vectb(3))
      vectc(3) = (vecta(1) * vectb(2)) -  (vecta(2) *  vectb(1))

      return
      end

************************************************************************
 
      subroutine  normalize(vector)

      integer n
      real vector(3), lengthvect, thelength

      thelength = lengthvect(vector)
      do 400 n = 1,3
        vector(n) = vector(n) / thelength
400   continue

      return
      end

************************************************************************
