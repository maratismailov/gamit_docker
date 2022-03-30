CTITLE CROSS_PROD
C
      SUBROUTINE cross_prod(VECA,VECB,vecc, SANG)

      implicit none 
 
C
C     ROUTINE TO COMPUTE THE SIN OF THE ANGLE BETWEEN THE BASELINES
C     FROM THE CROSS PRODUCT OF THE VECTORS.
c
c     Form:
c
c          C = A*B
c
C
c Variables
c ---------
c veca  -- the first vector in the cross product
c vecb  -- the second vector in the cross product
c vecc  -- the cross product of a and b
c
c sang  -- sine of the angle between the vectors
c
      real*8 VECA(3),VECB(3),SANG
C
      real*8 VECC(3)
C
      real*8 DOT3
C
C**** COMPUTE CROSS PRODUCT VECTOR
      VECC(1) = VECA(2)*VECB(3) - VECA(3)*VECB(2)
      VECC(2) = VECA(3)*VECB(1) - VECA(1)*VECB(3)
      VECC(3) = VECA(1)*VECB(2) - VECA(2)*VECB(1)
C
C**** NOW COMPUTE SINE OF ANGLE: sang commented out since it is not
c     needed
      SANG = SQRT( DOT3(VECC,VECC)/(DOT3(VECA,VECA)*DOT3(VECB,VECB)) )
C
      RETURN
      end
 
