       subroutine tran_part(x2,deriv)
C
C     Determine 7 parameter transformation partials for satellite
C     position coordinates given in system 1 and system 2. The
C     estimated parameters give the transformation from system 2
C     to system 1:  X1 = T + (1+SCALE)*R*X2 where T is a vector
C     of translation parameters, scale is a scalar, and R is the
C     linearized rotation matrix for rotations about the 1,2, & 3 axes.
C     The units of X and T are kilometers;the angles R1, R2, R2 are in arcs
C
C     S McClusky September 30, 1994
C
      implicit none
                                     
      include '../includes/dimpar.h'
      include '../includes/orbits.h'

      integer i,k
      real*8 x2(6),t(3),r(3),scale,deriv(3,7),convds
c
c     Radians to arcsec
      CONVDS= DATAN(1.D0)*4.D0/180.D0/3600.D0
C
C       Initial values of parameters
C
      DO I=1,3
        T(I)= 0.
        R(I)= 0. 
      enddo
      SCALE= 0.
C
C      compute derivatives
C
      K= 0
c
      DO I=1,3
        DO  K=1,7
          DERIV(I,K)= 0.D0
        enddo

C       Parameters 1-3 are translation - T(I)
          IF( I.EQ.1 ) DERIV(i,1)= 1.D0
          IF( I.EQ.2 ) DERIV(i,2)= 1.D0
          IF( I.EQ.3 ) DERIV(i,3)= 1.D0

C       Parameters 4-6 are orientation - R(I)
          IF( I.EQ.1 ) DERIV(i,5)= -X2(3) * (1.D0 + SCALE) * CONVDS
          IF( I.EQ.1 ) DERIV(i,6)=  X2(2) * (1.D0 + SCALE) * CONVDS
          IF( I.EQ.2 ) DERIV(i,4)=  X2(3) * (1.D0 + SCALE) * CONVDS
          IF( I.EQ.2 ) DERIV(i,6)= -X2(1) * (1.D0 + SCALE) * CONVDS
          IF( I.EQ.3 ) DERIV(i,4)= -X2(2) * (1.D0 + SCALE) * CONVDS
          IF( I.EQ.3 ) DERIV(i,5)=  X2(1) * (1.D0 + SCALE) * CONVDS

C       Parameter 7 is scale
          IF( I.EQ.1 )
     1    DERIV(i,7)= X2(1) + (   R(3)*X2(2) - R(2)*X2(3) ) * CONVDS
          IF( I.EQ.2 )
     1    DERIV(i,7)= X2(2) + ( - R(3)*X2(1) + R(1)*X2(3) ) * CONVDS
          IF( I.EQ.3 )
     1     DERIV(i,7)= X2(3) + (   R(2)*X2(1) - R(1)*X2(2) ) * CONVDS

      enddo

      RETURN
      END
