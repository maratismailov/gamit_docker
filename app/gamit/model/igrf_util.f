C Subroutines from IGRF that are used in all Versions

      SUBROUTINE DMDDEC (I,M,X)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      real*8 DE, EM, X
      integer*4 i,m
      
      DE = I
      EM = M
      IF (I.LT.0) EM = -EM
      X = DE + EM/60.0
      RETURN
      END
C 
      SUBROUTINE DDECDM (X,I,M)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      real*8 X, SIG, DR, T
      integer*4 i,m, isig

      SIG = SIGN(1.1d0,X)
      DR = ABS(X)
      I = INT(DR)
      T = I
      M = NINT(60.*(DR - T))
      IF (M.EQ.60) THEN
       M = 0
       I = I + 1
      ENDIF
      ISIG = INT(SIG)
      IF (I.NE.0) THEN
       I = I * ISIG
      ELSE
       IF (M.NE.0) M = M * ISIG
      ENDIF
      RETURN
      END
