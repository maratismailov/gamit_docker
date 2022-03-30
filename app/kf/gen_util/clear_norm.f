ctitle
 
      subroutine clear_norm(norm_eq, sol_vec, np)

      implicit none 
 
c
c     routine to clear the normal equations and the b vector before
c     acummulaing these values.
c
c Variables
c ---------
c norm_eq  -- real*8 square array containing the normal equations
c sol_vec -- real*8 vector containing the left hand side of the
c        equations
c np  -- the dimensioned size of the matrix and vector
c
c Local variables
c ---------------
c isize -- the numberof elements in the normal equations
c
*                                 ! the SOLVK parameter file
      include '../includes/kalman_param.h'
c
      integer*4 np 
c
      real*8 norm_eq(np,np), sol_vec(np)
c
* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters
      integer*8 isize

      data I8 / 1 /

c.... compute the size of normal equations array
      isize = (I8*np)*np
c
c.... use the vis scalar multiply routine to clear matrix
      call dwint8(0.d0, norm_eq,1, isize)
c
c.... clear the sol_vec the same way
C     call dwsmy(0.d0, sol_vec,1, sol_vec,1, np)
      call dwint(0.d0, sol_vec,1, np)
c
      return
      end
c.........................................................................
