CTITLE ROTATE_SOL

      subroutine rotate_sol( R1, sol )

*     Routine to apply the R1 solutiuon matrix to solution vector
*
* PASSED VARIABLES

*   R1(3,3)     - The rotation matrice for this site
*   sol(3)      - The solution vector (adjustment)

      real*8 R1(3,3), sol(3)

* LOCAL VARIABLES

*   i,j     - Loop couners

      integer*4 i,j

*   R1_sol(3)   - Product of R1*sol


      real*8 R1_sol(3)

****  Loop over each element
      do i = 1,3
          R1_sol(i) = 0.d0
          do j = 1,3
              R1_sol(i) = R1_sol(i) + R1(i,j)*sol(j)
          end do
      end do

****  Now return the new values
      do i = 1,3
          sol(i) = R1_sol(i)
      end do

***** Thats all
      return
      end


