      program test_phs

      implicit none 

      real*8 elev, dphs_l1, dphs_l2

      do elev = 5, 90, 1.0
         call phase_mod( elev, dphs_l1, dphs_l2)
         write(*,100) elev, dphs_l1, dphs_l2
 100     format(3F15.4)
      end do

      end

