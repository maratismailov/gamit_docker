      subroutine phase_mod( elev, dphs_l1, dphs_l2 )

      implicit none 

*     Routine to computed phase center model based on elevation angle

      include '../includes/const_param.h'

* PASSED VARIABLE

* elev  - Elevation angle (deg)
* dphs_l1 - Change in phase at L1 (mm)
* dphs_l2 - Change in phase at L2 (mm)

      real*8 elev, dphs_l1, dphs_l2

* LOCAL VARIABLES

* phs_el, phs_1l, phs_l2 - Elevation, phase corrections in mm for
*     L1 and L2 for rogue antenas only.
* err  - Estimated error in interpolation

      real*8 phs_el(17), phs_l1(17), phs_l2(17), err_l1, err_l2

* nn  - Position in table for interpolating from

      integer*4 nn
 
      data phs_el  /90.d0, 85.d0, 80.d0, 75.d0, 70.d0, 65.d0, 
     .              60.d0, 55.d0, 50.d0, 45.d0, 40.d0,
     .              35.d0, 30.d0, 25.d0, 20.d0, 15.d0, 10.d0 /
      data phs_L1  / 0.d0,  0.d0, -1.d0, -3.d0, -4.d0, -6.d0,
     .              -8.d0,-10.d0,-12.d0,-13.d0,-14.d0,
     .             -14.d0,-13.d0,-12.d0,-11.d0, -9.d0, -6.d0 / 
      data phs_L2  / 0.d0,  0.d0, -1.d0, -2.d0, -2.d0, -4.d0,
     .              -5.d0, -6.d0, -8.d0, -9.d0, -9.d0,
     .             -10.d0,-10.d0, -9.d0, -9.d0, -8.d0, -5.d0 /

****  Use the polyint routine for interpolation.
C     call polint(phs_el, phs_l1, 17, elev, dphs_l1, err_l1)
C     call polint(phs_el, phs_l2, 17, elev, dphs_l2, err_l2)
      if( elev.lt.0.d0 ) then
          write(*,110) elev
 110      format(' WARNING: Satellite at elevation ',f8.2,' deg')
      else if( elev.lt.10.d0 ) then
          dphs_l1 = phs_l1(17)
          dphs_l2 = phs_l2(17)
      else if ( elev.gt.90.d0 ) then
          write(*,120) elev
 120      format(' WARNING: Satellite at elevation ',f8.2,' deg')
      else
          nn = int((90 - elev)/5.d0) + 1
          if( nn.gt.15 ) nn = 15
          call polint(phs_el(nn), phs_l1(nn),  3, 
     .                elev, dphs_l1, err_l1)
          call polint(phs_el(nn), phs_l2(nn),  3, 
     .                elev, dphs_l2, err_l2)
      end if
        
          
****  Thats all
      return
      end
