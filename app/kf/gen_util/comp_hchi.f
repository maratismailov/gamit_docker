CTITLE COMP_HCHI
 
      subroutine comp_hchi( chi, work, org, pmu_parts, num_used)

      implicit none 
 
 
*     Routine to compute the horizonal Chi**2 given the current origin
*     value (3 translations, 3 rotations)
 
 
*   i,j,k,l     - Loop counters
*   num_used    - Number of sites used in origin
 
      integer*4 k,l, num_used
 
*   chi         - Value of chi**2 computed
*   horizontal  - horizontal displacement
*   org(6)      - current translations and rotations
*   pmu_parts(3,3,1)    - Polar motion/UT1 partials for XYZ for each
*               - site
*   pos(3)      - Computed new position adjustment after origin
*               - applied
*   radial      - Radial displacement at each site
*   total       - Total displacemement at site
*   work(3,2,1) - Working array.  Contains adjustmemt and radial
*               - partials
 
      real*8 chi, horizontal, org(6), pmu_parts(3,3,1), pos(3), radial,
     .    total, work(3,2,1)
 
 
      chi = 0.d0
 
*     Now sum up horizontal displacements for this choice
*     of origin
 
      do k = 1, num_used
 
          pos(1) = work(1,2,k) - org(1)
          pos(2) = work(2,2,k) - org(2)
          pos(3) = work(3,2,k) - org(3)
 
*         Now do the rotation
*                         ! Loop over x,y and UT1 values
          do l = 1,3
              pos(1) = pos(1) - pmu_parts(l,1,k)*org(3+l)
              pos(2) = pos(2) - pmu_parts(l,2,k)*org(3+l)
              pos(3) = pos(3) - pmu_parts(l,3,k)*org(3+l)
          end do
 
          radial = work(1,1,k)*pos(1) + work(2,1,k)*pos(2) +
     .             work(3,1,k)*pos(3)
 
          total =  pos(1)**2 + pos(2)**2 + pos(3)**2
          horizontal = abs(total - radial**2)
          chi = chi + horizontal
      end do
 
***** Thats all
      return
      end
 
