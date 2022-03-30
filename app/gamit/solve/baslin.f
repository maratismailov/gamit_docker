      real*8 function baslin(pos1,pos2)
c
c     J.L. Davis 870305
c
c     Routine to calculate the baseline length between two sites,
c     the 3-D cartesian coordinates for which are stored in POS1
c     and POS2.  The units for the baseline length are the same as
c     the units for POS1 and POS2, which must be identical.
c
c     Input variables
c     ---------------
c     POS1, POS2      3-D cartesian site coordinates
c
      real*8 pos1(3), pos2(3)
c
c     Internal variables
c     ------------------
c     i               Loop index
c
      integer i
c
c.... Initialize for summation
      baslin = 0
c
c.... Loop over X-Y-Z
      do 100 i = 1, 3
c
c....     Add the contribution for this coordinate
          baslin = baslin + ( pos2(i) - pos1(i) ) ** 2
c
  100 continue
c
c.... Take the square root to give length
      baslin = sqrt(baslin)
c
      end
