CTITLE COMP_PTIDE

      subroutine comp_ptide( pos_xyz, pmx, pmy, dXYZ, dNEU )

*     Routine to compute the pole tide correction for globk.
*     Based IERS 96 Conventions.

      include '../includes/const_param.h'

*     The inputs are the station postion in Cartesian XYZ coordinates,
*     the computes of pole position PMX and PMY (normal convention for
*     x and y pole position (mas).

*     The outputs are dXYZ cartesain coordinates and dNEU globk 
*     convention North, East and Up.

* PASSED variables
*-----------------
* pos_yxz(3)  -- Cartesian coordinates (m)
* pmx, pmy    -- Pole position (conventional X and Y pole position,
*                units mas)
* dXYZ(3)     -- Contribution to XYZ coordinates (m)
* dNEU(3)     -- Contribution to NEU coordinates (m)

      real*8 pos_xyz(3), pmx, pmy, dXYZ(3), dNEU(3) 

* LOCAL variables
* ---------------
* colat, long  -- Colatitude and longitude (radians)
* rot_mat(3,3) -- Rotation matrix from XYZ to NEU or visa versa depending
*                 on which call to rotate_geod
* pos_geod(3)  -- Geodetic colat, long and height (rad, rad, m)

      real*8 colat, long, rot_mat(3,3), pos_geod(3)

***** Compute the NEU contributions.  Convert the XYZ positions to 
*     Lat, long and height

      call XYZ_to_GEOD(rot_mat, pos_xyz, pos_geod)

      colat = pos_geod(1)
      long  = pos_geod(2)
  
*     Compute dNEU, convert mas to arc-seconds for pmx and pmy
* MOD TAH 140406: Replaced coefficients with more accurate calcualtion
*     (omega^2 * r_e^2)/2 * h_2/g   meters/radians  U
*     (omega^2 * r_e^2)   * l_2/g   meters/radians  NE
*     with omega = 7.292115e-5  rad/s
*          r_e   = 6371000  m  (approx mean radius of GPS sites) 
*          h_2 l_2 = 0.6207 0.0836  (IERS 2010 standards)
*          g     = 9.7803   m/sec^3
* echo 7.29211e-05 6371000 9.780318458 0.6207 0.0836 | awk '{printf(" %10.2f %10.2f \n",($1^2*$2^2)/(2*$3)*$4*(3.1415926535/180)/3.600,($1^2*$2^2)/($3)*$5*(3.1415926535/180)/3.600)}'
*      33.20       8.94  mm/mas
* IERS 2003 32 9 mm/mas; IERS2010 33 9 mm/mas (note 1 mm radial difference
*     between the two standards.
* MOD TAH 200220: d-3 is convert mas->arc sec: Result is meters
* Sign on dNEU(1) is for North; IERS stanard is dtheta which is South
      dNEU(1) =  0.00894d-3*cos(2*colat)*(pmx*cos(long)-pmy*sin(long))
      dNEU(2) =  0.00894d-3*cos(colat)  *(pmx*sin(long)+pmy*cos(long))
      dNEU(3) = -0.03320d-3*sin(2*colat)*(pmx*cos(long)-pmy*sin(long))

*     Now convert dNEU to dXYZ 
      call rotate_geod(dNEU, dXYZ, 'NEU', 'XYZ', pos_xyz, pos_geod,
     .                 rot_mat)

****  Thats all
      return
      end

