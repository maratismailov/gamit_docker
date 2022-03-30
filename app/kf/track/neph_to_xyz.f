CTITLE NEPH_TO_XYZ
 
      subroutine neph_to_xyz(t, pn, sys)

      implicit none
 
*     Routine to compute XYZ coordinates of satellite at time t for
*     sattellite number i.  SYS is 'I' or 'E' for inertial or Earth
*     fixed. (From Rinex file)
 
      include 'track_com.h'
 
*   t       - Time for computation (day number)
 
 
      real*8 t
 
*   i       - Satellite index number
*   pn      - Satellite PRN number 
 
      integer*4 i, pn
 
*   sys     - System for results.
 
 
      character*(*) sys
 
* LOCAL VARIABLES
 
*   j           - Loop counter in Kepler's equation
 
 
      integer*4 j
 
*   gm      - GM
*   eom     - Earth rotation rate (rads/sec)
*   a       - Semimajor axis
*   n0      - Mean motion
*   tk      - Time of epoch from toe (seconds)
*   n       - COrrected mean motion
*   mk      - Mean anomaly
*   ek      - Eccentric anomaly
*   vk      - true anomaly
*   sinvk, cosvk    - Sin and cos of true anomaly
*   pk      - argument of latitude
*   duk, drk, dik   - Coorections to arg of lat, radius and
*           - inclinations
*   uk      - Argumenr of latitude
*   rk      - radius at time tk
*   ik      - Inclination at tk
*   xpk, ypk    - Inplane coordiantes
*   omk     - Longitude of the asecdning node
*   rot_mat(3,3)    - Rotation matrix from XYZ to NEU
 
 
 
 
      real*8 gm, eom, a, n0, tk, n, mk, ek, vk, sinvk, cosvk, pk,
     .    duk, drk, dik, uk, rk, ik, xpk, ypk, omk, rot_mat(3,3) 

*     Find the index number given the PRN
      i = 0
      do j = 1, max_sat
         if( prn_sp3(j).eq.pn ) i = j
      end do

      gm = 3.986005d14
      eom = 7.2921151467d-5
 
      a = art(i)*art(i)
      n0 = sqrt(gm/a**3) 

      tk = (t-toe_jd(i))*86400.0d0
      n = n0 + dn(i)
      mk = m0(i) + n*tk
 
****  Solve Keplers equation
      ek = mk
      do j = 1, 10
          ek = mk + ecc(i)*sin(ek)
      end do
 
***** Get the true anomaly
      sinvk = sqrt(1-ecc(i)**2)*sin(ek)/(1 - ecc(i)*cos(ek))
      cosvk = (cos(ek)-ecc(i))/(1-ecc(i)*cos(ek))
 
      vk = atan2(sinvk, cosvk)
 
*     Argument of latitude
      pk = vk + w(i)
 
*     Correction terms
      duk = cus(i)*sin(2*pk) +cuc(i)*cos(2*pk)
      drk = crs(i)*sin(2*pk) +crc(i)*cos(2*pk)
      dik = cis(i)*sin(2*pk) +cic(i)*cos(2*pk)
 
      uk = pk + duk
      rk = a*(1-ecc(i)*cos(ek)) + drk
      ik = i0(i) + dik + idt(i)*tk
 
*     Get inplane coordinates
      xpk = rk*cos(uk)
      ypk = rk*sin(uk)
 
*     Compute the longitude of the ascending node
      omk = om0(i) + omd(i)*tk
 
*     If we are in Earth fixed frame account for rotation of Earth
      if( sys(1:1).eq.'E' .or. sys(1:1).eq.'e') then
          omk = omk - eom*(tk+toe(i))
      end if
 
*     Get three_d coordinates
      svs_xyz(1,pn) = xpk*cos(omk) - ypk*sin(omk)*cos(ik)
      svs_xyz(2,pn) = xpk*sin(omk) + ypk*cos(omk)*cos(ik)
      svs_xyz(3,pn) = ypk*sin(ik)
 
*     Now convert to latitude and longitude
      call xyz_to_geod( rot_mat, svs_xyz(1,pn), svs_loc(1,pn))

****  Thats all
      return
      end
 
 
 
