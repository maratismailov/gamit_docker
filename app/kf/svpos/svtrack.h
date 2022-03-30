 
*     Common block for svtrack progrma
 
*  max_sat       - Maximum number of satellites
 
      integer*4 max_sat
 
      parameter ( max_sat = 100 )
 
*  toe_jd(max_sat)  - Ephemeris times for each satelite (as jd)
*   toe(max_sat)        - Ephemeris time in seconds from start
*                       - GPS week.
*  svs_xyz(3,max_sat)    - XYZ of satellite
*  svs_loc(3,max_sat)    - Phi, Lambda, and height of sattellite
*  svs_ang(3,max_sat)    - Zenith angle, Azimuth and range to sat.
*  start, stop, step     - Start, stop and step for position
*  aode(max_sat)         - Age of ephemeris
*  af0(max_sat), af1(max_sat), af2(max_sat)     - Clock corrections
*  crs(max_sat), crc(max_sat)    - Radial corrections
*  dn(max_sat)   - Change to mean motion
*  m0(max_sat)   -Mean anomaly
*  cuc(max_sat), cus(max_sat)    - Correction for latitude
*  ecc(max_sat)          - Eccentricity
*  art(max_sat)          - root semimajor axis
*  cic(max_sat), cis(max_sat)    - Inclination corrections
*  om0(max_sat)          - longitude of node
*  i0(max_sat)           - Inclination
*  w(max_sat)            - argument of perigee
*  omd(max_sat)          - Omega dot
*  idt(max_sat)          - Inclination rate
*  cflg12(max_sat), pflg12(max_sat)      - Flags
*  weekno                - Week number
*  svacc, svhealth, tgd - Miscellaneous
*  aodc(max_sat)         - Age of clokc
*  site_xyz(3)           - Site XYZ
*  site_loc(3)           - lat, long of site
*  loc_rot(3,3)          - Rotation matrix from XYZ to NEU
 
      real*8 toe_jd(max_sat), toe(max_sat), svs_xyz(3,max_sat),
     .    svs_loc(3,max_sat), svs_ang(3,max_sat),
     .    start, stop, step, aode(max_sat),
     .    af0(max_sat), af1(max_sat), af2(max_sat),
     .    crs(max_sat), crc(max_sat), dn(max_sat), m0(max_sat),
     .    cuc(max_sat), cus(max_sat), ecc(max_sat), art(max_sat),
     .    cic(max_sat), cis(max_sat), om0(max_sat), i0(max_sat),
     .    w(max_sat), omd(max_sat), idt(max_sat),
     .    cflg12(max_sat), pflg12(max_sat), weekno,
     .    svacc, svhealth, tgd, aodc(max_sat), site_xyz(3),
     .    site_loc(3), loc_rot(3,3), lat_deg, long_deg
 
*  num_sat       - Number of satelites found
*  prn(max_sat)  - Prn nubers of sateliites (0 initially)
 
      integer*4 num_sat, prn(max_sat)
 
*  nav_file      - NAme of Rinex navigation file

      character*256 nav_file

* site_name      - Name of the site

      character*8 site_name

*---------------------------------------------------------------

      common / svtrack_com /  toe_jd, toe, svs_xyz,
     .    svs_loc, svs_ang, start, stop, step, aode, af0, af1, af2,
     .    crs, crc, dn, m0, cuc, cus, ecc, art,
     .    cic, cis, om0, i0, w, omd, idt,
     .    cflg12, pflg12, weekno,
     .    svacc, svhealth, tgd, aodc, site_xyz,
     .    site_loc,  loc_rot, lat_deg, long_deg,
     .    num_sat, prn,  nav_file, site_name
 
