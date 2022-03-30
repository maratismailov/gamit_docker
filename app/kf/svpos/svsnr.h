 
*     Common block for svtrack progrma
 
*  max_sat       - Maximum number of satellites
 
      integer*4 max_sat, max_data_types
 
      parameter ( max_sat = 100 )
* MOD TAH 151021:Updated code to handle large number of data types
      parameter ( max_data_types = 36 ) ! FIXED: If > 18 code needs mod

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
*  P1(max_sat), P2(max_sat), L1(max_sat), L2(max_sat)
*                        - Range and phase measurements (m)
*  S1(max_sat), S2(max_Sat) - SNR values at L! and L2
*  data_epoch            - Epoch of the data (jd)
*  omc(max_sat)          - Oobs minus C (m) 
 
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
     .    site_loc(3), loc_rot(3,3), lat_deg, long_deg,
     .    P1o(max_sat), P2o(max_sat), L1o(max_sat), L2o(max_sat),
     .    P1c(max_sat), P2c(max_sat), L1c(max_sat), L2c(max_sat),
     .    data_epoch, omc(max_sat), s1(max_sat), s2(max_sat)
 
*  num_sat       - Number of satelites found
*  prn(max_sat)  - Prn nubers of sateliites (0 initially)
*  num_chan      - NUmber of channels observed at this epoch
*  num_used      - Actual number of channels used.
*  chan(max_sat) - Channels observed at this epoch (PRN numbers)
*  num_data_types - NUmber of data types 
*  debug_start    = Start epoch for debug
*  debug_end      - End epoch for debug
*  code_type(max_sat) - Indicates if P-code or Xcorrelation (Set to
*                   0 for P-code L1, 1 for L2 X-correlation).  BAsed
*                   on C1 or P1 range measurement.
*  rcv_type       - Coded value for receiver type
*                 - 0  -- Linear SNR values
*                 - 1  -- Natural log SNR values
 
      integer*4 num_sat, prn(max_sat), num_chan,  chan(max_sat),
     .          nchi, num_used, num_data_types, debug_start, debug_end,
     .          code_type(max_sat), rcv_type

*  omc_OK(max_sat) - true if omc is OK else set false.
      logical omc_OK(max_sat)
 
*  nav_file      - NAme of Rinex navigation file
*  data_file     - Name of Rinex data file

      character*256 nav_file, data_file

* site_name      - Name of the site
* data_types(max_data_types)  - Type of data in the channels
* out_type       - Sets output coordinate type (XYZ or NEU)

      character*8 site_name, data_types(max_data_types)
      character*4 out_type

*---------------------------------------------------------------

      common / svtrack_com /  toe_jd, toe, svs_xyz,
     .    svs_loc, svs_ang, start, stop, step, aode, af0, af1, af2,
     .    crs, crc, dn, m0, cuc, cus, ecc, art,
     .    cic, cis, om0, i0, w, omd, idt,
     .    cflg12, pflg12, weekno,
     .    svacc, svhealth, tgd, aodc, site_xyz,
     .    site_loc,  loc_rot, lat_deg, long_deg,
     .    P1o, P2o, L1o, L2o, P1c, P2c, L1c, L2c,
     .    data_epoch,  omc, S1, S2,
     .    num_sat, prn,  num_chan,  chan, nchi, num_used, 
     .    num_data_types, debug_start, debug_end, code_type, 
     .    rcv_type, omc_ok,
     .    nav_file, data_file, site_name, data_types, out_type
 
