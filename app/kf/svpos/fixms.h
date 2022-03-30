 
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
*  P1(max_sat), P2(max_sat), L1(max_sat), L2(max_sat)
*                        - Range and phase measurements (m)
*  data_epoch            - Epoch of the data (jd)
*  data_sec              - Seconds tag for the data epoch
*  prev_site_clk         - Running value of clock.  May be entered
*                          through runstring.
*  init_dms              - Initial millisecond offset of clovk
*  omc(max_sat)          - Oobs minus C (m) 
*  apart(4,max_sat)      - Partial derivatives
*  kgain(4,max_sat),temp_gain(4,max_sat) - Kalman Gain and
*                         a temp intermeridate comp.
*  cov_parm(4,4), sol_vec(4) - Covariance matrix and solution
*                         vector for the parameters
*  wn(4)                  - Process noise model per step
*  acat_mat(max_sat,max_sat) - Covariance matrix of data plus
*                          predicted parameter noise
*  dx(4)                 - Adjustmenst to the parameters
*  out_spacing  - Spacing for output values (secs)
 
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
     .    C1o(max_sat), C2o(max_sat),
     .    P1c(max_sat), P2c(max_sat), L1c(max_sat), L2c(max_sat),
     .    data_epoch, data_sec, prev_site_clk, init_dms,
     .    omc_cl1(max_sat), omc_cl2(max_sat), omc_pl1(max_sat),
     .    omc_pl2(max_sat)
 
*  num_sat       - Number of satelites found
*  prn(max_sat)  - Prn nubers of sateliites (0 initially)
*  num_chan      - NUmber of channels observed at this epoch
*  num_used      - Actual number of channels used.
*  chan(max_sat) - Channels observed at this epoch (PRN numbers)
*  num_data_types - NUmber of data types 
*  debug_start    = Start epoch for debug
*  debug_end      - End epoch for debug
 
      integer*4 num_sat, prn(max_sat), num_chan,  chan(max_sat),
     .          nchi, num_used, num_data_types, debug_start, debug_end,
     .          p1f(max_sat), p2f(max_sat), c1f(max_sat), c2f(max_sat),
     .          l1f(max_sat), l2f(max_sat),
     .          p1i(max_sat), p2i(max_sat), c1i(max_sat), c2i(max_sat),
     .          l1i(max_sat), l2i(max_sat), plf, date_obs(5)

*  omc_OK(max_sat) - true if omc is OK else set false.
*  fix_phase       - Set true if option passed to adjust the phase values
*                    for millisecond jumps.  (This converts files back to
*                    standard).
      logical omc_OK(max_sat), fix_phase
 
*  nav_file      - NAme of Rinex navigation file
*  data_file     - Name of Rinex data file

      character*256 nav_file, data_file

* site_name      - Name of the site
* data_types(5)  - Type of data in the channels
* out_type       - Sets output coordinate type (XYZ or NEU)

      character*8 site_name, data_types(5)
      character*4 out_type

*---------------------------------------------------------------

      common / svtrack_com /  toe_jd, toe, svs_xyz,
     .    svs_loc, svs_ang, start, stop, step, aode, af0, af1, af2,
     .    crs, crc, dn, m0, cuc, cus, ecc, art,
     .    cic, cis, om0, i0, w, omd, idt,
     .    cflg12, pflg12, weekno,
     .    svacc, svhealth, tgd, aodc, site_xyz,
     .    site_loc,  loc_rot, lat_deg, long_deg,
     .    P1o, P2o, L1o, L2o, c1o, c2o, P1c, P2c, L1c, L2c,
     .    data_epoch,  data_sec, prev_site_clk, init_dms,
     .    omc_pl1, omc_pl2, omc_cl1, omc_cl2,
     .    num_sat, prn,  num_chan,  chan, nchi, num_used, 
     .    num_data_types, debug_start, debug_end, 
     .    P1f, P2f, L1f, L2f, c1f, c2f,
     .    P1i, P2i, L1i, L2i, c1i, c2i, plf, date_obs,
     .    omc_ok, fix_phase,
     .    nav_file, data_file, site_name, data_types, out_type
 
