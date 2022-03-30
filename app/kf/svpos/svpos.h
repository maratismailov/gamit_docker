 
*     Common block for svtrack progrma
 
*  max_sat       - Maximum number of satellites
*  max_ent       - Maximum entries per satellites
*  max_data_types  -- Maximum number of data types

      integer*4 max_sat, max_ent, max_data_types
 
      parameter ( max_sat = 100 )
      parameter ( max_ent = 100 )
* MOD TAH 151021:Updated code to handle large number of data types
      parameter ( max_data_types = 36 ) ! FIXED: If > 18 code needs mod
 
*  toe_jd(max_sat,max_ent)  - Ephemeris times for each satelite (as jd)
*  toe(max_sat,max_ent)     - Ephemeris time in seconds from start
*                       - GPS week.
*  svs_xyz(3,max_sat)    - XYZ of satellite
*  svs_loc(3,max_sat)    - Phi, Lambda, and height of sattellite
*  svs_ang(3,max_sat)    - Zenith angle, Azimuth and range to sat.
*  start, stop, step     - Start, stop and step for position
*  aode(max_sat,max_ent)         - Age of ephemeris
*  af0(max_sat,max_ent), af1(max_sat,max_ent), af2(max_sat,max_ent)     - Clock corrections
*  crs(max_sat,max_ent), crc(max_sat,max_ent)    - Radial corrections
*  dn(max_sat,max_ent)   - Change to mean motion
*  m0(max_sat,max_ent)   -Mean anomaly
*  cuc(max_sat,max_ent), cus(max_sat,max_ent)    - Correction for latitude
*  ecc(max_sat,max_ent)          - Eccentricity
*  art(max_sat,max_ent)          - root semimajor axis
*  cic(max_sat,max_ent), cis(max_sat,max_ent)    - Inclination corrections
*  om0(max_sat,max_ent)          - longitude of node
*  i0(max_sat,max_ent)           - Inclination
*  w(max_sat,max_ent)            - argument of perigee
*  omd(max_sat,max_ent)          - Omega dot
*  idt(max_sat,max_ent)          - Inclination rate
*  cflg12(max_sat,max_ent), pflg12(max_sat,max_ent)      - Flags
*  weekno                - Week number
*  svacc, svhealth, tgd - Miscellaneous
*  aodc(max_sat,max_ent)         - Age of clokc
*  site_xyz(3)           - Site XYZ
*  site_loc(3)           - lat, long of site
*  loc_rot(3,3)          - Rotation matrix from XYZ to NEU
*  P1(max_sat), P2(max_sat), L1(max_sat), L2(max_sat)
*                        - Range and phase measurements (m)
*  data_epoch            - Epoch of the data (jd)
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
 
      real*8 toe_jd(max_sat,max_ent), toe(max_sat,max_ent), 
     .    svs_xyz(3,max_sat),
     .    svs_loc(3,max_sat), svs_ang(3,max_sat),
     .    start, stop, step, aode(max_sat,max_ent),
     .    af0(max_sat,max_ent), af1(max_sat,max_ent), 
     .    af2(max_sat,max_ent),
     .    crs(max_sat,max_ent), crc(max_sat,max_ent), 
     .    dn(max_sat,max_ent), m0(max_sat,max_ent),
     .    cuc(max_sat,max_ent), cus(max_sat,max_ent), 
     .    ecc(max_sat,max_ent), art(max_sat,max_ent),
     .    cic(max_sat,max_ent), cis(max_sat,max_ent), 
     .    om0(max_sat,max_ent), i0(max_sat,max_ent),
     .    w(max_sat,max_ent), omd(max_sat,max_ent), 
     .    idt(max_sat,max_ent),
     .    cflg12(max_sat,max_ent), pflg12(max_sat,max_ent), weekno,
     .    svacc, svhealth, tgd, aodc(max_sat,max_ent), site_xyz(3),
     .    site_loc(3), loc_rot(3,3), lat_deg, long_deg,
     .    P1o(max_sat), P2o(max_sat), L1o(max_sat), L2o(max_sat),
     .    P1c(max_sat), P2c(max_sat), L1c(max_sat), L2c(max_sat),
     .    data_epoch, 
     .    omc(max_sat), apart(4,max_sat), kgain(4,max_sat),
     .    temp_gain(4,max_sat), cov_parm(4,4), sol_vec(4),
     .    wn(4), acat_mat(max_sat,max_sat), dx(4), data_noise,
     .    chi, out_spacing

      real*4 data_rxver  ! Version of data rinex file
 
*  num_sat       - Number of satelites found
*  prn(max_sat)  - Prn nubers of sateliites (0 initially)
*  num_ent(max_sat) - number of entries per satellite
*  num_chan      - NUmber of channels observed at this epoch
*  num_used      - Actual number of channels used.
*  chan(max_sat) - Channels observed at this epoch (PRN numbers)
*  num_data_types - NUmber of data types 
*  debug_start    = Start epoch for debug
*  debug_end      - End epoch for debug
*  proc_start, proc_end - Processing start and end
 
      integer*4 num_sat, prn(max_sat),num_ent(max_sat), num_chan, 
     .          chan(max_sat), nchi, num_used, num_data_types, 
     .          debug_start, debug_end, proc_start, proc_end

*  omc_OK(max_sat) - true if omc is OK else set false.
*  dump - Set true to dump data (DUMP added to out_type)
      logical omc_OK(max_sat), dump
 
*  nav_file      - NAme of Rinex navigation file
*  data_file     - Name of Rinex data file
*  dump_file     - File with range and phase data

      character*256 nav_file, data_file, dump_file

* site_name      - Name of the site
* data_types(5)  - Type of data in the channels
* out_type       - Sets output coordinate type (XYZ or NEU):
* MOD TAH 130330: Added DUMP to out_type to dump data

      character*8 site_name 
      character*2 data_types(max_data_types)
      character*8 out_type

      character*8 anal_typ   ! String with analysis type (LC,L1,L2)

*---------------------------------------------------------------

      common / svtrack_com /  toe_jd, toe, svs_xyz,
     .    svs_loc, svs_ang, start, stop, step, aode, af0, af1, af2,
     .    crs, crc, dn, m0, cuc, cus, ecc, art,
     .    cic, cis, om0, i0, w, omd, idt,
     .    cflg12, pflg12, weekno,
     .    svacc, svhealth, tgd, aodc, site_xyz,
     .    site_loc,  loc_rot, lat_deg, long_deg,
     .    P1o, P2o, L1o, L2o, P1c, P2c, L1c, L2c,
     .    data_epoch,  omc, apart, kgain, temp_gain, cov_parm, sol_vec,
     .    wn, acat_mat, dx, data_noise, chi, out_spacing,
     .    data_rxver, 
     .    num_sat, prn, num_ent, num_chan,  chan, nchi, num_used, 
     .    num_data_types, debug_start, debug_end, omc_ok, dump,
     .    proc_start, proc_end, 
     .    nav_file, data_file, site_name, data_types, dump_file, 
     .    out_type, anal_typ
 
