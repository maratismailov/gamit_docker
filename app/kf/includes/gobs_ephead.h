 
*     This is the declarations for the EPHEAD Block of the Gobs files
*     This block contains information common to all observations at
*     this epoch, and information common to all satellites from a
*     specific receiver.

*     These records are split into two parts.  The main ephead
*     block and then a series of records which contain the rcv
*     dependent information.  There are only as many of these
*     latter records as needed for the total number of sites
*     (This procedure is used to save space in the disk files.)
*     The user never sees the buffered recorded since that are
*     automatically read and written in the standard read and
*     write routines.
 
 
*   ephead  - This is the first word of the epoch header block
 
      integer*4 ephead
 
*   medium_flag(max_grcv) - Meteorological data status at
*           - this epoch.
*           - Bit  Meaning
*           -   1  Pressure Bad
*           -   2  Temperature Bad
*           -   3  Rel. Humidity bad
*           -   4  WVR data bad.
*   num_obs_this_ep  - Number of observations at this epoch
*   kd_flag(max_gkd_sites)    - Kinematic flag for the kinematic
*           - sites.  The bit meanings are:
*           - Bit  Meaning
*           -   1  Site in motion do not process
 
      integer*4 medium_flag(max_grcv), num_obs_this_ep,
     .    kd_flag(max_gkd_sites)
 
*   pressure(max_grcv)    - Surface pressure (mbar)
*   temp_C(max_grcv)  - Surface temperature (C)
*   rel_hum(max_grcv) - Relative humidity (0-1)
*   usr_hyd_ep(max_grcv)   - User hydrostatic delay by
*                       - site
*   usr_wet_ep(max_grcv)   - User hydrostatic delay by
*                       - site
*   est_ion_mod_ep(max_grcv)   - Estimate of the L1 ionosphere
*                         model value  at this epoch (cycles)
*   est_ion_los_ep(max_gsvs,max_grcv) - Estimate of the LOS
*                         delay (cycles)
*   est_atm_ep(max_grcv)  - Estimate of the atmospheris delay
*                           correction (cycles)
*   est_ataz_ep(2,max_grcv) - Estimate of the azimithal asymetrry
*                           correction (cycles)
*   kd_offset_L1(3,max_gkd_sites) - L1 anntenna position
*                       - fof the kinematic stations.  These are
*                       - NEU (m).  NOTE: Change of order from GAMIT
*   kd_offset_L2(3,max_gkd_sites) - L2 anntenna position
*                       - fof the kinematic stations
 
      real*4 pressure(max_grcv), temp_C(max_grcv),
     .    rel_hum(max_grcv), usr_hyd_ep(max_grcv),
     .    usr_wet_ep(max_grcv), est_ion_mod_ep(max_grcv), 
     .    est_ion_los_ep(max_gsvs,max_grcv), est_atm_ep(max_grcv),
     .    est_ataz_ep(2,max_grcv),  kd_offset_L1(3,max_gkd_sites),
     .    kd_offset_L2(3,max_gkd_sites)
 
*   epoch_mjd - Modified Julian date plus fractional days for the
*           - nominal epoch of theses observations.  (Actuall
*           - times of measurements may differ by several seconds
*           - from this time.
*   epoch_sod - Epoch seconds of day (used for greater accuracy in the
*               epoch.
*   rcv_clk(2,max_grcv)   - Receiver clock offset and rate
*           - at this epoch (seconds and seconds/second)
*   svs_clk(2,max_gsvs) - Satellite clock offset and rates at
*           - this epoch (seconds and seconds/second)
*   kd_site_pos(3,max_gkd_sites)    - Positions of the kinematic
*           - sites. (Cartesian XYZ (m)) (Cartesian XYZ (m))
*   kd_site_vel(3,max_gkd_sites)    - Velocity for kinematic/dynamic
*             sites (m/sec)
*   kd_site_acc(3,max_gkd_sites)    - Acceleration for kinematic/dynamic
*             sites (m/sec**2)
 
 
      real*8 epoch_mjd, epoch_sod, rcv_clk(2,max_grcv), 
     .    svs_clk(2,max_gsvs), kd_site_pos(3,max_gkd_sites),
     .    kd_site_vel(3,max_gkd_sites), kd_site_acc(3,max_gkd_sites)
 
*     End of block information:
 
*   last_ephead_word    - Last ephem header word (used to get
*           - size of block)
*   last_eph_rcv        - Last word of the receiver dependend
*             block.
*   last_eph_gkd        - Last word of the kinematic/dynamic site
*             block.
*   dummy_ephead(GO_RECL)   - Dummy words to ensure that go_recl word
*           - boundary falls in here
 
      integer*4 last_ephead_word, last_eph_rcv, last_eph_gkd,
     .          dummy_ephead(GO_RECL)

*   eph_rcv_buff(eph_buff_len) - A buffer to contain the rcv
*             dependent information.  (There is one of these
*             records per site.)

      integer*4 eph_rcv_buff(eph_buff_len)
 
*...................................................................
 
      common / gobs_ephead / ephead, num_obs_this_ep,
     .    epoch_mjd, epoch_sod, svs_clk,
     .    last_ephead_word,
     .    dummy_ephead

*     These values are copied to a site dependent buffer
      common / gobs_ephead_rcv / rcv_clk, medium_flag,
     .    pressure, temp_C, rel_hum, usr_hyd_ep,
     .    usr_wet_ep, est_ion_mod_ep, est_ion_los_ep, 
     .    est_atm_ep, est_ataz_ep, last_eph_rcv

*     These values are kinematic/dynamic site dependent
      common / gobs_ephead_gkd / kd_flag, kd_offset_L1, kd_offset_L2,
     .    kd_site_pos, kd_site_vel, kd_site_acc, last_eph_gkd 
 
      common / gobs_eph_buff / eph_rcv_buff
