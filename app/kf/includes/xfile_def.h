*
*     This file contains the definition of the quanities in the
*     xk-files.  It is used with the read_xf and write_xf subroutines 
*     in the Ghandlers library.

*     Gang Chen, July 14, 1996
*  ----------------------------------------------------
*     Maximum of some  values used in xk-file

*   xf_maxtxt   - Maximum lines of text at the beginning of the
*               - xfile
*   xf_maxsat   - Maximum number of satellites allowed in xfile
*   xf_ncvsit   - Maximum number of channels allowed in a
*               - receiver
*   xf_maxprm   - Maxiumum number of parameters.  The clean
*               - routine set this max number of solve parameters
*               - We set to the actual max possible in a x-file
*   xf_maxlab   - Maximum number of labels for parameter strings
*   xf_maxdat   - Maximum number of data types.
*   xf_maxcsb   - Maximum number of bias parameters saved in
*               - x-file header
*   xf_maxext   - Maximum number of extra variables in type 2
*               - record
*   xf_maxsav   - Maximum number of save variables in type 4
*               - record
*   xf_maxspr   - Maximum number of spare varibales in type 5
*               - record
*   xf_maxclk - Maximum number of terms in the clock poynomials
 
      integer*4 xf_maxtxt, xf_maxsat, xf_ncvsit, 
     .    xf_maxdat, xf_maxcsb, xf_maxext, xf_maxsav, xf_maxspr,
     .    xf_maxclk
 
* Parmeter values valid on 920812:
* Increased xf_maxsat to 32 931007:
* MOD TAH 940516:
* Increased xf_maxprm by 6 to account for the plar motion/UT1 offset
*         rate parameters
* Increaed maxlab to 22 for new partials. (This still seems to be
*         one more than needed).
* MOD TAH 950629:
* Increased number of partials per satellite from 9 to 15 for the
* new parameters (once-per-rev etc.)
* MOD TAH 090114: Increased xf_maxdat from 7 to 15 to allow for
*     new receivers.
* MOD TAH 150211: Increased xf_maxdat from 15 to 27 to allow for
*     Galileo observabales.
 
      parameter ( xf_maxtxt =  200 )
 
* MOD TAH 180329: Increased to 100 to handle GNSS
      parameter ( xf_maxsat =  100 )
      parameter ( xf_ncvsit =   20 )
      parameter ( xf_maxdat  =  27 )
      parameter ( xf_maxcsb  = 200 )
      parameter ( xf_maxext  =   2 )
      parameter ( xf_maxsav  = xf_maxsat )
      parameter ( xf_maxspr  =   2 )
      parameter ( xf_maxclk =  4 )

*     Define the file channel
      integer*4 xf_lu, obs_lu

      common /xf_lu_com/ xf_lu,obs_lu
*-------------------------------------------------------------------
*     Define the xk file name
  
*   xf_xfile   - output X-file name for this xfile 
*   xf_rfile   - Rinex -file name for this xfile 
*   xf_tfile    - T-file name for this xfile
*   xf_jfile    - J-file name fpr this xfile
*   xf_rxfile   - input X-file name for this xfile
      character*256 xf_xfile, xf_rfile,xf_tfile, xf_jfile, 
     .   xf_rxfile   
      common /xf_header/xf_xfile, xf_rfile, xf_tfile, xf_jfile,
     .        xf_rxfile
*-------------------------------------------------------------------  
* Define each of the types of records in a xk-file.
 
*-------------------------------------------------------------------
*     RECORD AREA TYPE 1:     xk-file maker's information

*     First set of the paramaters.  These should be updated when the
*     xk-file is updated.

*   xf_nversn   - Software version number multiplied by 100.
*   xf_ntext    - Number of text records at top of xk-file
 
      integer*4 xf_nversn, xf_ntext
 
*   xf_text(xf_maxtxt)  - Header lines
 
      character*256 xf_text(xf_maxtxt)

 
      common / xf_type1 / xf_nversn, xf_ntext, xf_text
 
*...................................................................
*     TYPE 2: Header series data block
 
*   xf_siteinfo  -  X position information of site 
*   xf_pos(3)     -  initial position of site (XYZ)
*   xf_neu(3)    -  initial position of site  Latitude, Longitude, Height(m)
*                    when NEU, SWU 
*   xf_posflg   - FLAG of site position  (NEU) (SWU) (XYZ)
   
      real*8  xf_pos(3),xf_neu(3)
      character xf_posflg*3, xf_siteinfo*80
      
      
*   xf_srttime(5)   - Start epoch (yr mon day hr min)
*   xf_srtsec      - Start epoch seconds 
*   xf_endtime(5)   - End epoch (yr mon day hr min)
*   xf_endsec      - End epoch seconds 
*   xf_srtdoy(2)   - Start epoch (yr days)
*   xf_srtsect      - Start epoch seconds in day
*   xf_start_mjd, xf_end_mjd -- Start and stop MJD for the xfile.

      integer*4 xf_srttime(5),xf_endtime(5),xf_srtdoy(2)
      real*8   xf_srtsec, xf_endsec, xf_srtsect, 
     .         xf_start_mjd, xf_end_mjd

*   xf_offarp(3)    - Up, North, East from mark to antena ref point (m)
*   xf_offsL1(3)    - Up, North, East L1 antenna offset (m)
*   xf_offsL2(3)    - Up, North, East L2 antenna offset (m)
*   xf_elvcut   - Solve elevation cutof (degrees)
*   xf_clock (xf_maxclk) - Site clock polynomial (sec,sec/sec,..)
*   xf_serwgt   - Weight given to this series (Nominal = 1.0) ! no used
*   xf_extra(xf_maxext) - Extra real*8 values in Type 2 header
 
      real*8  xf_offarp(3), xf_offsL1(3), xf_offsL2(3),
     .    xf_elvcut, xf_clock (xf_maxclk), xf_te, xf_tr,
     .    xf_serwgt, xf_extra(xf_maxext)
 

*   xf_nsat     - Number of satellites in this xfile
*   xf_ischan(xf_maxsat)    - channel numbers of the satellites.
*   xf_prn(xf_maxsat)    - PRN numbers of the satellites.
*   xf_ndat     - Number of data types.
* MOD TAH 180322: Killed D1 and D2 since never used
*   xf_dattyp(xf_maxdat,7)    - Data types saved in the x-file:
*              for each GNSS: 1 GPS  , 2 Glonass, 3 Galileio, 
*                             4 Beidou  (5-7 not used)
*               -   1  - L1 carrier phase (cycles)
*               -   2  - L2 Carrier phase (cycles)
*               -   3  - P1 Pseudo range (m)
*               -   4  - P2 Pseudo range (m) Pcode
*               -   5  - C1 Pseudo range (m) C/A code
*               -   6  - C2 Pseudo range (m) L2C code
*               -   0  - unknown
*   xf_mdat     -  Number of recording column in xf_obsv
*                  (L1,L2,P1,P2,C1,D1,D2 in order, zero if no obs)
*
*   xf_lambda(xf_maxsat,xf_maxdat)  - Wavelength factors for
*               - each of the datatypes on each of the
*               - satellites: Code is
*               -  0  - No data
*               -  1  - unambiguous, undoubled values
*               - -1  - ambiguous, undoubled values
*               -  2  - unambiguous, undoubled values
*               - -2  - ambiguous, undoubled values
*   xf_nepoch   - Number of epochs
*   xf_inter    - Obervation interval (secs)
*   xf_ircint   - Original receiver sampling interval before
*               - decimation
*   xf_typtime    - Type of time:
*               - 1  - UTC
*               - 2  - GPST
*               - ?  - Others unknown
*   xf_isessn   - Session number
*   xf_iy, xf_im, xf_id, xf_ihr, xf_min - Start epocj ymdhms ! not used
*   xf_ietide    - Bitmapped tide model description. ! not used
*                 Bit    Model
*                   1    solid earth tides
*                   2    K1 frequency dependant earth tide
*                   3    Pole tide
*                   4    Ocean Tide
*   xf_isptide  - Bitmapped short period EOR Model flag ! not used
*                 Bit    Model
*                   1    UT1
*                   2    pole
*   xf_nclock   - Number of terms in the clock polynomial
*   xf_jde      - (PEP) julian day of ephemeris 
*   xf_jdr      - (PEP) julian day of Earth rotation parameters! not used
*   xf_avlmet   - Availability of met data:  Binary coded
*               - BIT   Meaning
*               -  1 - Pressure available
*               -  2 - Temperature available
*               -  3 - Relative humidity available
*               -  4 - WVR delay available
*   xf_nslip    - Extra number of bias parameter flags
*   xf_nextra   - Number of extra real*8 values
 
      integer*4 xf_npart, xf_nsat, xf_ischan(xf_maxsat), xf_ndat,
     .    xf_dattyp(xf_maxdat,7), xf_lambda(xf_maxsat,xf_maxdat),
     .    xf_nepoch,  xf_isessn,
     .    xf_ietide, xf_isptide, xf_nclock, xf_jde, xf_jdr,
     .    xf_avlmet, xf_nslip, xf_nextra, xf_mdat, xf_prn(xf_maxsat)
      
      real*8  xf_ircint, xf_inter
      
      character   xf_typtime*4       
   
*   xf_islip(xf_maxcsb) - Epoch numbers of the extra bias flags
*   xf_slpst(xf_maxcsb)  - Satellite numbers for the extra bias
*               - flags (Given as PRN?)
 
      integer*2 xf_islip(xf_maxcsb), xf_slpst(xf_maxcsb)
 
*   xf_skd      - experiment type: S - static, K - kinematic,
*               - D - Dynamic ! not used
 
      character*1 xf_skd

*   xf_swver    - Version of the receiver software.
 
      real*4 xf_swver
 
*   xf_rcvrsw   - Receiver software type (e.g. ROG, TRM, GES) *6 ??
 
      character*3 xf_rcvrsw
      
*   xf_antmod   - Antenna phase center (NONE,ELEV,AZEL)

      character*4 xf_antmod
      
*   xf_frame    - Inertial ref. system (B1950 or J2000) ! no used
*   xf_precmod  - Precession model (IAU68 or IAU76)! no used
*   xf_nutmod   - Nutation model (IAU80) ! no used
*   xf_gravmod  - Gravity model for ARC ! no used
*   xf_sprmod   - Radiation-pressure model ! no used

      character*5 xf_frame, xf_precmod, xf_nutmod, xf_gravmod,
     .            xf_sprmod
     
      
*   xf_rctype   - Receiver type 
*   xf_rcvnum   - Receiver serial number
*   xf_anttyp   - Antenna type  
*   xf_antnum   - Antenna serial number 

      character*20 xf_rctype, xf_rcvnum, xf_anttyp, xf_antnum
      
*   xf_sitnam   - site name 
*   xf_makernam   -  maker name
*   xf_sitnam_l Long version of the site name (limited to
*               - 16 characters in most other places in GAMIT)
*   xf_openam   - RINEX Operator name/AGENCY
      character xf_sitnam*32, xf_openam*40, xf_makernam*32
 
      common / xf_type2_CH /xf_posflg,xf_typtime,
     .    xf_sitnam,xf_openam,xf_makernam,
     .    xf_frame, xf_precmod, xf_nutmod,
     .    xf_gravmod, xf_sprmod,  xf_skd, xf_rcvrsw, xf_antmod ,     
     .    xf_rctype, xf_rcvnum, xf_anttyp, xf_antnum,xf_siteinfo
     
      common / xf_type2_I4 /xf_srttime,xf_endtime,xf_srtdoy,
     .    xf_npart, xf_nsat, xf_ischan, xf_ndat,
     .    xf_dattyp, xf_lambda,  xf_nepoch,  xf_isessn,
     .    xf_ietide, xf_isptide, xf_nclock, xf_jde, xf_jdr,
     .    xf_avlmet, xf_nslip, xf_nextra, xf_mdat, xf_prn

      common / xf_type2_I2 / xf_islip, xf_slpst

      common / xf_type2_R8 / xf_pos,xf_neu,
     .    xf_srtsec, xf_endsec, xf_srtsect,
     .    xf_offarp, xf_offsL1, xf_offsL2, 
     .    xf_elvcut,xf_clock, xf_te, xf_tr,
     .    xf_serwgt, xf_extra,
     .    xf_inter, xf_ircint,
     .    xf_start_mjd, xf_end_mjd
     
      common / xf_type2_R4 /  xf_swver
*...................................................................
*     TYPE 3 Parameter block  ! not used at all  
 
*   xf_preval(xf_maxprm)    - Apriori values of the parameters used
*               - in the model run.  There are
*               - Element Parameter
*               -   1     geocentric Latitude (rad)
*               -   2     Geocentric Longitude (rad)
*               -   3     Radius (km)
*               -   4     Station clock offset (sec) -- This
*               -           value is not actually added to
*               -           Observed minus computed values
*               -   5     Station clock rate (sec/sec)
*               -   6     Station clock acceleration ((s/s)/sec)
*               -   7     Station atmospheric delay (Contains
*               -           zero currently)
*               - Looping over each satellite (n=0 to xf_nsat-1)
*               -   8+9*n IC X (km)
*               -   9+9*n IC Y (km)
*               -  10+9*n IC Z (km)
*               -  11+9*n IC Xdot (km/sec)
*               -  12+9*n IC Ydot (km/sec)
*               -  13+9*n IC Zdot (km/sec)
*               -  14+9*n Radiation Pressure 1 (0-1)
*               -  16+9*n Radiation Pressure 2 (km/sec)
*               -  17+9*n Radiation Pressure 3 (km/sec)
 
c      real*8 xf_preval(xf_maxprm)
 
*   xf_nlabel   - Number of parameter lables
*   xf_nparam   - Number of parameter whose apriori values
*               - are saved.
*   xf_islot(xf_maxlab) - Numbers giving the parameter types.
*               - (Definititions are defined in xf_rlabel) 
*   xf_idms(xf_maxlab)  - If this value is 1 then parameter
*               - should be converted to deg, min, sec; otherwize
*               - leave as is.          !  not used
 
c      integer*4 xf_nlabel, xf_nparam, xf_islot(xf_maxlab),
c     .    xf_idms(xf_maxlab)
 
*   xf_rlabel(xf_maxlab)    - These are the lables for the paramters
*               - The current definition, along with the
*               - parameter code number (xf_islot) are:
*               - Islot Label
*               -     1 Site latitude
*               -    11 Site longitude
*               -    21 Site height
*               -    31 Site atmosphere
*               -    41 Site clock epoch
*               -   101 Orbital Element 1
*               -   111 Orbital Element 2
*               -   121 Orbital Element 3
*               -   131 Orbital Element 4
*               -   141 Orbital Element 5
*               -   151 Orbital Element 6
*               -   161 Radiation Pressure 1
*               -   171 Radiation Pressure 2
*               -   181 Radiation Pressure 3
*               -   191 Along track
*               - These are not used in the SOLVG program.
 
c      character*20 xf_rlabel(xf_maxlab)
 
 
c      common / xf_type3 / xf_preval, xf_nlabel, xf_nparam, xf_islot,
c     .    xf_idms, xf_rlabel
 
*...................................................................
*     TYPE 4 Time tag block (one per epoch)
 
*   xf_jd       - JD Epoch of this measurement 
*   xf_mjd       - MJD Epoch of this measurement 
*   xf_sod      - Epoch of this measurement (UTC seconds of
*               - day) -- value recorded by ground station
*               - clock.  (Nominal time effectively)
*   xf_rclock   - Receiver clock offset (sec)
*   xf_rclk_s   - Sigma of Receiver clock offset (sec)
*   xf_save(xf_maxsav)  - Saved real*8 values

*   xf_antaz    - Antenna azimuth (degrees)
*   xf_kxyz(3)    - Kinematic site XYZ postion (m) 
*   xf_kvol(3)    - Kinematic velocity of site (m/s)
*   xf_kacc(3)    - Kinematic acceleration of site (m/s^2)
*   xf_kneu(3)    - Kinematic site latitude, longitude,radious (rad,rads, m) 
*   xf_L1ant(3)   - L1 antenna offset NEU (m) 
*   xf_L2ant(3)   - L2 antenna offset NEU (m) 

      real*8 xf_sod, xf_rclock,xf_rclk_s, xf_save(xf_maxsav),
     .   xf_kxyz(3), xf_kneu(3),xf_kvol(3), xf_kacc(3),
     .   xf_L1ant(3), xf_L2ant(3),  xf_antaz

      real*8   xf_jd, xf_mjd
     
*   xf_pres     - Pressure for this observation (mbar)
*   xf_temp     - Temperature (C)
*   xf_relhum   - Relative humidity (0-1)
 
      real*4 xf_pres, xf_temp, xf_relhum
 
*   xf_msat     - Number of satellites at this epoch
*   xf_iepoch   - Epoch number of this epoch
*   xf_iyr, xf_idoy - Year and day number of this epoch !! changed
*   xf_doy(2)   - Year and day number of this epoch
*   xf_okmet    - Validity of met at this time.  Met OK if bit
*               - set (see xf_avlmet)
*   xf_nsave    - Number of save variables
*   xf_kflag    - Flag for Kinematic x-file .  Meaning not
*               - clear.
 
      integer*4 xf_msat,  xf_iepoch,
     .    xf_doy(2), xf_okmet, xf_nsave, xf_kflag
 
*   xf_outflg    - written position option 'XYZ' or "NEU"
 
      character*4 xf_outflg

*   Options
    
*     xf_isposk - Position option
*		 0 - Static Record
*		 1 - Kinematic Record
*		 2 - Kinematic Record with velocity
*		 3 - Kinematic Record with velocity and acceleration

*     xf_isclk - Clock recording option
*     xf_ismet -  METEORIC recording option 
*     xf_isant -  ANTENNA recording option
 
      integer*4 xf_isposk,xf_isclk,xf_ismet,xf_isant
      common /xf_type41/ xf_isposk,xf_isclk,xf_ismet,xf_isant
      
      common / xf_type4 / xf_sod, xf_rclock,xf_rclk_s, xf_save,
     .    xf_kxyz, xf_kneu, xf_kvol,xf_kacc,
     .    xf_L1ant, xf_L2ant,xf_antaz, 
     .    xf_jd, xf_mjd,
     .    xf_pres, xf_temp, xf_relhum, xf_msat, 
     .    xf_iepoch, xf_doy, xf_okmet, xf_nsave, xf_kflag,
     .    xf_outflg
 
*...................................................................
*     TYPE 5  Data block for observation of each epoch  (in MPRN order)
*   xf_tau(xf_maxsat)     - total modeled delay (s)
*   xf_drate(xf_maxsat)    - total modeled delay rate (s/s)
 
*   xf_obsv(xf_maxdat,xf_maxsat)  
*               - Observed values (units depending of data type)
*               - Phases are in cycles at the appropriate freq.
*               - Ranges are in meters.
 
*   xf_omcs(xf_maxdat,xf_maxsat)  - Observed minus theoretical (units as
*               - above)
*               - All the OMC values are in cycles at the
*               - appropriate frequency.
*    xf_omcs is not used roight now, leave for future development 

*   xf_spare(xf_maxspr) - Spare real*8 values !!

 
      real*8  xf_obsv(xf_maxdat,xf_maxsat)
*      real*8 xf_tau(xf_maxsat), xf_omcs(xf_maxdat,xf_maxsat)
 
*   xf_elev(2,xf_maxsat)  - Elevation and rate of change (rad and rad/s)
*               - Later value not implemented yet.
*   xf_azimuth(2,xf_maxsat)   - Azimuth and rate of change (rad and
*               - rad/s) Later value not implemented yet.
*   xf_wvrdel(xf_maxsat)   - WVR delay (sec)
*   xf_atmdel(xf_maxsat)   - Atmospheric delay (sec)
*   xf_obswgt(xf_maxdat,xf_maxsat) - Observation weight
*   xf_ampL1(xf_maxsat)    - Signal amplitude at L1
*   xf_ampL2(xf_maxsat)    - Signal amplitude at L2
 
      real*4 xf_elev(2,xf_maxsat), xf_azimuth(2,xf_maxsat), 
     .       xf_wvrdel(xf_maxsat), xf_atmdel(xf_maxsat),
     .       xf_obswgt(xf_maxdat,xf_maxsat),
     .       xf_ampL1(xf_maxsat), xf_ampL2(xf_maxsat)

*   xf_data_flag - Data validity flag (AUTCLN, binary, see gobs_def.h)
*   xf_imet(xf_maxsat)  - Atmospheric delay flag
*                  0 -  no Atmospheric delay infor.
*                  1 -  zenith delay available only (m)
*                  2 -  WVR delay  available only (m)
*                  3 -  Both zenith and WVR delay  available

      integer*4 xf_data_flag, xf_imet(xf_maxsat)
 
*   xf_iprn(xf_maxsat)     - PRN number of this obervation
*   xf_okwvr(xf_maxsat)    - Indicates WVR data OK, -1 for none.
*   xf_ierfl(xf_maxsat)    - Error flag:
*               - Value meaning
*               -  (-11  Reweighted observation)
*               -   -2  Reweighted observation
*               -   -1  Unweighted data in Cview
*               -    0  Good observation
*               -    1  No observation
*               -    2  Deleted observation
*               -    3  Low amplitude observation
*               -    4  Low elevation angle observation
*               -    5  Not enough points for detection?
*               -    6  Outlier?
*               -   10  Bias flag needed at the epoch on this SV
*               -   98  Really OK (?)
*   xf_ndats(xf_maxsat)   - Number of data types at this epoch
*   xf_isnr(xf_maxdat,xf_maxsat)    - SNR ratio at this epoch
*   xf_ilck(xf_maxdat,xf_maxsat)    -  RINEX loss-of-lock indicator

 
      integer*2 xf_iprn(xf_maxsat), xf_okwvr(xf_maxsat),
     .    xf_ierfl(xf_maxsat), xf_ndats(xf_maxsat),
     .    xf_isnr(xf_maxdat,xf_maxsat),xf_ilck(xf_maxdat,xf_maxsat)
 

      common / xf_type5 /  xf_obsv, xf_elev, xf_azimuth,
     .    xf_wvrdel, xf_atmdel, 
     .    xf_ampL1, xf_ampL2, xf_obswgt, xf_data_flag,
     .    xf_iprn, xf_okwvr, xf_ierfl, xf_ndats, xf_isnr,
     .    xf_ilck,xf_imet
 
*---------------------------------------------------------------------
      character*156 xf_buffer
      common / xf_type6_com_CH /xf_buffer
