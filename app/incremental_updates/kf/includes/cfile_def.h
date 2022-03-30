*
*     This file contains the definition of the quanities in the
*     c-files.  It is used with the read_cf subroutines in the
*     Ghandlers library.

* MOD TAH 950622: Updated for new cfile format of 18 May 1995. Ver  930
* MOD TAH 980912: Updated for new cfile format of 15 Sep 1998. Ver  980
* MOD TAH 050206: Updated for new cfile format of 06 Feb 2005. Ver 1020
* MOD TAH 100902: Updated for new cfile format of 02 Sep 2010. Ver 1040
* MOD TAH 130119: Updated for new cfile format of 16 Jan 2013  Ver 1041
* MOD TAH 140327: Updated for new cfile format of 27 Mar 2014. Ver 1042
* MOD RWK 141209: Updated for new cfile format of  9 Dec 2014  Ver 1060
* MOD TAH 190521: Increased cf_maxlab to 33 for ECOMC model, maxorb to 22.
* MOD TAH 190711: Introduced copy of GLONASS frequencies for when c-file 
*                 is written and phase is remapped correctly.  Ver 1061
* MOD TAH 200126: Updated for svantdx being saved for L1/L2    Ver 1071

 
*     First set the paramaters.  These should be updated when the
*     cfile format is updated.
 
*   cf_maxtxt   - Maximum lines of text at the beginning of the
*               - cfile
*   cf_maxsat   - Maximum number of satellites allowed in cfile
*   cf_ncvsit   - Maximum number of channels allowed in a
*               - receiver
*   cf_maxprm   - Maxiumum number of parameters.  The clean
*               - routine set this max number of solve parameters
*               - We set to the actual max possible in a c-file
*   cf_maxlab   - Maximum number of labels for parameter strings
*   cf_maxdat   - Maximum number of data types.
*   cf_maxcsb   - Maximum number of bias parameters saved in
*               - c-file header
*   cf_maxext   - Maximum number of extra variables in type 2
*               - record
*   cf_maxsav   - Maximum number of save variables in type 4
*               - record
*   cf_maxspr   - Maximum number of spare varibales in type 5
*               - record
*   cf_maxclk - Maximum number of terms in the clock poynomials
*   cf_maxorb   - Maxium number of orbital partials by satellites
*                 (15 with Berne model on, 22 for ECOMC)
 
      integer*4 cf_maxtxt, cf_maxsat, cf_ncvsit, cf_maxprm, cf_maxlab,
     .    cf_maxdat, cf_maxcsb, cf_maxext, cf_maxsav, cf_maxspr,
     .    cf_maxclk, cf_maxorb 
 
* Parmeter values valid on 920812:
* Increased cf_maxsat to 32 931007:
* MOD TAH 940516:
* Increased cf_maxprm by 6 to account for the plar motion/UT1 offset
*         rate parameters
* Increaed maxlab to 22 for new partials. (This still seems to be
*         one more than needed).
* MOD TAH 950629:
* Increased number of partials per satellite from 9 to 15 for the
* new parameters (once-per-rev etc.)
*
* MOD TAH 970107:
* Added cf_maxorb to parameter list for max number of orbital partials.
* MOD TAH 980915: Increased number of orbit parameters to allow for 
*                 satellite offsets. (Changed from 15 to 18)
* MOD TAH 190524: Increased cf_maxorb from 18 to 22 for ECOMC
*     NOTE: Updates here should be made to ct_maxorb etc on c
 
      parameter ( cf_maxtxt =  50 )
* MOD TAH 200210: Increased from 32 to 35
* MOD TAH 201031: Increased from 35 to 45 for Beidou
      parameter ( cf_maxsat =  45 )
      parameter ( cf_ncvsit =  20 )
      parameter ( cf_maxorb =  22 ) 
      parameter ( cf_maxprm =  13 + cf_maxorb*cf_maxsat )
*                                     ! These are for site
*                                     ! position (3), clock offset,
*                                     ! rate and acc (3), atmosphere
*                                     ! and 15 orbital elments per SV
* MOD TAH 940516                      ! Added 6 for polar UT1 offset
*                                     ! and rate partials.
*                                     ! See cf_rlabel for definition.
* MOD TAH 95-812: Changed value to 26 to allow for 4 additional radiation
*                 parameters
* MOD TAH 980915: Added 3 to allow for satellite offsets.
* MOD TAH 190521: Added 4 to accomadte ECOMC model
      parameter ( cf_maxlab =  33 )
*                                     ! This values looks like a bug
*                                     ! (Cview routine use the number
*                                     ! of sites (ncvsit) for this
*                                     ! bound.  Actually there are
*                                     ! 15 parmeters types
* MOD TAH 940516                      ! Added 6 for polar motion/UT1
*                                     ! rate (to original 15, but this
*                                     ! is still one short of cfile
*                                     ! dimensioning.
      parameter ( cf_maxdat  =   5 )
      parameter ( cf_maxcsb  = 200 )
      parameter ( cf_maxext  =   2 )
      parameter ( cf_maxsav  = cf_maxsat )
      parameter ( cf_maxspr  =   2 )
      parameter ( cf_maxclk =  4 )
 
*-------------------------------------------------------------------
 
* Define each of the types of records in a c-file.
 
*     TYPE 1:

*   cf_nversn   - Version number. Same as model version multiplied
*                 by 100.
*                 Versions  930 autcln 3.07 and before
*                 Version   980 autcln 3.08 (both cfile versions can be
*                              read).
*                 Version 1020  autcln 3.22 (Reads all formats above)
*                 Version 1060  add GNSS variables
*   cf_ntext    - Number of text records at top of cfile
 
      integer*4 cf_nversn, cf_ntext
 
*   cf_text(cf_maxtxt)  - Header lines
 
      character*80 cf_text(cf_maxtxt)
 
*   cf_ntext    - Number of text records at top of cfile
*   cf_text - Header lines
 
      common / cf_type1 / cf_nversn, cf_ntext, cf_text
 
*...................................................................
*     TYPE 2: Header series data block
        
* MOD RWK 141209: Add frequencies
*   cf_fL1(cf_maxsat) - Higher L-band frequency
*   cf_fl2(cf_maxsat) - Lower L-band freuency
* MOD TAH 190711: Save copy of frequencies for GLONASS
*   sv_fL1(cf_maxsat) - Higher L-band frequency: Saved version
*   sv_fl2(cf_maxsat) - Lower L-band freuency: Saved version
*   cf_sec      - Start epoch seconds tag
*   cf_offarp(3)    - Up, North, East from mark to antena ref point (m)
*   cf_offsL1(3)    - Up, North, East L1 antenna offset (m)
*   cf_offsL2(3)    - Up, North, East L2 antenna offset (m)
C MOD TAH 200205: Added cf_antdaz as deviation from True N read from statoion.info
C                 with AntDAZ column.  Values from IGS log files
*   cf_antdaz   - Antenna deviation from True north
C MOD TAH 200126: Changed to cf_svantdx(3,2,maxsat) (added L1/L2 index) Ver 10.71
*   cf_svantdx(3,2,cf_maxsat) - Satellite antenna phase center offsets (m)
*   cf_elvcut   - Solve elevation cutof (degrees)
*   cf_clock (cf_maxclk) - Site clock polynomial (sec,sec/sec,..)
*   cf_te       - Seconds part of day for Ephemeris epoch
*   cf_tr       - Seconds part of day for Earth rotation
*   cf_UT1(2)   - UT1-AT and rate at start epoch  (sec, sec/day)
*   cf_xp(2)    - X-pole pos and rate at start (asec, sec/d)
*   cf_yp(2)    - Y-pole pos and rate at start (asec, asec/d)
*   cf_psi(2)   - Nutation in longitude and rate at start
*               - epoch (asec, asec/day. Rate not implemented)
*   cf_eps(2)   - Nutation in obliquity and rate at start
*               - epoch (asec, asec/day. Rate not implemented)
* MOD RWK 141209: Remove unused serwgt, overall weight to the session measurements
*   cf_atmlavg(3)    - Average N E U for atm loading          R*8
*   cf_hydrolavg(3)  - Average N E U for hydrolog. loading    R*8 
*   cf_extra(cf_maxext) - Extra real*8 values in Type 2 header

 
      real*8 cf_sec, cf_offarp(3), cf_offsL1(3), cf_offsL2(3),
     .    cf_antdaz, cf_svantdx(3,2,cf_maxsat), 
     .    cf_fL1(cf_maxsat),cf_fL2(cf_maxsat),
     .    sv_fL1(cf_maxsat),sv_fL2(cf_maxsat),
     .    cf_elvcut, cf_clock (cf_maxclk), cf_te, cf_tr,
     .    cf_UT1(2), cf_xp(2), cf_yp(2), cf_psi(2), cf_eps(2),
     .    cf_atmlavg(3), cf_hydrolavg(3), 
     .    cf_extra(cf_maxext)  
* MOD RWK 141209: Declaration needed for old-version code 
      real*8 cf_serwgt
 
*   cf_npart    - Number of partials stored in this cfile.
*               - (Could be 0, 5 or 15). Despite number here,
*               - cf_islot and cf_rlabel aways appear to be
*               - the same.  There are:
*               -  0 partials for a cfile with out partials.
*               -  5 partials when site pos, atmos, clock partials
*               - 15 partials when the 9 satellite orbit partials
*               -    are added.
* MOD TAH 980915: Added cf_norb for number of orbit partials
*   cf_norb     - Number of orbit partials (upto 18 for complete model).
*   cf_nsat     - Number of satellites in this cfile
*   cf_ischan(cf_maxsat)    - PRN numbers of the satellites.
*   cf_ndat     - Number of data types.
*   cf_dattyp(cf_maxdat)    - Data types saved in the c-file:
*               -   1  - L1 carrier phase (cycles)
*               -   2  - L2 Carrier phase (cycles)
*               -   3  - P1 Pseudo range (cycles @ L1) Pcode
*               -   4  - P2 Pseudo range (cycles @ L2) Pcode
*               -   5  - C1 Pseudo range (cycles @ L1) C/A code
*               -   6  - D1 L1 Doppler (Hz)
*               -   7  - D2 L2 Doppler (Hz)
*   cf_lambda(cf_maxsat,cf_maxdat)  - Wavelength factors for
*               - each of the datatypes on each of the
*               - satellites: Code is
*               -  0  - No data
*               -  1  - unambiguous, undoubled values
*               - -1  - ambiguous, undoubled values
*               -  2  - unambiguous, undoubled values
*               - -2  - ambiguous, undoubled values
*   cf_nepoch   - Number of epochs
*   cf_inter    - Obervation interval (secs)
*   cf_ircint   - Original receiver sampling interval before
*               - decimation
*   cf_mtime    - Type of time:
*               - 1  - UTC
*               - ?  - Others unknown
*   cf_isessn   - Session number
*   cf_iy, cf_im, cf_id, cf_ihr, cf_min - Start epocj ymdhms
*   cf_ietide    - Bitmapped tide model description.
*                 Bit    Model
*                   1    solid earth tides
*                   2    K1 frequency dependant earth tide
*                   3    Pole tide
*                   4    Ocean Tide
*                   5    Atm Tide (in conjunction with loading)
*   cf_isptide  - Bitmapped short period EOR Model flag
*                 Bit    Model
*                   1    UT1
*                   2    pole
*   cf_nclock   - Number of terms in the clock polynomial
*   cf_jde      - (PEP) julian day of ephemeris
*   cf_jdr      - (PEP) julian day of Earth rotation parameters
*   cf_avlmet   - Availability of met data:  Binary coded
*               - BIT   Meaning
*               -  1 - Pressure available
*               -  2 - Temperature available
*               -  3 - Relative humidity available
*               -  4 - WVR delay available (no longer supported)
*   cf_nslip    - Extra number of bias parameter flags
*   cf_nextra   - Number of extra real*8 values
*   cf_niextra  - Number extra I*4s (head.rec)
*   cf_ncextra  - Number of C*8 extras
*   cf_iextra(cf_maxext) - Extra I*4 variab'les

 
      integer*4 cf_npart, cf_norb, cf_nsat, 
     .    cf_ischan(cf_maxsat), cf_ndat,
     .    cf_dattyp(cf_maxdat), cf_lambda(cf_maxsat,cf_maxdat),
     .    cf_nepoch, cf_inter, cf_ircint, cf_mtime, cf_isessn,
     .    cf_iy, cf_im, cf_id, cf_ihr, cf_min, 
     .    cf_ietide, cf_isptide, cf_nclock, cf_jde, cf_jdr,
     .    cf_avlmet, cf_nslip, cf_nextra,
     .    cf_niextra, cf_ncextra ,cf_iextra(cf_maxext)
* MOD RWK 121416: declaration needed for old-version code
      integer*4 cf_iblk(cf_maxsat)
 
*   cf_islip(cf_maxcsb) - Epoch numbers of the extra bias flags
*   cf_slpst(cf_maxcsb)  - Satellite numbers for the extra bias
*               - flags (Given as PRN?)
 
      integer*2 cf_islip(cf_maxcsb), cf_slpst(cf_maxcsb)
 
*   cf_swver    - Version of the receiver software.
 
      real*4 cf_swver
 
*   cf_skd      - experiment type: S - static, K - kinematic,
*               - D - Dynamic
 
* MOD RWK 141128: Added cf_gnss for GNSS constellation
*   cf_gnss     - GNSS  ( G  R  C  E  J  I  )
 
      character*1 cf_skd,cf_gnss
 
*   cf_rcvrsw   - Receiver software type (e.g. ROG, TRM, GES)
 
      character*3 cf_rcvrsw
      
*   cf_antmod   - Antenna phase center (NONE,ELEV,AZEL)
*   cf_dryzen   - Hydostatic Zenith delay model          C*4  10.41
*   cf_wetzen   - Non-hydrostati Zenith delay model      C*4  10.41
*   cf_drymap   - Hydostatic mapping function            C*4
*   cf_wetmap   - Non-hydrostatis mapping function       C*4
*   cf_ionsrc   - 2nd order ion delay correction source  C*4  (10.42 ver)  
*   cf_svantmod(cf_maxsat) - SV antenna phase center model  C*4

      character*4 cf_antmod, cf_dryzen, cf_wetzen, cf_drymap, cf_wetmap,
     .            cf_ionsrc, cf_svantmod(cf_maxsat)
      
*   cf_frame    - Inertial ref. system (B1950 or J2000)
*   cf_precmod  - Precession model (IAU68 or IAU76)
*   cf_nutmod   - Nutation model (IAU80) 
*   cf_gravmod  - Gravity model for ARC 
*   cf_srpmod   - Radiation-pressure model 
*   cf_eradmod   - Earth-radiation model  (1042 version cfiles)
*   cf_antradmod - Antenna-radiation model (1042 version cfiles)

      character*5 cf_frame, cf_precmod, cf_nutmod, cf_gravmod,
     .            cf_srpmod, cf_eradmod, cf_antradmod

*  cf_magfield  - Name of magnetic field model for 2nd order ion (10.42 ver)
      character*6 cf_magfield

*  cf_speopmod  - Short-period EOP model (IERS92,IERS96) 
*  cf_etidemod  - E tide model (IERS92 IERS03)  	 
*  cf_otidemod  - Ocean tide model (OSO   NAO  )	    
*  cf_atmtide   - Atmospheric tides (ECMWF)		 
*  cf_atmlmod   - Atmospheric loading model & frame	  
*                 e.g., ECMWF CM
*  cf_hydrolmod - Hydrological loading model
*  cf_cextra (cf_maxext) - Spare C*8 space.

      character*8 cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide,
     .            cf_atmlmod, cf_hydrolmod, cf_cextra(cf_maxext)

*  cf_antmod_snx - Full name of antmod file (1040 version)
*  cf_svantmod_snx(cf_maxsat) - Full name of antmod file  (1040 version)

      character*10 cf_antmod_snx, cf_svantmod_snx(cf_maxsat)	    
     
*  cf_obfiln   - X-file name for this cfile
*  cf_tfiln    - T-file name for this cfile
*  cf_jfiln    - J-file name fpr this cfile
 
      character*16 cf_obfiln, cf_tfiln, cf_jfiln
      
*   cf_rctype   - Receiver type 
*   cf_rcvnum   - Receiver serial number
*   cf_anttyp   - Antenna type  
*   cf_antnum   - Antenna serial number
* MOD RWK 141209: Replace block number by the SV body description
*   cf_svantbody(cf_maxsat) - Body/antenna description ('BLOCK IIR, BEIDOU-2I , etc.)

      character*20 cf_rctype, cf_rcvnum, cf_anttyp, cf_antnum,
     .      cf_svantbody(cf_maxsat)
      
*   cf_sitnam   - Long version of the site name (limited to
*               - 16 characters in most other places in GAMIT)

 
      character*32 cf_sitnam
 
 
      common / cf_type2 / cf_sec, cf_offarp, cf_offsL1, cf_offsL2,
     .    cf_antdaz, cf_svantdx, cf_fL1, cf_fL2, sv_fL1, sv_fL2,
     .    cf_elvcut,cf_clock, cf_te, cf_tr,
     .    cf_UT1, cf_xp, cf_yp, cf_psi, cf_eps,
     .    cf_atmlavg, cf_hydrolavg, cf_extra,
     .    cf_npart, cf_norb, cf_nsat, cf_ischan, 
     .    cf_ndat, cf_dattyp, cf_lambda,
     .    cf_nepoch, cf_inter, cf_ircint, cf_mtime, cf_isessn,
     .    cf_iy, cf_im, cf_id, cf_ihr, cf_min, 
     .    cf_ietide, cf_isptide, cf_nclock, cf_jde, cf_jdr,
     .    cf_avlmet, cf_nslip, cf_nextra, cf_niextra, 
     .    cf_iextra, cf_ncextra, cf_islip, cf_slpst, 
     .    cf_swver, cf_antmod, cf_dryzen, cf_wetzen, cf_drymap, 
     .    cf_wetmap, cf_svantmod, cf_svantbody,
     .    cf_frame, cf_precmod, cf_nutmod, cf_gravmod, cf_srpmod,
     .    cf_eradmod, cf_antradmod, 
     .    cf_speopmod, cf_etidemod, cf_otidemod, cf_atmtide, cf_atmlmod,
     .    cf_hydrolmod, cf_cextra, 
     .    cf_antmod_snx, cf_svantmod_snx, 
     .    cf_skd, cf_rcvrsw, cf_obfiln, cf_tfiln, cf_jfiln, 
     .    cf_rctype, cf_rcvnum, cf_anttyp, cf_antnum, cf_sitnam,
     .    cf_ionsrc, cf_magfield, cf_gnss 
 
*...................................................................
*     TYPE 3 Parameter block
 
*   cf_preval(cf_maxprm)    - Apriori values of the parameters used
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
*               - Looping over each satellite (n=0 to cf_nsat-1)
*               -   8+9*n IC X (km)
*               -   9+9*n IC Y (km)
*               -  10+9*n IC Z (km)
*               -  11+9*n IC Xdot (km/sec)
*               -  12+9*n IC Ydot (km/sec)
*               -  13+9*n IC Zdot (km/sec)
*               -  14+9*n Radiation Pressure 1 (0-1)
*               -  16+9*n Radiation Pressure 2 (km/sec)
*               -  17+9*n Radiation Pressure 3 (km/sec)
 
      real*8 cf_preval(cf_maxprm)
 
*   cf_nlabel   - Number of parameter lables
*   cf_nparam   - Number of parameter whose apriori values
*               - are saved.
*   cf_islot(cf_maxlab) - Numbers giving the parameter types.
*               - (Definititions are defined in cf_rlabel)
*   cf_idms(cf_maxlab)  - If this value is 1 then parameter
*               - should be converted to deg, min, sec; otherwize
*               - leave as is.
 
      integer*4 cf_nlabel, cf_nparam, cf_islot(cf_maxlab),
     .    cf_idms(cf_maxlab)
 
*   cf_rlabel(cf_maxlab)    - These are the lables for the paramters
*               - The current definition, along with the
*               - parameter code number (cf_islot) are:
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
 
      character*20 cf_rlabel(cf_maxlab)
 
 
      common / cf_type3 / cf_preval, cf_nlabel, cf_nparam, cf_islot,
     .    cf_idms, cf_rlabel
 
*...................................................................
*     TYPE 4 Time tag block (one per epoch)

*   cf_sod      - Epoch of this measurement (UTC seconds of
*               - day) -- value recorded by ground station
*               - clock.  (Nominal time effectively)
*   cf_rclock   - Receiver clock offset (sec)
*   cf_zendel   - A priori zenith delay at epoch
*   cf_save(cf_maxsav)  - Saved real*8 values
* MOD RWK 141209: Remove the kinematic 'k' from site coordinate names
*   cf_latr     - Spherical latitude (rads)
*   cf_lonr     - Longitude (rads)
*   cf_radius   - Radius (km)
*   cf_L1Z      - L1 antenna offset up (m?)
*   cf_L1N      - L1 antenna offset North (m?)
*   cf_L1E      - L1 antenna offset East (m?)
*   cf_L2Z      - L1 antenna offset up (m?)
*   cf_L2N      - L1 antenna offset North (m?)
*   cf_L2E      - L1 antenna offset East (m?)
*   cf_antaz    - Antenna azimuth (degrees): Not used apparly
 
      real*8 cf_sod, cf_rclock, cf_zendel, cf_save(cf_maxsav),
     .    cf_latr_sph, cf_lonr, cf_radius, 
     .    cf_L1Z, cf_L1N, cf_L1E, cf_L2Z, cf_L2N, cf_L2E,
     .    cf_antaz
* MOD RWK 141209: Declaration needed for old-version code
      real*8 cf_klatr, cf_klonr, cf_kradk 
 
*   cf_pres     - Pressure for this observation (mbar)
*   cf_temp     - Temperature (C)
*   cf_relhum   - Relative humidity (0-1)
*   cf_atmlod(3) - Atmospheric pressure loading (NEU meters)
 
      real*4 cf_pres, cf_temp, cf_relhum, cf_atmlod(3) 
 
*   cf_msat     - Number of satellites at this epoch
*   cf_mprn(cf_maxsat)  - PRN numbers of the satellites at
*               - this epoch
*   cf_iepoch   - Epoch number of this epoch
*   cf_iyr, cf_idoy - Year and day number of this epoch
*   cf_okmet    - Validity of met at this time.  Met OK if bit
*               - set (see cf_avlmet)
*   cf_nsave    - Number of save variables
* MOD RWK 141209: Remove unused kflag
 
      integer*4 cf_msat, cf_mprn(cf_maxsat), cf_iepoch,
     .    cf_iyr, cf_idoy, cf_okmet, cf_nsave
* MOD RWK 141209: Declaration needed for old-version code
      integer*4 cf_kflag 
       
* MOD RWK 141209: Rename ksite to the new name
*   cf_sitecd    - Site occupation code for kinematic
 
      character*4 cf_sitecd
* MOD RWk 141209: Declaration needed for old-version code
      character*4 cf_ksite 
 
 
      common / cf_type4 / cf_sod, cf_rclock, cf_save, 
     .    cf_latr_sph, cf_lonr, cf_radius, 
     .    cf_L1Z, cf_L1N, cf_L1E, cf_L2Z, cf_L2N, cf_L2E, cf_antaz, 
     .    cf_zendel, cf_pres, cf_temp, cf_relhum,cf_atmlod, cf_msat,
     .    cf_mprn, cf_iepoch, cf_iyr, cf_idoy, cf_okmet, cf_nsave, 
     .    cf_sitecd, cf_kflag
 
*...................................................................
*     TYPE 5  Data block for each satellite (in MPRN order)

*   cf_svcepc   - Epoch offset for satellite (s)
*   cf_svcl1    - Non-offset polynomial contribution to satellite
*                 clock (this part is removed from the phase and
*                 needs to be added back to get full satellite clock
*                 error (L1 cycles). 
*   cf_tau      - total modeled delay (s)
*   cf_drate    - total modeled delay rate (s/s)
*   cf_obsv(cf_maxdat)  - Observed values (units depending on
*               - data type)
*               - Phases are in cycles at the appropriate freq.
*               - Ranges are in meters.
 
*   cf_omcs(cf_maxdat)  - Observed minus theoretical (units as
*               - above)
*               - All the OMC values are in cycles at the
*               - appropriate frequency.
*   cf_spare(cf_maxspr) - Spare real*8 values
*   cf_tmpart(cf_maxlab)    - Partial derivatives.  These are in the
*               - order of ISLOT.  For full solution,
*               - Entry  Type
*               -     1  dL1/dlat (cycles/rad)
*               -     2  dL1/dlong (cycles/rad)
*               -     3  dL1/radius (cycles/km)
*               -     4  dL1/atm del (cycles/cm)
*               -     5  dL1/clock offset (cycles/sec) (Almost
*               -           equals the frequency (except for
*               -          epoch offset affects
*               -  6-14  dL1/Orbital elements (cycles/km,
*               -           cycles/(km/sec), and cylces/unitless)
*               -    15  dL1/along track (NOT IN DATA FILES)
 
      real*8 cf_svcepc, cf_svcL1,
     .    cf_tau, cf_drate, cf_obsv(cf_maxdat), cf_omcs(cf_maxdat),
     .    cf_spare(cf_maxspr), cf_tmpart(cf_maxlab)
 
*   cf_elev(2)  - Elevation and rate of change (rad and rad/s)
*               - Later value not implemented yet.
*   cf_azimuth(2)   - Azimuth and rate of change (rad and
*               - rad/s) Later value not implemented yet.
* MOD RWK 141209: Add SV nadir angle
*   cf_nadang   - Nadir angle of SV antenna, radians 
* MOD RWK 141209: Remove unused WVR delay 
*   cf_atmdel   - Atmospheric delay (sec)
* MOD RWK 141209: Remove unsued observation weight obswgt)
*   cf_ampL1    - Signal amplitude at L1
*   cf_ampL2    - Signal amplitude at L2
 
      real*4 cf_elev(2), cf_azimuth(2), cf_nadang,  cf_atmdel,
     .    cf_ampL1, cf_ampL2 
* MODE RWK 141209: Declaaration needed for old-version code
      real*4 cf_wvrdel, cf_obswgt(cf_maxdat) 

*   cf_data_flag - Data validity flag (AUTCLN, binary, see gobs_def.h)

      integer*4 cf_data_flag
 
*   cf_iprn     - PRN number of this obervation
* MOD RWK 141209: Remove unused WVR flag 
*   cf_ierfl    - Error flag:
*               - Value meaning
*               -  -11  Reweighted observation
*               -   -1  Unweighted data in Cview
*               -    0  Good observation
*               -    1  No observation
*               -    2  Deleted observation
*               -    3  Low amplitude observation
*               -    4  Low elevation angle observation
*               -    5  Not enough points for detection?
*               -    6  Outlier?
*               -    7  Missing antenna model
*               -    8  SV missing from T-file.
*               -   10  Bias flag needed at the epoch on this SV
*               -   98  Really OK (?)
*   cf_ndats    - Number of data types at this epoch
*   cf_isnr(cf_maxdat)    - SNR ratio at this epoch
*   cf_nspare   - Number of spare values
*   cf_nparts   - NUmber of partial derivatives
 
      integer*2 cf_iprn, cf_ierfl, cf_ndats,
     .    cf_isnr(cf_maxdat), cf_nspare, cf_nparts
* MOD RWk 141209: Declaration needed for old-version code
      integer*2 cf_okwvr 
 

      common / cf_type5 / cf_svcepc, cf_svcL1, 
     .    cf_tau, cf_drate, cf_obsv, cf_omcs, cf_spare,
     .    cf_tmpart, cf_elev, cf_azimuth, cf_atmdel, cf_nadang,
     .    cf_ampL1, cf_ampL2, cf_data_flag,
     .    cf_iprn, cf_ierfl, cf_ndats, cf_isnr, cf_nspare,
     .    cf_nparts
 
*---------------------------------------------------------------------
 
