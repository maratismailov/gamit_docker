      SUBROUTINE MVERSN(VERS)

      CHARACTER*40 VERS
      CHARACTER*10 GETMAC,MACHIN
      character*45 libver
      integer*4 ierr, nblen
C.....
      MACHIN = GETMAC(1)
      WRITE (VERS,5, iostat=ierr)  MACHIN(1:nblen(machin))
    5 format ('10.72 2021/07/30 21:56 UTC (',a,')')   
      if( ierr.ne.0 ) call report_stat('WARNING','MVERSN','MODEL',
     .                                  vers,'IOSTAT ERR',ierr)
c    5 format ('%I%of %E% %U% (',a,')')
c     The above line contains magic SCCS keywords which will be
c     automatically updated by delta and get.
c     Please do not change them.

c     get library version
      CALL LVERSN(libver)
                  
C VERSION V.002   LAST LSI VERSION
C VERSION V3.10   IBM PC AT VERSION - 1/26/87
C                 FOUND ERROR IN OUTPUT IN DEFINITION OF FREQTR,
C                  MOVED TO MAIN PROGRAM - 1/24/87
C                 ADD COMMON ATMAX2 TO MAIN AND ATMRED AND
C                  MODIFICATIONS TO ATMRED - 1/26/87
C                 NOTE : USES OLD PRECESSION, PRCES.FOR
C                  NEW PRECESSION IN PRECES.IAU
C                 NOTE : RICK SAYS ANTENNA OFFSETS NOT DONE PROPERLY
C                  IN CORTYP FOR ANYTHING OTHER THAN L1 OFFSET
C                 MODIFIED COMPAR, OBSRED, OUTPUT & ITERAT
C                  C'S ARE MULTIPLIED BY A FACTOR OF 2 FOR MACROMETER
C                  DATA SINCE THE OBSERVABLE IS HALF-CYCLES,
C                  O'S ARE ONLY MODIFIED WHEN A FULL-CYCLE RECEIVER
C                  (E.G. TI-4100) IS MIXED WITH A HALF-CYCLE
C                  RECEIVER (MACROMETER), THUS C'S ARE ALWAYS
C                  MULTILPIED BY 2 FOR HALF-CYCLE INSTRUMENTS
C                  SO THERE IS NO PROBLEM ON ITERATIONS
C                 MINOR MODIFICATIONS IN DAINT,ARGMNT,UT1TID -1/29/87
C                 REMOVED DOUBLE THE OBSERVED PHASE OPTION FROM
C                  ROUTINE COMPAR - 2/12/87
C                 ADDED ANTENNA OFFSETS TO S-FILE, ROUTINES SSPOOL,
C                  POSDEF,REDURT,ITERAT - 2/16/87
C                 ADDED TROPOSPHERE AND WVR MODELS - SEE DOCUMENTATION,
C                  IMPLEMENTED BY JIM DAVIS AND YEHUDA BOCK - 2/16/87
C                  ROUTINES COMPAR,OUTPUT,ATMDEL,WTHHGT, ADDED
C                  DCHAO,WCHAO,WVRRED,ZENTMP,CFA,FFUN,MARINI,SAASZD,
C                  WPRESS
C                 MODIFIED S-FILE FORMAT TO INCLUDE SEPARATE RECORDS
C                  FOR THE ANTENNA OFFSETS - 2/17/87
C VERSION 3.11    SENT TO AERO - 2/18/87
C                 REMOVED TABS FROM CONVRT - 2/24/87
C                 REMOVED TABS FROM ATMRED - 3/3/87
C                 REMOVED ANTENNA HEIGHT FROM WTHHGT - 3/3/87
C                 LINK GDETIC.WGS TO CHANGE TO WGS72 ELLIPSOID - 3/6/87
C VERSION 3.12    SENT TO AERO - 3/6/87
C VERSION 3.20    THIS WILL INCORPORATE BOB'S CHANGES
C                 COMPILED PRCES.IAU AND SIDTIM.IAU - 3/12/87
C VERSION 3.21    RWK changes to simplify code and add the ability to
C                 perform calculations in an earth-fixed reference
C                 frame.  Begun 3/13/87. Continued 6/7/87
C                 Remove dimensions from all site quantities. Changes
C                 to 40 + routines
C VERSION 3.22    New C-File format.  Changes to COMPAR,SETUP,CHDRED
C                 CHDRED,WRTHED,OBSRED,CFOUT 6/24/87
c version 4.1     first version in sccs -kf 871102
c version 4.2     Thursday, January 28, 1988 kurt
c                 change UT1RED,POLRED,OPEN, so that tables
c                 can have general format, and so that they
c                 are assumed to be (links) in the local directory.
c                 Also add standard TAIUTC which reads LEAP.SEC
C version 4.3     Modify UT1RED, POLRED to read tables with
c                 only 2 lines of headers (unmodified from
c                 Harvard!)
c version 4.4     Thursday, February 11, 1988 Kurt
c                 Modify UT1RED, POLRED to include the magic
c                 multiplier.
c version 4.5     February 26, 1988 Yehuda
c                 Modify COMPAR for new error flag (=2)
c                 Modify COMPAR to fix bad logic in new error flag
c version 6.0     New release for the entire set of modules
C                  New error flags and logic -- MHM 880331
C                  OBSRED now determines low amplitude data (3)
C                  COMPAR now determines low elevation and epoch limited data (4)
C                 Satellite delete bug in SETUP - YB 3/6/88
C
C                 CHANGE NAME OF MODULE FROM 'COMPAR' TO 'MODEL'
c
c version 6.1     New GAMIT Release -- 880418 MHM
c version 6.2     Add earth tide effects on stations.  Changes to
c                  SITCOR, UT1TID.  New routine ETIDE.  RWK 4/29/88
c version 6.3     POLRED:  Fix format statement in POLE. head reader -- MHM 880501
c version 6.4     Add INCLUDE dimpar.fti and make variable dimensions
c                  for satellites, parameters, and observables
c                  Remove references to R files. Changes to MODEL
c                  CHDRED, OBSRED, PARTL, SATEL, SETUP, THDRED
c                  WRTHED.  Remove RHDRED.  Compile these and 2 others
c                  (CFOUT, XHDRED) with MAXSAT=7, MAXORB=9 ,MAXLAB=15 -- RWK 5/1/88
c                 Add computation and writing to C-file of pseudo-
c                  ranges.  Changes to MODEL, CFOUT, OBSRED, & SETUP
c                  Recompile these and 6 others (CHDRED, PARTL, SATEL
c                  THDRED, XHDRED, WRTHED) with MAXDAT=4 -- RWK 6/6/88
c                 Increase precision of time arguments for sidereal time by
c                  using JD, FRACT (or T seconds of day) instead of R*8 FJD.
c                  Reckon TRUN from IC epoch vice IC epoch minus one.  Changes
c                  in MODEL, SITCOR, ROTURT, SIDMAT, SIDTIM, & SATEL;
c                  add subroutine TIMINC from ARC -- RWK 6/13/88
c                 Copy changes to SCCS and COM directories -- RWK,MHM 880726
c version 6.5     Remove MODEL sat selection option (DRIVER also changed)
c                  Change to SETUP -- RWK 880727
c version 6.6     Correct bug in time epoch argument in MODEL -- RWK 880803
c version 6.7     Correct bug in dimension check in SETUP -- RWK 880816
c version 6.8     Allow different low-amp cutoff for TI and MAC II. Temporary
c                  fix using DATFAC as check in MODEL and OBSRED -- RWK 881010
c version 6.9     Fix bug in GDATUM to have WGS84 used for atm calculations
c                  when IDATUM=0 (geocentric coords) -- RWK 881101
c version 7.0    Fix argument calls (bug) for use of met tables.  Changes in
c                  in MODEL, ATMRED, WVRRED, ATMDEL.  Fix METMOD to read models
c                  when there is an input table -- RWK 881230
c                 Fix bug in GDATUM--used wrong (NAD27) 1/f with
c                  WGS84  a and SHFT --  RWK 890318
c version 7.1    Integer*4, Case Insensitive release -- MHM 890818
c                  Changes to ATMDEL,ATMRED,CFA,METMOD,SAASZD,SREAD,THDRED,WHDRED
c                  C-file READ/WRITE primitives added, dattyp,lambda replace
c                    ldual in CFOUT,CHDRED,OBSRED,SETUP,WRTHED,XHDRED,ZHDRED,OPEN
c                  Explicit L2 Phase center calculations added to CORTYP,SITCOR
c                  Replace DAINT with system function DINT in FUNCOF,NUTRED
c                  Clean logic, uninitialized variable in GDETIC,POLRED,PRCES,UT1RED
c                  All of the above, plus azimuth calculations added to MODEL
c version 7.2    Fix bug in failure to save character format array in POLRED,
c                  by SAVE-ing it.  At the same time, change all COMMON
c                  variables in POLRED, UT1RED, and NUTRED to SAVEs. -- rwk 891215
c                Read observation times explicitly from the X- or C-files, and
c                  compute receiver clock corrections epoch by epoch, reading a
c                  J-file.  Write the correction calculated from each satellite's
c                  SVDTs and ephemeris on the output C-file, (temporarily)
c                  using the SAVE vector.  Changes to MODEL,
c                  OPEN, SETUP, OBSRED, CFOUT, PARTL; new routine SCLOCK to
c                  compute SV clock corrections.  -- rwk 891226-900105.
c                Include errflg.fti from library, adding IER=7 as a receiver
c                  hardware flag (replacing IER=5, which had be used for
c                  "too few points" by AUTCLN).  Changes to MODEL -- rwk 900103
c                Write the MODEL version number, operator, and date to the
c                  C-file header.  Changes to MODEL, WRTHED, OPEN, SETUP. --rwk 900103
c                Fix bug in calculating clocks from polynomials:  the units
c                  scaling for the quadratic term was wrong for DRIVER/S-file,
c                  effectively causing the term to be ignore and injecting
c                  clock errors of tens of microseconds if a quadratic model
c                  was selected in DRIVER. Change to MODEL. -- rwk 900104.
c version 7.3    Fix bug in converting epoch times read from X-file in OBSRED.--rwk 900120.
c version 7.5    1) Back port from suns: OPEN, SETUP
C                2) Use library routines: PROPER and CROSS                  are in the library.
C                Kurt
c version 7.6    Change obsred to allow 0.3 seconds of slop. kurt,king 900309
c version 7.7    Yehuda's changes
c                Modify SREAD,SETUP&CORTYP to handle new station coordinate
c                 convention
c                Yehuda Bock 4/29/90 and K. Feigl to preserve backwards
c                 compatibility.
c version 7.8    Minor clean-ups.  K. Feigl 900510
c                MODEL: check RCLOCK for reasonability.
c                CORTYP: more comments for clarity.
c                OBSRED: clean up.
c                SREAD: clean up and debug.
c                UT1RED: better error messages
c version 7.9    Include Satellite clock offset (SVDT) in light-time calculation.
c version 7.10   Remove the change in 7.9.
c version 7.11   CLKACC now read from S-file in units of 1/second.
c                Modifications to SREAD and MODEL.
c                Old version of S-files no longer compatible.
c                The old S-file were (incorrectly) labelled ACCELERATION (1/DAY)
c                The units in the old S-file were actually day/s/s!  BAD! BAD!
c                The CLKTRM was never correct.
c                K. Feigl June 16, 1990.
c version 7.12   The following routines have been changed
c                so that dimpar.fti is in lowercase
c                 for SR10.2 compatibility:
c                CHDRED,OBSRED,PARTL,SATEL,SETUP,THDRED,WRTHED
c                Yehuda Bock September 12, 1990
c                Add recl=8000 for T-file open statement in OPEN
c                Copied Seiichi's changes to ATMRED, MIT version copied
c                 to ATMRED.OLD
c                Yehuda Bock September 14, 1990
c                Modify MODEL to improve W-file read efficiency.
c                Seiichi Shimada September 14, 1990
c                In OPEN, eliminate redundant Y/N questions and
c                 repair "jnone" option (no J-file)
c                In SETUP, move call to READJ and make contingent
c                 on IRCLK=2
c                In MODEL, bypass correction to computed phase
c                 due to satellite oscillators (second term), and
c                 ignore effect of satellite oscillators in the
c                 third term (compare to MODEL.OLD)
c                Yehuda Bock September 16, 1990
c version 7.13   In SETUP, call to READJ any time
c                Seiichi Shimada  Sept. 22, 1990
c                In SETUP, reinstate call to READJ contingent on
c                 IRCLK=2
c                Yehuda Bock, Sept. 23, 1990
c version 7.14   Kurt Oct 11, 1990 Homogenize Scripps and MIT:
c                All routines now have lower case include files.
c                Routines modified by Scripps:
c                  ATMRED: use Scripps version without changes.
c                  CHDRED: New C-file format supercedes Scripps version.
c                          lower case name for include file
c                  MODEL:  Read W-file only once per epoch as Seichi suggests.
c                          allow option to run without a J-file
c                  OBSRED: New C-file format supercedes Scripps version.
c                          lower case name for include file
c                  OPEN:   Use Scripps version without changes.
c                  SETUP:  allow option to run without a J-file
c                          lower case name for include file
c                  THDRED: lower case name for include file
c                  PARTL:  lower case name for include file
c                  WRTHED: lower case name for include file
c                Kurt Oct 19-23: little bugs:
c                  MODEL,CFOUT,OBSRED: obswgt should be dimensioned MAXSAT x MAXDAT
c                Kurt Oct 25:
c                  ETIDE:  return updated value of EVEC ! (10 cm effect!)
c                  CORTYP: remove extraneous variables
c                  OPEN:   remove extraneous variables
c                  POLRED: remove extraneous variables
c                Kurt Nov. 6:
c                  MODEL:  require J-file, KLOCK = 1,2 or 3
c                  SETUP:  require J-file, KLOCK = 1,2 or 3
c version 8.5    Kurt November 15
c                  correct initialization in MODEL
c Last MIT chang - up to here
c version 8.51   In SETUP only call READJ for option 3
c                In MODEL, skip second term computation for KLOCK=1,2
c                  Yehuda Bock November 20
c version 8.52   MAXNET=30, MAXSAT=18 in dimpar.fti
c                  Yehuda Bock 12/23/90
c version 8.53   MAXEPC=1500 in MAKEX.FTI
c                  Yehuda Bock 12/25/90
c version 8.54   In MODEL, modification in error flag check
c                  Yehuda Bock 12/28/90
c version 8.55   In MODEL & SETUP, use a J-file whenever it
c                  exists and for whatever receiver (Mini-Mac, too)
c                Yehuda Bock 12/30/90
c version 8.56   Clean up MODEL so that it will work with JNONE.
c                  This can produce very confusing results if differencing
c                  between 2 C-files, one produced with a J-file, one without.
c version 8.56   Add the offset between the SV phase center and the
c                  center of mass. Changes to MODEL.  king 910302
c                Modify the SV oscillator frequency terms for phase and
c                   PR to handle SA corrections in MODEL.  Kurt 910305.
c                Add printout of library version to MVERSN.  King 910307
c                Remove extraneous variables (IFRSTN etc) in OPEN  - king 910307
c                Read a cubic clock term in SETUP.  Kurt 910316
c version 8.7    SREAD: read a cubic clock term in SREAD.
c version 8.8    OPEN : use status='unknown' for P-file
c                       Yehuda Bock 4/2/91
c                MODEL: if RCLOCK is > 2000 s (truly obscene!),
c                       then set RCLOCK to zero, AND
c                       put the error flag to igunwt
c                       Kurt 910410
c version 8.9    SVANT: Attempt fix of SV offset for colinar sun, probe, earth
c                       King 910414
c version 8.10   Write S-file (L-file) sitnam, not C- or X-file sitnam on
c                   output C-file. Changes to SETUP, SREAD.  King 910510
c version 8.11   Fix problem in last change.   King 910510
c
c version 8.12   New cpu mathlib_sr10 Yehuda Bock, maxsat=15 910523
c                MODEL & OPEN:  Use 'Are you sure?' query to keep or
c                    delete the input C-file.  Bock/King 910701
c version 8.14   CHDRED: change output for antenna heights by one decimal place
c                Yehuda/Keith 910710
c version 8.15   Changes to add earth rotation values and their derivative to
c                    the C-file header. Also get met values in SETUP instead of
c                    MODEL.  Changes to MODEL, SETUP, SITCOR, UT1RED, POLRED,
c                    METMOD,WHDRED, ZHDRED.   King  911205
c version 8.16   Incorporate Scripps changes for kinematic option:
c                  Changed OBSRED to read new X-file epoch header
c                  In MODEL, get station coordinates from L-file
c                  Changed SETUP for X-file backwards compatibility
c                  Changed OPEN to read kinematic/static flag
c                  New subroutine LREAD reads an L-file
c                     yb & jfg 10/09/91
c                  Format changes in OPEN, fixed bug in MODEL
c                     yb 10/17/91
c                  Modify call to XHDRED in SETUP (read XTYPE parameter)
c                     yb 10/18/91  (changed XHDRED in library)
c                  MODEL: Always read coordinates from a l-file (for X/C file input)
c                     Add lowerc check for lat/lon flags    yb 11/27/91
c                  Add call to NUTTAB in SETUP to get nutation angles (not yet rates)
c                     Change 'phi' to 'psi' there and in WRTHED     King 911218
c                  Check UT1 type before adding tidal terms.  Changes to UT1RED
c                     SETUP, WRTHED, SITCOR.     King 911220.
c                  MODEL:  Avoid second call to CORTYP for static.   King 911221
c                  Replace SATEL with new version from orbits directory.
c                     Changes to MODEL, SETUP, THDRED. King 911221
c                     Call changed to library GSATEL.  King 920127
c                  Changes to use lfile., I-file and sited. file instead of the S-file
c                     for coordinates, receiver clock polynomial coefficients, and
c                     antenna offsets.  New routines READI and RSITED (both to be put in
c                     library eventually); changes to MODEL,SETUP, OPEN, and library LREAD.
c                  Remove TIMINC from Makefile (in library).
c                  Remove commons and extraneous INITPR from MODEL.
c                  Replace MATMPY (bug) with /orbits version moved to lib. King 920130
c                  Move NUTRED, POLRED, UT1RED, UT1TID, FUNARG, FUNCOF, SIDTIM
c                     to library.  Replace PRCES, NUTTAB, SITMAT with /orbits version
c                     moved to library (sense of matrices changed).   Changed Makefile.
c                     King 920131
c version 8.17     SETUP, OPEN: Read (in preference) station.info, sited., x-file for
c                     antenna offset information.   King 920324
c                  READI:  Check session number.   King 920324
c                  RSITED: Change yr,day,session tests to integer.   King 920324
c                  OBSRED: Fix check on X-file type.   King 920324
c version 8.18     SETUP:  Fix bug in call to NUTTAB.   King 920328
c                  RSITED: Fix bug in checking site name.  King 920328 / Dong 920331
c                  OBSRED, MODEL, SETUP:  Make xtype, kinflg checks consistent.  King 920331
c                  Recompile for modkin.fti change:  MODEL, OBSRED, SETUP, OPEN
c                       Makefile:  add modkin.fti to (OBJECTS) dependency.  King 920331
c                  SVANT:  Correct phase-center z-offset for BLK II SVs.    King 920423
c                  New C-file format:  Add SKD and IRCINT to READC2 and WRITC2 calls
c                       in CHDRED and WRTHED; replace KINFLG by SKD throughout; changes
c                       in MODEL, OBSRED, OPEN, SETUP.  King 920423
c                  Makefile:  Remove DSCRM, TTOPOS, WOBBLE  (no longer used).
c                  Explicit declarations: ATMDEL, ATMRED, CFA, CORTYP, GDATUM, GDETIC
c                       HEIGHT, IARRAY, MARINI, PARTL, SAASZD, SHFTOR, SREAD, WVRRED
c                       ZENTMP.
c                  Eliminate unused variables:  MODEL, OPEN, RSITED, SITCOR
c                       SVANT, WRTHED.    King 920424
c                  Add session number to X- and C-file.  Changes to SETUP, CHDRED,
c                       WRTHED.   King 920429
c version 9.1      Compiled at SIO by Bock 920503
c version 9.11     SETUP: backwards compatibility for IRCINT variable, remove old debug
c                       Bock 920512
c                  SETUP : fix call to readi BOCK 920514
c                  CHDRED : Spelling correction
c                  SETUP  : get year right after CHDRED
c                  WRTHED : correct output for IRCINT  Bock 920515
c version 9.12     SETUP  : slight change in logic for year.   King 920519
c version 9.13     SETUP  : set antenna height parameters and coordinates
c                           (found in modkin.fti) Bock 920525
c version 9.14     MODEL, OPEN : define iscrn in MODEL.   Bock/King 920608.
c         9.15     OPEN   : Fix cosmetic but politically important printout. King 920626
c         9.16     Makefile: Move THDRED and IARRAY to library.  King 920914.
c                  Remove unused variables:  ROTURT, SITCOR, SREAD.  King 920925
c         9.17     Makefile: Change ranlib order in Sun/HP version.  Feigl/King 921015
c                  AVCLCK : Change g to f format to avoid Sun compiler bug.  Fang/King 921019
c         9.18     SETUP : Add arguments to HISUB call.  King 930211.
c                  MODEL : Fix print statements for HP compatibility. Vigny/King 930224
c         9.19     MODEL : Don't reflag as low-el unweighted data.  Herring/King 930310
c         9.20     READI : Shorten printout line to accomodate DEC.  Oral/King 930413
c                  WRTHED: Increase format for # epochs.  Oral/King 930413
c         9.21     Mutiple zenith delay code:
c                   MODEL,SETUP,METMOD
c                  Bock 930408
c                  Undid zenith delay changes (comments left in place)
c                  Bock 930722
c         9.22     Remove zenith delay comments.   King 930818
c                  AVCLCK: Give BLK II's and BLK I's equal weight; print whenever
c                          a satellite is not used.  Herring (930823)/King 930914.
c                  MODEL, SETUP, WRTHED: Compute and save zenith delay value in
c                         parameter array on C-file header (move CORTYP call and
c                         add ATMDEL call in SETUP).  Remove check on elev = 0.
c                         in ATMDEL.   King 930914
c                  WRTHED, PARTL: Change units of atmospheric delay and partials
c                         from centimeters to meters.   King 930914
c                  AVCLCK: Actually implement 930914 changes (missed).
c                  ROTURT: Add sidereal time to arguments of (lib) SIDMAT.
c                          King/Morgan  931013
c         9.23a    (merge with 9.23 test version) SETUP: Fix elevation cutoff
c                  default from 14.84 to 15.0 deg.   King 931208
c         9.23     Updated Solid Earth Tide, Ocean Loading, Short period Earth
c                  orientation and Antenna phase centre Modeling.
c                  McClusky 931203
c                  ********New Earth Tide Modeling********
c                  MODEL: modified to update etide models, and to only call
c                         the sitcor routine once per receiver clock calc pass.
c                  SITCOR: modified to accept sidereal time, sidtm, from sidmat
c                          by way of roturt as it is needed in the computation
c                          of the earth tides done in etide.
c                          Also the logical unit number of the ocean tide file,
c                          IOTIDE,is passed so that Ocean tide file can be read.
c                  GDATUM: modified to read bitmapped tide model selection.
c                  OPEN: modified to open ocean loading file.
c                  ETIDE:modified and reorganised to compute IERS soild earth
c                        tides to 2nd order, frequency dependant K1 solid earth
c                        tide, pole tide, and, ocean loading dispalcements. All
c                        tide options selectable using bit mapped variable
c                       "itidety" scm.
c                  PTIDE: subroutine added to compute pole tide displacements
c                         pjm.
c                  OTIDE: subroutine added to compute ocean loading
c                         displacements pjm.
c                  SETUP: modified to read bitmapped tide setections.
c                  ROTATE_GEOD: subroutine from kf, rotates corrections in xyz
c                               to neu and vica versa.
c                  XYZ_TO_GEOD: subroutine from kf, converts xyz to geodetic
c                               ellipsiodal
c                  OCEARG:ADDED SUBROUTINE COMPUTES THE ANGULAR ARGUMENTS
c                         FOR SCHWIDERSKI COMPUTATION OF 11 OCEAN TIDES.
c                         IT WAS COPIED FROM kf/gen_util.
c                  TEST: subroutine added to convert lower case string to upper.
c                  DXYZTOG: subroutine by PJM to convert corrections in xyz
c                           to geodetic ellipsoidal values. not currently used.
c                  ********New Short Period Earth Orientation Modeling********
c                  ROTURT: subroutine modified to apply short period pole and
c                          ut1 corrections calculated in new routine sd_comp.f
c                          by (Herring, kf).
c                  SD_COMP: new subroutine from kf library added to model
c                           directory, for calculation of diurnal and
c                           semidiurnal pole and ut1.
c                  GDATUM : modified to read bitmapped code for turning on/off
c                           short period (diurnal-semidiurnal) earth orientation
c                           corrections.
c                  SETUP/MODEL: modified so that variable "iut1pol" passed to
c                               roturt bitmapped switch for control of short
c                               period corrections.
c                  ********New Antenna Phase Centre Modeling********
c                  OPEN: modified to open antenna phase centre model table
c                        called antmod.dat, unit 32.
c                  SETUP: modified to pass antenna code out to model for use in
c                         reading the antmod.dat table.
c                  MODEL: modified to apply the antenna phase error delay
c                         corrections calculated in new routine phasecc.f to
c                         modelled delay.
c                  PHASECC: new routine to calculate delay corrections cuased by
c                           antenna phase center varioations. Called from model.
c                  LINEAR: new routine to do linear interpolation of a vector.
c                  BILIN: new routine to do bilinear interpolation of a table.
c                  ********Modified P-FILE Output********
c                  OPEN: modified header output to pfile.
c                  CHDRED: modified header output to pfile.
c                  MEDMOD: modified output to pfile.
c                  READI: modified output to pfile.
c                  WRTHED: modified output to pfile.
c version 9.24     OPEN: Don't stop if no ocean tide or antenna phase center table.
c                        Temporary workaround.   McClusky/King 931227
c         9.94a   (MIT)
c                  ********Reorganisation of Model********
c                  MODEL: modified code taken out of main program and put
c                         into 5 new subroutines.
c                  AZ_EL: new routine compute azimuth and elevation of sat.
c                  EPOCH_CLK: new routine to compute receiver clock corrn each
c                             epoch from pseudoranges (if klock=3).
c                  KIN_COORD: new routine to get new lfile coords each epoch
c                             epoch if kinematic data being processed.
c                  PHSMOD: new routine to compute modeled carried beat phase.
c                  POLY_CLK: new routine to compute receiver clock corrns
c                            using polynomial coefficients.
c                        McClusky 931229
c version 9.25     MODEL,OPEN,GDATUM,SETUP,METMOD:
c                  Remove interactive run option, eliminate extraneous screen output
c                  SETUP: Do not call lib routine IMENU
c                  OPEN: Do not call lib routine LASK
c                        Bock 940108
c version 9.26     Merged SIO 9.25 and MIT 9.24a changes.   King 940111.
c                  MODEL, PHSMOD: Correct bug in call sequence; removed unused
c                        variables.   McClusky/King 940113
c version 9.27     PHASECC: Make Rogue and TurboRogue both use Rogue phase-center
c                        table.   McClusky/King 940127
c                  MODEL,AVCLCK, EPOCH_CLK, SETUP, SVANT: Change lblk1 from function
c                      (removed from library) to array filled by call to NSNPRN in
c                      order to get the right block for reused PRN numbers.  King 940124/940203
c                  SETUP, Makefile: Remove SREAD.  King 940203
c                  AVCLCK: Minor change in printout.  King 940204
c                  PHASECC: Comment extraneous printout.  McClusky/King 940321
c                  PHASECC: Comment more extraneous printout.  King 940404/940420
c                  PHASECC: Fix subscript in initialization.  McClusky/King 940506
c version 9.28     Changes to incorporate earth-rotation partials.  Simon McClusky Apr-May 94
c                    Change /lib/dimpar.fti for MAXLAB 20 --> 22.
c                    ROTURT: Add call to eopart.f where earth orientation parameter partials
c                            with respect to inertial site coords are computed. McClusky 940423
c                    WRTHED: Change to write out new prevals and rlabels for earth orientation
c                            parameters. McClusky 940423
c                    MODEL:  Add ut1part and polepart arrays to declarations, and to argument
c                            list of partl subroutine. McClusky 940423
c                    PARTL: Add code to calculate earth orientation partials wrt phase obs
c                    ETIDE: Argument list of roturt calls changed. McClusky 940423
c                    SITCOR: Argument list of roturt calls changed, and polepart and ut1part
c                            arrays passed back to main model module. McClusky 940423
c                    EOPART, Makefile: New subroutine to compute earth orientation partials wrt
c                            inertial site coordinates. McClusky 940423
c                    SETUP: Change code to set space for 21 instead of 15 partials when orbit
c                           partials are prsent in tfile. ie 6 extra eop partials. McClusky 940423
c                  SETUP: Add calling argument to library RSTNFO.  King 940511.
c                  PHASECC: Zero out tables.   McClusky  940516
c                  SETUP:  Store cublic term from I-file for use by autcln.  McClusky 940516
c version 9.29     SITCOR, ETIDE, EOPART, Makefile:  Replace call to ROTURT by call to library
c                      routines PNROT, SROTAT, PNS, and ROTNSP.  King 940523
c                  SETUP: Fix bug leading to zero coordinates if station.info missing.  King 940606
c                  PARTL: Copy in updated EOP version (missed in v 9.28).   King 940606
c                  WRTHED: Add zenith delay to P-file; minor format changes.  King 940608
c                  Makefile: Shift ranlib to be executed only once; remove unused dxyz1g, hisu,
c                          and ttopos.   Herring/King 940616
c                  Makefile: Remove TIMINC (now in /lib) and WOBBLE (no longer used).   King 940706
c                  Remove unused routines from directories:  PRCES, SIDMAT, TIMINC, UT1RED
c                          WOBBLE.   King 940706
c version 9.30     MODEL: Remove wetvar=0 initialization in epoch loop.  Fang/Bock/King 940705/940708
c version 9.31     MODEL, SETUP, SITCOR, EOPART:  Error in the calculation of pole rate and ut1 rate
c                         partials fixed by calculating partial with respect to observation span
c                         midpoint, not begining of the day.  McClusky 940715
c                  PHSMOD: Minor changes to debug print.  King 940805
c                  Convert to using GPST rather than UTC as internal time argument.  Changes to
c                         MODEL, ETIDE, OBSRED, OTIDE, READI, READJ(lib), SETUP, SITCOR, WRTHED.
c                         King 940729.
c                  SETUP: Correct error in calculating midpoint for nutation values.  King 940726
c                  MODEL: Add format to stop statement.  Bock 940810
c                  OBSRED:  Loosen time-check criterion from 0.3s to 1.0s to allow for bad
c                         Trimble clocks.   King 940815
c                  SETUP: Fix problem with matching I-file day.  King 940817
c                  READI: Fix bug with UTC I-files.  King 940828
c                  WRTHED: Fix format for P-file print.  King 940914
c                  READI: Fix bug with any kind of I-files.  King 940914
c version 9.32     OPEN, SETUP: No longer require an I-file.   King 940916
c version 9.33     MODEL, PHSMOD: Add General Relativistic time delay for Earth.  King 941116
c                  SETUP: Save elevation cutoff and antenna model for C-file in 'extra'.
c                  OPEN: Fix writing of run time.   King 941119
c version 9.34     Add transmitter / receiver antenna orientation dependent phase corrections
c                     MODEL, Makefile: Call new routine DIPOLE_COMP.
c                     SITCOR: modified to pass out xrot ( terrestrial -> inertial rotation matrix)
c                     SVANT: modified to pass out the sunhat vector ( sun wrt earth )
c                     DIPOLE_COMP: routine added to compute transmitter / receiver antenna
c                           orientation dependent phase corrections.   McClusky  941202
c                     PHSMOD:  Add corrections.  McClusky  941202
c                     ROTATE_GEOD:  Change dimensions from (1) to (*) to compile under debug.
c                     DIPOLE_COMP: Trial change of sign.  King/McClusky  941202
c                     DIPOLE_COMP: Change sign back.  King 941202
c                     DIPOLE_COMP, MODEL: Save previous phase for each satellite; change argument
c                           list.   King 941205
c                     DIPOLE_COMP: Try reversed sign again with GIG.  King 9401205
c                  PHSMOD:  Remove General Relativity term for test.   King 941205
c                  PHSMOD:  Remove dipole term for test.   King 941206
c                  Update MIT /active without General Relativity and the antenna dipole
c                     terms both commented out in PHSMOD for now.   King 941221
c                  OTIDE: Remove unused variable.   King 941221
c                  GDATUM, READI: Format changes to satisfy XL/RISC compiler.
c                  SETUP:  Remove initialization of unused variable.  King 940103
c version 9.35     SETUP:  Add rad param names to THDRED call.
c                  OPEN, PHASECC, OTIDE: If tables missing, set unit numbers=0 and stop.  King 950127
c                  WRTHED: Fix P-file print of 'extra'.   King 950127
c                  ***************** yaw model *******************************************
c version 9.36     MODEL: Add code to allow calculation of Satellite YAW angle. PT/SM 950316
c                  SOLRED: New routine to compute precise Sun position (from arc routine) PT/SM 950316
c                  LUNRED: New routine to compute precise Moon position (from arc routine) PT/SM 950316
c                  EPHDRD: Open and check luntab and soltab PT/SM 950316
c                  ETIDE,SITCOR: modified to pass in precise moon/sun position PT/SM 950316
c                  SVANT: pass sat unit vectors through from model/gps_coordinates PT 950322
c                  GPS_COORDINATES:  New routine, computes sat. unit vectors corrected for yaw attitude PT/SM 950322
c                  YAW_CONTROL_CRUDE: New routine, comuptes yaw angle for eclipsing/noon turn satellites PT 950322
c                  INITIAL_YAW_ANGLE: New routine, computes initial yaw ang. for sat starting in eclipse PT 950322
c                  DIPOLE_COMP: pass sat unit vectors through from model/gps_coordinates PT 950322
c                  SHADOW1: new routine to compute shadow factor on satellite ( from arc shadow.f). SM 950316
c                  OBSMOD:  phsmod renamed to obsmod, dipoldel and phctdel removed from range calc PT 950322
c                  PHSMOD: dipole term (+ive sign) and relativity turned back on PT 950322
c                  PHSMOD: remove dipole term from range calculation  PT 950323
c
c                  READ_YAW: New routine to read a yaw file   PT 950323
c                  GET_ECL_POS: New routine to compute sun, sat positions at start of all eclipses PT 950331
c                  GET_YAW_RATE: New routine to get correct yaw rate parameter for given sat eclipse PT 950324
c                  ECLIPSE_RISE: New routine to compute yaw angle of sat which rises in eclipse PT
c                  GSP_COORDINATES: allow -ive yaw rates  PT
c                  READ_YAW: read added flag to indicate eclipse yaw rate vs other  PT
c                  OPEN: added unit number for yaw file (=35) PT
c                  SETUP: open yaw file, or assign iyaw = 0 if no file name given  PT
c                  ETIDE: pass sun, moon coords into routine from model (not turned on yet) PT
c                  NORMALISE: New routine to normalise 3 x 1 vector  PT 950403
c                  GPS_COORDINATES: vectors normalised by calling subroutine  PT 950403
c                  ECLIPSE_RISE: vectors normalised by calling subroutine  PT 950403
c                  SETUP: replace nsnprn with svnav_read (allows yaw info in svnav.dat) PT950403
c                  READ_YAW: store eclipses in range 1hr before start - end of obs. span PT950404
c                  GPS_COORDS: renamed from gps_coordinates PT950404
c                  YAW_CTRL: renamed from yaw_control_crude PT950404
c                  SETUP: hour,min added to  call to svnav_read PT950405
c                  ECLIPSE_RISE,YAW_CTRL: allow 5 types of biases for yaw bias PT950406
c                  READ_YAW: stop if date in yaw file is invalid PT950407
c                  ETIDE: Use sun, moon coords passed in from model PT950411
c                  SETUP: arguments for precmod,nutmod,gravmod passed from THDRED, and
c                         back to MODEL PT950419
c                  MODEL: determine whether IAU76 or AENE63 precession to be used, pass to sitcor PT950419
c                  SITCOR: pass which precession model to use to PNROT and ETIDE PT950419
c                  ETIDE: pass which precession model to use to ROTSNP PT950419
c                  GET_ECL_POS: extra arguments in THDRED call PT950419
c                  READ_YAW: return a blank yfilename if error upon opening file PT950419
c                  MODEL: turn off yaw modelling if error occurs opening yfile PT950419
c                  SITCOR: passed IAU76 to srotat PT950420
c                  ECLIPSE_RISE: fixed bug when sat. rises during eclipse PT950420
c                  SETUP: Correct message in SUICID all.  King 950425
c                  MODEL,READ_YAW,ECLIPSE_RISE,GET_ECL_POS,GET_YAW_RATE,GPS_COORDS: change
c                     the size of yawmax from maxecl*3 to maxecl*5  PT950429
c                  MODEL,READ_YAW,ECLIPSE_RISE,GET_ECL_POS,GET_YAW_RATE,GPS_COORDS: change
c                     dimensioning of yaw arrays PT950501
c                  OPN_YAW: new routine to open yaw file and read number of entries PT95051
c                  makefile: add OPN_YAW, update routines calling dimpar.fti  PT950501
c                  MODEL,OPN_YAW: fix bug for when yaw file is not found PT950502
c                  OPN_YAW: changed error handling for bad read of yaw file PT950504
c                  GPS_COORDS: allow yaw rates down to 0.001 deg/sec PT950504
c                  GPS_COORDS, MODEL: compute partials inert. X,Y wrt yaw rate PT950508
c                  MODEL,GET_YAW_POS,SETUP,SITCOR,ETIDE: pass inertial frame and precession
c                      model to routines to compute inertial rotations PT950509
c                  YAW_CTRL: commence full yaw rate when satellite in less than 20% sun PT950510
c                  READ_YAW: fix duplicate statement numbers.  King 950512
c                  READ_YAW: stop if yaw file doesn't cover obs span.  Tregoning 950512
c                  MODEL: print frame, precesison model of tfile to screen  Tregoning 950512
c                  MODEL, GET_ECL_POS: Do yaw dimensioning locally, not in dimpar.fti.   King 950512
c                  READ_YAW: change stop to warning if yaw file doesn't cover obs span. Tregoning 950515
c version 9.36     WRTHED,SETUP,PARTL: Modified to handle 15 parameter orbital model. McClusky  950515
c                  SETUP, WRTHED: Generalize assignment of parameter labels and islots to
c                      handle 5, 20, or 26 partials.  King 950517
c                  Makefile:  Remove SD_COMP (in library).  King 950517
c                  GET_ECL_POS, SETUP: Add argument to THDRED call.  King 950517
c version 9.37     New C-file format:  MODEL, SETUP,CHDRED, OBSRED, WRTHED, CFOUT.
c                      King 950523
c                  OPEN:  Correct message for missing solar and lunar tables.  King 950529
c                  SETUP, WRTHED: Add antenna model to arguments passed.   King 950531
c                  Makefile:  Make all *.a files the same name.  King 950608
c                  WRTHED: Change names of radiation-pressure models to match ARC.  King 950609/950612
c                  GET_ECL_POS, MODEL, PARTL, SHADOW1, SVANT: Remove parameter statements for
c                     maxytp,maxyt2, now in dimpar.fti.   King 950612
c
c version 9.38     MODEL, EPHDRD:  Check reference frame of lunar/solar ephemeris.  King 950616
c                  EPOCH_CLK: use rate for satellite clock Tregoning 950619
c                  MODEL: Moved some things around  Tregoning 950619
c                  EPHDRD: Fix declaration of frames.  King 950619
c                  **Fix bug with elevation computation:  Remove old AZ_EL, add
c                     new AZ_ELEV, change Makefile.  Herring/McClusky/Tregoning  950728
c                  WRTHED: Fixed bug in format statement (line 324) Tregoning 950630
c version 9.39     SVANT:  Fix bug from 950322 (dx(3) contribution zero because
c                      line greater than 72 characters).  King 950714
c                  GPS_COORDS:  Minor changes to clarify IF structure.  King 950714
c                  MODEL, AZ_ELEV  : remove unused variables.   King 950714
c                  SITCOR:  Make implicit none and declare three variables.  King 950714
c version 9.40     MODEL, SETUP, CHDRED, WRTHED:  Add offset of antenna reference point
c                     from mark to C-file record 2.  King 950717
c                  SETUP: Assign names for receiver and antenna based on software version
c                     and antenna code.  (Temporary--put in XHDRED or on X-file eventually.)
c                     King 950718
c                  MODEL: Minor change in order and formating of messages.  King 950718
c version 9.41     READ_YAW, SOLRED: Remove duplicate declaration.  King 960719
c version 9.42     GDATUM: Fix overrun comment line and extra comma in format.  Sanli/King 950719
c                  OTIDE: Comment debug statements to satisfy DEC compiler.  Sanli/King 950719
c                  WRTHED: Add comment about DEC compiler warning.  Sanli/King 950719
c                  CHDRED, ETIDE, GDATUM, MODEL, SETUP, SITCOR, WRTHED:  Change names of
c                     itide/iut1pol to ietide/isptide to match SOLVE and h-file.  King 950720
c                  MODEL: Fixed bug when there are no eclipsing sats Tregoning 950720
c version 9.43     OBSRED: Fix mismatch calling arguments (isnr/amps) for READC5.  950725
c                  READ_YAW: Remove warning about tfile arc lengths. McClusky 950726
c                  SETUP: Extended warning to catch tfile arcs to short for new yaw model. McClusky 950726
c                  READ_YAW: Remove yaw warnings if yaw file doesn't cover obs span.  Tregoning 950726
c                      Merge McClusky and Tregoning changes.  King 950727
c                  GET_ECL_POS: Removed print statement  Tregoning 950726
c                  MODEL,GET_YAW_RATE: Stop if requested yaw rate time is more than 6 hrs since
c                     previous entry in eclipse_yaw array  Tregoning 950726
c                  MODEL: use correct pointer to yaw arrays (isat not ichan)  Tregoning 950726
c                  READI:  Change stop for missing site or day to warning.  King 950727
c                  GET_ECL_POS, LUNRED: Remove redundant dimensioning.  Sanli/King 950728
c                  MODEL,PHASECC: write to print file if no ant model found in
c                     antmod.dat   Tregoning 950727
c                  MODEL: Call GET_YAW_RATE only if SV is eclipsing. Tregoning 950728
c                  GPS_COORDS: Compute the sc_sun_dist.  Tregoning  950728
c                  GET_ECL_POS: Remove extraneous epoch print.  Tregoning/King 950728
c                  CFA: Avoid singularity at zenith (temporary fix).  Sanli/King 950801/950804
c                  CHDRED, OBSRED, SETUP: Remove extra commas and fix call list.  McClusky/King 950803
c                  GET_YAW_RATE:  Widen interval for testing for missing value in yaw
c                      file from 6 hr to 6 hr 10 min. Tregoning  950803
c                       Widen further to 6 hr 20 min. Fang/Tregoning/King 950804
c         9.44     SETUP: Fix mismatched declaration for islip,islpst.  Sanli/King 950807
c         9.45     AZ_ELEV: Make azimuth 0 to 360 vice -180 to 180.  King 950905
c         9.46     MODEL, SETUP:  Don't flag data in MODEL as low-elevation since this cannot
c                      be done consistently if CVIEW has unweighted points but left others
c                      as low-elevation.  We now rely on the MODEL clock checks and AUTCLN
c                      to detect bad data, and we apply low-el cutoffs only in AUTCLN and
c                      SOLVE.    King 951002
c                  OBSRED:  Don't reset low-el and '98' values to 0; allows use of
c                      old editing.   King 951003.
c         9.47     SETUP, CORTYP, KIN_COORD:  Store horizontal antenna offsets properly;
c                      clean up variable names.  King 951017
c                  SETUP:  Add antenna and receiver ids to lib/xhdred call.   King 951019
c                  Makefile:  Remove ROTATE_GEOD, XYZ_TO_GEOD to /lib.   King 951019
c         9.48     SETUP: Modified calling argument list of hisub, and added antenna
c                       (VPC) descriptor to the p-file. McClusky 951019
c                  OPEN:  Hard stop if no antenna phase model table.  McClusky  951018
c                  PHASECC: New format of antmod.dat table. New routine to read
c                       antmod.dat table called read_antmod added to /lib. McClusky 951019
c                  WRTHED: Cosmetic change to p-file output format. McClusky 951019
c                  LINEAR: Modified interpolation array dimensions. McClusky 951019
c                  BILIN: Modified interpolation array dimensions. McClusky 951019
c                  MODEL: Added iant to the argument list of setup.f . McClusky 951019
c                  Remove LINEAR and BILIN to /lib. Change Makefile.  McClusky  951019
c                  SETUP:  Change name of routine to read rcvant.dat.  King 951024
c         9.49     WRTHED: Expand format for antenna offsets by 1 space.  King 951027
c         9.50     SETUP: Get rcvr software and version from station.info; check against
c                       X-file header value.   McClusky/King 951109
c                  OPEN: Add message and suicid if batch file is short.  King 9501109
c                  AVCLCK, EPOCH_CLK, MODEL, POLY_CLK: Change warnings to standard form and
c                      write to P-file as well as screen.  King 951109
c         9.51     SETUP: Modified setup so that full station name written ion the cfile is
c                      always read off the station.info not the cfile or xfile header.  McClusky  951111
c                  WRHTED: Modified to write 16 char full station name onto cfile header.  McClusky  951111
c                  CHDRED: Modified to read 16 char full station name from cfile.  McClusky  951111
c                  OPEN: More explicit message if (newly required) antmod.dat missing.  King 951111
c                  EPOCH_CLK: Add subroutine name to warning for easier grep'ing.  King 951114
c                  EPOCH_CLK: Modified enormous clock error warning format. McClusky 951114
c         9.52     Makefile: Variable ranlib for compatibility with Solaris 2.  Fang/King 951208
c         9.53     SETUP: Compute start time for call to lib/rstnfo (station.info). King 960122
c                  YAW_CTRL: Fix bug in sign of yaw direction.  Tregoning 960201
c                  OPN_YAW: determine if yfile contains coords.  Tregoning 960202
c                  READ_YAW: read sat attitude from yfile if available.  Tregoning 960202
c                  GET_ECL_POS: various changes for new yaw approach.  Tregoning 960202
c                  SHAD_STEP: New routine to determine exact time and attitude that satellite
c                      goes into < 20% sunlight  Tregoning 960202
c                  OPEN: added comment about tmp yaw file unit number  Tregoning 960205
c                  YAW_ATTIT: New routine to compute unit x,y vectors of sat.  Tregoning 960205
c                  WRITE_YAW: New routine to update yaw file.  Tregoning 960205
c                  GPS_COORDS: Move large chunks of code into yaw_attit.  Tregoning 960205
c                  ECLIPSE_RISE: Simplify code now that yaw attitude is passed in.  Tregoning 960205
c                  MODEL: various bits to allow the above changes.  Tregoning 960202
c          9.54    Correct bugs in reading met parameters and WVR values from tables:
c                  MODEL, ATMDEL ,ZHDRED.  Walpersdorf/King 960215
c          9.55    Added report_stat calls to lots of routines: open, otide, model, cfout, avclck
c                    obsred, rsited, read_yaw, poly_clk, phasecc, atmdel, wvrred, atmred, gps_coords
c                    get_yaw_rate, epoch_clk, opn_yaw, write_yaw, ephdrd, setup, gdetic, height, wrthed
c                    zhdred, metmod, gdatum, readi, model. McClusky 960220,960223
c                  OPEN:  Remove extra comma.  Fang/King 960226
c                  OPEN:  Comment out for now the warning about the ocean tide table.  King 960302
c                  READ_YAW: Fix format for new yaw-file variables.  King 960304
c          9.56    MODEL: Fixed format of report_stat error message.  King 960326
c                  SETUP: Fix argument passed to report_stat.  King 960327
c          9.57    SITCOR: Fix erroneous calculation of sidereal time at leap-sec boundary.  King 960517
c                  MODEL, ATMDEL, SETUP, Makefile and new routines NMFH2P1, NMFW2: Add Niell
c                        mapping function.  King 960517
c                  MODEL, CFOUT: Change calling argument to avoid extra time variables.  King 960517
c                  MODEL, OPEN: Write start msg to status file; write yaw msg to warning file.  King 960529
c                  MODEL: Add report_stat calls for yaw and frame messages.  Tregoning 960531
c          9.58    MODEL, OPEN, SETUP:  Move all screen print to report_stat calls; remove old commented
c                      code and clean up comments and upper/lower case.  King 960531
c                  MODEL:  Clean up comments and printout of yaw code.  Tregoning/King 960604
c          9.59    LUNRED, SOLRED: Correct misspelling of 'iendfm', 'iendfs' (probably no
c                      effect on computations except speed, but not sure).  King 960711
c                  WVRRED: Declare 'message'.  King 960711  
c          9.60    Changes for new kf/gamit structure and Makefiles from unimake.  King 960718
c                        All routines:  removed trailing blanks and replaced lib/*.fti 
c                        includes by ../includes/*.h.
c          9.61    SETUP: Modified to write VPC table information into P-file. McClusky 960808 
c                  NMFH2P1:  Avoid changing value of variable 'doy' in southern hemisphere calculations.
c                        Niell/King 960814   
c                  OTIDE: Declare variable explicitly.  King 960814
c                  SETUP: Modified to fix VPC table information in P-file. McClusky 960815 
c          9.62    KIN_COORD: Get coords from X-file for Dynamic, L-file for Static & Kinematic. 
c                        Chen/King 960822
c          9.63    MODEL: Stop if GAMIT.fatal exists.  King 960912
c          9.64    WRITE_YAW:  Fix dimensioning of ecl_pos array.  Tregoning/King 961018
c                  GDATUM: Correct WGS84 flattening value.  Feigl/Murray/King 961030
c          9.65    YAW_ATTIT: Improve header comments   Tregoning 961113
c                  READ_YAW: Improve report_stat message   Tregoning 961113
c                  YAW_CTRL: Fix bug when yaw angles are near zero  Tregoning 961113
c                  MODEL: Fixed bug in eclipsing satellite summary. McClusky 961113
c          9.66    PHASECC:  Add more aliases.  King 970110
c          9.67    Changes to initialize and declare all variables for G77 compiler:  ATMDEL, CFA,  
c                       ECLIPSE_RISE, ETIDE, GPS_COORDS, LUNRED, MODEL, NMFH2P1, NMFW2, OCEARG, 
c                       OTIDE, PTIDE, SETUP, SOLRED, YAW_CTRL.  McClusky 970110  
c                  PHASECC: Use new library subroutine for alias checks.  King 970111  
c                  SETUP: Make 'swver' variable names consistent with station.info; fix 
c                       format for message to report_stat.   King 970111/970117/970122 
c          9.68    METMOD:  Fix report_stat call for missing met info.  King 970123 
c                  METMOD, OPEN, READI, SETUP:  Minor changes in print formats.  King 970127
c                  WRTHED: Fix print of nominal zenith delay.  King 970203  
c          9.69    SETUP: Change report_stat message and add comment suggesting a reason for 
c                       a static/kinematic mismatch between the X- and batch files.  King 970207
c                  Makefile, MODEL, PHASECC, SETUP, UPDATE_OFFS (new): Allow update of antenna 
c                       offsets during the session.   King 970214 970217
c                  SETUP: Correct typo from 970110 or 970122 causing bogus 'NOOL' model 
c                       instead of 'NONE'.  McClusky/King 970217  
c                  MODEL: pass extra arguments to get_ecl_pos (g77 bug)  Tregoning 970312
c                  GET_ECL_POS: initialise variables to zero, change argument list  Tregoning 970312
c                  MODEL, OPEN, SETUP, GET_ECL_POS:  Changes to pass T-file 'frame' rather than 'efixed' 
c                     from lib/thdred.  Force stop in OPEN if model input file says earth-fixed
c                     since the code is untested.  Tregoning 970321/King 970325 
c                  SETUP:  Add warning if antenna phase center model requested but antmod.dat
c                     table has 'NONE'.  King 970411 
c                  MODEL: Print out status only every 100 epochs, not 50.  King 970429
c                  SETUP: Fix swver check for case of input C-file.  King 970502   
c                  SETUP, UPDATE_OFFS: Round times to even minute for station.info.  King 970506  
c                  SETUP, PHASECC: Restructure P-file record and lib/read_antmod calling arguments 
c                     for better tracking of requested and found phase-center models.  McClusky/King 970509
c                  MODEL, SETUP: Remove obsolete elevation-cutoff code.  King 970509  
c                  METMOD: Minor changes to print formats.  King 970509
c                  SETUP: Fixed bug in P-file antanna model output format. McClusky 970514
c                  OBSRED: Check for bogus X-file by comparing epoch numbers.  King 970519 
c                  SETUP: Fixed bug from 970506 change forday calculation for reading station.info.  King 970528   
c          9.70    SETUP: Modify antenna ht printout to give values from station.info.  King 970604 
c          9.71    SVANT: Add (provisional) antenna offsets for Block IIR satellites.  King 970806 
c                  AVCLK, EPOCH_CLK, MODEL, SETUP, SVANT:  Replace logical lblk1 array by integer 
c                           iblksv array.  King 970806     
c                  SETUP: Remove unused 'iblock'; change message and fatal unmatching block #.  King 970808 
c                  MODEL, AVCLCK, EPOCH_CLK, SETUP, SVANT:  Fix bug causing svant failure if t-file has more 
c                          SVs than X-file; rename isvblk variable consistently.  King 970811  
c          9.72    Makefile, MODEL,GPS_COORDS,YAW_ATTIT,GET_ECL_POS: modifications for block_IIr satellites.
c                    BLOCK_IIR,GET_DUSK,CHECK_MODE: new routines for IIr satellite yaw computations  Tregoning 970822
c                  GET_YAW_RATE: comment out check on yfile times. This is temporary until I fix up  
c                           the code of get_ecl_pos etc to modify the yfile to allow "D" type y-bias
c                           (being D to indicate it is a dusk time)      Tregoning 970822     
c                  MODEL,SVANT: Fixed bug in using wrong index in iblksv array    Tregoning 970826
c                  GET_DUSK,SHAD_STEP,YAW_ATTIT: Add argument to call of yaw_attit  Tregoning 970829
c                  GET_ECL_POS,SHAD_STEP: pass block type to shad_step   Tregoning 970829         
c                  MODEL, SVANT: merged ANU and MIT changes.  King 970829
c                  YAW_CTRL: Fix bug in recovering nominal attitude. Set attitude to be nominal
c                      immediately after eclisping satellite enters > 20% sun.   Tregoning/McClusky  970905
c                  SETUP: Fix determination of partials from 'nintrs' for case of velocities on 
c                      T-file.  King 970915
c                  GPS_COORDS: Fix format of precise/crude warning, but then comment out since it 
c                      provides no real help to the user.  King 970920
c                  GPS_COORDS: Fix bug in definition of time variable   Tregoning 970922
c                  BLOCK_IIR:  Remove two debug statements.  King 970923
c          9.73    Clean-up:  remove SCLOCK (code put into EPOCH_CLK); remove old links
c                     in MIT active_hp to PHSMOD (now OBSMOD) , SD_COMP (now in /lib), 
c                     and WRITEE.  King 971020 
c                  READI: Set cubic term to zero if I-file entry missing.  Herring/King 971104
c                  NMFH2P1: Remove duplicate doy calculations that leads to height error of 3-9 mm for
c                     southern hemisphere latitudes of 30-75 deg, respectively.  Morgan/Niell/King 971105 
c                  MODEL: Remove from RCLOCK warning (erroneous) statement about use of polynomial.  King 971205 
c                  MODEL: Put eclipse summary into .status file, not screen (already in p-file).  King 980106
c          9.74    Changes for moving yaw code out of model--Tregoning 971201/King980106
c                   Move to lib: EPHDRD, NORMALISE, SHADOW1, SOLRED, TIMDIF, YAW_ATTIT
c                   Move to orbits for yawtab:  BLOCK_IIR, CHKMODE, GET_DUSK, GET_ECL_POS, GET_YAW_RATE,
c                      OPN_YAW, READ_YAW, WRITE_YAW, YAW_CTRL
c                   Remove altogether: ECLIPSE_RISE, INITIAL_YAW_ANGLE, SHAD_STEP 
c                   New:  SATATT
c                   Mods:  Makefile, GPS_COORDS, MODEL
c                   MODEL, GPS_COORDS: Changing calling arguments to use 6-vector to conform
c                      with changes to orbits/calc_attit and lib/get_att.  King 980113  
c                   MODEL: Remove requirement of a yaw file; change y-file header read.  King 980113  
c                   MODEL, SETUP: Make yfile name and unit name consistent with other files, change
c                        yfile to c*16; move y-file open and read code to setup; remove unused code
c                        and variables for e_time calculation.  King 980114
c                   SETUP:  Fix bug in checking y-file header.  King 980121
c                   MODEL, GPS_COORDS:  Set yaw to ideal yaw and don't call satatt if yaw 
c                         turned off.  King 980121 
c                   MODEL, GPS_COORDS, SVANT: Define a consistent 6-vector for SV coords.  King 980121
c          9.75     SETUP:  Add version number to y-file header.  King 980122    
c                   SATATT: Fix bug in interpolating yaw table.  Tregoning 980227   
c                   METMOD, OPEN, SETUP: Fix bugs in call to report_stat.  King 980227
c          9.76     SVANT: Added other possible block4 satellite z-offsets -0.63 (JPL value) 0.00 (COD value)
c                   Currently the original 1.72 value is being used . McClusky 05/22/98 
c                   GPS_COORDS: Debug format change. McClusky 05/22/98  
c          9.77     SVANT: Changed block4 (IIR) offset to 0.00 (COD EST) for testing. McClusky 05/27/98      
c                   SVANT: Changed block4 (IIR) offset to -0.63 (JPL EST) for testing. McClusky 05/29/98 
c                   SVANT: Changed block4 (IIR) offset to 1.2053 (COD IGS mail #1653 ) for testing. McClusky 05/29/98 
c          9.78     SVANT, SETUP, MODEL:  SV ant offsets read from svnav.dat by lib/svnav_read, no
c                      longer calculated in code.  King 980604
c          9.79     SETUP, WRTHED: Add BERN2 radiation-pressure model.  King 980708
c                   SETUP: Remove extraneous commas in format (IRIX on SGI warning).  Morgan/King 980901
c                   WRTHED: Change assignment of C-file islot values to be same as first M-file value. King 980904
c                   MODEL, CORTYP, ETIDE, KIN_COORD, PARTL, SETUP, SITCOR, UPDATE_OFFS: Change names of station
c                       partial matrices from 'parmt0/parmat' to 'sitepart0/sitepart'. King 980904
c                   MODEL, PARTL, SETUP:  Change calling sequences to include # orbit partials; change
c                       logic in PARTL to use this value. King 980904
c                   SETUP, UPDATE_OFFS:  Fix mismatched dimensions for parmt0 (apparently  
c                       harmless).   King 980905
c          9.80     MODEL, PARTL, SETUP, SVANT, WRTHED: Add SV antenna offset partials and a prior values
c                       to the C-file.   King 980905   
c                   SETUP, WRTHED: Change name of # of orbital ics+parameters from 'nics' to 'norbprm'
c                       for clarity, add 'norbpart to represent the # of orbital partials (same as
c                       'norbprm' if partials are integrated.  King 980908
c                   WRTHED: Remove clock rate and accel from preval vector; revamp and simplify 
c                       logic for assigning labels, slots, and a priori values.  King 980908  
c                   PARTL, WRTHED: Make estimation units of SV antenna offsets meters, not km.  King 980909
c                   CHDRED, WRTHED: Add 'norb' to calling arguments for lib/readc2.  King 980911 
c                   WRTHED: Change parameter names for greater consistency.  King 980916
c                   MODEL, CFOUT, OBSMOD, OBSRED: Save SV clock epoch and L1 phase correction for AUTCLN.  King 981002
c                   MODEL: Remove invalid declaration for 'sitepart'.  King 981007
c                   WRTHED: define nparam = 5 for case where orbits are not estimated  Tregoning 981028
c          9.81     MODEL, CORSIM(new), CORTYP, OBSRED, OPEN, SETUP, SIMRED(new), Makefile:  Create a simulation mode via 
c                     an input S-file instead of X-file.  King 981127/981210  
c                   OPEN, SETUP: Rename some variables and make OPEN implicit none.  King 981204  
c                   OBSRED:  Add missing (and mismatched) ampl1,ampl2 to calling sequence.  King 981208
c                   MODEL, SITCOR, ETIDE:  Clean up tide computations to allow only inertial input; retain
c                      E-fixed and NEU debug. King 981223
c                   PTIDE: Fix bug in order of vector.  Dong 981019/King 981228
c                   Remove implicit statement to satisfy SGI compiler; add implicit none (all
c                      routines previously compiled ok with explicit compiler switch:  
c                      ATMRED, CFOUT  OTIDE PARTL, PTIDE, WVRRED, ZENTMP.  Morgan/King 981231
c                   GDATUM: Fix message for report_stat call.  King 990129
c          9.82     GDATUM: Force use of new diurnal ut1/pole model.  King 990320   
c                   WRTHED: Correct bug in checking # orbit parameters.  King 990320
c                   CORTYP: Do not use comma as continuation character 
c                      (IBM and f90 restriction).  Fang/King 990324   
c                   SETUP: Add FATAL catch when intervals of obs file and yfilet are different   Tregoning 990505
c  version 9.83     OPEN, SETUP: Remove option of reading sited file.  King 990728   
c                   Changes for 4-digit years:  ATMRED, CHDRED, OBSRED, READI, SETUP, 
c                       SIMRED, UPDATE_OFFS, WRHDRED.  King 990728 
c                   OPEN: Add trap for input/output c-files names the same.  King 990802
c                   SIMRED: Add argument to lib/rsesfo call (for MAKEXP).  King 990817      
c                   OBSRED:  Correct Y2K problem msat=0 but valid date.  King 990903
c  version 9.84     OBSMOD:  Compute the SV clock terms written to the C-file (and used 
c                      by AUTCLN) at the actual observation epoch.    Herring 991006
c                   MODEL: Initialize tshad_old (crashes with Linux/g77 on Alpha).  McClusky/King 991101
c                   OPEN: Initialize efixed (crashes with Linux/g77 on Alpha).  McClusky/King 991101
c                   KIN_COORD: Always read from X-file in kinematic and dynamic modes.  Bock 991120
c                   SETUP: Use standard, not input name for antenna in P-file offset print.  King 991130
c  version 9.85     SIMRED:  Remove extra comma in read statement.  Fang/King 000215
c                   SETUP: Add missing argument to READJ call.  Fang/King 000215
c                   ETIDE, MODEL, OPEN, SETUP, SITCOR, Makefile; remove OTIDE:
c                     Read a new ocean tide (u-file) table to compute tides.  King 000224
c  version 9.86     MODEL: Add relativity correction for SV clocks.  Herring 000312/000317
c                   OBSMOD: Fix definition of svcepc to be measurement time.  Herring 000719
c                   ETIDE: Remove debug.  King/Herring 000829  
c                   ETIDE: Correct JD argument.  Herring/King 000830/000927
c  version 9.87     OPEN:  Temporarily use I-file instead of L-file to get U-file name.  King 001220
c                   PHASECC: Save 'sign' and 'antmod_out' variables to avoid loss on second call.  King 010112
c                   CFOUT: Fix declaration of 'ampl1', 'ampl2' to real*4.  King 010112
c                   OBSRED: Initialize to zero quantities not on X-file.  King 010112
c  version 9.88     MODEL: Remove unused shadow code; remove YAW_CONTROL_CRUDE.  King 010216
c                   MODEL, SITCOR: By-pass tide calculations for ipass=1 (saves 10% of
c                      run time, changes rcvr clock by 0.1 ns, no change in obs).  King 010216
c                   SETUP:  Add receiver info from station.info to P-file.   King 010223     
c  version 9.89     OPEN: Temporary trap for naming of u-file from i-file (need to fix later).  King 010719 
c                   OCEARG: Fix argument errors for Q1 and Mf tides.  Matsumoto/Shimada/King 020124    
c  version 9.90     SETUP, UPDATE_OFFS: Allow new-style station.info; change 'swvers' from R*8
c                      to R*4 to match rest of GAMIT.  King 020312 
c                   SETUP: Modify message re bad line in batch file.  King 020318 
c  version 9.91     ETIDE, MODEL, OCEARG, OCEARG2, SETUP, SITCOR, Makefile:  Add NAO ocean tide model.  
c                       Matsumoto, Shimada, King 020617  
c  version 9.92     METMOD (comments only), WHDRED.  King 020708     
c                   MODEL, OPEN, SETUP: Allow use of a GLOBK apr file as well as an l-file; use
c                     velocities to propagate position if movement > 1 mm/day.  King 020807  
c  version 9.93     Makefile, SETUP: Replace model/gdatum by lib/read_gdatum.f. King 021002 
c                   SETUP: Remove debug.  King 021104  
c  version 9.94     WDHRED: Add iostat check to read statements.  King 021204 
c                   ATMRED: Rewrite to clean up logic and prevent NaN if table overrun; replace call to 
c                            'check_y2k' with call to 'fix_y2k' to eliminate warnings. King 021204
c                   WHDRED: Remove line feed in last format statement. King 021205   
c  version 9.95     MODEL, SETUP, UPDATE_OFFS: Cleaner check of station.info changes with new-style file. King 030107 
c                   MODEL, SETUP: Fix bug in logic for old-style station.info.   King 030113 
c                   MODEL, UPDATE_OFFS:  Update antenna type when updating offsets. King 030115  
c                   MODEL, PHASECC:  Reset initial call variable when new station.info entry read.  King 030117
c  version 9.96     SETUP, UPDATE_OFFS: Remove obsolete 'icall' from call to hisub.  King 030418
c                   MODEL, SETUP, PHASECC: Replace call of read_antmod by call of get_antpcv; change table header
c                     in p-file to indicate PCVs now ordered by zenith angle, not elevation angle.  King 030418  
c                   SETUP: Read SV PCV model from input file, reorder printing of SV info to p-file, and
c                     get SV PCV info from antmod.dat.  King 030424  
c                   MODEL, PHASECC, SETUP: Implement corrections to SV antenna phase center.  King 030502
c                   Makefile, MODEL,  CFOUT, CORTYP, CORSIM, SETUP, OBSRED, UPDATE_COORDS (new), WRHED:  
c                     Major revisions to use sb update_coords to do station.info, velocity, and eq/rename 
c                     updates.  Remove sb update_offs.  Rationalize units so that modkin.h krad and kpos, and 
c                     gdetic.dat semi and shft are kept in m (consistent values on files (and with antenna 
c                     phase-center offsets) and converted to km locally for model calculations.  Revisions to 
c                     includes/modkin.h to eliminate kradk, add kpos(3) and group variables more logically.   King 030509
c  version 9.97     MODEL, PHASECC, SETUP: Set 'newant' in model.f, not lib/get_antpcv.f; fix bug in passing svantmod to phasecc.f   King 030521 
c                   PHASECC:  Fix bug in assigning L2 PCV correction.  King 030521   
c                   SETUP: Change label of PCV corrections in P-file back to 'Elevation angle'. King 030521
c                   ATMDEL, METMOD, WHDRED: Resplace function upperc by subroutine uppers. King 031027 
c                   DIPOLE_COMP: Shorten in-line comment to less than 72 col (ifc warning). King 031027
c  version 9.98     MODEL, OPEN : Allow reading of the name of a local scratch file to speed up processing by avoiding NFS.  McClusky 031215
c                   UPDATE_COORDS: Reset 'newant' flag with each call so that coordinate- and phase-center update is
c                      avoided (slows down MODEL dramatically).  King 040115 
c                   MODEL: Simplify calculation of nadir angle for SV PCV corrections.  King 040119
c  version 9.99     CORTYP: Add missing division by cos(lat) in applying longitude offset.  King 040318
c  version 10.01    Makefile, MODEL, ATMDEL, ETIDE, IMFW2 (new), IMFH1P0 (new), LININTER (new), OPEN (comments only),
c                      RDGRID (new), SETUP, SITCOR, UPDATE_COORDS, VMF (new). Changes to apply atmospheric 
c                      loading and use the IMF, VMF mapping functions.  Parameter maxamtlod added to includes/dimpar.h.
c                      Tregoning 040123; merged by King 040328
c                   ATMRED: Fix bug in reading W-file.  King 040526
c  version 10.02    ATMDEL, MODEL, OPEN, SETUP, Makefile, remove WVRRED and ZHDRED:  Appropriate z-file for 
c                      met print file (no longer allow WVR input, never implemented).  Modify batch file 
c                      logic (backward compatible) for w- and  z- files (see also fixdrv/fversn.f).  King 040610/040525 
c                   METMOD, WHDRED: Tweek P-file printout for atmospheric models.  King 040626
c                   SETUP: Fix bug in new met print code.  King 040626
c  version 10.03    SETUP: Add azimuth-dependent PCVs to p-file.  King 040902
c                   PHASECC: Skip PCV evaluation if elevation less than zero. King 040903
c  version 10.04    OPEN, SETUP: Add 'pcncod' argument to read_rcvant call and get the DCB 
c                     corrections from dcb.dat.  King 040930/041011   
c                   MODEL, OBSMOD:  Apply DCB correction.  King 041011 
c                   MODEL, SETUP, OPEN: Skip DCB corrections if file missing or not needed. King 041111        
c  version 10.05    MODEL, SETUP:  Apply Ashtech codeless clock correction if not changed in RINEX file.   King 0412109
c  version 10.06    CHDRED, SETUP, XHDRED, WRTHED: Pass more models to SOLVE; revise c-file format; remove
c                      'extra' variables from xhdred call. King 050201
c                   OPEN: Don't convert file names to uppercase for print. King 050203   
c                   SETUP: Add code to read the header of the atm loading file.  Tregoning/King 050208
c                   PHASECC, SETUP: Make selection by PRN rather than Block #.  King 050209
c                   MODEL, SETUP, SVANT, WRTHED: Keep units of svantdx vector in meters so that it's 
c                     consistent on the c-file  with the parameter vector and the ground antenna offsets. King 050209/0214
c version 10.07     MODEL, ETIDE, IERS2003_ETIDE (new), SETUP, SITCOR: Add IERS 2003 solid-Earth tide model,
c                     pass the model name from the input batch file.  Watson/King 050216
c                   CHDRED: fix compile bug - added missing ")"  Tregoning 050218    
c                   ETIDE, SETUP: Correct misnamed old solid-E tide model to IERS92 (not IERS96). Watson/King 050218
c version 10.08     MODEL, CFOUT: Store atm loading values in 'save' vector of C-file Record 4.  King 050224 
c                   IMFH1P0, IMFW2, OPEN, RDGRID:  Remove unused variables.
c version 10.09     SETUP: Correct file label 'svant.dat' to 'svnav.dat' in SV PCV print to p-file. Herring 050315
c                   MODEL, ETIDE, SETUP, WRTHED: Use c*3 'otide_source' instead of c*8 'otidemod' for file 
c                     formats (OSO or NAO); use 'otidemod' for actual model name written to C-file.  King 050413 
c                   LININTERP: Fix fatal error caused by small roundoff problem.  Tregoning 050415
c                   MODEL: Reset 'newant' to .false. after updating coordinates and phase-center model so as
c                      as avoid rereading antmod.dat at every epoch thereafer (slows down MODEL factor of 50).  King 050425
c version 10.10     PHASECC, SETUP, UPDATE_COORDS:  Add radome to read_rcvant and hisub calls.  King 050618
c                   PHASECC, SETUP: Expand size of 'svantcod' to 8 characters  to accomodate BLK_IIR-A and BLK_IIR-B. King 050618
c                   SETUP: Allow iblk=5 for Block IIR-A SVs.  King 050622    
c                   SETUP: Write warning to p-file if antenna+radome combination not found. King 050623
c version 10.11     ETIDE, PTIDE, SETUP: If ietide bit 5 (16) set, remove the mean pole position of 2000 when
c                     computing the pole tide.  Herring/King 050627
c version 10.12    PHASECC, SETUP, UPDATE_COORDS: Add 'warnings' argument to hisub and get_antpcv calls.  King 050719
c version 10.13    SETUP, PHASECC: Add radome subsitution and minimum elevation of PCV tables to arguments
c                    for /lib/get_antpcv.  King 050819      
c                  MODEL, SETUP OPEN (+ includes/errflg.h): Flag observations at elevation angles less than the
c                     minimum for th PCV model.  King 050906
c                  ETIDE, MODEL, SETUP, SITCOR: Add atm tides read from the u-file.  Tregoning 050718/King 050906
c                  ATMDEL, GMF, MODEL, SETUP, Makefile: Added Global Mapping Function GMF   Tregoning 050831
c                  ETIDE: Fix failure to initialize atm tidal loading variable.  King 051005
c version 10.14    MODEL, SETUP, UPDATE_COORDS: Rearrange logic to call hisub2 rather than hisub.  King 051006
c                  SETUP: Put substituted radome type ('NONE') into 'anttyp' variable for c-file header
c                      to be read in SOLVE for h-file header.  King 051019
c                  PHASECC: Add comment for Block IIR-M.  King 051021
c                  SETUP: Keep original radome name in 'anttyp' for c- and h-file headers, but use
c                      '-' in column 16 to indicate substitution, '+' if no substitution, and 
c                      blank if radome is NONE or UNKN.  King 051026  
c                  SETUP: Allow iblk=6 for BLOCK IIR-M SVs.  King 051027
c                  SETUP: Add report_stat warning for missing E-tide model.  King 051116
c version 10.15    UPDATE_COORDS: Fix warning for too-late start of station.info file.  McClusky/King 060130
c                  SETUP: replace NMF with GMF as fallback mapping function
c                       when no values for VMF/VMF1  Tregoning 060207  
c                  ATMDEL, SETUP, MODEL, VMF1, Makefile: added VMF1 mapping function Tregoning/King 060227 
c                  OPEN:  Fix casefold problem to allow MODEL to read input c-files.  Herring 060317 
c version 10.16    MODEL, SETUP, METMOD, OPEN, VMF1_HT, RDGRID, ATMDEL, READ_METRNX (new), GPT (new): 
c                      Read VMF1 from a global grid, read atm params from met rienx files and, if no file present,
c                       default to using values from the GPT model   Tregoning  060615
c                  MODEL: Remove 100-epoch print and add total # epochs to 'Normal stop' message. King 060718
c                  OPEN: Change messages for 'w-file' to 'met file'.  King 060718
c                  ETIDE, MODEL, OPEN, SETUP, SITCOR:  Change name of u-file to reflect multiple (non-otide) contents.  King 060719
c                  MODEL, ATMDEL, ETIDE, OPEN, READU(new), SETUP, SITCOR, modkin.h:  Generalize code 
c                     for u-file values; put all read values into common; remove METMOD.    King 060720
c                  GMF, IMFH1P0, IMFW2, VMF1: Make calling arguments consistent. King 060722
c                  modkin.h: Rename include to model.h and rename commons within to reflect broader usage. King 060722
c                  SIMRED: Change calling arguments for rsesfo (see lib/lversn.f). King 060815   
c                  SETUP: Open y-file as 'old', not 'unknown', so proper message when missing.  King 060823  
c                  MODEL, ETIDE, SITCOR: Move atm loading interpolation to sitcor; standardize debug print.  King 060825
c                  MODEL, ATMDEL, CFOUT, SETUP: Remove WVR variables; set default in cfout.  King 060828 
c                  MODEL, ATMDEL, SETUP; Move atmospheric delay variables from calling arguments to
c                     new common /metcom/ in model.h.   King 060828
c                  IMFH1P0: Change 'd' to 'cd' for debug lines for gFortran.  McClusky 060828
c                  ATMDEL: Fix report_stat messages.  Shimada 060914  
c                  SITCOR: Fix round-off problem in interpolating loading.  Tregoning/King 060925
c version 10.17    ETIDE, Makefile: Move ocearg.f and ocearg2.f to /lib.  King 061017
c                  ETIDE, SETUP:  Get corrections to center-of-mass offset of ocean tides.  King 061031
c                  READU: Increase line size for NAO99 OTL.  Shimada 061108
c                  WRTHED: Remove the 'E' or 'M' designation in otidemod since CE always corrected
c                    to CM.   Herring/King 061111 
c                  MVERSN: Suppress screen print (status message enough). King 061112
c version 10.18    SETUP:  Fix 30-SV limit in format statement. King 061128
c                  SETUP: Remove incorrect check for reasonableness of autcln cutoff angle.  King 061130
c                  LININTERP: Trap the case when the last observation epoch exactly equals the last
c                    u-file epoch.  Tregoning 061204
c version 10.19    SIMRED: Change calling argument for rsesfo (see /lib/lversn.f) King 061207      
c                  ATMDEL: Fix bug in when to allocate ufile ZHD values  Tregoning 061215
c                  READU: Fix bug of failure to initialize 'end_site'.  Herring 061222
c                  ATMDEL: Fix bug in matrix dimensions for interpolating 'map' values from u-file;
c                    restructure logic to better handle combinations of met values . King 061229
c                  Makefile: Remove GPT, FFUN, SAASZD, WPRESS, added to /lib for /grdtab.  Tregoning/King 061229
c                  OPEN: Write existence of u-file to GAMIT.status.  King 070109 
c                  ATMDEL, LININTERP, SETUP, OPEN, CHECK_MET (new) GET_MET_SOURCE (new), Makefile,  model.h: 
c                     Rework atmospheric delay calculations to better account for missing met values.  King 070110
c                  MODEL, ETIDE, SITCOR: Make atml values for lininterp single precision to match
c                      atm delay changes of 070110.  King 070116
c                  OPEN:   Open and assign units for eq_rename file.  King 070116
c                  ATMDEL: Correct temperature (deg C) value at altitude for c-file record (no effect
c                     on zenith delay calculation).  King 070123
c                  Makefile: Move wpress.f back from /lib (no longer needed by grdtab).  King 070124
c                  OPEN, UPDATE_COORDS: Change filename from 'eq_rename.dat' to simply 'eq_rename'. King 070125
c                  OPEN, SETUP: Adjust p-file printout.  King 070201
c version 10.20    UPDATE_COORDS: Fix format for p-file message with new coordinates. King 070206
c                  MODEL, OBSRED: Add site name to three fatal messages.  McClusky/King 070326
c                  VMF1: Fix typo caught by Solaris compiler.  King 070327
c                  OBSRED: Changed valid eopch test to allow up to 2 sec clock drift. McClusky 070402
c                  MODEL: Issue warning for missing DCB only once per SV. McClusky/King 070403
c                  ATMDEL, READU:  Add a 9th value for VMF1 map list files.  King 070405
c                  MVERSN: Remove print unit number from 'lversn' call.  King 070416
c                  SETUP, PTIDE: Remove unneeded 'iscrn' variable from calls to 'polred' and  'ut1red';
c                     removed unneeded 'jd' 'fract' from call to 'ut1tid'.  King 070416
c                  MVERSN: Remove call to 'clrscr'.  King 070418
c version 10.21    IONDEL, MODEL, OBSMOD, READ_IONEX, SETUP, Makefile: Add ability to read an IONEX 
c                    file and compute the 2nd and 3rd order ionopheric corections.  King 070515 
c                  ATMDEL: Fix bug in time argument for GMF and VMF1.  Smythe/Tregoning/King 070517
c version 10.22    SITCOR, MODEL : Change name of frame rotation matrix from xrot to pnsmat throughout
c                    and allow reuse for tides and ionspheric delay as well as the dipole component
c                    of the phase.   King 070725 
c                  IONDEL, MODEL, OPEN, SETUP:  Allow use of an IONEX (f-) file for computing the
c                    2nd & 3rd order ionopsheric delay.  King 070727
c                  EPOCH_CLK, ETIDE, OBSMOD, READI, SETUP, SITCOR: Remove unused 'iscrn' from calls 
c                    to /lib routines nuttab, pnrot, readj, rotsnp, srotat.  King 070906
c                  ROT_GFILE: Fix badly constructed do loop.  King 070906
c                  AVCLCK, AZ_ELEV, CHDRED, CHECK_MET, DIPOLE_COMP, EOPARTL, ETIDE, GPS_COORDS
c                     IERS2003_ETIDE, MVERSN, OBSMOD, OBSRED, PARTL, PHASECC, POLY_CLK, PTIDE, 
c                     SETUP, SIMRED, SITCOR, SVANT : Remove unused calling arguments. King 070910
c                  Makefile: Remove CFA, MARINI, ZENTMP (no longer used). King 070910  
c                  EOPART, OBSRED, OPEN, WHDRED: Comment out unused labeled statement.  King 070910
c version 10.23    SITCOR: Special code to handle end-of-year boundary in u-file.  King 071023/071102/071116
c                  SETUP:  Read nutation model from the table and compare with the t-file header.  King 071221/071227
c version 10.24    MODEL, IONDEL, OBSMOD, OPEN, READ_IONEX; Makefile; new routines CALL_MAG, 
c                    CALLPPT, ION_CORR, MAG, PIERCE_PT, SLANT_TEC, VEC_XYZ; model.h: Modifications for 2nd 
c                    and 3rd order ionosphere. E. Petrie/R.King 080117
c                  CFOUT: Fix report_stat call.  King 080118  
c                  DIPOLE_COMP: Fix problem with dacos(dot(x,y)) when dot(x.y) = 1.0d0. Bug in gfortran 
c                     compiler with optimization (ok without optimization).  King 080122 
c                  IONDEL: Fix bug in calculating interpolation indices.  Petrie 080124   
c                  MAG: Make explicit none and declare all variables. Herring 080207   
c                  IONDEL: Substitute gamit 'transp' for intrinsic 'transpose, not universally available. King 080208
c version 10.25    ETIDE: Fix bug in calling arguments for ptide.  King 080214
c                  IONDEL: Reduce the number of unused variables.  King 080214
c                  Makefile: Remove unused get_ion_source.  King 080219  
c                  GET_MET_SOURCE: Fix bug causing fatal when RINEX file is bad. Petrie/King 080219
c                  SETUP, WRTHED: Add warning and blank out otide model name if grid or list file present but
c                      not applied. King 080220                                              
c                  ETIDE: Add to p-file 
c                  SETUP: Add application of pole tide, and mean removed, to p-file.  Herring/King 080306
c version 10.26    MODEL, SETUP, includes/errflg.h:  Remove fatal and then skip modeling if an SV
c                      is on the x-file but missing from the t-file.  King 080423      
c                  UPDATE_COORDS: Fix bug in updating position after station.info change.  King 080502  
c                  PHASE_ECC, SETUP, UPDATE_COORDS: Fix bug in recalculating PCV model after station.info
c                    change; add 'debug' to calling arguments for /lib/get_antpcv.  King 080505
c version 10.27    MODEL, SETUP, WRTHED: Use span midpoint for coordinates as well as EOP.  King 080508 
c                  OPEN, SETUP, UPDATE_COORDS: Remove calls to oldstyle hisub and rstnfo; rename 
c                     hisub2, rstnfo2 to hisub,rstnfo.  King080509
c                  SETUP, UPDATE_COORDS:  Calculate E-fixed coordinates only once; if a rename during span,
c                   mark the data from the shorter segment bad.  Herring/King/McClusky  080509             
c                  OBSRED:  Trap bad x-file line.  Herring 080510
c                  OPEN: Print opening of eq_rename file to GAMIT.status.  King 080513  
c                  MODEL: Change report_stat call for too many iterations from 'status' to 'fatal'. McCluskey 080515
c                  MODEL, UPDATE_COORDS: Fix bug in getting PCVs if change mid-session.  King 080522
c                  SETUP: Change p-file print for antenna offset to get original and DHARP in HISTORY file. King 081124
c                  SETUP: Put subsitution of NMF for CFA mapping functions here to avoid many warnings in atmdel. King 081219
c version 10.28    IONDEL, CALL_MAG, CALLPPT: Set thin layer height at 450km in IONDEL and pass this value through. Petrie 081211
c                  IONDEL: Set up so VTEC interpolation will work for whole day with CODE IONEX files with fewer than 13maps. 	
c                      Fix 2 minor VTEC interpolation bugs in IONDEL. Petrie 081211
c                  SLANT_TEC: Use appropriate mapping function for CODE IONEX files according to date. Petrie 081211    
c                  SETUP: Remove check on obsolete CFA mapping function. King 081218
c                  READ_IONEX: Fix bug in format for 'EPOCH OF LAST MAP'; fix subroutine name in report_stat calls.. King 090506  
c                  Makefile, IONDEL, IONO_CORR, MAG_DIPOLE, Makefile: Allow use of dipole field to compare with IGRF.  Petrie, King 090509
c                  CALL_MAG, CALLPPT, IONDEL, IONO_CORR, MAG_DIPOLE, PIERCE_PT, SLANT_TEC:  Cleanup of 
c                      comments and obsolete code. Petrie 090522
c version 10.29    SETUP: Fix time argument in calling GPT.  Fund/King/Tregoning 090630.
c                  OPEN: Comment out debug on magfield.  King 090728
c                  SIMRED: Fix bug in reading number of data types.  King 090827
c                  SETUP:  Add trap for unreasonable geodetic height.  King 090910
c                  SETUP: Change default stop from '2100 0' to '2100 1' to avoid trap in julday. King 090922
c                  UPDATE_COORDS: Put update from station.info message into warnings as well as status.  King 091120
c                  SETUP, PHASECC: Add site name to call of /lib/get_antpcv to be used in error messages. King 100111
c version 10.30    MODEL, OPEN: Assign unique names to the scratch c-files for running in a cluster environment. McClusky 100127
c version 10.31    ATMDEL, CFOUT, CHDRED, CORTYP, IONDEL, MODEL, OBSRED, SETUP, SIMRED, 
c                     UPDATE_COORDS, WRTHED : Consolidate site-code variables into model.h common.  King 100207
c                  Makefile, GET_ANTINFO, MODEL, SETUP, UPDATE_COORDS: New routine to extract and print antenna PCV
c                    information if changed during a session.  King 100207 
c version 10.32    Makefile, ../includes/model.h, MODEL, UPDATE_COORDS, SETUP, SITCOR, READ_METRNX
c                    READU, OPEN, OBSRED, GET_MET_SOURCE, ETIDE, ATMDEL, CFOUT, CORSIM, SIMRED, READ_IONEX
c                    IONDEL, GET_ANTIFNO (new):  New subroutine to get and print updated station.info values
c                      mid-session, added station.info and unit numbers to model.h commons, remove more kinematic
c                      code.  King 100209                                                                        
c                  OBSRED: Set all weights to 1.0 (not just those for valid data).  King 100209 
c                  IONDEL, IONO_CORR, SLANT_TEC: Correct mapping function for STEC (affects pre doy 252 2001).  Petrie 091027/King 100219
c                  IONO_CORR: Interpolate Nmax with VTEC rather than STEC, prevent negative interpolated Nmax for 
c                      very low TEC.  Petrie 091027/King100219
c                  OBSMOD: Correct the addition of the higher order ionospheric terms to the pseudoranges.  Petrie 091027/King 100219
c                  GET_ANTINFO: Fix bug in failure to close autcln.cmd, resulting in an infinite loop.  King 100331
c                  UPDATE_COORDS: Round to even second, not even minute; remove unused 'icall'.   King 100331
c version 10.33    SETUP: Allow iblk=7 for BLOCK II-F SVs.  King 100529
c                  SIMRED: Add format example (comments only). King 100817
c version 10.34    MODEL, Makefile, CFOUT, CHDRED, GET_ANTINFO, GET_ANTPCV, GET_SVANTPCV, OBSRED, PHASECC, WRTHED:  
c                     Add 10-character antenna model names and epoch-by-epoch atmosopheric loading 
c                     variables the C-file; move get_antpcv and get_svantpcv from lib to model; 
c                     add variables to common.  King 100906
c                  READU: Fix bug in component (UNE/NEU) order for atmospheric tidal loading, King/Tregoning 100828
c                   SETUP UPDATE COORDS: Add rcvers, rcvrsn, antsn to rstnfo call. 
c version 10.35    Makefile, MODEL, IONDEL,  CALL_MAG, IONO_CORR, MAG11 (new), OPEN: Change model 
c                    request for higher order ionosphere from 'dipole' or not, to 'magfield'; alter 
c                   units.  Araszkiewicz/Petrie/King 101007
c                  SETUP: Add receiver and antenna serial numbers to p-file.  King 101104
c version 10.36    GET_ANTINFO: Fix bug in logic for setting +,- or blank for radome.  King 101110/101111
c                  GET_ANTPCV: Fix bug in assigning L1 L2 offsets. King 101111
c                  GET_SVANTPCV: Fix bug in computing interpolation interval.  King 101111  
c                  PHASECC: Remove unused variable; rearrange debug print.  King 101111
c                  MAG, MAG11:  Use explicit declarations. McClusky  King 101112
c                  READ_IONEX: Fix minor typo caught by Intel compiler.  McClusky 101112 
c                  SETUP: Correct p-file print name of atmospheric tide model.  Shimada/King 101201
c version 10.37    WRTHED: Allow UCLR1 radiation-pressure model. 110125  
c version 10.38    PTIDE: Replace the IERS 1996 convention for the mean pole with the IERS 2010 convention. King 110524
c version 10.39    SETUP: Set correctly the SP EOP tide model variable for the c-file and h-file (ok in p-file 
c                    and model applied if specified in the sestbl.  King 110725
c version 10.40    READ_METRNX: Read value free-format (though illegal under RINEX v2) to allow reading
c                    improperly written met files.   King 111214 
c                  CHDRED: Correct read statements for 101110 format changes (missed in model).  Herring 111215
c                  SETUP: Fix naming of SP EOP tide model. King 111227
c version 10.41    GET_ANTPCV: Allow 0.5-degree ANTEX files. Moore/McClusky/King 120104
c         10.42    GET_ANTINFO, GET_ANTPCV, SETUP, UPDATE_COORDS, model.h: Put anttyp into common and
c                     use it rather than antcod exclusively in antenna model code. Retain antcod only
c                     for station.info and hisub.  King 120726
c                  MODEL, OPEN, READ_BATCH(new), IONDEL, OBSRED, SETUP, model.h, Makefile: Separate the reading 
c                      of the batch file from the opening of files and remove some unncessary variable.  King 120726
c                  READ_BATCH: Add a fatal check for skd='S' (only static currently supported)
c                  GET_ANTPCV, MODEL, OPEN, dimpar.h, :  Read a site-dependent PCV model   Moore/King 120726 
c                  MODEL, OPEN, READ_BATCH, SETUP: Fix logic in passing/checking ion models. King 120731
c                  MODEL, SETUP: Fix bug in passing atmflg.  King 120731
c                  GET_ANTINFO, GET_ANTPCV, PHASECC, UPDATE_COORDS, model.h:  Put newant in common, remove redundancies, 
c                       and fix bugs in updating PCVs.  King 120801
c                  SETUP: Add i/o error number to yaw-file open failure. King 120801
c                  SATATT: Add commented-out debug for yaw values.  Petrie 120906
c version 10.43    ATMDEL, CFA, OPEN, READ_BATCH, SETUP, WPRESS: Add GPT2 tropo model. King 121211/121219
c version 10.44    CHDRED, WRTHED, SETUP: Add dryzen and wetzen to C-file header. King 130116
c                  SETUP, OPEN: Fix bugs for simulations.  King 130213
c                  OPEN: Fix bug in opening z-file.  McClusky 130314
c version 10.45    SETUP: Add eradmod and antradmod to thdred call. King 130329
c                  READ_BATCH: Fix failure to initialize scratch_dir. Herring/King 130419
c version 10.46    SETUP: Fix bug introduced with version 10.43, erroneously overriding the VMF1 met
c                    models with GPT if GPT among the met_source options; if gpt.grid is present, both
c                    GPT hydrostatic zenith delays and mapping functions were used, if not, VMF1 mapping functions 
c                    but  GPT zenith delays were used.  McClusky 130923
c                  IONDEL: Replace null strings by blanks in report_stat calls. King 131011       
c                  READU: Fix number of arguments in report_stat call. King 131011
c                  MODEL, OPEN: Increase the PID in the c-file name to 6 digits for large machiens. Herring 131206
c version 10.47    GPS_COORDS, MODEL, SATATT, SETUP: Changes for new (10.51) version of the y-file. King 140125
c version 10.48    CHDRED, SETUP, WRTHED: Add E-radiation and antenna-radiation models to the c-file header. King 140327
c                  CHDRED, SETUP, WRTHED: Add ionsrc and magfield to the c-file header. King 140401
c version 10.49    PTIDE: Updated coefficients to be constistent with IERS 2010 standards 
c                  (Radial coeeficient change from 32->33.20 mm/mas). Herring 140406
c version 10.50    MODEL: Fixed a bug causing o minus c (omc) for ichan 31 always getting initialzed to 0.0. McClusky 140410  
c version 10.51    MODEL, Makefile: Move satatt.f to gamit/lib since used by ARC for the UCLR1 and UCLR2 
c                     radiation-pressure models; remove 'iepoch' in satatt debug call.  King 140416  
c                  READU: Correct format for reading ATMOSMAP values from a VMF1 list file. King 150121 
c version 10.52     Makefile, IONDEL,  CALL_MAG, MAG12 (new), OPEN: Introduction of IGRF12. Shimada 150123  
c version 10.53    MODEL, SETUP, CHDRED, WRTHED: Changes for non-GPS GNSS. King 140123
c                  MODEL, SITCOR,   remove obsolete code for E-fixed computations.   
c                  Major mods to put variables in common, change C-file and svnav.dat format, and to 
c                    support GNSS: MODEL, Makefile, AVCLCK, CFOUT, CHDRED, CORTYP, EPOCH_CLK, ETIDE, GET_ANTINFO
c                    GET_ATMDEL (renamed from ATMDEL), GET_SVANTPCV, IONDEL OBSMOD, OBSRED, OPEN, PHASECC
c                    READ_BATCH, SETUP, SIMRED, SITCOR, SVBODY_COORDS, UPDATE_COORDS, WRTHED.  King 150123 
c                  PHASECC: Remove superfluous definition of pi.  Herring 150502 
c version 10.54    SETUP, GET_DCB2(new), Makefile.:  Modify for v2 dcb.tab; add start/stop times to svnav_read call.  King 150528
c                  UPDATE_COORDS: Change function itimdif to integer*8.  King 150528
c                  GETDCB2: Add fatal if DCB not found, and warning if extrapolating.  King 150603
c                  CHDRED, OBSRED, SETUP: Correct problems with calling arguments and times with 
c                    c-file input. King 150612
c                  GET_DCB2, SETUP: Use SVN not PRN in reading new-format dcb.dat.  King 150612 
c                  SETUP: Correct printout of DCBs when using old-format file. King 150617
c                  MODEL: Correct bug in invoking block number for yaw.  King 150924
c version 10.55    SETUP, SVBODY_COORDS: Replace block # with SV body type. King 151124
c                  GET_SVANTPCV, PHASECC, SETUP, : Remove extraneous jd0 argument in call to get_antinfo; 
c                    add jd to arguments for get_svantpcv.  King 151228
c                  GET_ANTPCV, GET_SVANTPCV: Add 'gnss' to read_antex call so to use GPS calibrations 
c                    for ground antennas if none available for other GNSS.  King 151229
c                  SETUP: Fix report_stat message for missing SV PCVs.  King 160108
c                  OBSRED: Remove kinematic/dynamic reads of x-file.  King 160420
c version 10.56    GET_ATMDEL, GET_MET_SOURCE, READU, SETUP, model.h:  Clean up the setting of sources for
c                    P, T, H and the echo to the p-file.   King 160426
c                  GET_DCB2: Add SVN to warning message about extrapolated values. Morgan/King 160503
c                  READ_METRNX: Replace free-format read with RINEX-standard F7.1 to guard against
c                     -9999.9 values with no preceding space.  King 160612 
c                  SETUP: Set IRNSS frequencies and change fatal to warning. King 160822
c                  DIPOLE_COMP: Add initialization of integer phase (no effect with most operating 
c                    systems). King 160901
c version 10.57    ETIDE: Add pole-position to returned values from lib/rotsnp King 170412 
c version 10.58    SVBODY_COORDS: Remove the sign reversal for IIR body coordinates to be consistent
c                     with IGS standards and the new orbits/kouba_gps.f routine.  King 170505
c                  MODEL, SATATT, SVBODY_COORDS: Pass the ievent variable to key orbit-normal 
c                     attitude.  King 170706       
c                  GET_ANTINFO, GET_ANTPCV, GET_SVANTPCV, SETUP : Set and pass the two ANTEX frequency 
c                     codes to be used. 
c                  GET_DCB2: Skip warning if not GPS.  King 170531 
c                  SETUP: Change Beidou L1/C1 designations to L2/C2 to reflect change in nomenclature
c                     with RINEX 3.02.  King 170602 
c                  SETUP: Allow single-frequency processing for GPS.  King 170814 
c version 10.59    SETUP: Move code to get Glonass frequencies read. Herring 180319
c                  SETUP: Add printout of tidal waves to OTL model line.  King 180319
c version 10.60    MODEL, ETIDE, OPEN, SITCOR, SVBODY_COORDS, Makefile: Add reading of planetary 
c                     ephemeris file; add arc.h commons; remove SOLRED and LUNRED in favor of /lib 
c                     versions; make Moon vectors dimension 6 for consistency.   King 180320
c                  MODEL, CHDRED, EPOCH_CLK, ETIDE, GET_ANTPCV, GET_ANTINFO, GET_MET_SOURCE, OPEN, 
c                     READ_BATCH, SETUP, UPDATE_COORDS: Add new common includes/units.h to provide 
c                     unit numbers removed from includes/model.h. King 180320
c                  MODEL, CHDRED, OBSRED, READ_BATCH, SETUP, UPDATE_COORDS, WRTHED: Set frame from 
c                     t-file header; change obsolete 'skd' to 'gnss' in the batch file and calls to 
c                     lib/xhdred; remove'skd' from includes/model.h.  King 180322 
c                  SETUP, UPDATE_COORDS: Remove obsolete 'skd' from calls to lib/rstnfo. King 180322
c                  MODEL: Call /lib/evrtcf to get Everett coefficients for solred and lunred. King 18-324
c                  SETUP: Move code again to have the frequencies set for getting SV antenna offsets. King 180324
c                  OPEN: Open the nutation file only if it is needed (no 'nbody' for Sun and Moon). King 180328
c                  MODEL, OPEN. SETUP: Replace data/neweph/ with check on the existence of 'nbody' to
c                     invoke the new ephemeris code.  King 180329
c version 10.61    WRTHED: Fix label and # parameters for BERN2/ECOM2.   Herring/King 180405
c                  SETUP: Remove careless line causing storage overwrite. King 180507 
c version 10.62    GET_ANTPCV: Fix logic in trapping missing ANTEX frequency. King 180531
c                  READU: Read the ocean tide waves from the u-file.  King 180531
c                  SETUP: Change traps for wrong frequencies from fatals to warnings. King 180606
c version 10.63    GET_SVANTPCV: Allow for L1/L2 PCO differences (are so far only for Galileo
c                    and small).  King 180716
c                  GET_ANTPCV, GET_SVANTPCV, includes/model.h: New calling arguments for lib/linear.f 
c                     (see lib/lversn.f) to use floating rather than integer intervals. King 180720
c                  GET_ANTINFO: Restore p-file printout of elevation-dependentt PCVs.  King 180720
c                  UPDATE_COORDS:  Add missing 'debug' to call of get_antinfo.  Gegout/King 180907 
c                  SETUP: Fix trap on frequency so that it works for L1-only.  Fabris/King 181121
c version 10.64    ATMRED,  AVCLCK, CHDRED, EPOCH_CLK, GET_ANTINFO, GET_MET_SOURCE, OPEN, READI
c                       UPDATE_COORDS, WRTHED, ZHDRED: Change (back) 'luprt' to 'iprnt' for
c                       consistency with model.h. King 190425
c                  WRTHED, dimpar.h: Increase slots for 4 additional ECOM2 radiation paramters. King 190425
c                  SETUP: Set npart and norbpart correctly for ECOMC (19 orbit parameters). King 190523 
c                  SETUP: Added antpwr to svnav_read call Herring 190702.     
c version 10.65    MODEL, CHDRED, GET_SVANTPCV, PHASECC, SETUP, SVANT, WRTHED, model.h:  Use separate 
c                   PCOs and PCVs for SV antenna offsets if available. King 190826 
c version 10.66    SITCOR, SETUP: Allow use of new IAU sofa inertial RF models (IAU0A, IAU0C and IAU06). 
c                  McClusky 190801. Merged Herrinng 190918.
c version 10.67    ETIDE, IERS2010_ETIDE.  Added IERS10 Solid earth tide model.
c                  McClusky/Herring 200123.
c                  MAG, MAG11, MAG12, MAG13; IONDEL (prints) ICRF13 added for magnetic field past 2020.  
c                  Added igrf_util for common routines Herring 200124
c                  MODEL, WRTHED. svantdx changed to (2,3,nsat) for L1/L2.In wrthed preval is LC valuue. 
c                  MODEL.h to added antdaz to common.
c                  SETUP, UPDATE_COORdS, GET_ANTPCV, WRTHD, CHDRED, PHASECC, GET_ANTINFO: All updated
c                    to change PCO and PCV corrections based on antenna azimuth antdaz. Herring 200205.
c                  SETUP, ETIDE, PTIDE, READ_BATCH: Updated for IERS20 mean pole.  Called call to
c                    mean_pole in libraries/comlib; select model starting with bit 7 (IERS20), bit 5 (IERS10)
c                    bit 3 ZERO mean,  Improved documenation in model.h  Herring 200219.
c version 10.68    SETUP: Allowed both L6/C6 (preferred) or L7/C7 (if no L6/C6) to be processed 
c                    together.  There may be ambiguity resolution issued with mixed data, Changed nutation 
c                    default to IAU0A to be consisent with fixdrv. Herring 200223/200303.
c version 10.69    MDMAKE: Updated decoding of Earth Rotation sestbl. entry. TAH 200505.
c                  SETUP: Updated documentation for Desai and Sibois model. Added new lower
c                    frequencies for G,E and C systems TAH 200512.
c                  GET_SVANTPCV, GET_ANTPCV: Updated for 3 atxfrq choices and 4 characters. TAH 200512.
c                  UPDATE_COORDS: Tested antenna changes to see if messages should be printed 
c                    about changes. TAH 200527
c                  SETUP: Fixed setting alternative C06 frequency for Beidou when lfreq 7 option used TAH 200617.
c version 10.70    SETUP: Added capability to read and use GPT3 model from gamit/lib/. King 20201216/Floyd 20201218.
c version 10.71:   SETUP: Multiplied height by 1000 to convert to m in call to gpt3. Floyd 20210308
c version 10.72:   SETUP: Added L5 (B2a; B2b = L7) lower frequency choice. Floyd 20210729

      RETURN
      END


