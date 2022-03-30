Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego, 1994.  All rights reserved.
      SUBROUTINE OVERSN(VERSION)

      implicit none

C     Update history for Orbit utilities BCTOT, NGSTOT, TTONGS, TTOICS, TTOASC
C     R. King   26 December 1991

      CHARACTER*128 VERS
      character*45 libver
      CHARACTER*120 VERSION
      CHARACTER*10 GETMAC,MACHIN

      integer*4 nblen
      integer*4 ierr

      MACHIN = GETMAC(1)

c     get the orbit utilities version
      WRITE (VERS,5,iostat=ierr) MACHIN(1:nblen(machin))
    5 format ('10.21 2021/02/02 20:37 UTC (',a,')')
      if( ierr.ne.0 ) call report_stat('WARNING','OVERSN','ORBFIT',
     .                                  vers,'IOSTAT ERR',ierr)

c     get library version
      CALL LVERSN(libver)
      
      version = ' ' 
      WRITE(version,10) trim(vers),trim(libver)
   10 FORMAT('ver. ',a,' Library ver. ',a)

C History:  Comments through version 8.5 apply to BCTOT only

C Version 1.00  In BCTOT, T-files will have NAVSTAR numbers - this is a
C                temporary measure until all modules are PRN based
C               YB 1/9/88
C               Another temporary measure - for full B-file (several messages
C                per satellite), use only first message until fit is in place.
C               In previous B-file, only the last message was saved.
C               YB 1/10/88
C               Bug in TMAKE and FILENM - 1/13/88 YB
C               Sent to AERO - 1/13/88 YB
C Version 1.10
C               FILENM.NEW and BREAD.NEW read new E-file format from MAKEX
C                YB 1/15/88
C               Bug in TIMESE - 1/16/88 YB
C Version 1.20  Added new FICA Block 9 E-files in BREAD, added new routine
C                RDFICA - 1/18/88
C
c version 2.2   Wednesday, January 20, 1988   8:23:19 pm (EST)
C                delta s.dimpar.fti
C                 8 sats, 9 ICs
C                2.2
C                delta s.srotat.ftn
C                 PJD is now an integer*4
C                2.2
C                delta s.cpyrgt.ftn
C                 Fix spacing.
C                2.2
C                delta s.trot.ftn
C                 General clean up of problems due to type conversion.
C                2.2
C                delta s.taiutc.ftn
C                 Complete re-write, making all JDs integer*4.
C                2.2
C                delta s.secmin.ftn
C                 Fix problems due to type conversion.
C                2.2
C                delta s.thdprt.ftn
C                 Fix problems due to type conversion.
C                2.2
C                delta s.bctot.ftn
C                 Fix problems due to type conversion.
C                2.2
C                delta s.at_subs.ftn
C                 Hilite is now just a dumb print.
C                2.2
C                delta s.openb.ftn
C                 Record length = 2048 on apollos.
C                2.1
C                delta s.iedit.ftn
C                 This is the apollo version.
C                2.2
c                delta s.brdxyz.ftn
c                 Fix problems due to type conversion.
c                2.2
c                delta s.tmake.ftn
c                 Fix problems due to type conversion.
c                2.3
c                delta s.timecn.ftn
c                 Correct Yehuda's mangled conversion to integer*2
c                2.2
c                delta s.rdfica.ftn
c                 I/O status on apollos must be integer*4
c                2.2
c                delta s.bread.ftn
c                  I/O status on apollos must be integer*4
c                2.2
c                delta s.filenm.ftn
c                 Standard tables on Apollos
c                2.2
c
c    version 2.3
c    Thursday, January 28, 1988   2:29:50 pm (EST)
c                change to read the format statements from the
c                tables:
c
c                delta s.polred.ftn
c                 Read a general polar motion table.
c                2.2
c                 delta s.ut1red.ftn
c                 Read a general UT1 table.
c                2.2
c                delta s.taiutc.ftn
c                 Read LEAP.SEC instead of //BULLEN/DATA/GPS/TABLES/LEAP.SEC
c                2.3
c                delta s.openb.ftn
c                  Open tables in working directory.
c                2.2
c   version 2.4  January 29, 1988
c               change UT1RED, POLRED to read integer values
c               from tables.
c   version 2.5  Thursday, February 11, 1988
c                delta s.ut1red.ftn
c                 Read magic multiplier (RMULT) from file header, too.
c                2.4
c                delta s.polred.ftn
c                 Read magic multiplier (RMULT) from file header, too.
c                2.4
c                delta s.bctot.ftn
c                 Correct comparison of number of sats in X- and T- files.
c                2.3
c   version 2.6  Now makes an earth-fixed and inertial G-file, new
c                 routine GMAKE and modified TROT and BCTOT - YB 4/18/88
c
c   version 6.1  New GAMIT release -- 880417 MHM
c           6.2  Temporary delete of GMAKE calls in BCTOT, until
c                 velocity calculations are corrected -- 880418 MHM
c                Velocity calculations corrected in BRDXYZ
c                Call GMAKE again in BCTOT -- 880418 YB
C                In TMAKE, change default number of Ic's to 8 in T-file header
c                In BCTOT, initialize radiation parameter slots
c                Comment out call to GMAKE until velocity components checked
c                YB 4/19/88
c   version 6.4  POLRED:  Correct format statement for POLE. header read
c                 MHM 880501
c   version 6.5  DIMPAR:  Change MAXSAT to 7 and MAXICS to 9 to match all
c                 other GAMIT modules -- rwk 880805
c   version 6.6  Set NICS=MAXICS in TMAKE for writing on header of
c                 T-file -- rwk 880815.
c   version 6.7  Add interpolation scheme to get G-files.  Use PRN
c                 numbers vice NS numbers.  Change DIMPAR to match
c                 other modules (MAXORB vice MAXICS).  Changes to BCTOT,
c                 GMAKET, SATEL, THDPRT, TMAKE, TIMESE, DIMPAR -- rwk 880819
c   version 6.8  Restructure to clean up code and add a check on E-, X-file
c                 satellite correspondences.  Remove unnecessary I*2s.
c                 Remove TIMECN in favor of libratry functions. Allow lowercase
c                 filenames. Changes to BCTOT, TMAKE, BREAD, BRDXYZ, SATEL,
C                 TIMESE, TINPUT, XHDRED.
C                 New routine CHKSVS.  Removed TIMECN, HMS, FJLDAY, FJDTIM-- rwk 890821
c   version 7.1  Write I*4 T-file, and add more explicit declarations.
c                 Changes to BCTOT, TMAKE, TROT. -- rwk 880822
c   version 7.2  Begin to port to sun. Use universal BREAD including subfr1.
c                 Kurt 900505
c   version 7.3  Correct BRDXYZ to use same units as BREAD.
c   version 7.4  SUBFR1 now dimensioned to 8 elements.
c   version 7.5  Change ROTCRD from (1) to (*) kurt 900828
c   version 8.1  correct satname formats in BCTOT
c                correct XHDRED call in TINPUT.
c                There is a bug remaining if the T-file (E) exists.
c                Kurt 901114
c   version 8.3  Fix IEDIT bug:  Add nbatch to GMAKET.
c                king and mhm 910104

c   version 8.5  Scripps changes Nov 90 - May 91
c                Incorporate batch and interactive modes
c                BCTOT, TINPUT, GMAKET
c                In FILENM, Do not care if T-file already exists, overwrite it
c                  also correct UNKNOWN in FILENM (was UKNOWN),
c                  lowercase E for earth-fixed T-file name
c                In GMAKET, open G-file as unknown
c                In GMAKET, character*1 DUMMY
c                Yehuda Bock 20 Nov. '90
c                New Makefile, -cpu mathlib_sr10, maxsit=maxnet=37, maxsit=15
c                    Yehuda Bock 5/24/91
c   version 8.6  Combine in orbits directory with NGSTOT and TTONGS:
c                    Changes to RTOT, OPENB, FILENM changed and renamed to FILNMB,
c                    slightly new versions from NGSTOT directory of NUTRED,
c                    UT1RED, POLRED, ROTCRD, ROTSNP, SROTAT, and TIMINC.
c                    Replace SECMIN by library SECDIF in BCTOT and TMAKE.
c                    Changes to BCTOT and TMAKE for array order and T-file unit
c                    Changes to BRDXYZ, BCTOT and TMAKE to check for bad s-major axis
c                 Create OVERSN viz BVERSN to include changes to all routines
c                    Add call to BCTOT,TTONGS, NGSTOT, TTOASC, TTOICS (not TTOT)
c                 Modify TTONGS to use 5 fewer epochs on each end in order to
c                    maintain the (11-pt) T-file interpolation accuracy.
c                    Modifications to TTONGS, WSP1HD, SDTRIT.   King 911226
c                 Change SDELT to be consistent with MODEL in GSATEL, TTONGS,
c                    and TTOICS.    King   911226
c                 Clean up GSATEL a bit (first step in a longer process).  Fix
c                    a bug leading to the use of undefined YY values; Prohibit
c                    interpolation within 5 epochs of either end of the T-file.
c                    Add # epochs to call of GSATEL in TTOG, TTOICS, & TTONGS.
c                 Incorporate TTOG, removing call to OPNCHK and changing the
c                    calls to GSATEL.  Change the Makefile.  King 911228
c                 New call to library DS2HSS for leap-year problem; changes
c                    to WSP1HD, TTOICS, TMAKE, THDRIT, SDTRIT, BCTOT.  Bock 920102
c                 Fix UTC-GPST bug in TTONGS; changes also in SDTRIT.  King 920104
c                 Call new library UT1RED to include option of tidal terms: SROTAT.
c                 New calling sequence for XHDRED:  BCTOT, TINPUT.     King 920123
c                 Remove IEDIT calls in TIMESE and GMAKET; replace SATEL by
c                     GSATEL in GMAKET.    King 921025
c                 Makefile:  Remove GSATEL, now in library.   King 920127
c                 Remove UT1RED, UT1TID, POLRED, NUTRED, NUTTAB, FUNARG, FUNCOF, SIDMAT
c                     SIDTIM, and PRCES to library; reconcile calling routines:
c                     TROT, ROTSNP, PNROT, SROTAT    King 920131
c                 Fix bug in /lib/taiutc and /lib/nutred causing infinite loop
c                     in BCTOT and TTONGS.    King 92/2/19
c version 9.1     BCTOT : Changed hardwired PRN check > 37 to be consistent
c                   with reade in library -- this is the largest PRN number
c                   planned - I believe
c                     Yehuda Bock 920503
c version 9.11    TINPUT : add missing variable isessn to xhdred call
c                 Remove atsub,edtsub, gettim,getdat from directory
c                  and Makefile (all in lib)
c                     Yehuda Bock 920507
c                 BCTOT:
c                  initialization of nbrd by ibrd-1 instead of ibrd
c                  check for PRN = 0 moved to proper spot
c                  fix e-file record picking for no records past midpoint
c                 Bock/Genrich 920510
c                 Remove dversn from Makefile
c                 Bock 920513
c version 9.12    NGSTOT : Fix branching bug.
c                 TIMECN : Fix typo.
c                 TTOT :  Remove machine dependence.
c                 Remove unused variables:  BCTOT, BRDXYZ, FILNM, GETICS, GMAKET,
c                        OPENB, PNROT, ROTSNP, SDTRIT, SROTAT, THDRED, THDRIT,
c                        THDPRT, TIMECN, TIMESE, TINPUT, TTOASC, TTONGS,
c                        TROT, TTOG, TTOASC, TTOICS, TTOT, WSP1HD
c                 Makefile:  Remove TABLES. King 920623
c                 BCTOT: Add calling argument for library READE.  King 920624

c version 9.13    GMAKET:  Remove dummy variable from Sun version.
c                 IDINT4:  Remove from directory (not used).
c                 ROTCRD:  Change (1) to (*) on MIT Apollo and Sun.
c                 ROTSNP, TIMECN, TINPUT:  Remove extraneous declaration.
c
c version 9.14    BCTOT : Fix file open problem with earth-fixed G-file
c                         Fix IC logic problem when only 1 entry per SV.
c                           King 920923
c version 9.15    TINPUT : Round fraction of day to nearest whole number
c                           of seconds. Bock 930115
c version 9.16    TINPUT : Add antcod to XHDRED call.   King 930212
c                 Makefile:  Removed unused RHDRED.   King 930212
c version 9.17    BCTOT, SORTBC, orbits.fti: Common block for ephem to handle
c                         large values of maxbrd, remove unnecessary array t.
c                 INDEXX: Change (1) to (*).    Murray/King 930222
c version 9.18    TTONGS : Add an option to select either SP1 or SP3 format
c                 WSP3HD : newly added to TTONGS module for SP3 header
c                 SDTRIT : modified so that it can optionally write out
c                          either SP1 or SP3 records
c                 NGSTOT : Can use SP3 format, as well as SP1
c                 RSP3HD : newly added to NGSTOT module for SP3 header
c                 TDTRIT : modified to handle SP3 format
c                 GETICS : generalized to accomodate SP3 format
c                 NGSOPN : let user know about SP3 format
c                 Makefile: Add RSP3HD, WSP3HD.
c                 Fang/Bock 930407
c                 WSP3HD : Add additional space in SP3 comment lines
c                 Bock/Fang 930804
c version 9.19    Add SIO orbit-differencing code to MIT Apollo and stdrel.
c                 Makefile, New routines: ORBDIF, INTERP, ORBRMS, RHDRED, RDSPDT
c                 ORBDIF, ORBRMS: Change to generic intrinsics.
c                 Remove ARGMNT from directory and Makefile (in lib). King 930825
c version 9.20    TTONGS: Get day-of-year from T-file header, not filename.
c                       Remove tabs.   King 930902
c
c version 9.21    0RBDIF: ORBRMS subroutine changed to write maximum orbit
c                   differences in the header of the plot-files, used later
c                   for scaling axis in the plotorb script file.
c                   The conversion of X/Y/Z orbit differences to radial/along
c                   /cross track components was simplified with errors in the
c                   cross track differences in the old version corrected.
c                     McClusky 930908.
c                 SROTAT:  Add sidereal time to arguments of (lib) SIDMAT.
c                     King  931013
c                 NGSTOT: Fix fatal bug when X-file input.  Add/remove debug
c                     from TDTRIT.   King 931013
c                 THDPRT: Add year to calling arguments for NSNPRN.  King 931018
c                 Increase formats for 32 satellites: RSP3HD, WSP1HD, WSP3HD.
c                     King 931018
c                 GETICS: Fail more gracefully with bad file; remove unused
c                         variable.
c                 NGSTOT: Fix computation of ic epoch.   King 931025
c version 9.22    ORBDIF: Modified to read t as well as ngs files.
c                         Modified the interpolation of files from interp.f
c                         to the gamit library interpolator gsatel.f.
c                 ORBRMS: Modified logic, orb arrays passed in now contain
c                         position and velocity at matched epochs of orbit_iref,
c                         and orbit_i2nd, NOT orbit_ref, orbit_i2nd, orbit_iref-
c                         orbit_iref(1ms off) and orbit_iref-orbit_i2nd(1ms off)
c                 TTEMP:  New routine to read NGS files and temporarily writes
c                         them in tfile format for later interpolation by
c                         gsatel.f.
c                 RDSPDT and INTERP no longer called.
c                        McClusky 940122
c                 ORBDIF: Fix declarations.  Add HMS2JD to Makefile.  King 940122
c                 ORBITS.FTI:  Increase maxbrd from 500 to 750.  King 940207
c version 9.23    Makefile: Add dependencies explicitly.  Fang/King 940215
c                 WSP3HD:  Reduce comment line to 60 characters.  Fang 940412
c                 Makefile: Change HP compiler flags from +e +E1 to +U77.  Fang/King 940506
c version 9.24    Modify ROTSNP, SROTAT routines and move them, with PNROT, to /lib.
c                   Changes also to Makefile and TROT.   King 940520
c                 Makefile: Shift ranlib to be executed only once.  Herring/King 940624
c                 ORBDIF: Modified to fix bug occuring when comparing orbits with different tabular
c                         intervals. Problem arose when ref_orbit had an earlier starting time than
c                         the 2nd_orbit.   McClusky/King 940705
c                 Makefile: Remove DAYNUM, DOT, SECSMD, SECSUM, TIMINC (in lib) ;
c                           GETRHD, GETRSV, NGSPRG (not used).  King 940707
c                 Remove unused ARGMENT, ATSUB, EDTSUB, GETRHD.     King 940707
c version 9.25    TTONGS, ORBDIF, TTOG, Makefile: Replace THDRED with library version
c                      which now returns GPST rather than UTC.
c                 SDTRIT, TTONGS:  Remove unused variable from calling argument.
c                 NGSTOT, TTONGS, THDRIT, Makefile:  Use THDRIT for TTONGS as well
c                      as NGSTOT; remove THDPRT.    King 940714
c                 BCTOT, ORBDIF, TTOT, WSP3HD: Declare undeclared variables.  King 940714
c                 WSP1HD, WSP3HD: Correct screen info (GPST not UTC).  King 940714
c                 NGSTOT, TTONGS, TTOT, THDRIT, TDTRIT: Pass nintrs.  King 940715
c                 TTONGS: Convert input keywords to uppercase.   King 940715
c                 NGSOPN: Add diagnostic warnings for incorrect sp3 filename.  King 940718
c                 GETICS: Fix reading of SP3 file.   King 940718
c                 Declare undeclared variables:  CHECK, CHKSVS, FJDTIM, FJLDAY, GEOXYZ,
c                      GMAKE, INTERP, LJUST
c                 Write G-files for GPST, and consolidate routines:  Makefile, BCTOT,
c                      GMAKE, GMAKET (removed), NGSTOT, TTOG.   King 940719
c                 RSP1HD, RSP3HD: Fix time designation on header print.  King 940719
c                 GMAKE, BCTOT, NGSTOT, TTOG:  Add program name to call.  King 940720
c                 TTONGS, TROT: Fix minor bugs.  King 940805
c                 BCTOT, TMAKE : Fix bug in writing broadcast g-file.   King 940811
c                 TINPUT, TIMESE: Change default T-file interval to 15 minutes and
c                      question header to 'GPST'.   King 940811
c                 GMAKE: Blank header field to avoid writing nulls.  940815
c                 BCTOT: Fix problem with IC epoch.  940815
c                 TINPUT, GMAKE: Continue to fix problem with IC epoch.  940817
c version 9.26    ORBDIF: cleaned up initial comments!!!!!! McClusky 940907
c                 ORBRMS: modified code to allow plot files containing satellite ground tracks and
c                 orbit difference vectors in terms of NEU to be created. McClusky 940907
c                 ROTATE_GEOD, XYZ_TO_GEOD, DVDOT, DVSWP: routines from herring KF directories
c                 added to orbits directory. ###[Should be added to LIB directory as also used in
c                 MODEL]###. McClusky 940907
c                 New module for estimating orbital differences: new ORBFIT, INVER2, EST_FIT.
c                      Change Makefile.   McClusky 940913
c version 9.27    ORBDIF:  Allow simultaneous solution for global parameters.
c                     New subroutines: ELEM, EOPART (modified from model), GPGT (from solve)
c                     INVER2 (from solve), PLT_POSTFIT, READ_INPUT, TERR_PAR, TRANPART, TMERGE.
c                     Changes to ORBRMS, Makefile.  McClusky 9400921/941105
c                 TTONGS:  Minor format change to make ANSI standard.  King 941205
c                 TTOASC: Minor change to accomodate new 9 parameter srp model tfile...
c version 9.28    ORBDIF: THDRED argument list modified (added icsnam). McClusky 941222
c                 ORBFIT: THDRED argument list modified (added icsnam). McClusky 941222
c                 TROT: THDRED and THDRIT argument list modified (added nics,icsnam). McClusky 941222
c                 TTONGS: THDRED argument list modified (added icsnam). McClusky 941222
c                 GMAKE: THDRED argument list modified (added icsnam), also simplified
c                        writing of gfile header. McClusky 941222
c                 THDRIT: Argument list modified (added nics,icsnam), code modified to write orbital
c                         parameter list into comt(2) line of tfile header. McClusky 941222
c                 BCTOT: Put in default orbital parameter list as data statement, modified
c                        TMAKE argument list (added icsnam). McClusky 941222
c                 TMAKE: Modified argument list, (added nics,icsnam), THDRED argument
c                        list modified (added nics,icsnam). McClusky 941222
c                 TMERGE: Modified argument list of THDRED and THDRIT (added nics,icsnam). McClusky 941222
c                 NGSTOT: Modified argument list of THDRIT (added nics,icsnam). McClusky 941222
c                 BCTOT, WSP3HD, TTOT, TTONGS: Minor changes to make ANSI standard.  King 950103
c                 NGSOPN: Fix bug in format for reading GPS week and day.  Banks/King 950103
c                 BCTOT, NGSTOT: Use fixed-dimension array for parameter names.  King 950106
c                 TTOG: T-file open duplicates sb gamke; remove.   Fang/King 950126
c version 9.29    Makefile:  Remove unused TIMECN (also delete from directories).  King 950210
c                 TMERGE: Modified to print headers of merged tfiles. McClusky 950316
c                 EST_FIT:  Fix minor (HP) compile errors.  King 950330
c                 TTOASC: If nics=8, set nics=9 to accomodate NOAA t-files.  King 950330
c                 Makefile:  Move INDEXX to /lib.  McClusky/King 950401
c                 TTONGS: Compiled using Lieske precession computation PT 950413
c                 TROT: add precession,nutation,gravity model arguments in thdred PT 950415
c                       pass true/false into rotsnp based on precession model used PT 950415
c                       pass prec,nut,grav models as arguments to thdrit PT 950415
c                 TTONGS: add prec,nutation, gravity model arguments in thdred call PT950415
c                 THDRIT: write prec,nut,grav models on tfile header PT 950415
c                 TMERGE: add prec,nut,grav arguments to thdred,thdrit calls PT 950415
c                 NGSTOT: hardwire prec (IAU76), nut (IAU80), grav (blank) models into t header PT950415
c                 TMAKE: hardwire prec (IAU76), nut (IAU80), grav (blank) models into t header PT950415
c                 ORBDIF: Correct call to THDRED to match 950415 update.  King 950418
c                 ORBRMS: Fix bug in skipping plot.  King 950420
c                 TERR_PAR: Add IAU76 logical to SROTAT call.  King 950420
c                 ORBRMS: change plot range calculation for rad/along/cross Tregoning 950425
c                 ORBRMS: allow an extension on plot files Tregoning 950425
c version 9.30    Makefile: Reorganised so that each module has its own library PT 950430
c                 ORBFIT: Major reorganisation of original code. Now orbfit can estimate satellite
c                         orbital differences for all satellites simultaneously, allowing for both
c                         IC adjustment and numerous terrestrial and inertial transformations. McClusky 950502
c                 EST_FIT: Modified to estimate all SV/GLOBAL parameters simultaneously. McClusky 950502
c                 READ_INPUT: Modified to read fixed/free parameters, and excluded satellites. McClusky 950502
c                 PLT_POSTFIT: Modified to handle new postfit residual array format. McClusky 950502
c                 TRAN_PART: Cleaned up code and comments (old name TRANPART). McClusky 950502
c                 GLOB_PART: Cleaned up code and comments (old name TERR_PAR). McClusky 950502
c                 ORBFIT,GLOB_PART: fix up calls to inert.-terr. routines Tregoning 950510
c                 TTONGS,NGSTOT,TROT,GMAKE,TMAKE,THDRIT: changes to allow B1950/J2000 and
c                       IAU68/IAU76  Tregoning 950510
c                 ORBFIT: Output the inertial frame, prec models in both tfiles PT950511
c                 EST_FIT: Fix bug in calculation of right hand side of normal equations. McClusky 950511
c                 PLT_POSTFIT: Fixed bug in creation of plotfile name McClusky. 950511
c                 NGSTOT: Fixed bug - declared "frame" variable  Tregoning 950512
c                 EST_FIT: Fix another bug in calculation of right hand side of normal equations. McClusky 950512
c version 9.31    GMAKE, NGSTOT, ORBDIF, ORBFIT, THDRIT, TMAKE,  TMERGE, TROT, TTONGS: Add radiation
c                   pressure model argument to THDRIT and lib/THDRED calls.  King 950517
c                 THDRIT: removed debug  Tregoning 950524
c                 GLOB_PART: Remove unused variable.  King 950525
c version 9.32    GMAKE, ORBDIF, ORBFIT, SDTRIT, TTOICS:  remove parameter statement for
c                    maxytp,maxyt2, now in dimpar.fti.   King 950612
c                 GMAKE:  Fix format for G-file header.   King 950619
c                 TTOT, Makefile: added earth <-> inertial tfile conversion  Tregoning 950620
c                 ORBFIT,EST_FIT: output estimated IC's to an svs file  Tregoning 950620
c                 ORBFIT: Allow "0" for no plot files   Tregoning 950620
c                 ORBFIT: Fixed bug in format statement. McClusky 950630
c                 PLT_POSTFIT, EST_FIT:  Write .rms file even if no plots.  Tregoning 950714
c version 9.33    ORBFIT, TTOT:  Remove duplicate declaration.  King 950719
c                 OPENB: Add status to open statement (required by DEC).  King 950719
c                 SDTRIT: Add branch to currently inaccessible messages for unexpected
c                     end-of-file.  Sanli/King 950719
c                 THDRIT: Make debug dump comments to satisfy DEC compiler. Sanli/King 950719
c                 BCTOT: Fix logic in setting initial PRN #.  Sanli/King 950719
c                 TTOT:  Stop if number ICs missing from T-file header.  McClusky/King 950719
c                 ORBFIT, READ_INPUT: Add undefined 'iscrn' to call list.  Sanli/King 950728
c version 9.34    GTOG: remove unused variable from argument list  Tregoning 950807
c                 OPEN_EOP,GTOG: pass unit numbers as arguments   Tregoning 950807
c version 9.35    TINPUT: Expand the span for the t-file from broadcast by an additional
c                     hour to accomodate the possibility of a yaw file.  King 950926
c                 TINPUT: Allow batch input with no X-file.   King 951010
c                 OPENB: Name print file 'bctot.prt' vice 'trotat.prt'.  King 951011
c                 New BCTOT using closest epoch:  Mods to BCTOT, BRDXYZ, Makefile;
c                   new routine CLOSEST_EPOCH; remove TMAKE.  King 951011
c                 TTOASC: Add times to output.   King 951011
c                 BCTOT: Now writes out an arc input file. McClusky 951012
c                 ARCINP: New Subroutine to write out arc input file. McClusky 951012
c                 GMAKE: Modified to pass gfile name out of routine. McClusky 951012
c                 NGSTOT: modified calling argument to gmake. McClusky 951012
c                 TTOG: modified calling argument to gmake. McClusky 951012
c                 BCTOT: Modified to write correct IC names on gfile, and
c                   read info off the sestbl. if available. McClusky 951012
c                 ORBFIT: Fixed bug in writing out gfile header. McClusky 951013
c version 9.36    TINPUT: Add receiver and antenna ids to lib/xhdred call.  King 951019
c                 Makefile: Remove ROTATE_GEOD, XYZ_TO_GEOD, DVDOT, DVSWP to /lib.  King 951019
c                 BCTOT:  Remove duplicate declaration.  King 951024
c version 9.37    BCTOT:  Remove duplicate declaration.  King 951103
c version 9.38    Makefile: Variable ranlib for compatibility with Solaris 2, and
c                     propagate 951019 changes into stdrel.  Fang/King 951208
c                 ARCINP: Fix compile bug for DEC  Tregoning 951214
c                 BRDXYZ: Fix compile bug for DEC  Tregoning 951214
c version 9.39    Added report_stat routine to: arcinp,bctot,brdxyz,check,chksvs,
c                 closeb.f,elem.f,eopart,est_fit,filnmb,getics,glob_part,gmake,
c                 gtog,hms2jd,ngsopn,ngstot,open_eop,openb,orbfit,oversn,plt_postfit,
c                 read_input,rotcrd,rsp1hd,rsp3hd,sdtrit,tdtrit,thdrit,tmerge,
c                 tran_part,trot,ttoasc,ttog,ttoics,ttongs,ttot,wsp1hd,wsp3hd. McClusky 960229
c                 Makefile: Add oversn.o for programs gtog and tmerge.  McClusky 960229
c version 9.40    TTOASC: Fix HP runtime bug writing IC's to dump file. Variable text.  McClusky 960312
c                 NGSTOT: Fixed bug. Nominal direct radiation pressure parameter now
c                 written out as 1.d0 not 0.d0. McClusky 960312
c version 9.41    NGSTOT:  Write T-file headers with 7-pt reduction for interpolation to avoid
c                    ungraceful abort in MODEL.  King  960319
c                 GETICS, NGSTOT:  Remove unused calling arguments.  King 960319
c                 TIMCHK (not currently used):  Change abort to use report_stat.  King 960319
c                 TDTRIT: Remove conversion by undefined GPS-UTC (not needed).   King 960319
c                 NGSTOT: Remove 7-pt reduction--MODEL/setup checks for this already. King 960327
c                 TTOASC: Modified dimension of data record array allowing velocity type
c                   T-file to be sucessfully read. McClusky 960424
c version 9.42    GLOB_PART: Fix erroneous calculation of sidereal time at leap-sec boundary.  King 960517
c                 Makefile: Fix typo in fti dependencies for NGSTOTLIB.  Tregonning/King 960529
c                 BCTOT: Remove redundant report_stat call.  King 960618
c version 9.43    ORBFIT: Fix bug in rounding seconds part of epoch in output G-file.  King 960709 
c version 9.44    Changes for new kf/gamit structure and Makefiles from unimake.  King 960718
c                     All routines:  removed trailing blanks and replaced lib/*.fti 
c                     includes by ../includes/*.h.  
c                 SPANUP:  Declare calling argument explicitly.  King 960815
c version 9.45    GMAKE: Fix typo in question about new IC time. McClusky 960916
c version 9.46    BCCHECK, Makefile:  New routine to check BC ephemeris elements.  King 961104
c                 BRDXYZ:  Change report_stat call from 'BCTOT' to 'ORBITS'.   King 961105 
c                 WRITEE, Makefile:  Add routine to write out a new nav file.  King 961107
c                 BRDXYZ:  Correct report_stat call.  King 961107
c                 BCTOT: Change call to (lib) READE for BCCHECK quantities.  King 961107
c version 9.47    G77 compiler changes:  
c                 PLT_POSTFIT, ORBRMS: Removed call to intrinsic function datan2d, as G77 does
c                    not support intrinsic trig functions whose argumants are in degrees. McClusky 970113
c                 FJLDAY, FJDTIM, ORBFIT, GMAKE, OPENB, ARCINP, ELEM, BCCECK, TTOG,I NVER2: Initialised 
c                   previously uninitialised variables, removed unused variables, and explicitly
c                   typed all unexplicitly typed variables. McClusky 970113   
c version 9.48    NGSTOT: modified to make IC time 12:00 if centre of sp3 file is within 1hr
c                   of midday. Also fiex bug in initialization of satics(9,*). McClusky 970317 
c                 NGSTOT: Trap and omit from t/g file any satellite on the sp3 file whose ICs
c                   do not satisfy reasonableness for semi-major axis (to catch 0.0 entries).  King 970320
c                 ORBFIT: avoided round off problem in epoch calculation (not fixed). McClusky 970317
c                 TTOT: fix bug in logic of call to thdred.  Tregoning 970320
c                 ARCINP, GMAKE, ORBDIF, ORBFIT, TMERGE, THDRIT, BCTOT, NGSTOT, TROT, TTONGS,
c                   TTOT: Replace efixed variable with better use of frame variable.  Tregoning 970321 
c                 NGSTOT, RSP1HD, RSP3HD: Fix bug in counting SVs, add check of maxsat.  King 970327
c                 THDRIT: Fix position of 'frame' in T-file header to match arc/wrthed.  King 970508
c                 BCTOT, OPENB, TROT : Use different print file name for OPENB calls from BCTOT
c                    and from TROT. King 970515/970516 
c                 INVER2: Fix program name in report_stat call.   King 970515
c                 NGSTOT: Add warning if days mismatched between SP3 and T files.  King 970516
c version 9.49    TTOICS: Fix declaration of buffer.  Fang 970805
c                 BCTOT, TIMESE, TINPUT: If no X-file, read session.info for start/stop time.   King 971028
c                 ARCINP: Fix bug in setting end-day.  King 980102  
c version 9.50    Add new program YAWTAB.  Tregoning 971207 / King 980106
c                   Moved from /model and modified:  BLOCK_IIR, CHKMODE, GET_DUSK, GET_ECL_POS, GET_YAW_RATE,
c                        OPN_YAW, READ_YAW, WRITE_YAW, YAW_CTRL.
c                   New: CALC_ATTIT, GET_ATT, YAWTAB 
c                   Mods: Makefile. 
c                 YAWTAB: Convert input to command-line arguments.  King 980106
c                 GET_ATT: Fix missing byte number for declaration of yaw_bias_tmp.  King 980107
c                 Add maxyaw to ../includes/orbits.h.  Use instead of 'yaw_entries' to dimension
c                   arrays in BLOCK_IIR, CALC_ATTIT, GET_ATT, GET_ECL_POS, GET_YAW, OPN_YAW,
c                   READ_YAW, WRITE_YAW, YAWTAB.  King 980107 
c                 READ_YAW: Initialize ecl_num.  King 980112
c                 Makefile, YTOASC:  Add program.  King 980112 
c                 YAWTAB:  Add T-file name, interval, and # epochs to y-file header.  King 980113
c                 YAWTAB, OPN_YAW, READ_YAW: Change t- and y-file names to character*16 for consistency with 
c                    rest of software.  King 980114  
c                 GET_ATTIT:  Remove vlight and twopi from calling sequence for lib/shadow1.  King 980121 
c                 YAWTAB, YTOASC: Add version number to y-file header.  King 980122 
c                 TINPUT:  Remove duplication declaration.  Morgan/King 980210
c                 YAWTAB, YTOASC:  Remove duplicate declaration.  Morgan/King 980223
c version 9.51    YAWTAB: Use two epoch intervals - one for calcs and one for output yaw table. Tregoning 980227
c                 BCTOT, GETICS, ORBFIT, TROT, TTOASC, TTOT, TTOICS: Fix bugs in report_stat call.   King 980227
c                 BLOCK_IIR:  Fix number of arguments in call to chkmod to check return to
c                     nominal mode.  King 980227 
c                 GET_ATT, YAWTAB:  Remove unused variables. 
c         9.52    GET_ECL_POS:  Write dusk times to GAMIT.warning; remove overwriting of y-file. Tregoning 980312
c                 GET_DUSK: Fix calling sequence for yaw_attit.  Tregoning  980312 
c                    (Note added 980703: This was missed in MIT version, not fixed until 980703)
c                 ARCINP, BCCHECK, BCTOT, BLOCKIIR, BRDXYZ, CHKSVS, ELEM, EST_FIT, GET_ECL_POS, 
c                    GET_YAW_RATE, GETICS, GLOB_PART, GMAKE, GTOG, INVER2, NGSOPN, NGSTOT, OPEN_EOP,
c                    OPENB, OPN_YAW, ORBFIT, PLT_POSTFIT, READ_YAW,  RSP1HD, RSP3HD, SDTRIT, TDTRIT,
c                    TMERGE, TROT, TTOASC, TTOICS, TTOT, WRITE_YAW, YAWTAB, YAW_CTRL,  YTOASC: Changes to
c                    implement use of binary tabulated yaw file; change report_stat calls to use program name for 
c                    module and orbits/[subroutine] for source.    King 980313
c         9.53    YAWTAB:  Correct bug when int < 30 and equal to X-file interval.  Tregoning/King 980401  
c         9.54    BLOCK_IIR: Modified to do calculations fro Block IIR satellites in a RHS; matches
c                      changes in /lib/yaw_attit.  McClusky 980521
c                 CALC_ATTIT: Debug format change. McClusky 980521
c         9.55    YAWTAB: Add SV antenna offsets to svnav_read call.  King 980615
c         9.56    GET_DUSK: Fix calling sequence for yaw_attit (rwk missed at 980312)
c         9.57    BCTOT: Add BERN2 radiation pressure model.  King 980708
c                 ORBFIT,EST_FIT:  Use g/t-file parameter names for ICs and rad-press.  King 980709`
c                 BCTOT: Add BERN1 radiation pressure model.  King 980714
c                 ORBFIT: Get correct rad-parm names from tfiles.  King 980716
c                 PLT_POSTFIT:  Remove extraneous comma in format (SGI IRIX warning).  Morgan/King 980901
c                 GPGT: Remove implicit statement to satisfy SGI compiler.  Morgan/King 981231  
c         9.58    BLOCK_IIR1: comments added to old routine block_IIr Tregoning 990201
c                 CALC_ATTIT: new logic to handle 2 modes of IIr behaviour Tregoning 990201
c                 ELEM: Fix report_stat call.  King 990205
c         9.59    ORBFIT: Modified to allow estimation of terrestrial rotation rate parameters. McClusky 990219
c                 READ_INPUT: Now reads rotation rate fixed/free entry from the orbfit.cmd file.
c                             Default: t_rat 0 0 0 (fixed not estimated). McClusky 990219 
c                 EST_FIT: Now allows estimation of terrestrial rotation rate. McClusky 990219 
c                 BCTOT,BRDXYZ: Pass PRN number to BRDXYZ. Bock 990223
c                 ORBFIT: Call topens to open t-files (to get machine-architecture check). King 990305
c                 NGSTOT: Clear arrays to avoid bad ICs with PC Solaris.  Fang 990308 
c         9.60    TROT: Turn on short-period UT1 and pole terms, to be used by BCTOT, NGSTOT,
c                   and TTONGS.  King 990310
c                 TROT: Change short-period UT1/pole model to 'Ray'.  King 990319
c                 BCTOT: Do not use comma as continuation character (IBM and f90 restriction).  Fang/King 990324
c         9.61    YAWTAB,GET_ATT: initialise some variables to zero (HP problems) Tregoning 990614
c         9.62    Use library functions for JD conversions: Remove FDJTIM, FJLTIM, HMS2JD, MDYJUL, change Makefile,
c                   GMAKE, ORBFIT, RSP1HD, RSP3HED, TDHRIT:  Use library routines for Julday.  King 990728
c                 Changes for 4-digit years: BCCHECK, BCTOT, GETICS, GMAKE, NGSTOT, RSP1HD, RSP3HD, READ_YAW
c                   SDITRIT,THDRIT, TIMESE, TTOASC, TTOICS, WRITEE, WRITE_YAW, WSP1HD, WSP3HD.  King 990728
c                 ORBFIT:  Fix naming of .rms, .fit, and plt_ files when user leaves off decimal in name.  King 990805
c                 TTONGS: Make question about output file name more truthful.  King 990805
c                 TINPUT: Add argument to lib/rsesfo call (for MAKEXP).  King 990817
c                 ORBFIT: Fix bug in Y2K code.  McClusky/King 990819
c                 ORBFIT: Initialize interpolation variables (abort on Alpha/Linux).  King 991029
c         9.63    ARCINP: Fix Y2K bug.  Fang/King 000103
c                 GMAKE:  Fix Y2K bug.  McClusky/King 000103 
c         9.64    SPANUP: Fix G77 intolerance of 'idnint' . McClusky 000104
c                 NGSTOT: Fix bug when more than one SV has bad entries.  Sugimoto/King 000202
c                 SPANUP: Re-fix g77 bug, using 'int' (truncation, like 'idnint') not 'nint' 
c                    (rounding).  Voci/King/McClusky 000306 
c         9.65    TTONGS, RDSVCLK (new), Makefile:  Read a clock file generated by autcln
c                    and write values on to the sp3 file.  King 000330      
c                 WSP3HD: Fix format for 4-digit GPS week.  King 000330.
c                 TROT: Convert screen print to report_stat.  King 000330
c                 TTONGS, RDSVCLK, WSP3HD:  Fix declarations.  King 000401
c                 BCCHECK: Fix initialization of nsat.  Voci/King 000503
c                 TTONGS: Fix help line for command-line arguments.  King 000905
c         9.66    NGSOPN: Trap sp3 file name greater than 16 characters.  King 001010
c                 GET_DUSK:  Fix time argument in all to solred.  Tregoning/King 001011
c                 SDTRIT: Fix default clock-value overflow.  Fang 001012
c         9.67    ORBFIT:  Major changes to use with multiple T-files.  King 001220
c                   Changes to: Makefile, ORBFIT READ_INPUT, PLT_POSTFIT
c                   Remove:  EST_FIT,GLOB_PART
c                   Add:  KEPLER, NORMINC, NORM_SOLVE, PARTL, WRITE_G, WRITE_SUMMARY, XYZ2RAC
c                   New include:  orbfit.h
c                 NGSTOT, GETICS, NGSOPN,TDTRIT, RSP1HD, RSP3HD: Allow command-line arguments,
c                      no screen print.  King 001220
c                 ORBFIT, PARTL, WRITE_G:  Correct duplicate variable declarations.  King 001226
c                 PLT_POSTFIT:  Remove unused variable.  King 001226
c                 BCCHECK, TDTRIT: Fix uninitialized variables for HP compiler.  King 010103   
c                 NGSTOT: Fix bug in shifting ICs when a bogus SV removed.  King 010103 
c         9.68    YAWTAB, GET_ATT:  Initialize and pass two variables to avoid their being undefined
c                   (and potentially non-zero on the HP)--'ji0' to avoid gsatel miscalculation at tabular
c                   boundaries, 'yaw_bias' to avoid discontinuity at beginning of eclipse.  King 010110
c                 YAW_CTRL: Save 'last_shad' to avoid step and use of first yaw at beginning of
c                   eclipse.  King 010111
c                 YAW_CTRL:  Change position of initialization to avoid undefined 'last_time_t'
c                    (no apparent effect but dangerous).  King 010111
c                 GET_ECL_POS: Fix logic for Block IIR SVs to avoid calculating dusk time 
c                     (gsatel fatal) when there is none within the arc; remove unnecessary Y2k
c                     fix after call to jul2hma.  King 010115
c         9.69    BCTOT, CHKSVS, SPANUP, TIMESE, TINPUT: Add command-line option; clean up code, fix three bugs.  King 010122
c                 GMAKE: Allow overwriting of old g-file in batch mode.  King 010123
c                 NGSTOT: Remove extraneous debug print of 'ians'.  King 010129
c                 OVERSN, BCCHECK, BCTOT, GTOG, NGSTOT, ORBDIF, ORBFIT, TMERGE, TTOASC, TTOG,
c                     TTOICS, TTOT, YAWTAB:  Increase size of 'version'.  King 010202
c                 WSP3HD: Modify first-line format to add the format version and position code
c                     to conform to IGS standards. Fang 010215
c                 WRITE_SUMMARY: Write to ORBITS.status the overall fit of sp3 file.  King 010219  
c                 ORBFIT, orbfit.h: Increase, but also trap, maximum epochs on combined T-file. King 010505
c                 RSP1HD:  Add 'batch' to calling arguments (mismatch with ngstot.  Tregoning/King 0100518 
c                 ORBFIT: To avoid DEC/OSF1 compiler warning, comment out trap of maxtfil and maxorb too 
c                    large for unit assignments.  Tregoning/King 010518   
c                 ORBFIT, PARTL:  Replace undefined variables in calculation of tdiff for eop partials.  King 010518
c                 BCCHECK, WRITE_SUMMARY:  Fix non-portable hollerith continuation.  Tregoning/King 010518  
c                 ORBFIT.H, ORBFIT:  Make substitution of orbfit.h for maxepc promised at 010505 but overlooked; 
c                    change comments in orbit.f. Tregoning/King 010613
c                 TINPUT: Skip callling spanup if input dates all zero (interactive mode, no xfile or session.info) King 010618
c                 BCTOT: Fix dimension of 'satprm' from maxsat to 6; assign precmod and nutmod in get_models. King 010619
c                 ORBFIT: Use rcpar rather than getarg to avoid problem with HP UX-11.  Herring/King 010718  
c                 ORBDIF:  Remove undefined 'satname' and unused 'satnam'. Morgan/King 010731  
c                 TTOICS:  Fix undefined 'iyear'.  Morgan/Kng 010731 
c                 SDTRIT:  Fix undefined variable in error message.  Morgan/King 010731  
c         9.70    READ_INPUT:  Fix uninitialized 'nsvs_excl'.  Charade/King 020627  
c         9.71    RSP1HD: Fix report_stat calls to reflect correct subroutine.  King 021206
c                 NGSTOT, GETICS, RSP3HD, TDTRIT:  Changes for sp3-c format; fix check of
c                   files name vs date.  King 021206 
c         9.72    SDTRIT, TTONGS:  Increase length of string for clock file name; correct units on
c                   clock output.  Herring 030210  
c                 RDSVCLK: Allow two entries (value and sigma) for clocks.  Herring 030318
c         9.73    TTONGS, SDTRIT: Allow command-line or interactive input of start, stop times. King 030814
c                 ARCINP: Remove superfluous character beyond col 72 (ifc warning). King 031027  
c                 ORBDIF: Fix bug in setting frame from SP3 header.  King 031027 
c                 TMERGE: Fix length of 'status' variable.  King 031027
c         9.74    TDTRIT: Fix bug in reading the coordinates line of an SP3-C file.  Munekane/King 040129
c         9.75    ORBFIT, READ_INPUT, WRITE_SUMMARY, ORBFIT.H: Now automatically edits bad satellites
c                        if the max_fit_tol entry in the command file is used. McClusky 041201
c         9.76    ORBFIT, ORBRMS:  Change a fatal to warning for different SP# reference frames;
c                    change a format statement.  Herring 041218
c                 WRITE_G: Generalized code for writing svs_file for GLOBK to include write entries for
c                    for all possible radiation-pressure parameters.  Herring/King 041220
c                 READ_YAW: Add test on 'M' (for lunar eclipse), though no difference in behavior
c                    coded.  (Temporary change?) King 050103
c                 WRITE_G, WRITE_SUMMARY: Remove debug and add blank lines in rms file. Herring 050104 
c         9.77    TINPUT: Remove 'extra' variable from call to lib/xhdred. King 050129
c         9.78    YAWTAB, CALC_ATTIT, GET_ECL_POS, WRITE_YAW:  Treat iblk=5 (as well as 4) as a 
c                    Block IIR satellite.  King 050707
c                 YAWTAB, CALC_ATTIT, GET_ECL_POS: Allow iblk=6 for Block IIR-M SVs.  King 051027
c         9.79    ORBFIT, orbfit.h: Increase maximum epochs by 6 to accommodate extra epochs now used for
c                    integration.  King 060131
c                 TDTRIT: Fix bugs in reading SP3-C files.  Murray/King 060202
c                 ORBDIF, RDSPHD : Fix error in reading PRN list.  Herring  060202
c         9.80    TINPUT: Change arguments for rsesfo (see lib/lversn.f) King 060815
c                 TDTRIT: Fix data type error reading SP3-C pos/clock uncertainty fields.  McClusky 060829
c         9.81    NGSTOT, RSP3HD, TDTRIT, TTONGS, SDTRIT, WSP3HD:  Apply center-of-mass corrections for ocean tides. King 061030
c                 TTONGS, WSP3HD, SDTRIT: Fix header for CMC and allow NONE.  King 061105
c         9.82    ORBFIT, PARTL(comment only) , READ_INPUT:   Changes to allow > 30 SVs.  King 061127
c                 RSP3HD, TTONGS, WSP3HD: Change SP3-C format for line 22 as per IGS msg.  King 061129 
c         9.83    TINPUT: Change arguments for rsesfo (see /lib/lversn.f).  King 061207
c                 OVERSN: Remove call to clrscr; remove print unit number from 'lversn' call.  King 070416
c                 GTOG: Remove unneed arguments from call to 'rot_gfile'.  King 070416      
c         9.84    TDTRIT: Fix improper declaration of 'sigx'.  Balandin/King 070618 
c                 ORBFIT: Use a user-specified rather than generically named command file.  McClusky 070713
c         9.85    GTOG, PARTL,TROT: Remove unused 'iscrn' from calls to /lib routines pnrot, rotsnp, and srotat.  King 070906
c                 BCCHECK, BCTOT, CLOSEB, GTOG, NGSTOT, ORBDIF, ORBFIT, OVERSN, TMERGE, TTOASC, 
c                    TTOG,TTOICS, TTONGS, TTOT, WSP3HD:   Remove unused screen & print units in calling sequence.  King 070910  
c                 BCTOT, EOPART, FILNMB, NGSOPN, OPENB, OPEN_EOP, ORBDIF, ROTCRD, TMERGE, TRAN_PART: Remove unused label.   King 070910    
c                 BRDXYZ, BCTOT, EOPART, GET_DUSK, GET_ECL_POS, GMAKE, NORMINC, ORBFIT, TDTRIT, TIMESE, 
c                    TINPUT, TTOG, YAWTAB: Remove unused calling argument.  King 070910   
c                 TDTRIT: Add dummy print statement for unused calling arguements.  King 070910        
c                 TTOT: Fix format for reading selection.  King 070910  
c                 Makefile, WRITE_YAW: Remove unused routine.  King 070910
c                 PARTL: Fix bug in calling argument for eopart.  King 070928   
c                 WSP3HD: Remove tab to satisfy gfortran. Herring 071103     
c         9.86    GMAKE, GTOG, NGSTOT, ORBFIT, WRITE_G, orbfit.h : Accommodate nutation model and 
c                    gravity model on g-file.  King 071221/071227   
c                 YAWTAB, READ_YAW, GET_ECL_POS, GET_DUSK:  Don't calculate dusk time for 
c                    Block IIR if a lunar eclipse. .  King 080129        
c                 YAWTAB, READ_YAW, OPN_YAW:  Read beta angle from y-file.  King 080130      
c         9.87    YAWTAB: Fix calling arguments for oversn; remove unneeded units common. King 080212
c                 NGSOPN: Change logic for getting GPS week and day, to work with igu orbits. King 080321
c                 WRITE_G: Remove colon from end of message.  King 080321
c                 BCTOT:  Fix bug in reading brdc file (getting only one SV); get nutation model from 
c                   the file header.  King 080428        
c         9.88    TINPUT: Change antenna offset arguments for /lib/xhdred.  King 080509  
c                 ORBDIF TTEMP: Create unique temporary t-file names. McClusky 080514
c                 WRITE_SUMMARY: Add a trap on zero change in chi2; scale; change units of velocity
c                   print from m/s to mm/s.  Herring 090227
c         9.89    TTONGS, WSP3HD:  Change header entries from sp3-b to sp3-c.  King 090908  
c                 WSP3HD, SDTRIT: Fix bugs in format.  King 090925    
c         9.90    SP3TOT, SCANSP3, Makefile: New version of NGSTOT to handle bad SVs.  King 091015
c         9.91    YAWTAB, CALC_ATTIT, GET_ECL_POS: Allow iblk=7 for Block IIF SVs.  King 100529    
c                 WRITE_YAW: Write time and attitude for BLock IIR-M (missed in     update)
c                    and II-F.  King 100529
c                 BCTOT: Allow UCLR1 radiation-pressure model.  King 110124
c                 RSP3HD, SDTRIT, WSP3HD:  Fix bug in test for use of OTL mode, in
c                   writing the PCV comment in in the SP3 header, and in assigning
c                   a clock value when no clock file is available. King 110416  
c         9.92   NGSTOT, GETICS, RSP3HD, TDTRIT: Allow sp3 file to contain non-GPS satellites,
c                   but skip reading of these. King 120317 
c                RSP3HD: Correct bug in last change when SV id is blank rather than 'G'.  King 120322
c                GTOG: Use lib/read_gfile rather than (now obsolete) /lib/ghdred. King 120711
c         9.93   ARCINP, BCTOT, GET_ECL_POS, NGSTOT, ORBDIF, ORBFIT, THDRIT, TMERGE, 
c                   TROT, TTOT, WRITE_G, YAWTAB, orbfit.h:  Add Earth radiation and antenna 
c                   radiation to reading and writing g-file and t-file headers. King 130329
c                TTNOGS: Fix bug in calling thdred.  King 130426
c         9.94   GTOG, NGSTOT, ORBRMS, TTOT: Fix missing field separate for $ in format statements to 
c                   satisfy the IBM AIX compiler. King 131011
c         9.95   BCTOT, ARCINP, GET_MODELS:  Add earth and antenna radiation models. King 131217
c                SORTBC, TRAN_PART: Add dimpar.h to get maxsat for maxyawsv in orbits.h.   King 131231 
c         9.96   YAWTAB, KOUBA_YAW (new), SVSUN_ANGLES (new), YTOASC, Makefile: implement the 
c                  Kouba yaw model, replacing the old Bar-Sever model; remove obsolete 
c                  routines block_IIr, calc_attit, chkmode, get_att, get_dusk, get_ecl_pos, 
c                  get_yaw_rate, opn_yaw, read_yaw, write_yaw, yaw_ctrl.  King 140125
c                YAWTAB, KOUBA_YAW, Makefile; remove obsolete get_yaw; remove ../includes/yawtab.h.  King 140218 
c                YAWTAB: Issue a warning of a mis-matched bias only once.  Herring/King 140327
c         9.97   Changes for GNSS:  ARCINP, BCCHECK, BCTOT, GMAKE, NGSOPN, NGSTOT, ORBDIF, ORBFIT, READ_INPUT, RSP1HD, 
c                   RSP3HD, SDTRIT, THDRIT, TINPUT, TMERGE, TROT, TTOASC, TTONGS, WRITE_G, WSP3HD, YAWTAB. King 150108
c                TROT: Change to use IERS2010 rather than Ray model for short-period UT1/pole. King 140325
c                TTONGS: Fix bug in calling arguments for thdred.  Herring/King 150502
c         9.98   NGSTOT, YAWTAB: Add start/stop times to svnav_read call.  King 105020
c                BCTOT: More mods needed to handle GNSS correctly.  King 150615
c         9.99   RSP3HD, Makefile: Move rsp3hd to /lib for use by /makex.  King 150730
c                NGSTOT: Fix format for writing Beidou satnam; allow for SV altitudes to 50,000 km.  King 150922 
c        10.00   GHDRED, WRITE_G: Remove obsolete SRP model parameter names, change BERN2 to new model. King 151021 
c                NGSTOT, TROT: Enable the reading of short-period EOP from sestbl., instead of original hardwire code. Wang 151022
c        10.01   NGSTOT, GETICS, TDTRIT: Allow reading of any (single) set of GNSS SVs.  King 151103
c        10.02   KOUBA_YAW, YAWTAB, YTOASC: Replace GPS block #s by SV body type. King 151112
c        10.03   Makefile, BCCHECK, BCTOT: Remove bccheck (incompatible with current GNSS code
c                 and rarely used; match arguments for modified lib/reade.   King 151119   
c                YAWTAB: Fix bugs in writing y-file and yawtab.out. King 151124
c        10.04   BCTOT, GMAKE, TTOG,  Makefile (remove TINPUT and TIMESE).  King 151210
c                KOUBA_YAW: Add code for Beidou. King 160105
c                YAWTAB: Change check on coded yaw to include Beidou.  King 160107.
c                BCTOT: Check iflag from lib/reade for bad record. King 160218
c                NGSTOT, GETICS, TDTRIT: Clean up SV selection from SP3 file.  King 160602
c                BCTOT, ARCINP, BRDXYZ, OPENB: Add screen print to help find errors; simplify time 
c                   computations. King 160804
c                CLOSEST_EPOCH: Change obsolete 'suicid' call to report_stat fatal. King 160805
c                TROT: Fix file name of program in report_stat warnings. King 160812.
c                TTOASC: Trap nsat=0 before creating print formats. King 160812
c                NGSTOT: Convert gnss_sel to uppercase.  King 160812
c                YTOASC; Fix bug in reading current y-file version. 160822
c                YAWTAB: Set a place-holding yaw=0 to allow testing of IRNSS. King 160822
c         10.05  ORBRMS: Fix format for writing maxdneu.  King 170411 
c         10.06  SP3TOT, TROT: Add pole-position in call to lib/rotsnp.  King 170412 
c         10.07  YAWTAB: Fix '0' PRN echo in status message. King 170505
c                YAWTAB, KOUBA_YAW: Pass y-bias flag rather than value and compute the
c                    value within the yaw routine. King 170505
c                YAWTAB and new routines KOUBA_GPS, KOUBA_GLONASS, KOUBA_BEIDOU, KOUBA_GALILEO,
c                    MU_ANGLE, Makefile:  Use separate routines for each GNSS; remove the sign reversal
c                    for GPS IIRs; new code to handle beta angles < 0.07 (GPS, Glonass, Galileo);      
c                    code Galileo from Kouba and web notes, not eclips_May17.  King 170522
c                ORBFIT, WRITE_SUMMARY: Allow no estimation (can now replace ORBDIF). King 171007
c                NGSTOT: Fix failure when the sp3 file header has otlmod NONE.  PFang/King 171107
c         10.08  ORBDIF, WRITE_SUMMARY: Added 3-D RMS to outputs, increased significant digits 
c                    in summary outputs and added MEAN output in ORBDIF.  Removed blank lines in 
c                    ORBDIF.   Added runstring option to ORBDIF. Herring 180115.
c                NGSTOT: Allow sp3d.  King 180206
c                YAWTAB: Added 'isat' argument to call of kouba_glonass.  Herring 180319 
c         10.09  GHDRED, GTOG: Changed start time variables from jd0,t0 to jdb,tb to be consistent
c                    with includes/arc.h and other modules. King 180320       
c                YAWTAB: Add new common includes/units.h to provide unit numbers removed from 
c                    includes/arc.h. King 180320
c                BCTOT, CLOSEB, NGSTOT, OPENB, OPEN_EOP, ORBFIT: Open the nutation table only if needed. King 180329
c                YAWTAB: Remove data/neweph/ in favor of check on 'nbody' present. King 180329 
c                YAWTAB: Fix bugs: missing fjd argument in call to ephdrd, missing call to evrtcf; missing isat 
c                    argument in call to kouba_glonass.  Herring/King 180404 
c                BCTOT: Fix format for setting SV name.  King 180417
c                YAWTAB: Fix undefined fjd when using old luni-solar code. King 180425 
c                YAWTAB: Correct typo in comment on new version. King 180518
c          10.10 GHDRED: Allow ECOM1 and ECOM2 as SRP models.  King 181108
c                BCTOT: Fix nutation model assignment to IAU00 for MHB_2000 case.  King 181108
c          10.11 KOUBA_GPS: Temporarily make GPS Block IIIA behave like Block IIF. King 190109
c          10.12 GHDRED, ORBFIT READ_INPUT    : Allow 13-parameter ECOM2 radiation-pressure model.  King 190425 
c                NGSTOT: Change initial SRP model from 'SPHRC' to 'ECOM1'.  King 190430
c                GHDRED, WRITE_G: Increase radiation parameters to 19. King 190524
c                TTOASC: Fix IC format to work for 13 SRPs.  King 190611
c                ORBRMS: Changed rms output format to be consistent with orbfit Herring 190625.
c                NGSTOT/YAWTAB/BCTOT: Added antpwr(/x) to svnav_read call Herring 190702.
c          10.13 Several routines listed below are modified to allow use of SOFA based IAU2000A and IAU2006
c                interial reference frames. McClusky 190801
c                BCTOT: Allow inertial frame  to be read from sestbl. Default is new IAU2000A frame
c                GTOT: Only a small bug fix in writing out precmod
c                SP3TOT: Changed default inertial frame to SOFA IAU0A (IAU2000A). 
c                        We should fix thisroutine to read sestbl at some point?
c                NGSTOT: Allow inertial reference frame to be read from sestbl. Default is new IAU2000A frame.
c                TROT: Inertial reference frame default is new IAU2000A frame. (Should fix to read from sestbl.)
c                PARTL: Added call to rotsnp_sofa to get prec/nut/sidtm vales for new IAU SOFA referecne frames
c                EOPART: Just debug changes. 
c          10.14 YTOASC: Add number of SVs to ascii printout.  King 190927
c                Makefile, Y2ORBEX, QUATERNION, WRITE_ORBEX: New routine to convert the y-file yaw angles 
c                       into quaternions and write these to an ORBEX file.  King/Herring 190927 
c                NGSTOT: Set the precession model default to IAU76 if there is no sestbl. King 191001
c          10.15 READ_INPUT/WRITE_SUMMARY: Changed output to single list when max_fit_tol exceeded.
c                        Fixed reading ofsvs_exclude.dat (McClusky/Herring 200121). 
c                YAWTAB: Replaced '32' with '200' in print statement for SVs. King 200203
c                KOUBA_GPS: Turned off debug being turned on.  Herring 200211.
c                NGSTOT: Set J2000 default precession to IAU0A. Herring 200304.
c          10.16 ORBITS.H Increased to 9 sp3 files to be fit Herring 200331.
c                ORBFIT: Modified to handle up to 9 days of data as individual days (+-1.25 hrs missed each
c                       day or a continuous 9 day SP3 file. Herring 200331.
c                PLT_POSTFIT: Added MJD to output residuals so time in known Herring 200331.
c                ORBRMS: Added MJD to output residuals Herring 200331.
c                ORBDIF: Passed refepoch (PEPJD) into orbrms so that MJD can be output. Herring 200331.
c          10.17 YAWTAB.F: Added one day to date for call to svnav_read because the start of the t-
c                          file may occur before the start time of a PRN record. MAF (2020-04-14, MIT)
c          10.18 GET_SVCOORDS.F, PARTL.F, SP3TOT.F,TROT.F: Replaced arbitrary setting of iut1pol bit mapped 
c                        diurnal/semidiurnal models with value read from sestbl. TAH 200505
c                READ_SP3, WRSP1HD: Updated 32I to 50I to allow for 35 Beidou satellites TAH 200618
c                WRITE_SUMMARY: Increased bad list format to allow for more than 10 bad satellites TAH 2007006.
c                ORBDIF: Deleted temporary t-files on exit; cleaned up debug.  TAH 200714. Modified to
c                        only when SP3 files input. TAH 200810.
c                ORBDIF: Initialized rmsnam to stop null in output when not specificed. TAH 201020.
c                ORBRMS: Initialized ext to stop null in output when not specificed. TAH 201020.
c                NGSTOT: Added explicit ends to antbody string in generating satnam to avoid errors with 
c                        some compilers for long antbody strings. TAH 201031.
c                KOUBA_GPS: Cleaned yrate code indenting and variable use.  TAH 201111.
c          10.19 SDTRIT,RDSVCK: Added GNSS to call and reading so correct clock added from multi-gnss
c                        clock file generate with merge_igs_clk (see sh_mergesp3 also). Added EOF
c                        line at the end of SP3 file, TAH 201211.
c          10.20 GETICS: Updated 2020 to 2100 as end date. TAH 210101
c          10.21 TTONGS, WSP3HD, SDTRIT: Updated to pass sampling rate to allow (integer) higher sampling
c                        rate than t-file (sample rate 3 output 5 minute sp3 from 15 minute t-file)
c                        TAH 210122.
c                ORBDIF: Removed deleting 2nd t-file if t-file rather than sp3 file used. TAH 210201
c                YAWTAB: Updated to write enough epochs to handle 5-minute t-files. TAH 210202
c
      RETURN
      END

