      SUBROUTINE LVERSN(VERS)
      CHARACTER*45 VERS
      CHARACTER*10 GETMAC,machin
      INTEGER NBLEN
C.....
      machin = getmac(1)
      vers = ' '
      WRITE (VERS,5) machin(1:nblen(machin))
    5 format ('11.47 of 2021/07/30 21:48 UTC (',a,')')

C Version 8.1       First version with LVERSN documentation
c                   Comments from Kurt's updates 90/11/12:
c                   IMENU: Allow only 5 mistakes
c                   NSNPRN: Scripps update
c                   PICKFN: SR 9.7 version
c                   READC2, READC4, READC5, WRITC2, WRITC4, WRITC5: New C-file format
c                   READJ: Read the whole thing into memory
c                   TOPENS: Larger dimension
c                   XHDRED: Read receiver software version
c Version 8.2       RRXHED and RRINEX - Modify for RINEX Ver 2 - rwk 91/01/22
c version 8.3       RRXHED: correct read of character*40 buff40.
c                   IEDIT:  make batch work for 120 characters.
c version 8.4       LBLK1:  Fix bug confusing PRNs 2 & 3 - rwk 91/03/02
c version 8.5       TAIUTC: updated for Leap second Dec 31, 1990 Kurt&King
c version 8.6       GETARG: dummy routine added for Apollos to
c                            perform the sun library function
c version 8.7       READE:  do not read a comment line as data
c                           Kurt 91/04/08
c version 8.8       READJ:  Suppress warning of duplicate records - rwk 91/4/9
c                   DIMPAR: Increase MAXDAT to 5 to handle NOBTYP on RINEX
c version 8.9       RRXHED: Allow lower-case keywords in RINEX header - rwk 91/5/24
c version 8.10      ERRFLG: Add flag for good L1, bad L2 to use on short lines
c                            - yb (by rwk) 91/6/30
c version 8.11      ICLARG: New function to read command-line arguments
c                            - Kurt 91/7/2
c version 8.12      WRITC5: Set azmdot and atmdel to zero until MODEL
c                           actually calculates them.
c                           - Kurt 91/7/18 for Yehuda.
c version 8.13      XHDRED: Read antenna height to 1/10 mm.
c                           Should be backwards compatible.
c                           - Kurt 91/7/18 for Yehuda.
c version 8.14      TAIUTC: Update stop date to 31 Dec 1991
c version 8.15      THDRED and SOD2JD **to be added from NGSTOT - King 91/11/8
c                   TIMINC added from NGSTOT
c                   READC4, WRITC4: mods for kinematic
c                   XHDRED: Read another significant place for radius (0.1 mm)
c                            - yb 91/08/25
c                   Modify TAIUTC to read leap.sec file  Yb 91/11/08
c                   Reinstate TAIUTC that does not need table   YB 91/11/16
c                   READC2, WRITC2: Change 'phi' to 'psi' for nutation in longitude.
c                            - king 91/12/18.
c                   New routine UPPER1 to gradually replace UPPERC - King  91/12/30
c                   Mods to LREAD for old-style L-file and to use UPPER1, UPPERC
c                            - king 91/12/31
c                   PICKFN: Minor change to add spacing to messages.  King 91/12/31
c                   DS2HMS, TIMCON, and new LEAPYR to handle leap-year
c                           problems - Bock and King 92/1/2
c                   RRINEX: Pass loss-of-lock indicator - Oral and King 92/1/23
c                   Add  UT1RED, UT1TID, FUNARG, FUNCOF from model/sccs - King 92/1/23
c                   TIMCON: Fix sow to sod bug when no UTC/GPS conversion - King 92/1/23
c                   GSATEL added from orbits/sccs; in MODEL call to SATEL replaced
c                          by call to GSATEL (model/sccs/satel was identical to
c                          orbits/sccs/gsatel after  92/12/96 mods).  King 92/1/27
c                   TAIUTC: Replace with table-reading version that keeps the leap.sec
c                          file open and reuses values when possible.  King 92/1/19
c                   Add POLRED, NUTRED, and SIDTIM from /model/sccs.
c                   Add modified PRCES, SIDMAT, and NUTTAB from /orbits/sccs.
c                   Add TRANSP, MATMPY, PNS, SNP, and ROTMAT from /orbits/sccs.
c                           King 92/1/31
c                   Add HISUB from /makexp/sccs, changing calling arguments to R*8.  King 92/2/25.
c                   NUTRED:  Add END= branches for both header and data reads. King 92/2/14
c                   NUTRED:  Fix bug in backspace logic.  King 92/2/20
c                   TAIUTC:  Fix bug in saving values.   King 92/2/20
c                   Add RSTNFO (new routine for station.info) to library.  King /92/2/20
c                   These changes copied from MIT to SIO 92/2/26 (?)
c version 8.16      NUTRED:  Fix finally bug in backspace logic.  King 92/2/27
c                   RSTNFO:  Read in format line before description line from station.info
c                            Bock 92/2/28
c                   RSTNFO:  Read in extra record from station.info for type=2, icall=3
c                            Fix year,day checking logic
c                            Bock 92/03/01
c                   XHDRED:  Add SKD variable (Static/Kinematic/Dynamic)
c                            Bock 92/03/06
c version 8.17      RSTNFO:  If session number = 0, set = 1.   King 920324
c                   WWRITE:  Remove CWRITE call in Sun version. Now same as Apollo. tah/rwk
c                   XHDRED:  Set SKD to 'S' for old-format case.  King 920331
c                   Recompile for modkin.fti change:  LREAD, READC4, WRITC4
c                            Add modkin dependency to Makefile.    King 920331
c                   GSATEL:  Fix bug in check of large values of Y-vector.  King 920402
c                   RSTNFO:  Increase size of format to 90 characters.  King 920414
c                   HISUB :  Add station and date for error message.  King 920415
c     include '../includes/dimpar.h' MAXNET from MAXSIT, change MAXSIT to MAXNET in
c                                setting MAXPRM.  Reorder and add comments.  King 920423
c                   READC2, WRITC2:  Add SKD and IRCINT to header record.  King 920423
c                   Remove unused variables identified by Sun compiler:
c                          DLINPK  GSATEL  HISUB  PRCES  READE  RRINEX  TAIUTC
c                          TIMCON  WDMS    WTIME  WRXHED RHDRED   King  920424
c                   Add RSESFO (from fixdrv) to lib. Change Makefile
c                   READC2, WRITC2, XHDRED:  Add session number to header.  King 920429.
c                   RSTNFO: Fix bug in line read for version 2.   King 920429
c version 9.1       READE:   Update idumb to check PRN numbers > 24
c                            (I've used 37 which I believe is the maximum PRN number)
c     include '../includes/dimpar.h'AT=18, MAXSIT=37=MAXNET
c                   Makefile : remove cview.fti (not used)
c                            Bock 92/05/03
c version 9.11      RSESFO : If SN = 0 => SN = 1
c                            Bock 92/05/07
c         9.12      Fixed bug in IERRCM - if statement for igl2bd was missing
c                   Shimon 92/5/23
c                   READC4,WRITC4 : change debug D to comment C
c                   Bock 92/05/25
c version 9.13      TIMCON : Remove INT4 for Sun compatibility.   King 92/05/08
c                   RSTNFO : Fix bug in reading daily entries.    Murray 92/06/03
c                            Read values only if column 1 blank.  King 92/06/04
c                   LREAD : Hard stop if site not found in l-file Bock 92/06/07
c                   RSTNFO: Fix read of headers to agree with 92/06/04 change.
c                              Bock/King 92/06/08
c                   HISUB : Hard stop if height measurement type not found in station.info.
c                              King  92/06/11
c version 9.14      MAPAMP: Add Rogue and Ashtech mappings.  Oral/King 92/06/22
c                   READE : Add calling argument to make more robust.   King 92/06/24
c                   HISUB : Correct bug DHTGP-->SLTGP for TRMSLD.  King 92/06/25
c                   RRINEX: Allow for COMMENT lines within data file.  King 92/06/25
c version 9.15      MONDAY: Fix bug with retained leap-year argument.  Fang/King 92/07/06
c                   RSTNFO: Fix bug when station missing and last line blank. King 92/07/06
c version 9.16      LASK : Change Apollo to make Sun/Apollo compatible.
c                   MAPAMP: Update Sun with missed 920622 Apollo changes.
c                   Makefile:  Remove IEEE from Apollo (already off Sun)
c                   Clean up 'CD==>' 'c' differences : PNS, READC1, READM1, READM2
c                          READM3, TIMINC, UT1TID, WRITC1, WRITC5, WRITM1, WRITM2
c                          Bock/King 920905
c                   Makefile:  Add THDRED, IARRAY, CLKERA
c                   RSTNFO : Fix bug with reading comments, and remove obsolete format.
c     include '../includes/dimpar.h'X.FTI  : Increase MAXTXT, MAXLIN to 100.
c                   READE : Remove unused ISAT from call.
c                   XHDRED: Add stop for NCHAN > MAXSAT.
c                              King 92/09/21
c version 9.17
c                   HISUB : Add slant height for Ashtech antenna.  Bock  92/10/01
c                           Add SERCEL antenna.   Feigl/King  92/10/01
c                   RRXHED: Rename 'C2' to 'P2' to handle serial P-code Trimble
c                           change to RRINEX made and then undone)  Bock/Calais 10/02/92
c                   RRINEX: Handle nobtyp>5. Feigl/King  92/10/05

c versinclude '../includes/dimpar.h'rease MAXSAT to 20.   Bock/King 92/10/15
c                   Explicitly declare variables in 26 subroutines:
c                      ARGMNT  ATSUB   CDATE   CONVERT   CROSS   DAYJUL   DAYNUM
c                      DBLAS1  DLINPK  DOT     GETDIR    GSATEL  IDOY     JULDAY
c                      LJUST   MENEDT  NSNPRN  ORDSIT    PROPER  RAND5    RJUST
c                      SNP     TRANSP  UPPERS  UPNAM1    UT1TID   --King 92/10/16
c version 9.19      RRINEX: Check if C1 = 0 or P1 = 0 in choosing correct L1 presudorange
c                           (important for AS data) Bock/Fang 92/10/26
c                   HISUB: Add UNAVCO standard for measuring slant heights for Trimble SST antennas
c                          Bock 92/10/27
c version 9.20      HISUB: Add DMCHRG as an alias for ROGSNR to mean the Dorne-Margolin
c                          antenna with JPL choke ring.   King 92/10/28
c                   WRITM1, WRITM2:  Remove debug print statements.   King 92/11/3
c                   HISUB: Added Rogue information in comment form
c                          Add SLBCR for Dorn Margolin slant height to bottom of choke ring
c                          Bock 92/10/29; 92/11/5.
c version 9.21      READE: Avoid error if CR at end of header.   King 92/11/10
c                   TAIUTC: Better message when file missing.    King 92/11/10
c                   Makefile:  Removed CWRITE and CEDIT for HP version in MIT stdrel
c                          (not used, already out at SIO?)   King 92/11/18

c version 9.22      CLKERA:  Make K-file site-name check only a warning.  King 92/12/17
c                   HISUB :  Change DHTGP to SLTGP for TRMSLD.   Oral/King 92/12/22
c                   READJ :  Change tdiff warning to avoid g-format compiler bug.
c                               Fang/King  92/12/23
c                   HISUB :  Add TRMSSE (same as TRMSST) Bock 92/12/29
c                   MOPENS:  Increase record length Bock 93/01/11
c                   XHDRED:  Update to 24 satellites Bock 93/01/11
c version 9.23      Consolidate SUN & Apollo versions according to RWK memo 1/19/93
c                   Changes include: upperc INCLUDE => lowerc include
c                                          INCLUDE => include
c                                    D = > CD
c                   Routines affected: LBIAS,LGOOD,LMARG,READC3,READC5,WRITC2,WRITC3
c                                      (all above changes copied from Apollo to SUN)
c                                      ERRCOD,LREAD,READC4,RFICA,RRXHED,THDRED,WRITC4,WRITM3
c                                      (all above changes made on Apollo)
c                   Remove ieee.f from SIO Apollo directory
c                   Copy from Apollo to SUN: CLKERA, READJ (change include statement)
c                   Bock 93/01/21
c version 9.24      Add SETTYP (from /makex) to /lib.  Change Makefile.  Feigl/King 93/1/15
c         9.25      WRXHED:  Correct minor glitches in header.   Beavan/King 930129
c         9.26      Consolidate Apollo and SUN versions
c                   Changed on Apollo:
c                   COPENS(still machine specific),HISUB,IERRCM,READM1,READM2,READM3
c                   RSESFO,WRXHED
c                   Changed on SIO SUN:
c                   IERRCM include '../includes/dimpar.h'it=32, maxcfl=maxsit
c                   CVIEW.FTI: ncvsit = 32 Bock 93/02/10
c                   Changed HISUB and COPENS on MIT Apollo to make identical to
c                       SIO changes (all are cosmetic).

c version 9.27
c                   XHDRED: Add antenna code to return.   King 930212.
c                   HISUB:  Add  RINEX height to return and two-option call.  King 930215
c         9.28      HISUB:  Add WM-102 antenna type.   King 930310
c                   Makefile, SUB_CHAR : Add new routine to substitute character strings,
c                        called initially by ctox.   Herring/King 930311
c                        New HP-compatible version in /stdrel.  Fang/King 930311
c                   WRITC5: Write out the atmospheric delay and azimuth rate
c                       (the latter still dummy in MODEL).   King  930312
c                   Makefile:  restore CLKERA to /stdrel version.
c          9.29     TAIUTC: Fix message on failed open.  King 930316
c          9.30     HISUB:  Add SLBGN (bottom inside of notch) for Trimble.
c                        Agnew/King  930401
c          9.31     Add subroutine DOYGWK for conversion between GPS week and day
c                    to Day of year and year (both ways) - used by FIXDRV (NGMAKE)
c                   Fang 930405
c                   RSTNFO: Convert rcvcod to uppercase.   King 930406
c                   DOYGWK: Declare variables explicitly.  King 930407
c                   Makefile:  Remove leftover testprograms from Sun Makefiles and
c                       Sun, Apollo directories.  King  930413.
c          9.32     ORBSIT: Update version from FIXDRV Bock/Fang 930419
c          9.33     SETTYP:  Write only 3 data types for Trmble SSP.  Oral/King 930502
c                      include '../includes/dimpar.h' MAXATM and change MAXPRM.  Bock
c                   HISUB: Added IGS site/antenna information as documentation.
c                       Bock 930816
c                   RSTNFO: Remove debug, Genrich fix of kinematic code
c                   Bock 930822
c          9.35     RSTNFO: Read CR correctly.  Oral/King 930910
c                   RSNTFO: Fix branch error.   King 930918/930921
c                   Add for new solve: BLANK, GETCMD, COUNT_ARG, LIFT_ARG,
c                       MATCH_NAME, STRJST. Makefile.  Dong/King 931002
c                   Add KBIT for bit-mapped variables. Makefile.  King 931013
c                   SIDMAT: Add sidereal time to output arguments, for tides.
c                        King/Morgan  931013
c                   NSNPRN, THDRED: Provisional modification for duplicate PRN #s.
c                        Should be made more robust eventually.  King   931018
c                   WRXHED, XHDRED.  Increase formats for 32 satellites.  King 931018
c                   Remove unused variables: GETMAC (Sun), IEDIT, LOWERA,
c                        LOWERC, MENEDT (Sun), PICKFN, UPPERA, UPPERC, WRITC3,
c                        WRXNAV.  King 931018
c                   RSESFO: Correct fatal bug when session.info has #sats < MAXSAT.
c                        King 931019
c                   XHDRED, THDRED: Change print line to match changes in MODEL and
c                        lessen confusion in P-file.   McClusky/King  931220/931223
c                   NSNPRN: Increase dimensions to 50 SVs.  Vigny/King  931220
c                   READJ: Remove comment on non-existent variable svdt and rename
c                        variable dt to svdt.   King 931222
c          9.36    HISUB: Correct TurboRogue choke-ring width.  Heights measured to
c                        bottom of pre-amp and input DHPAB were ok; heights measured
c                        to bottom of choke ring and input as DHPAB or DHBCR or SLBCR
c                        placed the phase center 3 cm too high.
c          9.37    WRXHED: Fix RINEX header.  Bock/Fang 940102
c                  RSESFO: Fix problem of maximum number of satellites in session.info file
c                          (same problem as in rsceno in MAKEX).  Bock 940104
c                  RSESFO: Fix problem with old format, streamline, and declare
c                        variables.  Add makex.fti to Makefile dependency.  King 940107
c                  PROPER: Update Copyright notice.  King 940107
c          9.38    HISUB, SETTYP: Add Leica receivers.  Morgan/King 940111
c          9.39    NSNPRN, THDRED : Convert NSNPRN  to subroutine, return block # and
c                        mass, generalize for duplicated PRN numbers (use Block 2 only
c                        if year > 1992)..  Change also calls in makex/wxhead and
c                        solve/uporb, and add to arc/redsat.   King 940119
c                  NSNPRN: Clean up a bit, debug.  Bock 940121
c          9.40    Remove logical function LBLK1 from library (calls in /model replaced
c                        by NSNPRN) in order to avoid mistakes when PRN #s reused.
c                        Change Makefile.   King 940122
c          9.41    RRINEX: Fix bug causing loss of (significant) data for a receiver
c                        (e.g. Leica) that has a non-zero value for the tenth-microsecond
c                        column of the time tag.  McClusky/King 940126
c                  RSESFO: Fix bug in reading sats with old format.  King 940207
c          9.42    HISUB: Fix incorrect Rogue antenna values (3 mm), add ROGAOA antenna
c                        code (Dorne Margolin B - Allen Osborne design.  Bock 940222
c          9.43    RRXHED: Initialize RINEX header variables to avoid getting leftover
c                        from previous file when missing in new file.  King 940322
c          9.44    HISUB: Change offsets for 4000-series TI antennas (L1 -6 mm, L2 -14 mm);
c                         Change offsets for Ashtech antennas (L1 +33 mm, L2 + 13 mm).
c                           McCluskey/King 940326
c                  Makefile:  Correct syntax causing a problem in the HP (Sun ok either
c                         way.  Fang/King 940405
c          9.45    SETTYP:  Minor bug in warning for P-code Ashtech with FICA input.
c          9.46    Changes to unify Sun and HP routines - Fang/King 940428
c                     CHKERR, FERROR, GETDAT:  Use SIO HP routines.
c                     GETMAC, Makefile:  Use SIO getmac.c instead of MIT getmac.f.
c          9.47    SETTYP: Deal with SSE RINEX files with C1 & C2
c                   observables (these are really P1 & P2).  Bock 940429
c          9.48    RRINEX: Handle a new antenna offset from RINEX 2 files, and pass
c                    the values out. (Called by /makex makex.f & rxscan, and by
c                    clean/rinex.f.   King 940503
c                  GETMAC:  Redimension buffer to match Sun uname parameter list.  Fang/King 940503
c                  Makefile: Change HP compiler switch from +e +E1 to +U77.   Fang/King 940503
c                  RSTNFO: Change message when entry not found.  King 940506
c          9.49    RSTNFO: Add calling argument - changes needed in /fixdrv, /makexp, /makex,
c                     /model, /utils, /hi.   King 940511.
c                  DIMPAR.FTI:  Increase MAXLAB to 22 to account for earth-rotation partials;
c                               Set MAXDM using variable dimensions.   McClusky/King 940511
c          9.50    RSTNFO: Rework logic for robustness and to read times for static
c                       observations.    King 940517
c                  RRINEX: Remove debug print.   King 940517
c                  RRINEX: Declare antenna variables real*8.   King 940518
c          9.51    Add earth-rotation rotines for common use of /arc, /model, and /orbits.
c                     Changes to Makefile; add SROTAT (modified), ROTSNP (modified), and
c                     PNROT from /orbits and SD_COMP from /model.   King 940523
c          9.52    HISUB: Add offsets for FRPA-2 antenna.  King 940601
c                  RSTNFO: Fix bug with multiple sessions.  King 940601
c                  HISUB: Add new Ashtech antenna with larger ground plane.  Bock 940602
c                  SROTAT: Declare kbit logical (HP complains).  King 940613
c                  GETDIR: Change 'ls' to '\ls' to avoid aliases (works on HP, neutral
c                       on Sun for now).   King 940620
c                  Makefile:  Move ranlib to OPTLIB dependencies to save time.  King 940620
c                  RSTNFO: Add warning for 0. slant height.   King 940621
c                  SECSMD, Makefile: Add from /orbits.   King 940707
c          9.53    UTC2GPS, Makefile: Add modified version from /fixdrv.   King 940712/940722
c                  XHDRED:  Set (existing) flag for UTC/GPS time.  King 940713/940722
c                  THDRED:  Convert T-file times to GPST if in UTC.  Changes needed for
c                           compatibility in /fixdrv, /model, and /orbits (which still
c                           has its own thdred).   King 940713
c                  ROTSNP, PNROT: Change variable name TDTUTC to TDTGPST and make units
c                           seconds, not fraction-of-day (changes made also in /arc (SBFN, TIMRED),
c                           /model (ETIDE, OTIDE, SETUP), and /orbits (TROT).   King 940715
c                  SROTAT:  Add comments that calling argument is UTC not GPST.
c                  PROPER:  Update copyright notice.   King 940715
c                  SROTAT, ROTSNP:  Add xpole, ypole to SROTAT calling argument for use in
c                           Earth-rotation partials in /model (only other call).  McClusky/King 940718
c                  THDRED:  Modify screen print.  King 940719/940721
c                  Makefile, RUNTIM: New routine to return the date and time in a character
c                           string.   King 940719
c                  READJ:  Convert to using GPST internally.   King 940720
c                  LVERSN: Blank version to avoid non-ascii characters.   King 940728
c                  CLKERA: Assume input is GPST (change also in /utils/plotk).  King 940728
c          9.54    HISUB:  Add Ashtech direct ht to top of ground plane.  McClusky/King 940809
c                  HISUB,MAPAMP,SETTYP: Add TOPCON GP-R1DP & GP-R1DY systems.  Bock 940810
c                  CLKERA: Increase maximum jumps from 100 to 240, allowing jumps every
c                           6 minutes for 24 hours.   King 940819
c                  XHDRED: Fix format for header print.  King 940914
c          9.55    WDMS: Added option +3 and -3 so that given a radians input a simple conversion
c                        to degree minutes seconds could be made. McClusky 940915
c                  HISUB: Added graphic for Ashtech large ground plane antenna.  Bock 941008
c                  POLRED, UT1RED:  Fix infinite loop under certain backspace conditions.
c                         (Change made 940930)  King 941105
c                  LEAPYR: Allow 2-digit as well as 4-digit input.  King 941114
c          9.56    THDRED: Fix print for expanded rad parameters.  McClusky 941222
c                  Changes to satisfy XL compiler on RISC machines:
c                     CHKERR, PROPER, RAND5  King 950102
c                  Remove EDTSUB IEDIT from directory and Makefile.  Used only for interactive menus,
c                     now obsolete.  King 950102
c                  GETUSR: Add call to getenv to handle background jobs.  King 950106
c                  THDRED: Check for number parameters > maxorb.  King 950106
c                  DIMPAR.FTI: Set MAXDAT = 7 (not 5) to handle extra observables in RINEX.
c                        Fang/Bock
c                  RRXHED: Fixed bugs in reading PGM line of rinex header, and number of
c                          comment lines. McClusky 950201
c                  RRXHED: Allowed blank lines in header for RINEX version 2.  King 950202
c                  RRINEX: Catch nobtyp = 0.   King 950202
c          9.57    TIMCON: Replace call to SUMDAY.  King 950210
c                  SECSMD, SUMDAY, Makefile: Remove (no longer used, and don't account for year
c                        boundary). King 950210
c          9.58    THDRED:  In order to accomodate NOAA T-files, don't allow nics=8 (hope there
c                        are no other T-files around with this value.  King 950330
c                  JUL2HMS, Makefile: New routine (combines several other time conversion routines). Tregoning 950331
c                  AMAG, Makefile:  New function to do vector magnitude.  Tregoning 950331
c                  INDEXX, Makefile: New routine from orbits that creates a sorted index of an array. McClusky 950331
c                     include '../includes/dimpar.h'ax number eclipses for yaw-rate and velocity-impulse parameters.
c                         Tregoning/King  950331
c                  THDRED:  Fix bug in converting from NS to PRN.   King 950401
c                  SVNAV_READ, Makefile:  New routine for reading svnav.dat with yaw rates; written for
c                        ARC by S. McCluskey 950330; moved to /lib by King 950403.  Remove old NSNPRN.
c                  THDRED:  Replace call to NSNPRN by call to SVNAV_READ.  King 950403
c                  SVNAV_READ: Stop if SV not found.  King 950405
c                  SVNAV_READ, THDRED: Add hrs,min to call.  McClusky/King 950405
c                  Makefile:  Update .fti dependencies for ORDSIT, GSATEL, AMAG3, RSESFO.  King 950406
c                  SVNAV_READ:  Stop if invalid date.  Tregoning/King 950407
c                  SVNAV_READ:  Stop if maxnsn exceeded.  King 950410
c                  THDRED: Correct times in SVNAV_READ all.  King 950410
c          9.59    IERRCM: Reordered returned flag when comparing two flags. LOEL now
c                          lower priority than unwt. McClusky 950412
c                  LLOEL, Makefile: New function to evaluate if variable passed is a low elevation flag. McClusky 950412
c                  ROTSNP: new precession argument - IAU76: true,  AENA 1963: False PT 950414
c                  PNROT: new precession argument  - IAU76: true,  AENA 1963: False PT 950414
c                  PRCES: replace data statement with new precession argument PT 950414
c                         fix if statement for IAU76 coding to compute prec. angles on first call PT 950413
c                  THDRED: determine precession model used from tfile header and pass it out PT 950415
c                  ROTSNP, SROTAT, SIDMAT, SIDTIM  : Fix sideral time to be consistent with precession
c                         constant (IAU 1968 or 1976).  King 950420
c                  INDEXX: Handle n=1 case properly.  King 950421
c                  READM3: Fixed an incorrect suicide statement from writm3 to readm3. McClusky 950421
c                  PRCES, SIDTM, SIDMAT:  Final values, correct logic, and clearer comments for IAU68
c                         and IAU76 expressions.  Tregoning and King 950424
c                  ROTSNP,PNROT,PRCES,SROTAT,SIDMAT,SIDTIM: Compute snp matrices for whatever
c                         inertial frame and/or precession model is requested  PT950409
c                  THDRED: Read inertial frame from tfile header PT950509
c                  RSTNFO: Added code to print some useful information when an error reading
c                          the station.info file is encountered. McClusky 950516
c                  READM1: Zero out m-file arrays before reading.  King 950516
c                  THDRED: Add argument for radiation pressure model.  King 950517/950518
c                  Makefile: Replace KBIT with complete BIT_UTIL from HP1000 library; Add new
c                           routine IER2DF to convert old GAMIT error flags to AUTCLN flags.  King 950518
c                  GETCMD:  Increase line size to 256 characters to accomodate 15 orbit parameters.  King 950518
c          9.60    New C- and M-file formats:  READC1, READC2, READC4, READC5, READM1, READM2
c                         WRITC1, WRITC2, WRITC4, WRITC5, WRITM1, WRITM2.  King 950520/950607
c                  HISUB: Add new Ashtech Dorne Margolin with choke ring antenna.  Bock 950527
c                  THDRED: Make sure radiation-pressure model keyword is uppercase.  King 950608
c                    include '../includes/dimpar.h'L:  Add T-file dimensions to dimpar.fti.  King 950612
c                  READC2, WRITC2:  Correct spelling of 'srpmod' in call list.  King 950613
c          9.62    READM1, READM2, WRITM1, WRITM2:  Add decimation factor to record 2.  King 950616
c                  PRCES: Allow J2000 for reference epoch.  King 950620
c                  PNROT: Fix comments.  King 950620
c                  HISUB: Add comment on Ashtech Dorne-Margolin antenna.  Bock 950620
c                  PRCES: added comments/explanation   Tregoning  950621
c                  PRINTMAT,Makefile: new routine to print matrix (for debugging)  Tregoning 950621
c                  THDRED: assign J2000 as default inertial frame for efixed Tregoning 950623
c          9.63    READC2, WRITC2:  additional mod to format to add antenna reference point
c                      offset from mark, for SINEX file.  King 950717
c                  PRINTMAT: Remove unused variable.  King 950717
c          9.64    Makefile: Add lines for DEC.  Sanli/King 950719
c                  ERRCOD: Remove unused (and undefined) 'ighdwr'.  Sanli/King 950719
c                  GETDAT: Add (commented) lines for DEC.  Sanli/King 950719
c                  GETDIR: Add status for open (needed by DEC), and correct bug in assignment
c                      of ftell(6) call.  Sanli/King 950719
c                  RAND5: Make function GAUSS double precision to match calls; replace 'jtime'
c                      with 'i3(6)' in RAND_$INIT.  Sanli/King 950719
c                  STRJST: Define function in length=0 case.  Sanli/King 950719
c                  RSTNFO: Remove unreachable statement.  Sanli/King 950719
c                  READC2, WRITC2: Change names of variables itide/iut1pol=> ietide/isptide
c                      for conformity with SOLVE/h-files.  King 950720
c                  Makefile:  Use different routines, bit_utils_sun.f and bit_utils_hp.f
c                      depending on machine.   King 950720
c          9.65    READC4, WRITC4: Remove duplicate declaration.  King 950721
c                  ORDSIT:  Change test on dimensions from 'maxnet' to 'ndim'. King 950725
c                  READJ:  Add channel and prn # to warning message.  McClusky/King 950726
c                  RRINEX: Add warning if # SVs on record exceeds dimensions.  King 950726
c                  BIT_UTIL_SUN: Remove redundant dimension statement.  King 950727
c                  IER2DF: Fix mistyped variable names (bug).  King 950727
c                  WRITM2: Fix word-count.  King 970727
c                  THDRED: fixed logic for defining frame/precession for efixed  Tregoning 950724
c                  Makefile, ROT_GFILE, GHDRED: moved from ARC to library, and modified
c                       for GTOG  Tregoning 950725
c                  RXREAD: moved to library from ARC  Tregoning 950725
c                  ROT_GFILE,GHDRED: removed common statements  Tregoning 950727
c                  LIFT_ARG: Correct comments.   King 950729
c                  ROT_GFILE:  Remove reference to arcout file.  King 950801
c                  COPENS:  Remove conversion of status to lowercase (problem on DEC,
c                       possibly with altering of character string.  Sanli/King 950804
c                  Makefile, SD_COMP: Remove block data portion of file to new routine
c                       SD_COMPBD so that it can be listed first in Makefile (DEC
c                       requirement).   Sanli/King 950804
c          9.66    GHDRED: Correct format statement and add radiation pressure model names to
c                       warning message (no effect on integration).  McClusky/King 950807
c                  PRCES: Add 'd0' to 0.5 in J2000 code.   King 950807
c                  ROT_GFILE: removed unused variable from argument list  Tregoning 950807
c                  SUICID: Remove question about core dump.   King 950814
c                  CLKERA: Fix omission of UTC-to-GPST conversion in computing residuals;
c                        require values spanning longer time for segment fit.  King 950817
c          9.67    POLRED: Remove extraneous slash in assignment of 'nrbmin'.  Morgan/King 950825
c                  READE: Fix branching into block IF.  Morgan/King 950830
c                  WRITC5: Assign previously undefined 'wvrdel4'.  Morgan/King 950830
c          9.68    RRXHED: Avoid infinite loop with empty RINEX file.  Herring/King 950912
c          9.69    Add RDSEST (from /fixdrv) to library; change Makefile.  King 951003
c                  WRXHED: Add 'END OF HEADER' to rinex header (for cview)  Tregoning 951002
c                  RDSEST: Avoid case sensitivity for sestbl. keywords (calling arguments
c                      and table entry.   King 951011
c                  READE:  Accept only PRNs <32, vice <37.  King 951011
c          9.70    RRXHED, WRXHED: Pass 20-char names, rather than integers for receiver and
c                      antenna serial numbers (change also makex/rhead,rxscan and
c                      utils/merge_rinex,argo2rx,xtorx.   King 951019
c                  XHDRED: Extract and pass out the 20-character receiver and antenna
c                      names and serial numbers (change also clean/readx; makex/makek;
c                      fixdrv/imakef,seschk,xcheck,xcsasts; orbits/tinput;
c                      utils/countx,plotk,xtorx).   King 951019
c                  Add ROTATE_GEOD, DVDOT, DVSWP, XYZ_TO_GEOD from /orbits and
c                      /model.  Change Makefile.   King 951019
c                  REC_ANT: New routine to read receiver/antenna translation table rec_ant.dat. McClusky 951019
c                  READ_ANTMOD: New routine to read information from antmod.dat table. McClusky 951019
c                  HISUB: Major changes this routine now uses ARP --> L1/L2 phase centre offsets
c                      read out of the antmod.dat table. McClusky 951019
c                  BILIN: Moved bilinear interpolation routine from model to lib. McClusky 951019
c                  LINEAR: Moved linear interpolationroutine from model to lib. McClusky 951019
c                  DIMPAR.FTI: Added maxaz and maxel parameters to dimpar.fti. McClusky 951019
c                  REC_ANT: Modified to read rec_ant.dat in free format. McClusky 951020
c                  READ_ANTMOD: Fixed a bug reading comment lines in the antmod tables. McClusky 951020
c                  READ_RCVANT (rename from REC_ANT), Makefile: Change name of table with rcvr/antenna
c                      full names to rcvant.dat.  Tregoning/McClusky/King 951024/25
c                  XHDRED:  Correct format in reading antenna serial #.  King 951024
c          9.71    RSTNFO:  Change warning to stop for slant height = 0.  King 951026
c                  HISUB:  Add to message on unfound measurement type; add comments about
c                      inconsistencies in old code for ASHXII SLBGP and TRMSLD SLTGP.  King 951027
c                      Fix declaration of 'switch'.  King 951102
c                  GETMAC:  Increase buffer size to suit Solaris 2.  Fang/King 951102
c          9.72    HISUB: Make 6 mm corrections in radius to notch value for TRMSST and fix
c                      documentation for Ashtech antennas.   Agnew/King 951108
c                  READC1: Write more explicit message if C-file empty.  King 951109
c                  HISUB:  Add and subtract aliases for DM and Ashtech antennas.  King 951109
c                  READ_RCVANT: Fix character comparison for rcvr match.  Tregoning 951109
c                  XHDRED: Initialize with blanks the rcvr/ant codes.  King 951109
c          9.73    READ_RCVANT:  Make message more explicit if file missing.  King 951111
c                  HISUB:  Disallow 'TOPCON' in favor or 'TOPP12' or 'TOPGD3' for antennas
c                     accompanying the GP-R1DP (Ashtech P12) and GP-R1DY (Ashtech Z-12)
c                     receivers, respectively.  King 951111
c          9.74    RSTNFO: Catch slope height near zero.  McClusky/King 951113
c                  RDSEST: Fix bug in (9.69) change to make entries case insensitive.  King 951114
c                  HISUB: Correct by 6 mm ARP-to-TGP distance for Ashtech L and P from D. Agnew
c                     measurements; correct omitted pre-amp height (35 mm) for Ashtech SLBCR
c                     measurement.  McClusky/King 951114
c          9.75    Makefile: Variable ranlib to make compatible with Solaris 2.  Fang/King 951207
c                  READ_ANTMOD: Removed free format read of antmod.dat (DEC problems)  Tregoning 951215
c                  READ_RCVANT: Removed internal free format read (DEC problems)   Tregoning  951215
c          9.76    Makefile, BAD_OPTION, CASEFOLD, READ_LINE: Globk routines added to library to handle
c                     free format character reads better.  Tregoning 960102
c                  HISUB:  Fix documentation (comments only) for Ashtech GD3.  Agnew/King 960118
c                  READ_ANTMOD: Fix one more DEC-compatible problem. Tregoning  960123
c                  READ_ANTMOD: Removed free format read for L2 as well!  Tregoning  960124
c                  READ_ANTMOD: Replace fixed format reads with read_line calls.   Tregoning 960129
c                  CLKERA: Allow a second pass through the k-file if the clock synch has failed
c                      first time through--may be a bad first epoch.   King/McClusky 960215
c                  READ_LINE: Change input array dimension from (1) to (*) to work with
c                      bounds-checking on.   McClusky/King 960216
c                  READE: Modified roution to filter out bogus ephemeris records, and to
c                         properly read rinex version2 navigation files. McClusky/King 960216
c          9.77    Makefile, CASEUNFOLD, REPORT_STAT, TRIMLEN:  Globk routines added to library to
c                         handle status, warning, error reporting  (temporarily change report_stat
c                         to call GAMIT gettim and getdat rather than Globk 'hpsystime').  King 960217
c                  REPORT_STAT, CHECK_ASCII, Makefile : Update with new version from libhp1000.   King 960219
c                  Add report_stat calls to 30 routines (through 'read_rcvant', alphabetically).  King 960221
c                  Makefile:  remove unused MENEDT.  Also remove from directory (already out of
c                      Makefile) unused DUMMY and LOWER.  King 960221
c                  HISUB: Extra comma in format statement.  Fang/King 960226
c                  Makefile, REPORT_STAT, SYSTIME:  Copy Sun version of SYSTIME from kf
c                      lib, removing HPSYSTIME, and copy from kf updated REPORT_STAT to handle
c                      removal of '/'s in module name so that lib routines can pass the
c                      calling module instead of 'LIB'.   King 960226
c          9.78    READE: changed call to get_prog_name --> rcpar inline with new report_stat.
c                  REPORT_STAT: fixed compile bug, missing definition of integer variable i. McClusky 960227
c                  READ_ANTMOD: Detect old-style antmod.dat tables and abort. McClusky/King 960229
c                  Changed 50 routines to incorporate calls to report_error.  King/McClusky 960301
c                  GHDRED: Fix format statement used to write message for report_error.  Fang/King 960302
c                  CLKERA: Remove duplicate declaration.    King 960303
c                  ORDSIT, READC1, READC2, READC3, READC4 : Move 'implicit none' from below 'dimpar.fti'.
c                  READC5, WRITC1, WRITC2, WRITC3, WRITC4, WRITC5: Add implicit none.   King 960303
c                  LREAD:  Fix print of site in report_stat call.  King 960318
c                  READM3: Check for array sizes before reading m-file header.  King 960322
c                  REPORT_STAT: Changed file names from [MODULE].* to [GAMIT].*.  King 960327
c                  READ_RCVANT:  Fix 'message' size in call to report_stat.  King 960327
c                  HISUB: Fix write format for 'message' variable use by report_stat. McClusky 960329
c                  RRINEX:  Remove extra comma in calls to report_stat.  King 960401
c                  REPORT_STAT: Change to give FIXDRV its own status files to avoid erasure in batch
c                      run.  King 960401
c                  ROT_GFILE: Fix write format for 'message' variable used by report_stat.  King 960401
c                  RRINEX: Fixed embedded subroutine RXERR, which contained incorrect use of the routine
c                     RCPAR. RCPAR was called as a subroutine instead of an integer function. McClusky 960415
c          9.79    Modified multiple routines to expand the "message" character variables from *80 --> *256:
c                  CLKERA,DBLAS2,MAPAMP,NUTRED,RDSEST,READC1,READM1,READM2,RFICA,ROTSNP,SVNAV_READ,TAIUTC,THDRED,
c                  TIMCON,UT1RED,XHDRED. McClusky 960429
c                  POLRED, UT1RED, NUTRED: Modified to correctly handle error messages for out of table
c                     bounds requests. McClusky 960429
c                  XHDRED: Change null to blank assignment of 'rctype' (caught by DEC compiler).  Bock/King 960508.
c          9.80    ROTSNP: Fix erroneous calculation of sidereal time at leap-sec boundary.  Fang/King 960517
c                  GHDRED: Fix conversion from GPST to UTC at leap-sec boundary.  King 960517
c                  HISUB: Correct SLPAC (corner of base) and add SLPAS (side of base) for Trimble
c                         SXD and SLD antennas; improve comments for DM antennas   King 960525
c                  ORDSIT: Move location of integer declarations above character (DEC compile error). Tregoning 960527
c                  RSTNFO, SVNAV_READ, TIMCON: Replace '' strings with ' ' in report_stat calls  Tregoning 960527
c                  WRITC5:  Fix order of variables 'ampl1, ampl2, obswgt' in writing.  Tregonning/King 960528
c                  REPORT_STAT:  Fix bug in 'CLEAR' call.  Tregonning/King 960528
c                  RSTNFO: Add site code to fatal message if no match; check for chronological order
c                     of entries in station.info.   King 960531
c          9.81    READJ: Fix error in write format of variable "message" when satellite was missing from
c                     the J-file. McClusky 960604
c          9.82    REPORT_STAT: Use module name for MAKEXP (as well as FIXDRV) files.  King 960610
c                  CLKERA, DBLAS2, FERROR, GETDIR, GSATEL, LASK, PICKFN, RDSEST, READC1:  Fix report_stat
c                      calls to pass prog_name instead of 'LIB'.   King 960610
c                  RSTNFO:  Get prog_name correctly for header entry point.  King 960610
c                  RSTNFO:  Bypass session time check if called by MAKEXP.  King 960614
c                  REPORT_STAT: Fix bug from mistyped 'ierr/jerr', which created a segmentation fault
c                       on a 'CLEAR' call with Solaris 2.  Herring/King 960617
c                  LREAD:  Add catch for illegal geodetic-coordinate entry.  King 960619
c          9.83    SUB_CHAR: Remove redundant TRIMLEN from this routine (DEC complains). Bock 960609
c                  READE: Change "lunit" to "lu" .  Bock 960609
c                  GHDRED: If rad-model name missing from g-file, assume spherical if
c                       nics=9 (as before) and Berne if nics=15.   Bock/King 960623/960710
c                  REPORT_STAT: Incorporate /kf version changes to avoid overwriting constants
c                       and set report filename logic to work for either GLOBK or GAMIT, either
c                       batch or interactive. King 960710 
c          9.84    Changes for new kf/gamit structure and Makefiles from unimake.  King 960718
c                        All routines:  removed trailing blanks and replaced lib/*.fti 
c                           includes by ../includes/*.h.  
c                        Removed CHECK_ASCII, CASEFOLD, CASEUNFOLD, GETMAC, SYSTIME, SUB_CHAR,
c                           and TRIMLEN, now included from libraries/libcom. 
c                  GHDRED: Fix bug in call to report_stat.  King 960724 
c                  READ_ANTMOD: Added new arrays to the list saved for later read_antmod calls. McClusky 960808 
c                  Declare variable explicitly:  ATSUB, CDATE, COPENS, DAYNUM, GENIV1, 
c                        GETDIR, LASK.  King 960814
c                  RAND5:  Replace intrinsic 'and' by machine-dependent 'cand' (in libraries/comlib). King 960815
c          9.85    WRITC4:  Remove 'stop' with kflag.ne.0.  Chen/King 960822
c          9.86    Makefile, LOWERA, UPPERA:  Remove these two routines, which are dangerous 
c                       and no longer called in GAMIT.  King 960829
c          9.87    RSESFO: Trap # SVs exceeding dimensions or session.info format.  King 961016
c                  Makefile: Removed RAND5 routines (save for later use).  King 961017
c                  READ_ANTMOD: Make all failures fatal.   Fang/King 961030
c                  READE: Set program for report_stat files to 'ORBITS' for BCTOT, BCCHECK.
c                     King 961105
c                  READE: Add time-of-transmission and header info to callling arguments
c                     to support writing out of new e-file.  King 961107
c                  RSESFO: After trap, set the format dimensions to 32, not the current nsat.  King 961107    
c                  CLKERA: Raise tolerance on synch to allow 0.5 ms in 120s.  King 961118
c                  CLKERA: Send non-synch message to report_stat as well as screen to allow
c                      use of messages to filter later solutions.  King 961125
c                  CLKERA: Modified to allow estimation clock coefficients for receiver clocks
c                          that have high allen variance's but drift rates less 0.4 ms/epoch. McClusky 961127 
c          9.88    HISUB: Add alias DMGASH for ASHDMG.  King 970110
c                  Makefile, ANT_ALIAS (new), READ_ANTMOD:  New routine to get antenna code aliases;
c                          convert READ_ANTMOD to lower case and add comments.  King 970111
c                  Remove RCPAR, READ_LINE (duplicate routines in libraries/comlib). Morgan/King 970111
c                  Changes to initialize and declare all variables for G77 compiler:
c                      BILIN, CLKERA, DBLAS1, DBLAS2, DLINPK, GEOXYZ, GETDIR, GHDRED, 
c                      LASK, LINEAR, LJUST, NUTRED, POLRED, PRCES, RJUST, ROT_GFILE,
c                      RSTNFO, SIDMAT, SIDTIM, SROTAT, STRJST, UT1RED, WFICA. McClusky/King  970111
c                  RSTNFO: Change 'rcvers' => 'swvers' to avoid confusion with RINEX variable.  970111
c                  XHDRED: Replace print by call to report_stat.  King 970113 
c                  READ_ANTMOD:  Remove redundant declaration; correct type in call to read_line.  
c                      King 970117/970122
c                  READM1, WRITM1:  Update version to cfmrg version 9.31 (gradient parameters).  King 970123
c          9.89    ANT_ALIAS, HISUB, READ_RCVANT:  Add alias call in hisub and read_rcvant; 
c                      fix alias names for DM antennas. King 970127  
c                  RSESFO:  Add trap for bogus format lines.  King 970129
c          9.90    LINEAR, BILIN: Fixed bug that allowed -ve elevation angles to corrupt 
c                        antmod.dat table interpolation. McClusky 960129
c                  READM1: Fix report_stat message for obsolete m-file.  Herring/King 970206
c                  RRXHED: Initialize start,stop times to 0; fix comment.  King 970210
c                  RRINEX: Fix comment.  King 960210
c                  Makefile, INDEXI: New routine to sort an integer array (analagous to indexx
c                      for real*8.   King 970211
c                  RSTNFO: Correct report_stat message.  King 970212 
c                  GETDIR:  Remove message when files not found (done in calling programs). King 970212
c                  SETTYP: Stop if ndat = 0.  King 970212
c                  RRINEX: Fix report_stat calls for time glitches.  King 970213  
c                  JULDAY: Fix report_stat calls.  King 970219
c                  WTIME:  Change label format from 'a20' to 'a' (needed for rxscan). McClusky/King 970221
c                  RRINEX: Fix format problem.  King 970303 
c                  FERROR, RRINEX: Fix length of filename.  King 970303
c          9.91    CONVRT: Removed unused variable DEGS. McClusky 970305  
c                  SD_COMP, SD_COMPBD: Modified so G77 would correctly find BLOCK DATA. McClusky 970307
c                  GETDIR: Mods to make LINUX version work. Commented out ftell and fseek calls. McClusky 970311
c                  DS2HMS: Mods to make LINUX version work. Problem with time format. McClusky 970311
c                  Makefile, SORT_STRING:  New routine for alphabetizing and removing duplicate
c                      entries; replaces ORDSIT (removed).   King 970307 
c                  WRXHED, WRXNAV: Make work with new makex/fic2rx.  King 970311/0313 
c                  SETTYP: For Trimble SST, force l2fact=2 if ndat=3 (bad RINEX translation) 
c                      or if swver=4.0 (new, artificial, for partial serial P-code tracking).
c                      King 970314   
c          9.92    Makefile, READD: Add version from /ctox and use also for /makex/xtorx. 
c                      King 970321  
c                  TOPENS: Stop if binary tfile is incompatible with computer 
c                      architecture (ie tfile is incorrect binary format).  Tregoning/McQueen 960204  
c                  GETDIR: Initialise istat = 0   Tregoning 970320/King 970325  
c                  THDRED: Replace variable efixed with new logic for variable frame. Tregoning970321
c                           (note: there is no longer an inertial frame defined for
c                            the EFIXD case - is this a problem?)  
c                  XHDRED: Initialised ischan array before reading header. McClusky 970314
c                  GETDIR: Reinstated the ftell and fseek calls for unit 6 . McClusky 970325
c                  GETTIM: Fix dimension of i3 (6 -> 3).  McClusky 970325  
c                  THDRED: Correct format bug in assigning frame. King 970326    
c                  TOPENS: Check if tfile exists before reading from it.  Tregoning/King 970336
c          9.93    CLKERA: Generalize code to allow easy shift among orders of polynomial, for
c                      both the no-jumps and jumps cases.  Leave default at linear (no jumps) and
c                      cubic (with jumps).  Making the jumps poly linear makes the fit more
c                      robust for initial gaps, and thus less likely to cause SOLVE numerical problems,
c                      but this not necessary with autcln use_phase_clocks option. King 970430/970508
c                  XHDRED: Fix format in maxtxt-exceeded message.  King 970501
c                  SETTYP: Fix bug in 970314 change: set lambda for L2=-2 (not +2).  King 970505
c                  HISUB:  Add Trimble choke ring.  King 970505 
c                  EVEN_MINUTE, Makefile:  New routine to round start or stop times to the even
c                      minute for comparison with session.info and station.info entries.  King 970506 
c                  THDRED: fix small bug in defining frame for case 'UNKWN'; fix a comment and
c                      a few format literals.  Tregoning/King 970421/970508
c                  READ_ANTMOD, HISUB: Restructure calling arguments for better tracking of requested
c                      and found phase-center models.   McClusky/King 970509
c                  READE: Fix message for bogus etoe.  King 970514    
c                  RRXHED: Initialize sampling interval to zero to avoid picking up last value.  King 970519
c                  CLKERA: Modified so that machine dependant round off is eliminated. McClusky 970521
c                  CLKERA: Make jump detection a <<change>> in residuals; add both fits to
c                      extended printout; add polynomial orders for jump and no-jump fits to
c                      calling arguments. King 970522
c                  HISUB:  Add ASHG3R for USCG radomed version of Geodetic III.  King 970529 
c                  THDRED: Don't issue a 'no inertial frame' warning for E-fixed files.  King 970529  
c          9.94    HISUB: Fix bug in use of ASHGD3 and ASHG3R.  King 970603      
c                  RRINEX, SVNAV_READ:  Avoid Hollerith line-splitting for picky DEC compiler.  King 970606 
c                  TOPENS:  Check for binary compatibility only if status=old in call.  McClusky/King 970606    
c                  TOPENS:  Minor changes in message. King 970612
c                  COPENS:  Check for binary compatibility.  King 970612
c          9.95    WTIME: Minor format change to output format. McClusky 970626 
c                  HISUB:  Add Leica SR299I, SR299E, SR299W but only for DHPAB (ARP).  King 970702   
c                  UT1RED: Remove slash in call to report_stat.  King 970707
c                  HISUB:  Add Leica height hook measurement.  King 970717
c                  HISUB: Fix report_stat message format error.  King 970718   
c                  RRXHED: Fix to ignore different wavelength factor for SV.  McClusky/King 970722  
c          9.96    HISUB: Add SCIGN antenna options for Ashtech.  Bock/Genrich 970728; King 970804
c                  SVNAV_READ: Compress Hollerith format to satisfy DEC.  King 970805  
c                  GETDAT: Declare IDATE external to satisfy DEC.  King 970805
c                  Makefile:  Remove lblk1.f (obsolete from 9.40).  King 970806
c                  SETTYP: Determine ndat solely from idat for RINEX-in, X-out.  King 970807  
c                  RRXHED: Read properly with warnings different wavelength factors for different 
c                      SVs (but still set all equal on output).  King 970811 
c                  XHDRED: Trap the rrxhed bug.   King 970811
c                  RXERR: Fix subroutine name in report_stat calls.  King 970811
c                  RRINEX: Fix formats in report_stat messages.   King 970817
c                  WRXHED: Write 4-digit yr for start, stop times.  Agnew/King 970829 
c          9.97    SVNAV_READ: Fix bug in format of read of svnav.dat (table also changed to
c                      avoid the problem even without the code change).  Tregoning/King 970901
c                  PICKFN:  Allow a response of '0' to return 'none', giving the option of
c                      having the calling program take alternative action.  King 970902 
c                  READC5:  Add 'data_flag' to debug print. King 970902 
c                  SETTYP:  Add documentation of Trimble firmware versions.  King 970902 
c                  WRXNAV: Correct assignment of clock epoch and write 8th line for
c                        RINEX 2.  King 970903
c                  READE: Convert warning to report_stat call.   King 970903
c                  GETDIR: Redirect 'No match' to tmp.getdir and then remove.  Herring/King 970905
c                  RSTNFO:  Skip blank lines to avoid getting the wrong entry.  Herring/King 970905
c                  RFICA: Initialize arrays before reading to get rid of junk in the case of short
c                        records (particular nc=14 but only 12 entries for blk101).  King 970917
c          9.98    READ_RCVANT:  Intialize 'answer' to avoid erroneously missing antenna type.  King 971008 
c          9.99    DAYJUL: Added 1.d-10 kluge to jd so that linux real -> int bug is avoided. McClusky 971027 
c                  SETTYP: Add warning for SST RINEX translated at SSE.  McClusky/King 971027
c                  DAYJUL: New routine. Now only integer arithmatic is used in calculations. McClusky 971028 
c                  CLKERA: Correct bug leading to double counting of jumps.  King 971205
c                  HISUB:  Increase TGP - ARP by 1 mm; GP radius by 10 mm.  King 971208
c                  READE:  Change FATAL to WARNING for bad RINEX nav records.  King 971209
c                  GETDIR: Add to comment about 'ls' info about Gnu compilers.  King 971216
c                  RSTNFO: Minor change in message format.  King 971216
c                  RRINEX: Replace references to ferror with rxerr.  King 971217 
c                  HISUB, ANT_ALIAS:  Add Leica choke ring, consolidate other Leicas.  King 971229
c                  SETTYP:  Add comment re Trimble serial P-code firmware.  King 971229 
c                  SETTYP: Fix bug in 971027 change.  King 971230  
c                  HISUB: Add distinction between outside and inside of holes for ASHGD3.  King 971231
c         10.01    Added from /model for yawtab:  EPHDRD, NORMALISE, SHADOW1, SOLRED, TIMDIF, YAW_ATTIT.
c                    YAW_ATTIT modified slightly.  Modified Makefile. Tregoning 971201 / King 980106 
c                  YAW_ATTIT:  Change S/C coord vector dimensions to match those in orbits/calc_attit.f.  King 980113  
c                  RRXHED: Minor mods to fix some bound checking problems. McClusky 980115 
c                  SHADOW1:  Use 6-vector for SV coordinates; remove vlight and twopi from calling 
c                           sequence.  King   980121    
c                  Makefile, JEL, JELF: Move /solve indexing routines to /lib in order to facilitate
c                       their use in clean/polyft.f in conjunction with HP compiler bug.  Herring/King 980121   
c         10.02    WTIME: Added new calling argument 'GPS' or 'UTC' to decide the printed output time. McClusky 980123 
c         10.03    GETDIR: Fix to force 'ls' under 'csh' to avoid problem executing under a chron 
c                        (which uses 'sh').  Herring 980204 
c                  GETCMD:  Skip comment lines and change search for 'end' (used only for FONDA, actually)
c                          to avoid quitting on 'end' in first 10 columns.  King 980217  
c                  GETDIR:  Restore the '\' before 'ls' to avoid problems with aliases.   (But there 
c                      may still be problems on some compilers that interpret `\` as an escape.  Herring/King 980224  
c                  RRXHED:  Fix call to report_stat.  King 980228   
c                  READD: Fix calling arugments for sort_string for multisession.  King 980303
c         10.04    INVER2:  Add from merging of /solve and /orbits version.   King 980313
c         10.05    YAW_ATTIT: Turned of the LHS calculation for Block 4 satellites. McClusky 980518
c         10.06    WTIME: No longer call lowers to convert the ostype variable. McClusky 980518 
c                  YAW_ATTIT: Switch Block 4 from left- to right-handed system (matches changes in
c                      /orbits/block_iir.  McClusky 980521
c         10.07    SVNAV_READ: SV ant offsets read from svnav.dat.  King 980603
c                  READE:  Add more error detection for short nav files.  King 980605
c                  HISUB:  Add Ashtech marine antennas.  Matt King / R. King  980612
c         10.08    GHDRED: Add CODE9801 (BERN2) and SRDYB (previously omitted by mistake) radiation 
c                     pressure models.  King 980701
c                  New routine CROSS_UNIT, Makefile: new routine to compute the unit vector
c                     from a cross product.  King   980702
c                  HISUB:  Correct typo in naming Ashtech marine antenna (MAR->MRA). M. King/R. King 980706
c                  New routine ANGVEC, Makefile:  new routine to compute the (signed) angle
c                     betweeen two vectors.  King 980707  
c                  GHDRED: Add BERN1 radiation pressure model.  King 980714
c                  XHDRED: Add checks for read errors.  King 980720
c                  Makefile:  Remove ljust and lunit (moved to libraries/comlib).  King 980721
c         10.09    GETCMD: Add icmd=0 option to indicate input command found but no value
c                     returned (ok in some cases).  King 980804
c                  XHDRED: Pick out correctly the rcvr serial # for FICA headers.  King 980806
c                  WRXHED: Removed extraneous comma (SGI complained).  Morgan/King 980901
c                  Makefile, FERROR:  Move to libraries/comlib since called by getdir_sun.   Morgan/King 980901
c                  XYZ_TO_GEOD: Change earth radius and flattening to WGS84 values to match
c                     XYZ_to_GEOD in kf/includes/const_parms.h    Matt King, R King 980901
c                  ANT_ALIAS: Update and match to RINEX standard the Leica antenna names.  King 980902
c         10.10    READC1, WRITC1, WRITC2, READC2:  Change version to 9.80 to reflect changed preval assignments 
c                    and addition of 'norb' parameter.  King 980907/980911
c                  XHDRED: Blank out ends of strings read from X-file to avoid propagating
c                    nulls into RINEX via XTORX.  King 980911
c                  READM1, READM3, WRITM1, WRITM3:  Add multiple gradients to M-file.  King 980914
c                  CLKERA: Fix format statement for message about < 3 pts sychned.  McClusky/King 980916
c                  XHDRED: Fix bug in reading receiver type.  King 980921
c                  SETTYP: Make sure all four columns are written on X-file.  King 980921
c                  RRXHED: Allow END OF HEADER to end header with RINEX 1 (strictly RINEX 2).  McClusky/King 980923
c                  READC5, WRITC5: Add SV clock parameters to C-file (still version 9.80).  King 981002
c                  READC5: Replace 'suicid' by 'report_stat'.  King 981003
c                  RSESFO: Allow 4-digit year; make implicit none.  King 981204
c                  Makefile, RAND5:  Add random-number generator for MODEL simulation.  King 981208
c                  TIMINC:  Replace with new version, more robust on Linux.  Fang/King 981230
c                  Remove implicit statement to statisfy SGI compiler, and add implicit none;
c                    all routines compiled ok previously with explicit compiler switch:
c                    ARGMNT, CONVRT, CROSS, DAYNUM, DOT, FUNARG, FUNCOF, GSATEL, NUTTAB,
c                    PNROT, PNS, ROTMAT, SNP, TIMINC, TRANSP:  Morgan/King 981231   
c                  READM1, WRITM1: Fix comments on rfname (x, not c).  King 990105
c                  HISUB: Add check for slant ht smaller than radius.  King 990122
c                  HISUB: Add DHBGP for Trimble SST, SSE antennas.  King 990122
c                  GHDRED:  Change suicid call to report_stat and fix message.  King 990123
c                  RAND5: Fix declaration conflicts (g77). King 990129
c                  TIMINC: Fix declarations for DEC compiler.  Fang/King 990213
c                  RAND5: Use /comlib 'cand' wrapper for 'and' functions.  Herring 990217
c                  XHDRED:  Require end-of-header 'END' to be uppercase to avoid detection
c                      from UNAVCO 'End of DB comment'.  van Domselaar / King 990218      
c                  XHDRED: Replace last fix with new one testing on number of characters
c                      for 'END' or 'end'.  Fang / King 990224 
c                  HISUB: Add new Trimble antenna code (TRMICG). Bock 990228
c                  TIMCON:  Fix GPS Week > 1000 problem.  Bock 990315
c         10.11    SROTAT: Call new dirunal pole/UT1 routine (comlib/ray.f).  King 990319
c                  RAND5: Replace instrinsic 'and' by wrapper 'cand' (in libraries/comlib).
c                    (not sure why lost from 960815 and 990217 changes).  Herring/King 990322 
c                  CLKERA, RSTNFO: Do not use comma as continuation character 
c                    (IBM and f90 restriction).  Fang/King 990324
c                  HISUB, ANT_ALIAS (comment only): Add Geotracer and UNAVCO Micropulse antennas.  
c                     King 990410
c                  RRXHED, RRINEX:  Add check on nobtyp>maxobt and trap duplicate L1 entries
c                     (Geotracer 2200 RINEX files with nobtyp=9).  King 990415 
c         10.12    ANT_ALIAS, HISUB: Changed Leica and TI antennas to match agreed SIO/MIT conventions
c                      with names similar to IGS standards.   King 990423
c                  IMENU: Add double-read kluge to avoid problem reading stdin after pickfn
c                      call.  Fang/King 990429
c                  READ_RCVANT:  Fix subroutine name in report_stat call. King 990429
c                  RRXHED, SETTYP : Clean up logic and warnings for L1-only observations.  King 990505
c                  IMENU: Remove double-read kluge and also non-return in format statements to
c                     fix once and for all the problem called by reading stdin after pickfn.  King 990513
c                  HISUB: Replace TRMSLD with TRSLMC (from new IGS table) and trap TRMSLD; add
c                     revisions of ASHL12 to be designated ATGEOB and ATGEOC.  King 990529
c         10.13    SETTYP: Fix test for ndat<4 to apply always, even for L1-only.  Murray/King 990628
c                  HISUB, ANT_ALIAS (rcvant.dat). Set std names for Leica internal antennas to
c                     LC299I/LC399I with aliases LC_299/LC_399, SR299I/SR399I.  King 990628
c         10.14    Changes to support 4-digit years:  New routines CHECK_Y2K, FIX_Y2K, new version of 
c                     JULDAY,CLKERA, DAYJUL, DOYGWK, GHDRED, HISUB, READE, READJ, RRINEX, RSESFO, RSTNFO, 
c                     SVNAV_READ, THDRED, TIMCON, WRXHED, WRXNAV, WTIME, XHDRED.  King 990727  
c                  ANT_ALIAS, HISUB (+antmod.dat): Add separate entries for Leica 300-series antennas (though same as 
c                     200-series)  Morgan/King 990810
c                  HISUB: Fix bug in SLBGP for Trimble 4000 SL (10877.10, round GP).  Morgan/King 990812
c                  ANT_ALIAS: Correct comment (only) for LC299I.  King 990813
c                  GHDRED: Remove (dummy) option to read times from X-file.  King 990816
c                  RSESFO: Convert to free-format reads; add argument to by-pass day check; 
c                      allow check on day only if yr and session = 0.  King 990818
c                  READ_RCVANT: Fix comment on size of ant/rcvr names.  King 990819 
c                  HISUB: Add ASHGG1, ATDMCB, ATDMCB, ATDMCE, variations on Ashtech DM antenna.  Bock 990810
c                  RRXHED:  Fix Y2K problem with RINEX Version 1 files.  Herring/King 990916  
c         10.15    RRINEX: Fix logic for handling errors.  King 991012
c                  RSESFO: Trap blank line at end and exit gracefully.   King 991012
c                  HISUB: Add more codes for Ashtech DM chokering antennae.  Bock 991027
c                  RRXHED: Fix year variable for TIME OF LAST OBS.  Herring 991116
c         10.16    HISUB: Add add all missing antennas now available on NOAA GRD web page.  King 991127  
c                  READ_ANTMOD: Allow continuing with warning if antenna missing.   King 991127/991130
c                  HISUB: Allow UNKNWN as antenna code; fix 'unknown antenna' fatal for 4-digit year.  King 991130
c         10.17    HISUB: Add RDGMRA, TRBMRA, ASHGG1, TRSTL1, and MPLL12; added DHBCR for LC_303 and LC_504 and 
c                       warning for DHPAB with LC_504.  King 991229
c                  GHDRED:  Fix 2YK bug in reading batch-file input headers. van Domselaar/King 991230 
c                  GDHRED:  Allow 4-digit years in IC epoch.   King 000103
c                  RSESFO: 'Line' variable too short.  Fang/King 000104  
c         10.18    HISUB: Fix format for year (i2 -> i4) for fatal message.  Herring/King 000126
c                  CLKERA, READJ, RFICA: Fix minor syntax errors in format statements.  Fang/King 000215
c                  READJ: Fix read for old format, add iostat checks.  King 000302
c         10.19    RRINEX: Fix dimension problem with RINEX variables obs, jssi, jlli (maxobt, not maxdat). King 000428
c                  HISUB: Add TRUSCG.  King 000505
c                  TOPENS: Remove problematic message about T-file incompatibility.  King 000525
c         10.20    RSTNFO: Stop if sitcod blank.  King 000808
c                  HISUB: Warning only if called by HFUPD.  King 000815
c                  UT1TID:  Make implicit none.  King 000816    
c                  TOPENS:  Define program name for report_stat call.  Herring/King 000912
c                  RRINEX: Fix logic in decoding event flag.  King 000925
c         10.21    HISUB:  Add Topcon 72112 Turbo SII antenna; group Trimble 14177 (ST L1) with
c                          14532 (SST L1/L2).  King 001003
c                  READD: Trap number of x-files greater than maxnet.  King 001102
c                  HISUB: Add GEOCHR antenna with L1PHC (and DHPAB) only.  King 001106
c                  READD: Remove extra comma in last report_stat call.  King 001204
c                  CLKERA: Fix dimension of dterm2 to match Linkpack routine.  King 001204
c         10.22    HISUB: Add Ashtech 701945C_M (ATDM1C).  King 001213
c                  GHDRED: Don't print rad-model warning to screen (in arcout).  King 001220
c         10.23    NUTRED: Fix format for 'out-of-range' message.  King 010123
c                  RSESFO: Clarify warning about input session number.  King 010123 
c                  RRINEX: Change logic so that C1 will not overwrite P1.  King 010215
c                  RRXHED, WRXHED: Changes for RINEX version 2.10:  version and interval 
c                      now real, extra decimal place for time of first observation.  King 010301
c         10.24    CLKERA: Write the FIXDRV clock terms only to fixdrv.out.  King 010308         
c                  RRINEX:  Pass out the SV-type character to distinguish GPS from GLONASS, etc. King 010313
c                  GPS_ONLY, Makefile: New routine to resort a RINEX data record to include only GPS.  King 010313
c                  RRINEX: Fix dimensions of issi/illi to (,2) to match calling programs.  King 010313
c                  RRINEX: Allow more than 12 SVs on epoch line.  King 010321     
c         10.25    HISUB: Corrected if statement for NO503R.  van Domselaar / King 010511
c                  RSTNFO: In kinematic mode, fix start/stop assignment, and return trkcod blank 
c                     to indicate eof.  King 010515
c                  RRXHED: Change format for reading sampling interval from f10.3 to f10.0; this 
c                     satisfies a quirk with (at least some) Linux/g77 systems, and works on others also.  King 010515  
c                  Makefile, RSTNFO2, WSTNFO: Routines to read/write new station.info format (RSTNFO2 to be 
c                     renamed RSTNFO after further testing).  King 010515  
c                  ICLARG: Change to call rcpar from /comlib rather than getarg.  Herring 010610
c                  CLKERA: Remove extraneous comma.  Herring 010610    
c                  RSTNFO: Set trkcod blank if EOF.  King 010613
c                  READE:  Set transmission time (trans_sow) from FICA FF(3).  King 010618
c                  GHDRED: Increase # epochs for interpolation from 11 to 12 to cover unusual yawtab case. 
c                      Matt King / R King 010618   
c                  ICLARG: Remove unused variable.  King 010619  
c                  GHDRED: Increase # epochs for interpolation from 12 to 13 to cover unusual yawtab case. 
c                      **NOTE:  This changes slightly the orbital interpolation in MODEL, presumably for the better.
c                      Matt King / R King 010717/010719 
c                  HISUB: Add four Ashtech, one Novatel, one Magellan, one Sokkia, two Trimble, and
c                      one Aeroantenna antennas.  King 010718   
c         10.26    CLKERA:  Remove outliers > 5 ms from k-file fit for i-file.  King 010727    
c                  RRXHED: Remove extraneous parenthesis.  Morgan/King 010731 
c                  RAND5:  Move initialization of 'initflag' to block data.  Morgan/King 010731
c                  RRINEX: Remove extraneous comma.  Matheussen/King 010731   
c                  CLKERA: Intialize quadratic and cubic terms in linear case to avoid inversion problem;
c                      make matrix inversion warning-only; fix missing argument in report_stat call.  Herring 010802   
c                  ANT_ALIAS:  Add alias for ATDMG2 (ATDMGG). King 010821   
c                  ANT_ALIAS:  Add alais for TRMZEP (TRZEPH) and TRMZGP (TRZEPG) to compensate for 
c                       inconsistent value in rcvant.dat with release 10.05.  King 010910 
c         10.27    HISUB: Add Sercel DSNNAP_002 (SRNAP2).  King 020117/020122      
c                  RRXHED:  Remove lone reference to FERROR.  King 020306   
c         10.28    Full implementation of new station.info format: Mods to RSTNFO, RSTNFO2, Makefile   
c                       new routine CHECK_OLDSTINFO.   King 020308 
c                  JULDAY: Change maximum year to 2100 in check (for reading station.info).  King 020311  
c                  RSTNFO, RSTNFO2, WSTNFO: Change 'swver' from R*8 to R*4 to match rest of GAMIT.  King 020319
c                  JELF:  Comment out code to use floating point index computation (no longer need for HP?).  Herring 020322
c         10.29    RRINEX: Fix format error for 'time runs backwards' warning.  McClusky 020417
c         10.30    LREAD, CHECK_CRD_FILE, Makefile:  Allow reading of a GLOBK apr file as well as an l-file; 
c                     changes also to /includes/modkin.h.  King 020807 
c                  JD_TO_DECYRS, JD_TO_YMDHMS, YMDHMS_TO_JD, Makefile: Routines from kf/gen_util, added to support apr-file code.  King 020807
c                  UPNAM2, UPNAM3, Makefile:  Move upnam2 from /fixdrv to /lib to be consistent with upnam1 even though
c                     upnam2 is not used anywhere else at the moment; add upnam3 for apr files.  King 020807
c                  DEGDMS, Makefile: Moved from /tform to support use in /makex.  King 020809    
c                  HISUB:  Correct inversion of SLLGP and SLHGP for Ashtech III antennas.  King 020920
c                  LINEAR:  Rewrite interpolation to avoid bizarre g77 (2.96) compiler bug.  Herring/King 020924    
c                  HISUB: Add SOKGD3 (700969A), assumed same as Ashtech 700718A (ASHGD3).  King 020925 
c                  HISUB: Add AeroAntenna AT2775_42 (AERA42).  King     
c         10.31    Makefile, READ_GDATUM:  New routine to read old or new style gdetic.dat; replaces
c                     model/gdatum and tform/gdatum, and called in SOLVE.  King 021002  
c                  Makefile, DMSDEG:  Move from /tform (to join degdms).  King021003    
c                  HISUB: Add slant height measurements for TRMZGP (Zephyr).  King 021007    
c                  UPNAM2: Use 'newchr' rather than 'upchr' which doesn't exist any more.  King 021022  
c         10.32    HISUB: Add slant height measurements for the Ashtech Geodetic IV w/ GP.  King 021202     
c                  RSTNFO2: Change the header keywords from 6-characters to longer, more descriptive labels.  King 021212
c                  RSESFO: Add trap for number of SVs exceeded dimensions.  Tregoning 021219  
c                  RRXHED, SETTYP: Trap dual-freq obs in RINEX header of single frequency data.  King 021220   
c                  WSTNFO: Revision to use 'human-readable' headers.  Bock/Scharber/Jamason/King 021220
c                  RSTNFO2, ITIMDIF(new), Makefile: More complex code to check whether a session or epoch is within
c                     thes station.info entries.  King 030107   
c                  RSTNFO2: Remove the check for duplicates; change open-ended stop to 9999 999 0 0 0. King 030109
c                  SETTYP:  Fix bug in detecting single-frequency Ashtech.  King 030113
c                  ITIMDIF: Trap open-ended dates (9999 999) and also integer overflow.  King 030114
c                  READ_ANTMOD: Change name of initial-call variable from 'first' to 'newant'.  King 030117
c                  READ_RCVANT: Add iostat check to all reads.  King 030205
c                  RSTNFO2: Fix length of rcvrsn string; trap bogus start/stop times.  King 030205  
c                  RSTNFO2: Fix logic in check_span to avoid warning if not in span.  McClusky/King 030212 
c         10.33    LREAD, CRD_FILE_TYPE: Rewrite to identify a new- (apr-) style L-file by reading it. King 030215  
c                  HISUB:  Add Topcon Hiper GS (TOPHIP).  King 030224    
c         10.34    RSTNFO2: Fix column count in reading to avoid extraneous comments.  King 030310 
c                  WSTNFO:  Fix dimension and character declaration for 'value'.  King 030331    
c                  RSTNFO2: Fix column count for reading comments.  King 030402  
c                  HISUB: Add DHBGP for ASHL12 (PGC< 1995).   King 030410
c                  SETTYP: Extend Trimble codeless overrides (wrong wavelength in header) to other receivers toc
c                     to handle TEQC protocol of always using L2 factor = 1 and overriding with LLI bit.  King 030417 
c         10.35    Makefile, GET_ANTPCV (new), READ_ANTEX (new), READ_ANTEX_HEAD, READ_ANTMOD (new version), HISUB:  Reorganize
c                    routines to read antenna phase center offsets and variations from a table.  Main routine now
c                    get_antpcv, which calls routines to read either ANTEX of GAMIT antmod.dat.  King 030418/030423.
c                  HISUB: Remove superfluous 'icall' and put dhpab calculation at bottom.  King 030418   
c                  Makefile, GET_SVANTPCV (new):  Add calculations of SV phase center values from antmod.dat.  King 030425
c                  GET_SVANTPCV: Fix bug for recognizing BLOCK IIR SVs.  King 030430
c                  GET_SVANTPCV: Replace elevation and azimuth by nadir angle.  King 030502
c                  LREAD: Relax requirement that apr-style file have site name in column 2.  McClusky/King 030506
c                  LINEAR:  Change warning from elevation angle to zenith angle.  King 030507
c         10.36    LREAD:  Calculate all forms of coordinates (dms, radians, xyz) for either spherical of Cartesian
c                     L-file; convert to km; save in includes/modkin.h (extensively modified and reordered).  King 030509
c                  READC4, WRITC4: Changed name of radius variable from 'kradk' to 'krad', reflecting change in
c                     includes/modkin.h (no change in units on C-file: km).  King 030509
c                  Makefile, DEGDMS, DMSDEG, RADDMS(new), DMSRAD(new): Replace old versions with new
c                      ones that take the sign character explicitly and do just one task. Remove 
c                      little-used CONVRT.  King 030510
c                  ANT_ALIAS: Add aliases for Topcon antennas.  King 030519
c                  GET_ANTPCV, GET_SVANTPCV: Don't reset 'newant' here (passed from calling routine);
c                     initialize svantmod array.  King 030521
c                  RADDMS: Fix bug in sign of angle in radians.  King 030521
c                  HISUB:  Correct dimension of LC_504 by 1-2 mm and remove warning message.   Casula/King 030527
c         10.37    SORT_STRING: Fix logic in removing duplicates.  King 030528
c                  HISUB, RSTNFO, RSTNFO2: Allow DHARP (as well as, temporarily, DHPAB).  King 030529
c                  DECYRS: Fix bug in JD interpretation (off by 0.5 days).   King 030604
c                  Makefile, DECYRS_TO_YDHMS (new).   King 030605
c                  DMSDEG: Fix bug in test of negative degrees.  King 030605
c                  RTNFO2: Fix infinite-loop possibility with station.info missing *SITE line. King 030606
c                  DEGDMS: Fix bug. King 030623                           
c                  ANT_ALIAS: Add ATDM1D and ATDM1E.  King 030623
c                  LREAD:  Fix bug in quotes around iostat for reading minutes of latitude. King 030711
c         10.38    GHDRED: Increase number of epochs padded on ends of integration from 13 to 14 to
c                    handle occasional yaw problem.  King 030721   
c                  HISUB: Restore height-hook offset for LC_504 to 0.350.  Casula/King  030729        
c                  XHDRED: Trap empty X-file on first read; prevent infinite loop reading comments 
c                      when file empty.  King 030805
c                  LREAD: Fix two minor compile bugs.  Herring 030815
c                  WSTNFO: Add first-column comment character to calling sequence and write statement. King 030815
c                  HISUB: Add Ashtech 700228.A with extended groundplane (ASHLEX).  King 020816
c                  RSTNFO: Add fix_y2k call for stop time.  King 030926
c                  GET_ANTPCV: Rewind iant if model not found to avoid ungraceful abort on next read.  King 031013
c                  HISUB: Add DHTGP for ASHLEX.  King 031013
c          10.39   GET_ANTPCV, READ_ANTEX, READ_ANTEX_HEAD: Add ANTEX version and SINEX code to ANTEX reads.  King 031015
c                  RRXHED: Add code to convert (illegal) lowercase observables to uppercase.  JiangWP/King 031016
c                  HISUB: For Leica antennas, distinguish between 200/300-series height-hook offsets (.345m) and 
c                      500-series offsets (.350m).  King 031021
c                  GET_SVANTPCV: Fix calling arguments for READ_ANTEX to include ANTEX version and SINEX code.  King 031021
c                  HISUB, RSTNFO: Trap unreasonable year before call to julday, to report site name. King 031025
c                  STNFO: Fix failure to skip comment lines in next-record-read (icall=2) case.  King 031025 
c                  HISUB: Issue warning on Leica antenna only once.  King 031025
c                  THDRED: Subsititute sb uppers for fn upperc avoid length mismatch (ifc comiler). King 031026
c                  YAW_ATTIT: Fix line > 72 character (comment only)  
c          10.40   GETDAT: Added fix_y2k call to make sure system times are returned in y2k format. McClusky 031022
c                  GRTDAT: Getdat is now calls /libraries/comlib/systime, not the idate routine directly. Mcclusky 031022
c                  GETUSR: Modified to exclude a usrname that is null not just blank. McClusky 031029
c                  WSTNFO: Fix bug in column width for 'Station Name'   King 031106
c                  READ_RCVANT: Add ability to open gg/tables directly without local link; close unit number first
c                    before opening.  Herring 031111
c          10.41   RSTNFO: Trap bad 'sitcod' entry.   King 031205
c                  HISUB: Add NavCom AN2004, AN2008, SF204G, and RT301S.  King 031217
c                  IC4ARRAY, Makefile:  Add function to file the index of a character*4 (e.g. site) token in 
c                     and array.   King 031231   
c                  HISUB: Add slant height measurements for the Ashtech Marine antennas.  McClusky/King 010406
c                  CLKERA: Increase tolerance for bad PR from 5 ms to 300 ms to allow a jump to be inserted.  McClusky/King 040115
c                  GET_ANTPCV, GET_SVANTPCV: Add 'sinex_code' to call to get_antex_head and use for getting model name
c                       for p- and h-files.  King 040120
c                  GET_SVANTPCV: Fix bugs in computing SV antenna PCV corrections.  King 040120
c                  HISUB: Add Novatel WAAS antenna with NOV600 element (NOV6WA). King 040123   
c                  LINEAR: Trap and set output = 0. any zenith angle greater than the table boundary. King 040202
c                  GET_SVANTPCV: Fix undefined 'nel' when SV PCV model not found on ANTEX file.  Tregoning/King 040205
c                  RDSEST: Don't convert output to uppercase when it's a Unix file name.  King 040223
c                  HISUB: Add SOK600 (DHARP only) to Novatel group (SOK600 same as NOV600). King 040317
c                  HISUB: Add Leica AX1202 (LC1202) (DHARP only).  King 040407    
c                  COPENS: Fix report_stat message format.  King 040524
c          10.42   CRD_FILE_TYPE: Fix identification of L-file type so that short files can be used.  King 040614
c                  LREAD: Trap extra line at bottom indicating eof.  King 040616
c          10.43   LREAD: Fix bug with embedded comments in apr-style L-file.  King 040727
c                  GSATEL: Fix check of end-time, too tight for orbdif w/ matching spans.  King 040806
c                  HISUB: For Trimble Zephyr, make the top of the plastic cover 3 mm above the bottom. King 040825
c          10.44   GET_ANTPCV, Makefile:  Replace bilin.f with Tregoning version from  utils/grdtab.f (but make real*8)
c                    for az/el  PCV interpolation (attach rather than making it a separate program); set output 
c                    antmod = 'NONE' rather than blank if input = 'NONE'.   King 040903
c                  READ_ANTEX: Fix bug in reading azimuth-dependent values.  Schulte/King 040903
c          10.45   READ_RCVANT, RSTNFO2, GET_ANTPCV: Add reading of the code for converting C1 and P2' for cross-correlating receivers. King 040930
c                  HISUB:  ASH701945A_M (ATDMRA), assumed identical to ASH701945B_M (ATDMRB). King 040930
c                  ANT_ALIAS: Add TRBTMA with alias TRBRMA as a result of IGSCB typo.  King 041001
c                  HISUB: Add ASH701945G_M (ATDM1G).  King 041013
c                  RSTNFO2: Add start/stop times to fatal trap for earlier/later.  King 041014
c          10.46   RSTNFO2: In check_span, flag 'early' if start 'le' rather than 'lt' time -60 to get warning
c                    rather than 'no overlap' (fatal).  May need more tweaking here.  King 041111
c                  LREAD: Fix site variable in report_stat call for missing site.  McClusky/King 041128
c                  HISUB: Add JNSMARANT_GGD (JNSMAR).  King 041222
c                  READ_ANTEX: Initialize jdstart, jdstop to zero so that they will not be used.  Herring 041227
c                  RSTNFO2: Change initialization of radome from NONE to UNKN.  Herring/King 041230
c                  PJDHMS, Makefile: Move from ARC.  King 050103 
c                  GET_ANTPCV: Change name of subroutine bilin to bilin8 (real*8) to avoid confusion with
c                      real*4 version in attached to utils/interp_atm.f.   King 050112
c          10.47   READC1, WRITC1, WRITC2, READC2:  Change C-file version to 10.2 
c                       (models and extra added to record 2.  King 050201    
c                  XHDRED: Remove 'extra' variables from calling arguments.  King 050201
c                  Remove obsolete getdir.f from directory (replaced in 1998 by routines in /comlib. King 050201
c                  RSNTFO2: Fix bug in that nlist was not saved along with first_call.  Tregoning 050201
c                  ANT_ALIAS: Add TRBRMA as alias to TRBTMA.  Shimada/King 050204
c                  GET_ANTPCV, GET_SVANT_PCV, READ_ANTEX: Add antenna SN (or PRN) for ANTEX. King 050208.  
c                  GET_SVANT_PCV: Add channel # and PRN # to calling arguments.  King 050209 
c                  READ_ANTEX_HEAD, READ_ANTEX: Allow versions up to 1.2.   King 050209
c                  GET_ANTPCV: Initialize unused serial number to blank.  King 050209
c                  WSTNFO: Truncate comments with last non-blank character.  King 050215 
c                  READ_RCVANT: Issue warning not fatal if calling program is mstinf or mstinf2.  King 050216
c                  HISUB: Fix misnamed TRBMRA to TRBMTA.  King 050216
c                  HISUB: Add TPSPG_A (TPSPGA) and TPSPG_A+GP (TPSPGG).  King 050311
c          10.48   HISUB: Add slant height for Topcon choke ring (TPSC3D or TPSC3R).  King  050314
c                  HISUB: Correct typo in ID for TPSPGA and TPSPGG.  King 050314
c                  HISUB, ANT_ALIAS: Correct typeo in TRBMRA/TRBMTA.  King 050315  
c                  HISUB: Add Topcon CR4 (TPSCR4 and, with snowcone radome TPSCC4). King 050315
c          10.49   GET_ANTCV: Set model to 'NONE' when not found on antmod.dat.  King 050426
c                  LREAD: Fix bug causing problem for MAKEX when duplicate entries for GLOBK-style l-files. Herring 050506
c          10.50   READ_ANTEX, READ_ANTEX_HEAD: Update to version 1.3.  King 050618
c                  READ_RCVANT, RSTNFO2, GET_ANTPCV, HISUB: Add radome to calling arguments and checks; issue warning if
c                    radome is missing.   King 050618
c                  GET_SVANTPCV: Expand 'svantcod' and block id to include IIR-A and IIR-B.  King 050618
c                  READ_ANTEX: Set radome = NONE if blank in ANTEX file.  King 050622
c                  GET_SVANTPCV:  Fix report_stat warning for missing SV antenna.  King 050622
c                  GET_ANTPCV: If antenna+radome not found, use antenna+no-radome model.  King 050622
c                  SVNAV_READ: Read a 'B' or 'b' after block number and change '4' to '5'. King 050623
c                  GET_SVANTPCV:  Suppress warning of missing SV antenna when relative model used. King 050623
c                  JD_TO_DECYRS: Correct spelling error in comment.  King 050627
c                  HISUB: Add antennas from NGS 04 table:  AERA41, AERA62, AER160, ASHRAS, ATPROM ATDMRD, 
c                    ATDMRE, THAZMX, LC399A, MPLCHR, MPLL1W, MPREGW, MPLWNW, NOV512, NO7022, NO7023,   
c                    SENS14, SEN14R, SENS49, SEN49R, SENS96, SEN96R, SOK702, SOKA11, SOKA12, SOKSTR,
c                    SPP8GP, PSHIG, TPSLIT, TPSLG2, TPSLG3, TPSLGG, TRMI13, TR4800, TR5800
c                  READ_RCVANT: Match only the antenna name (not radome) when converting.  King 050628
c                  ANT_ALIAS: Add JNS3CD and JNS3CR for TPS3CD and TPS3CR.  King 050708
c          10.51   GET_ANTPCV, HISUB, RSTNFO2: Add logical value 'warnings' to decode_values call; set = T for 
c                     GAMIT calls but may be F for calls by mstinf2 and hfupd.  Herring/King 050719  
c          10.52   RSTNFO2, READ_RCVANT: Do not set radome=NONE if unknown; set radomes from GAMIT 
c                     code or IGS name for a few cases.  King 050812/050819   
c          10.53   GET_ANTPCV: Don't change the original radome code from UNKN but use NONE to
c                     read an ANTEX file; change 'warning' to 'fatal' if still not found.  King 050818
c          10.54   GETCMD, MCHKEY, Makefile: Make mchkey a separate routine to be called from other
c                     library routines and modules (needed for fixdrv/somake).  King 050819
c                  GET_ANTPCV, HISUB: Change calling arguments to pass out radome subsitution flag
c                     and minimum elevation of the PCV table.  King 050819
c                  READ_RCVANT: Warn, not stop for hfupd (as well as mstinf2) calls.  King 050902
c          10.55   XHDRED: Change to allow reading of ARP instead of L1 L2 phase center offsets. King050928 
c                  Makefile, HISUB2: New routine to read hi.dat table instead of using internal values. King 050929
c                  HISUB: Round radius for TI4100 to nearest 0.1 mm for consistency with hi.dat table.  King 051001
c                  ANT_ALIAS: Make RDGMRA alias for ROGMRA (not clear why first code arose).  King 0501003
c                  SVNAV_READ: Fix bug causing iblk=5 to be assigned to the wrong table entry (problem
c                     only when the temporary '4b' scheme used in svnav.dat).  King 051019
c                  GET_SVANTPCV, SVNAV_READ:  Allow iblk=6 for Block IIR-M (also, temporarily, read '4M' or '4m') King 051021
c                  READ_ANTMOD: Remove redundant read of offsets; set zen2=90 for minelv check. King 051108
c          10.56   ANT_ALIAS: Add 'zero-centered' aliases (JNSZCT JNSZCR) for JNS versions of TPSC3D TRPC3R. King 051118
c                  GET_ANTPCV: Use different warning for radome mismatch when radome UNKN.  King 051129
c                  READ_ANTEX: Add space in warning; initialize SINEX code to blank.  King 051130
c                  DMRRAD: Fix trap for negative degrees.  Nercessian/King 051208
c                  RRXHED: Change trap for C2 meaning P2 to work only for old observations since 
c                     there is now a true C2.  King 060111
c                  GET_ANTPCV:  Initialize SINEX code to avoid undefined characters.  King 060130
c          10.57   READ_ANTEX: Remove gratuitous warning about missing radome.  King 060211
c                  RRXHED: Add special fix for column misalignment of 'nwave'.  King 060222
c                  READM1: Add traps for exceeding maxsit, maxsat, maxprm.  King 060222
c          10.58   RSTNFO2:  If makexp is the calling program, only a warning, not a fatal is invoked when
c                    a linked RINEX file has data outside the range of station.info entries.  King 060522.
c                  RRXHED: Change the read of the RINEX version number to work for both integer and real.
c                     King 060526
c                  HISUB2: Fix bug that allowed 0. for HI when htcod missing.  McClusky/King 060724
c          10.59   RRXHED:  Trap bogus first line.  King 060707
c                  RSTNFO2: Add iostat to read of start, stop times.  McClusky/King 060707
c          10.60   RSESFO: Changing arguments to return hr,min rather than GPS week, and assume
c                     that session.info is in GPST.  King 060815
c                  LREAD: Change name of include to from modkin.h to model.h.  King 060818 
c                  GETCMD: Added iostat to read of solve batch file. McClusky 060818
c                  COPENS, READC1: Move check for binary compatibility from copens to readc1.  King 060904
c          10.61   DLAPACK: Added lapack inversion routines. McClusky 060908
c                  MAKEFILE.GENERIC: addedd dlapack. McClusky 060908
c                  INVER2: modified to use lapack inversion routines. McClusky 060908
c                  LREAD: Add site name to report_stat message.  King 060912
c          10.62   OTLCMC, Makefile: New routine to read the center-of-mass correction for ocean tidal
c                     loading models from a (new) file oltcmc.dat.  King 061031
c                  GETCMD: Add iostat in reading file.  McClusky 061106
c                  TOPENS: Change logic for detecting binary incompatiblity so that it works for 
c                    gFortran. King 061106/061107
c                  OTLCMC: Add 'status=old' to open statement.  King 061107
c                  COUNT_ARG, NBLEN, BLANK, LIFT_ARG, STRJST : Fix problematic character declaration for gFortran.  McClusky 061108
c          10.63   DLAPACK, INVER2:  Use condition number (IER 140) instead of pivot element size (IER 130)
c                     to test ill-conditioning of matrix (fixes problem with dependent biases).  McClusky 061115c
c          10.64   RSESFO: Remove reset of session so that isessn=0 can be used to indicate matching of
c                    day only (not year or time).  King 061207
c          10.65   INVER2: modified tolerance on the RCOND check for dependent parameters. McClusky 061222
c                  Makefile, FFUN, GPT, WPRESS, SAASZD: Move routines from MODEL to support GRDTAB.  
c                      Tregoning/King 061228
c                  READ_ANTEX: Fix (normally harmless) column-comparison and 'int' rounding.  Herring/King 061229
c                  TOPENS: Remove unused variable.  King 061229
c                  LREAD: Add site name to report_stat message on multiple entries.  McClusky/King 070110
c                  READ_RCVANT: Add convert_antpcv to list of programs invoking warning, not fatal; rewind
c                     file after warning to allow multiple calls.  Herring/King 070115  
c                  Makefile, LREAD, RENAME_MATCH: Finish coding to read eq_rename.dat.  King 070116/070118
c                  FIX_Y2K: Fix bug when year >2099.  King 070117
c            10.66 Makefile; remove saaszd, ffun, wpress, no longer needed by grdtab.  King 070 124
c                  LREAD: Change filename from 'eq_rename.dat' to simply 'eq_rename.  King 070125
c                  Makefile: Set g77-specific FFLAGS to get optimization lowered to 1 to work around
c                    gcc/g77 3.4.x compiler bug with Lapack routines.  McCluskey/King 070130
c                  DLAPACK: Add a dummy print statement to subroutine dlamc1 as an expedient (temporary?) workaround 
c                    for an apparent compiler bug in gcc/g77 3.4.x under level 3 optimization  McCLusky/King 070131.
c            10.67 SVNAV_READ: Add iostat to open failure message. King 070208
c                  HISUB2: Fix alias for Ashtech choke rings (error is only 1 mm).  King 070221          
c                  RSTNFO2: Initialize and save a few more variables.   King 070305
c                  LREAD: Fix bugs in using eq_rename to update coordinates.  King 070307     ,
c                  INVER2: Changed bias condition matrix test from abs value to bias to bias ratio. McClusky 070319
c                  INVER2: Removed condition matrix test. Condition number - RCOND passed in argument list for interpretation
c                          by calling routines. McClusky/King 070319
c                  GPT: Change 'end subroutine' to 'return', then 'end' to satisfy Solaris compiler.  King 070326
c            10.68 Makefile: Remove sb atsub.f.  King 070416
c                  EPHDRD, GHDRED, GPS_ONLY, GSATEL, OCEARG, READD, READJ, SHADOW1, SVNAV_READ, THDRED, TAIUTC,
c                     UT1TID, XHDRED, YAW_ATTIT:  Remove or use unused variables  and labels to avoid compiler warning.  King 070416 
c                  HISUB: Add antenna code to warning message.  King 070416
c                  HISUB2: Temporarily add dummy statement to use input argument 'warnings'.  King 070416
c                  LVERSN: Remove unused arguement 'iun'. King 070416
c                  NUTRED, NUTTAB: Remove unused 'iscrn' from nutred call.  King 070416
c                  POLRED, UT1RED, SROTAT: Remove unused 'iscrn' from polred call.  King 070416
c                  ROT_GFILE: Change logic to avoid obsolete got; removed unused variables
c                     from calling sequence.  King 070416
c                  RRINEX: Remove unused arguments from calling sequence. King 070416  
c                  SROTAT: Remove unused 'jd' 'fract' from call to 'ut1tid' King 070416
c            10.69 RSTNFO2: Redetemine entry headers if 'mstinf2' is the calling program since there may be more
c                    than one station.info file; add site id and day to error message on decoding entries. King 070501
c            10.70 GET_ANTPCV: Fix bug in interpolating AZEL values.  Wang/King 070519
c                  SETTYP: Add code to override mistranslated L1-only data from TR4600 (swver=2.99).  King 070607
c                  SETTYP: Trap RINEX wave factor = 0, issue warning and set lambda = 1.  King 070711
c                  RSTNFO2: Make test for missing receiver or antenna smarter; uncomment 'warnings' set (bug).  King 070726
c                  READ_ANTEX: Fix minor bug in testing for blank serial number (translation case). King 070726
c                  NUTTAB, PNROT, READJ, ROT_GFILE, ROTSNP, SROTAT:  Remove unused 'iscrn'.  King 070906
c                  DBLAS2: Comment out unused format statement.  King 070906
c                  TAIUTC: Fix up format for reporting date past table.   King 070912
c                  NYDAYS, Makefile: Add function to get number of days in a year.  King 071023
c                  TIMINC: Add 1.d-15 to sod to avoid roundoff incrementing day.  King 071105
c                  GET_ANTPCV: Add temporary code to set IGS05_1451 to I2_I1451 so that the first
c                     two characters will be unique in the 4-character c-,h-file designation.  Herring/King 071106
c                  READ_RCVANT: Add missing close statement before return after receiver-not-found warning. King 071113
c                  RRXHED: Change highest allowed RINEX version from 2.10 to 2.11 (warning only).  King 071221
c            10.71 GHDRED: Add nutation model and gravity field to g-file.  King 071221/071227
c                  RRINEX: Fix format statement for report_stat messaage. McClusky 080130
c            10.72 Makefile: Remove xyz_to_geod in favor of version now in /libraries/comlib;
c                      add MMPLY from ../solve.  King 080307
c                  SETTYP: Fix test for incorrect L2 wavefactor to account for 0 (L1-only receiver). King 080331 
c                  GET_ANTPCV, HISUB: Add 'debug' to calling arguments for get_antpcv.  King 080505
c                  Makefile, CHECK_OLD_STNFO: Remove old-style hisub and rstnfo and and rename hisub2 to 
c                    hisub, and rstnfo2 to rstnfo; keep check_oldstnfo temporarily, with fatal message. King 080509
c                  LREAD: Add site name to error messages.  McClusky 080510
c                  LREAD: Fix logic for handling the first rename (a problem only when a break mid-session). King 080513
c                  GET_ANTPCV: Fix bug in getting PCVs for mid-session change of antmod.  King 080522
c                  LREAD: Fix problem of missing pre-rename site.  King 080522
c                  LREAD: Fix problem with two many renames.  Herring 080709
c                  RSTNFO: Change message about incompatible station.info format.  King 080717
c                  ANT_ALIAS: Make ATDMR2 equivalent to ATDMRB.  King 080930
c                  READ_ANTEX_HEAD; Trap too many comments.  King 081006
c            10.73 WSTINFO: Change declaration of 'line_comment' from C824 to *(*).  King 081020
c                  READ_RCVANT: Don't print warning if receiver or antenna name is '--------'.  King 081020
c                  RSTNFO: Increase size of end-of-line comment to 30 characters (for mstinf).  King 081031
c                  Makefile: Remove JD_TO_YMDHMS, JD_TO_DECYRS, YMDHMS_TO_JD, DECYRS_TO_YDHMS, DS2HMS, LEAPYR (in /comlib). 
c                  CHECK_OLDSTNFO: Uncomment 'old_stinf = .true.' to satisfy compiler (though not needed). King 081118
c                  RRXHED: Allow first line of RINEX to be lowercase (warn only). King 081128   
c                  DBLAS1: Mods to conform to f90.  Morgan 081203
c                  RRXHED, RRINEX: Mods to accommodate # observables = 10. King 081218
c                  SETTYP: Trap and fix bad header for Ashtech codeless when swver = 1.0.  King 090203
c                  RSTNFO: Add trap for typos leading to unreasonable year entry in station.info. King 090410
c            10.74 READDF: Change name from 'readd' to avoid conflict with /libraries/comlib/readd used by /kf.  King 090415
c                  RRINEX: Allow more than 10 observation types; fix format in debug print.   King 090514
c                  RRINEX: Fix bug in 090514 changes, causing problem if numobt < 6.  King 090520 
c                  RSTNFO: Fix bug in detecting '-----' in rctype. King 090529
c                  RSTNFO: Fix bug in calling arguments introduced on 090529. King 090602
c            10.75 HISUB: Change 'hisub2' to 'hisub' in report_stat calls. King 090826
c                  RSESFO: Fix typo iin error message. King 090826
c                  RSTNFO: Fix bug in distinguishing between 'Receiver Type' and 'Receiver SN', and 'Antenna Type'
c                    and 'Antenna SN'.   Charade/King 090915
c                  JULDAY, RSTNFO: Fix trap for unreasonable day.   Chen Jixiang/King 090917/090922
c                  LREAD:  Change default stop time to 2100 1 rather than 2100 0 to avoid possible problem in julday. King 090922
c                  RSTNFO: Fix error traps in reading header and also in reading values.  King 091104
c                  SETTYP: Fix logic for setting lambda when L1-only.  King 091228
c                  GET_ANTPCV: Add site code to calling argument for use in error messages.  King 100111
c            10.76 RSTNFO, WSTNFO: Remove obsolete 'trkcod' variable.  King 100207
c                  GET_ANTPCV: Clarify documention of output variables.  King 100207
c                  Makefile: Remove CHECK_OLD_STNFO (needed only for transition). King 100209
c                  LREAD: Add dimpar.h and change 'ksite' to 'sitecd' to accommodate changes in model.h.   King 100209
c                  ANT_ALIAS: Added SOKGD3 as alias for ASHGD3.  King 100326
c                  RSTNFO: Change span-match slop from 60s to 5s.  King 100331
c                  LINEAR: Limit warning message for out-of-bounds interpolation to once.  Herring/King 100406
c            10.77 RRINEX; Increase allowed SVs at an epoch from 24 to 48; initialize 'tform' variables to
c                    avoid binary characters in debug output.  King 100407
c            10.78 GET_SVANTPCV:  Allow iblk=7 for Block IIF  King 100529    
c            10.79 READC1, READC2, READC4, WRITC1, WRITC2, WRITC4:  Add the full 10-character ground and
c                   SV antenna model to record 2, and the atmospheric loading values to record 4.  King 100827
c                  RSTNFO: Add rcvers, rcvrsn, antsn to calling arguments and initilize blankInitialize to blank. King 100908
c                  READ_ANTEX, READ_ANTEX_HEAD:  Generalize for more than 2 frequencies; allow ANTEX version 1.4.  King 100928
c                  READ_ANTEX: Fix bug in matching antenna serial numbers.  Petrie/King 101026
c             10.80 RRINEX, SETTYP : Allow use of C2 if P2 not available.  King 101028
c                   READ_ANTEX: Fix bug in setting antenna serial number.  King 101104
c                   Makefile: Move get_antpcv, get_svantpcv back to /lib since used by kf/htoglb.  King 101110
c             10.81 RSTNFO: Fix case sensitivity for radome.  King 101111
c                   RRXHED: Fix minor typo caught by Intel compiler.  McClusky 101112
c             10.82 WRXHED: Allow more than 9 observables to be written in RINEX header (merge_rinex). King 101231
c             10.83 GHDRED: Allow new Univ College London radiation model. Petrie/King 110118
c             10.84 SROTAT: Allow IERS2010 model for short-period EOP.  Watson/Herring 110214
c             10.85 RRINEX: Allow C2 to be used if P2 zero.  King 110420
c             10.86 RRXHED, RRINEX: Fix bugs in reading RINEX files with more than 10 observables. Herring/King 110722
c             10.87 READ_ANTEX: Fix bug in setting stop time.  King 110818
c             10.88 DBLAS1, DBLAS2, DBLAS3 and DLAPACK updated to current versions of routines called by GAMIT/GLOBK programs. Floyd 111102
c             10.89 READ_ANTEXT: Allow for 0.5 degree ANTEX files. Moore/McClusky/King 120104.
c             10.90 RRXHED/RRINEX:  Allow reading 20 observables. Herring 120126/120202
c             10.91 Makefile, new routine READ_GFILE, remove GHDRED. King 120518
c                   ANT_ALIAS: Add equivlanece of the Altus APSAPS-3 and TI Asahi TIAPENSMT8883G antennas. King 120521 
c             10.92 Makefile, READ_GFILE, new routine ASSIGN_SRPNAME:  Separate out the assignment of radiation-pressure
c                    parameter names; add the time_type to the calling arguements of read_gfile.   King 120531
c             10.93 GPT2(new), Makefile: Replace GPT with GPT2 for a priori P,T,RH, and replace
c                     GMF with GPT2 for mapping functions.  King 121211/121219
c                   READ_GFILE: Change test on UTC/GPST to warning, not fatal. King 121219
c             10.94 READ_ANTEX: Correct bug causing bad L2 values whenever the 'FREQUENCY RMS' entires are
c                    present (most recent NGS entries).  King 130130
c             10.95 READC1, READC2, WRITC1, WRITC2: Add dryzen and wetzen to record 2 of c-file. King 130117
c                   ANT_ALIAS: Add original provision name for ASH802147_A. King 130305
c             10.96 THDRED: Add Earth radiation and antenna radiation model names to call. King 130320
c                   RRINEX: Fix bug in reading more than 36 SVs.  Herring 130424
c             10.97 SROTAT, new UTLIBR, Makefile: Add libration terms to UT1. King 130822
c                   POLRED, UT1RED: Fix 1-day interpolation error and bug preventing many-day integrations in ARC. Hayashi/King 130822        
c                   ANT_ALIAS: Add JAVTRIANT as an alias of JAVTRIANT_A as per IGS change. King 130822
c                   READ_GDATUM, READ_ANTEX; Correct # arguments in one call to report_stat.  King 131011
c                   RSTNFO: Correct typo in setting blank string. King 131011
c             10.98 ROUND6H, Makefile: Move this routine from being embedded in various grdtab source files. King 131025
c                   RSTNFO: Allow fewer than 5 '-----' in decoding swver. King 131105
c                   READ_ANTEX: Allow blank instead of '0' in frequency code on ANTEX file.  King 131121
c             10.99 YAW_ATTIT: Restore the sign reversal for the S/C body X-axis for Blk IIR SVs to be
c                     consistent with the USAF definition and the Kouba yaw model.  King 140108
c                   New routine BLOCK_NAME. King 140123
c             11.01 READC1, READC2, WRITC1,WRITC2: Add E-radiation and antenna-radiation models to the c-file 
c                    header; new C-file version is 10.42. King 140327
c                   READC2, WRITC2: Add ionsrc and magfield for 2nd-order ion corrections to the c-file
c                     headers; kept version 10.42.  King 140401
c             11.02 Makefile, SATATT: Move SATATT from /model to /lib (used by ARC for the UCLR1 and UCLR2
c                     radition-pressure models.  King 140416 
c             11.10 Major mods to change C-file and svnav.dat formats, and to support GNSS: Makefile, GET_GNSS,
c                     LREAD, READC1, READC2, READC4, READC5, READE, READ_GFILE, RINEX, RRXHED, RSESFO, SVNAV_READ, 
c                     THDRED, WRITC1, WRITC2, WRITC4, WRITC5, XHDRED.  King 141206   
c                     GPSBLOCK: New temporary routine to get GPS block # for radiation pressure and yaw. King 150109
c                   EPHDRD: Fix program name in report_stat calls.  King 150108 
c                   GPSBLOCK: Fix bug for IIF.  King 150316
c             11,11 SVNAV_READ, READ_GFILE: Add start/stop times to calling arguments for snav_read; 
c                   don't abort if called by DCBTAB, DCBTAB2, or ARC. King 150520/150527
c                   ITIMDIF, RSTNFO, SVNAV_READ: Make the itimdif function I*8.  King 150528
c             11.12 RSP3HD, Makefile: Move rsp3hd from /orbits for use by /makex. King 150730
c                   TIMCON: Correct format for conversion error in report_stat message. King 150731
c             11.13 ITIMDIF: Fix erroneous documentation of sign (no change in routine). King 150831
c             11.14 SVNAV_READ: Add GNSS code to error messages. King 150918
c                   GPSBLOCK: Fix erronesous program name in report_stat call. King 150922
c             11.15 GET_GNSS, RRXHED, RRINEX, SEL_OBTYP (new), mv SETTYP from /lib, Makefile: 
c                     Add the ability to read RINEX 3.  King 151014
c                   GET_GNSS: Remove debug print statements King 151020
c             11.16 ASSIGN_SRPNAMES: Remove old models and change parameters for BERN2.  King 151021
c                   SROTAT: Enable the Gipson short-period EOP model by Gipson. Wang 151020
c                   SD_COMP, Makefile: Move to libraries/comlib. Wang/King 151022
c                   RRINEX: Fix problem with no non-gnss SVs at RINEX epoch. King 151022
c                   RRINEX, SEL_OBTYP: Fix bug in order of observables (P1 over C1) and allow for C1 or C2 
c                     to be used if P1 or P2 is not available. King 151028
c             11.17 RRXHED: Change 'RCV CLOCK OFFS APPL' from fatal to warning. King 151106
c             11.18 RRINEX: Fix bug in selecting C2 as backup for P2.  King 151110
c             11.19 YAW_ATTIT: Replace GPS block numbers by sv body type. King 151112  
c                   LREAD: Add code to remove comments delineated by ! . Herring/King 151113
c             11.20 READE, READJ, RSP3HD: Mods to read all GNSS systems.  King 151125
c             11.21 READJ: Allow reading of pre-GNSS j-files.  King 151127
c                   SEL_OBTYP: Add report_stat call for unknown GNSS. King 151130
c                   READE: Add time-system conversion from UTC to GPTS for Glonass. King 151207
c                   RRINEX: Add code to restore leap second to Glonass pseudoranges. King 151208  
c                   READE: Get reference time from clock epoch.  King 151208 
c                   TIMCON: Add the UTC-to-GPST conversion. King 151209
c                   RSP3HD: Edit comments for selected SVs. King 151211
c                   WTIME: Fix error message for unrcognized conversion time. King 151211
c                   RSP3HD: Fix report_stat module names.  King 151215
c                   READE: Zero out unused bclock and bephem values and change sign of TauN. King 151222
c                   SEL_OBTYP: Add C1 as backup for P1 for Glonass.  King 151228
c                   READ_ANTEX: Match GNSS system as well as frequency.  King 151229
c                   RRINEX: Fix bug in reading RINEX 3 data values when requested SVs not first. King 160106
c                   RRXHED: Changing warning to test verions > 3.02 rather than 2.11.  King 160107
c                   RRXHED, SEL_OBTYP: Add code to change Beidou C2/L2 to RINEX standard C7/L7, and
c                     fix type is selecting the observation type. King 160107
c                   READ_ANTEX: Add BEIDOU totest on no radome.   King 160108
c                   SEL_OBTYP: Add L2/C2 as higher-frequency Beidou option (change with RINEX 3.03).  King 160108
c                   RRINEX: Remove the code adding the leap-second to GLonass pseudoranges. Herring/King 160113
c                   READE: Fix bug when a blank line at the end of the nav file. King 160118
c                   READE: Flag records with bogus year to be skipped by calling program. King 160219
c                   SEL_OBTYP: Add RINEX 2.12 observable codes.  King 160219
c                   MERGE_RINEX, WRXHED: Add arguments for gnss and rxtime. King 160323
c                   RSP3HD: Change 'line' length from 69 to 80 for sp3-d.  King 160408 
c                   READ_RCVANT: When converting from full name to GAMIT code, using 16, not 15 character
c                     to compare.  King 160419 
c                   XHDRED: Remove warning for short or missing expanded data line, now ok.  King 160420
c                   RSP3HD: Remove sp3gnss from calling arguments (never needed)..  King 160602 
c                   READE: Add IRNSS with temporary fix to incorrect week number   King 160805
c                   SEL_OBTYP: Add IRNS observables. King 160808
c                   RRXHED: Fix warning for RINEX version (allow up to 3.03). King 160817
c                   SEL_OBTYP: Give C1W precedence over C1C for GPS; switch order of frequencies for IRNSS. King 160819
c                   READ_ANTEX: Add IRNSS to if-statement checking whether an SV antenna.  King 160823
c             11.22 CLKERA: Add ability to read version 2 k-files.  King 160901 
c             11.23 ROTSNP, ROT_GFLE: Add pole position values to calling sequence. King 170412
c             11.24 SATATT, YAW_ATTIT: Fix calling argument and remove sign reversal for IIR body axes 
c                      (now IGS convention). King 170523 
c                   XHDRED: Read the observable types from the x-file.  King 170531
c                   SEL_OBTYP: Add L1X, C1X, L5X, and C5X to Galileo observable options. King 170531 
c                   READ_ANTEX: Refine logic for getting GNSS calibrations.  King 170601
c                   READE: Trap week number originating at W1356.  King 170602 
c                   SEL_OBTYP: Add L1X, C1X, L2X, C2X to Beidou observable options. King 170607
c                   WRXHED: Fix rxobtp to character*3.  King 170720
c                   EPHDRD: Change unit numbers for lunar and solar ephmerides to 34 and 35 in order
c                      to assign the ocean tide tables 33.  King 170925
c             11.25 RSP3HD: Mods to read sp3d.   King 180206 
c                   EPHDRD: Remove calling arguments in favor of commons in arc.h.  King 180319
c                   LUNRED, SOLRED, Makefile: Move from ARC. King 180319
c                   THDRED: Change stop time from jds,ts to jdf,tf to be consistent with 
c                      other modules.  King 180320 
c                   READE:  Fix clock sign error for Glonass. Herring 180320
c                   EPHDRD: Change start-time variable names here and in includes/arc.h from jd0,t0 
c                      to jdb,tb to be less ambiguous and consistent with other modules. King 180320 
c                   EPHDRD, LREAD, LUNRED, SOLRED: Add new common includes/units.h to provide unit numbers
c                    removed from includes/arc.h or includes/model.h King 180320
c                   XHDRED, THDRED: Change 'skd' to 'gnss' in calls.  King 180322 
c                   RSTNFO: Remove obsolete kinematic code.  King 180322 
c                   Makefile: Move EVRTCF from /arc to /lib.  King 180324 
c                   SEL_OBTYP: Add 'C2 ', and 'C2P' as backups for Glonass.  Herring/King 180326  
c                   MHB_2000, NUTTAB, Makefile: Add nutation routines from kf/gen_util/precess. King 180328
c                   EPHDRD: Pass fjd, not jdb,tb in tests for ephemeris span. King 180404
c                   READDF: Skip first two lines of d-file; remove multi-session options. King 180427
c                   EPHRED: Fix bugs in setting times.  King 180430
c                   SEL_OBTYP: Fix bug in selecting RINEX L1 phase and L2 pseudorange for Glonass. King 180507
c                   POLRED, UT1RED: initialize jd1 to 0 (probably not necessary). King 180507
c                   TIMCON: Change error-checking limit for GPS week from 2000 t0o 5000 and fix the
c                     error-message format.  Fang/King 180513/180517
c                   SEL_OBTYP: Fix another bug in setting Glonass pseudorange observable. King 180522
c                   SEL_OBTYP: Fix bug in setting GPS pseudorange observable. King 180531 
c             11.26 SEL_OBTYP: Add L2X and L7X to the Beidou frequency list. King 180716       
c                   EPHRED: Fix bad read of the nbody file. King 180716
c                   EPHTRP: Fix indexing problem that might cause a memory fault. Gegout/King 1808016
c                   LINEAR: New more general version (used currently for antmod interpolations)
c                       in MODEL.  King 170720
c             11.27 ASSIGN_SRPNAMES: Add ECOM1 and ECOM2 to check on valid names. King 181108
c             11.28 GPSBLOCK: Temporary kluge to make Block IIIA the same as Block IIF for 
c                     Earth radiation.   Herring/King 190129
c                   SEL_OBTYP: Remove debug print.  King 190226 
c             11.29 Increase maximum line length to 1024 and associated maximum
c                   observables allowed (also "MAXOBT" in gamit/includes/makex.h)
c                   to 63, to accommodate long records in RINEX 3 files. Floyd 190404
c             11.30 ASSIGN_SRPNAMES, READ_GFILE: Update checks for currently supported radiation-pressure models. King 190425 
c                   READDF, READM1, READM2, WRITM1, WRITM2, dimpar.h: Change all 'maxnet' to 'maxsit'  
c                     since multi-session no longer supported.  King 190502 
c                   READ_GFILE: Allow 19 orbital parameters. King 190524 
c                   CLKERA: Fix bug in extra read of 'clkoff'.  Gegout/King 190529
c             11.31 THDRED: Assign names for SRPs 16-19 from data statement of comt(3). King 190618
c                   SVNAV_READ: Cleaned up code and logic for reading igs_metadata.snx file linked to svnav.data,
c                     Added antpwr to return arguments.  Set 0 when svnav.dat read. Herring 190702.
c                   READ_SVSINEX: Finished routine and added to lib. Removed version imbedded in svnav_read.  Herring 19071.
c                   THDRED/READ_GFILE/: Added the missing arguments in the svnav_read call.  Herring 190701.
c                   THDRED: Fixed format for satellite ics output, fixed Orbital parameter types output Herring 190703.
c                   RRXHED: Skip reading of SV-dependent wavelength factors since not supported. King 190815
c                   SATATT: Change dimensions of 'yatt' from 'nsat' to 'maxsat' King 190826
c             11.32 New subdirectory sofa containing IAU SAFA package. Compiled lib called libsofa in lib directory McClusky 190801.
c                   CRS_TRS, ROTSNP_SOFA: New routines to use IAU SOFA McClusky 190801.
c                   ROTSNP modified to use new SOFA routines for inertial referecne frame McClusky 190801.
c                   TIME_TT, TIME_TT_SEC new routines for use with SOFA library McClusky 190801.
c                   GPT3_L new routine for reading GPT3 grid McClusky 190801.
c                   READE: Set bclock(4) for SBAS.  Wang/King  191202 
c                   RSP3HD: Rework the coding for reading an sp3d header in support > 102 SVs.  King 200202 
c             11.33 THDRD, READC1, READC2, WRITC1: WRITC2. Updated for L1/L2 Satellite PCO values in c-file.  
c                      Increased version to 10.71; svantdx(2,3,nsat) Herring 200126.
c                   READC5, WRITC5, WRITM1, RSTNFO, WSTNFO to addd antdaz antenna azimuth and update version Herring 200205
c                   READ_ANTEX, READ_ANTEX_HEAD: Removed fatal on number of frequencies; made header running more
c                     flexible. Herring 200209.
c                   SEL_OBTYP: Higher piority added for L6X, L6A, L6I for Beidou.  Herring/King 200223
c                   READ_SVSINEX: Fixed logic and error reporting when some entries (e.g. TX_POWER) are missing. Herring 200228.
c                   SVNAV_READ: Made subroutine name consistent in report_stat calls.  Herring 200228.
c             11.34 RSTNFO: Changed length of comment field in stinfo lines.  Cleaned up (illegal) use of LEN and passed strings. 
c                   Corresponding length mods made to kf/htoglb/mstinf.[fh] Herring 200317
c             11.35 GET_IUT1POL: New routine to get iut1pol value from sestbl. Called rather than arbitrary set values in 
c                     different routines Herring 200205
c                   ROTSNP_SOFA, SROTAT: Implemented new Desai and Sibois diurnal/semidiurnal ocean tide EOP model. TAH 200505
c                   HFEOP_DESAI: New module to implement the Desai and Sibois model. TAH 200505
c                  :MAKEFILE.GENERIC: Added hfeop_desai.f90 and get_iut1pol.f: routines. TAH 200505
c                   ROT_GFILE: Replaced setting iut1pol to zero with value read from sestbl. TAH 200505.
c             11.36 SEL_OBTYP: Removed 11.33 update and replaced with user selected lower frequency 
c                     for G5, C7, E6, E7 and F8. TAH 200511.
c                   READ_ANTEX: Added extra atxfrq to have primary and secondary F2 choice.  Increased
c                     length to 4-characters to allow selection to be marked.
c             11.37 RSTNFO: Updated to allow a gap between end and start times of up 6hrs.  This happens 
c                     with firmware updates where the gap can be 1-minute to a few hours.  There is probably
c                     no data in gap but GAMIT is still generating epochs duringthe gap (makex and maybe 
c                     model fatals). Cleaned up indenting around modified code.TAH 200524.
c                   SEL_OBTYP: Fixed initialization loop to 6 (from 4) for iobtypx; Changed secondary choice 
c                     to P2 for Glonass when C2 in header.  TAH 200526.
c                   READE: Trap and fix bad years in GLONASS entries TAH 200605.
c                   READ_SVSINEX: Updated fatal message when GLONASS frequency not found to report the 
c                     SV and PRN. TAH 200608. 
c                   RSP3HD: Updated 32I to 50I to allow for 35 Beidou satellites TAH 200618
c             11.38 MAKEXP: Assign sp3 name when running manually. Floyd 20200928 (originally 190121)
c             11.39 RRXHED: Fix warning for RINEX version (allow up to 3.04). Floyd 20201016
c             11.40 RSP3HD: Remove comma before output variable(s) on lines 172 and 193. Floyd 20201130
c             11.41 EPHRED: Remove comma before output variable(s) on lines 151, 155, 160, 164, 452, 467, 482, 497, 509 and 522;
c                   RSP3HD: Added SP3 versions a (or blank) and b to reading of header lines (3-7 and 8-12);
c                   XHDRED: Remove comma before output variable(s) on line 203. Floyd 20201201
c             11.42 GPT2: Removed "gpt_filnam" variable to accommodate moving the file opening
c                         and initial reading to the calling routine. King/Floyd 20201216
c                   GPT3: Added subroutine for GPT3 model from a combination of TU Wien routines gpt3_5.f90
c                         and gpt3_1.f90, with saasthyd.f and asknewet.f appended. King 20201216/Floyd 20201218
c             11.43 TIMCON: Updated 2020 to 2100 as end date. TAH 210101.
c                   RRINEX: Updated 2020 to 2100 as end date. TAH 210101
c             11.44 RRXHED: Updated warning for RINEX version (allow up to 3.05). Floyd 20210415
c             11.45 WRITM1: Updated version to 1072 for the change in slot numbers to accommodate 45 Beidou satellites TAH 210701.
c             11.46 SEL_OBTYP: Added L2Y phase and C2Y pseudorange observable for GPS. Floyd 20210708
c             11.47 SEL_OBTYP: Added L5/C5, L5X/C5X and L5P/C5P for BeiDou B2a. Floyd 20210730
c
      return                                       
      end
