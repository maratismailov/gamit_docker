Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.

      SUBROUTINE MVERSN(VERSN)

      implicit none

      integer nblen

      character*10 getmac,machin
      character*40 vers
      character*45 libver
      character*109 versn

c     get makex version
      machin = getmac(1)
      write (vers,5) machin(1:nblen(machin))
    5 format ('10.29 2020/12/16 03:13 UTC (',a,')')

c     get library version
      call lversn(libver)
      versn = ' '
      versn= vers//' Library ver. '//libver
C
C Version 1.14  Put MAKEX in SCCS and backdate history.

C Version 1.15  Get the version string correct - kurt 871103

C Version 1.16  Fix time conversion routines, and test for gaps - kurt 871108

C Version 1.17  Install super-polished time_converter - kurt 871118

C Version 1.18  Examine quality vectors to catch cycle slips - kurt 880111

C Version 1.19  Read only least significant 2 byes of quality vectors - kurt 880114

C Version 1.20  E-file format is now FICA block 9.  Use new TAIUTC with hardwired
c                 //bullen/data/gps/tables/leap.sec - kurt 880118

C Version 1.21  Read the next day's data file if needed to fill scenario - kurt 880217

C Version 1.22  Fix minor bugs in going over weeks - kurt 880221

C Version 1.23  Now write only error flags 0 or 5.  Handle depleted FICA files from
C                 NGS2FIC, according to "frequency plan" - durt 880410

C Version 1.24  Do not write sats no in scenario to K or E-files - kurt 880412

C Version 1.25  Do not write bogus epochs to K and E files.  "Bogus" includes epochs outside
C                 the scenario, too - kurt 880413

C Version 1.26  Do not allow the week number to go backwards or forwards by more than
C                 one week - kurt 880414

C Version 6.1   New GAMIT release 6.1 -- mhm 880417

C Version 6.2   Count the number of records, not just good ones.  Fix gap counting.
c                 Hard wire batch mode.  Clean up logic about good data.  Change name of
C                 makex.ins to makex.fti, and include it in more subroutines - kurt 880418

C Version 6.3   Keep sats in order specified in th scenario file - kurt 880423

C Version 6.4   Finally get the time-tag-checker working - kurt 880426

C Version 6.6   Time tags beyond the end of the scenario are no longer considered 'bad'.
C                 Make sure that eh ephemeris is included in the scenario before calculating
C                 a clock offset for the K-file - kurt 880427

C Version 6.7   Incorporate time tag improvements, and TI ROM data - kurt 880725

C Version 6.8   Allow turning off writing a clock file - kurt 880727

C Version 6.10  Modify DOFICA to read the ROM phase data correctly at end of the sampling
c                 interval, as does GESAR - kurt 880803

C Version 6.13  Handle new, weird CORE v 4.8, which samples differently.
C                 Print out the time tags of the first 5 records - kurt 880808

C Version 6.14  Handle weid time tags of CORE v 4.1, flagged as 4.11 - kurt 880825

C Version 6.15  When opening second file, adjust week numbr accordingly - kurt 880906


C Version 6.16  Read new FICA blocks, remove frequency plan bias by counting
C                 seconds from start of scenario, rather than start of week - kurt 881114

C Version 6.17  Changes to accommodate MACROMETER and variable-length GESAR FICA blocks.
C                 Split DOFICA into DOFICA, BLK6, BLK9, BLK401, BLK70, and BLK670 - king 881006
C               Increase channel dimensions in MAKEX.FTI from 4 to 8.
C               Change BATCH file to include software and version number, in MAKEX,
C                  WXHEAD, RBATCH.  Allow gaps in FICA file.   - king 881007
C               Changes to RFICA and WTIME - king 881103
C
c version 6.18  is missing due to damaged p.mversn.ftn
c
c
c Version 6.19  Read FICA block 80.  These blocks are the same as the block 70s
c                 developed at MIT, but never released.  Clynch at ARL gave us
c                 approval to use blkid=80.  MAKEX will continute to read the
c                 block 70s for compatibility, but of course this is dangerous!
c                 A warning will be issued.
c                 Changes to DOFICA, add routine blk80. - kurt 880218

c Version 6.20  Remove start time code from MAKEX to new routine SETTIM.
c               Remove declaration of ephemeris quantities from MAKEX (not used).
c               Rationalize use of variables RCVRSW, SWVER, and SWVERU.
c               Remove reading of header information from MAKEX to new routine RHEAD,
c                 and add call to new routine BLK1001 for MACROMETER II.
c                 Change MACROMETER II data routine from BLK670 to BLK1080,
c                 but keep old routine for now.  Changes to MAKEX, DOFICA, BLK101,
c                 BLK670; new routines SETTIM, RHEAD, BLK1001, BLK1080.
c               Increase MAXOBS to 900 in MAKEX.FTI to handle 30s data.
c               -    king 890220
c
c version 6.21  Allow a time tag with bad seconds but assumed good week number.
c                 This is a less conservative approach to evaluating data quality.
c                 It may fail if the week number rolls over within one file.
c                 kurt 890223
c
c Version 6.22  Increase epoch dimensions to 1100, satellite dimensions to 12.
c                 Changes to dimpar.fti and MAKEX.FTI
c               Add PRN 14 (=NS14) and PRN10 (=NS13) to PRN2NS - King 890322
c
c Version 6.23  Add provisional MINI-MAC blocks.  Set UTC/GPS time flag
c                  in calls to TIMCON.  Use lowercase comparisons. Changes
c                  to MAKEX, RSCENO, BLK1080, BLK1180, BLK670, DOFICA,
c                  RCOORD, RSITED, SETTIM, STNCRD, WCLOCK, WTIME.
c                  kurt 890505 (comments and delta's by king)
c               Implement new X-file format and new library calls.  Changes to
c                 MAKEX and WXHEAD; new routines SETTYP and CHGOBS; GETDAT added
c                 to bind file.
c               Remove write of uninitialized trakid to E-file
c                 in WORBIT.
c               Rationalize satellite dimensions and loop indices.  Change RCVRCH
c                 in MAKEX.FTI and dimensions to MAXCHN. Introduce new variable
c                 NCHAN into MAKEX to denote channels in the actual receiver.
c                 Remove SATS_IN_SKY from MAKEX.FTI and dimensions; use MAXSAT.
c                 Changes to MAKEX, RSCENO, RPHASR, PRN2NS, ISGOOD, INSCEN, BLK6,
c                 BLK1180, BLK401, BLK70, BLK670, BLK80, BLK1080, DOFICA.
c                 king 890727
c Version 6.24  Thursday, August 3, 1989   6:26:48 pm (EDT) kurt
c         minor modifications for dealing with minimac data
c         RPHASR: fix call to ferror
c         CHGOBS: Do not modify data here.  It has already been done.
c         OPENF,RCOORD, RBATCH:  Call standard ferror
c         SETTIM: Minimac code is now MIN.
c Version 6.24 Monday, August 7, 1989   8:19:00 pm (EDT) kurt
c         Minimac compatability.
c         STNCLK:  Sat clock correction already applied, so SVDT = 0
c         BLK1180: FICA block 1180s are defined to have Doppler phase,
c                  which is positive when Doppler phase is increasing
c                  (or when range is decreasing).
c         WXHEAD: Put in a space after 'on'.
c Version 6.25 Monday, August 7, 1989   8:45:44 pm (EDT) kurt
c        Minor clean up.
c
c Version 6.27 Sunday, August 13, 1989   7:33:38 pm (EDT) kurt
c        Handle GESAR version 1.6 properly.
c        changes to settim only.
c Version 6.28 didn't work
c Version 6.29 Tuesday, August 15, 1989   12:16:59 pm (EDT) kurt
c        Handle ROM pseudoranges correctly.
c        The pseudorange should be propagated to the same
c        epoch as the phase is observed.  This is 59.08 seconds. GPST
c        BLK401: check for cycle slips, propagate pseudorange.
c        STNCLK, CHGOBS and SETTIM: add ROM flag
c        Also look at 0.001 seconds for Minimac data.
c Version 6.31 Friday, August 18, 1989   8:45:31 pm (EDT) kurt
c        Fix BLK401 to correctly handle the fact that the ionosphere affects
c        the phase and group delay rates differently.
c
c        Create MAKEK from MAKEX August 29, 1989.  king/kurt
c
c version 6.32 Sunday, October 29, 1989   5:17:03 pm (EDT) kurt
c        Minor repairs to main, stnclk and worbit.
c
c         Created MAKEJ.    December 22, 1989.
c
c Version 7.1  Remove PHASER format and replace by RINEX.
c              Remove pseudorange file and replace by SVCLK file
c              Write time tags on X-file.
c              Add Trimble as receiver type.
c              Add receiver name to header record.
c              Replace WGS72 by WGS84 for geodetic coordinates in clock calculation.
c              Write a J-file from FICA input if requested.
c              Save sub-frame 1 (and some sub-frame 2) stuff in BLK9 and DOFICA
c              Changes to MAKEX, RHEAD, SETTYP, WXHEAD, DOFICA, BLK9;
c                 New routines MAPAMP, WSVCLK, RRINEX, RXERR, RRXHED;
c                 Remove CHGOBS and RPHASR.
c version 7.3 Yehuda Bock:
c              In MAKEX write out ephemeris to E-file regardless of pseudorange obs
c              YB 13 April '90
c              In MAKEX change option to "UNKNOWN" instead of "NEW" for
c              X,E & K files
c              In SETTIM, define NCHAN (=4) for ROM receiver
c              In ISGOOD add TMODE to parameter list (modify MAKEX
c               and FIC2RX in call to ISGOOD), FIC2RX not used, neither
c               is SETAMP which uses ISGOOD as well
c              Remove CHGOBS from Makefile
c              YB 16 April '90
c version 7.4 Kurt Feigl
c              Remove TMODE from variables passed to ISGOOD.
c              Reasons:
c              1) TMODE is only defined for TI receivers.
c              2) It is not in RINEX files.
c              3) The only test on it was: (tmode(j) .ge. 1 .or. tmode(j) .le. 10)
c                 which does not exclude anything.
c              The old routine is preserved as ISGOOD.OLD
c              Purge SETAMP and CHGOBS, per Yehuda's suggestion.
c version 7.5 Kurt Feigl
c              Use RINEX Signal Strength Indicators (0-9 scale) for amplitudes.
c version 7.9 Kurt Feigl
c              Incorporate Yehuda's changes for Right-handed coordinate convention.
c              WXHEAD:  change coordinates to have 'N','S','E','W' flags
c              RCOORD:  change coordinates to have 'N','S','E','W' flags
c                       backwards compatability with ' ' and '-' preserved
c                       improve Yehuda's code to handle southern hemisphere, too.
c              STNCRD:  change to Right-handed system, remove underscores.
c              MAKEK:   calls to XHDRED and GETEPH changed
c              XHDRED:  Eliminate in favor of library routine of same name,
c                       but with longer call list, now includes LATFLG and LONFLG
c              GET_EPHEM: change name to GETEPH.  Now calls library routine BREAD,
c                        to read E-file or RINEX N-file. change call to include subfr1
c version 7.11 Kurt Feigl
c              MAKEX:  RINEX ISSI 4 is now good, rather than low amplitude.
c version 7.12 Kurt Feigl
c              MAKEX    eliminate ighdwr (hardware) flag.
c              ISGOOD   examine RINEX ISSI flag to determine goodness.
c version 7.13 Peter Morgan
c              SETTIM   Minimac special version 1.50 for GOTEX Tsukuba data
c version 7.14 Kurt Feigl
c              MAKEX    RINEX ISSI = 3 is good, ISSI = 2 is low amplitude.
c version 7.15 Bob King  Saturday, 16 June 1990
c              MAKEX, DOFICA, BLK9, WRINAV - increase dimension of subfr to 8
c                  storing the HOW word for time of ephemeris record
c      MAKEK:  modify BREAD to use same units (meters) as MAKEX/BLK9
c              Change orbit file to an input rdorb file- MAKEJ, MAKEK, MAKEX.FTI,
c              CLOSEM, GETEPH -- rwk 900706
c              Remove special code for Mini-Macs: MAKEK, STNCLK
c              **NB: Mini-Mac X-file should be translated to RINEX standard
c              Fix bug for C/A PRs in GETPR -- rwk 900716
c version 7.18 Make Bob's changes - Kurt
c              Read and write RINEX-type E-files. Change all files references
c                for FORBIT etc. to FWRORB, and FEXTR2 to FRDORB.
c              Change header descriptions to be more specific re input RCVRSW
c              Add check of pseudorange reasonability to function ISGOOD.
c              Add printout summary of rejected (unreasonable) data points.
c       MAKEK: Change name of BREAD to READE in GETEPH. -- rwk 900731
c version 7.19 Bob King  - Tuesday, 31 July 1990
c              Remove ability to write E- and J-files.  Assume they're there.
c              Correct Mini-Mac pseudorange to RINEX standard (add back SVDT).
c              Change name of BREAD to READE in GETEPH.
c              Changes to MAKEX, RBATCH, SETTIM

c version 7.20 Bob King - Thursday 11 October 1990
c              Add RCVRSW, SWVER to X-file header (to identify Mini-Macs)
c               Changes to MAKEX and WXHEAD
c              Yehuda Bock - Tuesday, 11 September 1990
c              Modify SETTIM to deal with Trimble v4.10
c              Modify MAKEX and RBATCH to read session number from input
c               makex.batch file, e.g.,
c               (A4,A2,A1,1X,A3,2X,A3,1X,F4.2)
c               0010901.127  TRM 4.1
c              Yehuda Bock - Friday, 28 September 1990
c              Reorder some code so that two character year is output
c               to screen (MAKEX)
c              Change RSCENO & MAKEX to read session number from
c               sceno file
c              Check Trimble v4.11 (not v4.10) in SETTIM
c              Yehuda Bock - Sunday, 30 September 1990
c              Add 4.10 version for Trimble receivers (SETTIM)
c              Yehuda Bock - Wednesday, 3 October 1990
c              Modified SETTIM to accept Ashtech data
c               for the time being assumed a bogus v1.00
c              Yehuda Bock - Monday, 8 October
c              Modified RSITED & MAKEX to accept session number from
c               makex.sited file
c              Modify RSITED & RBATCH for two character year
c              Output bug in RSCENO
c              Yehuda Bock - Wednesday, 10 October
c              If ISSI = 1, set flag to low amplitude instead of deleting. -kurt 901023
c version 8.2  Release 8, handle Trimbles.  King and Feigl 901112
c              MAKEX: K-file interval = 120 s for Trimbles.
c              Accept issi=1,2 for Trimble (Temporary patch).
c
c version 8.3  Add FICA Blks 1201 and 1280 for CIGNET Trimble
c              Changes to RHEAD, DOFICA
c              New routines BLK1101, BLK1201, BLK1280
c              Bob King - 19 Dec 90
c              Add FICA Blks 1301 and 1380 for CIGNET Rogue
c              Add Rogue software versions
c              Changes to RHEAD, DOFICA, SETTYP, SETTIM
c              New routines BLK1301, BLK1380
c              Bob King - 29 Dec 90
c version 8.5  Add call to library LVERSN in MAKEX  rwk 91/01/19

c version 8.6  Add TRM version 4.30 kurt 91/03/21
c              MAKEJ Option for Selective Availability; new JREADC.  kurt 91/3/21

c version 8.7  MAKEJ : In open statement for J-file, use status='unknown'
c        -8.8    in all cases, & remove extraneous code
c              SETTIM : add Rogues and Mini-Mac v1.64
c              SETTIM : add Trimble 3.12,3.22,3.28,4.30
c              SETTIM : Assume that all Mini-Mac RINEX data passed
c               through Werner's conversions & have correct time tags
c              SETTIM : for GESAR, set NCHAN = 4 (was set to 6??)
c              ISGOOD, MAKEX : add new error flag (igl2bd=8), comment out
c               in MAKEX for now (had a bug in MODEL which must be sorted out)
c              STNCLK : fix print statement
c              SETTYP : reorder C1 ahead of P1 check for RINEX, this is for
c               Rogue DSN data which has five data types (C1 and P1), still
c               need to make changes to properly read DSN data
c              MAKEX : modify to J-file need not be read
c              MAKEX : add epochok check for k-file, to avoid those
c               annoying duplicate records in k-files
c              MAKEFILE : add errflg.fti
c              YB 4/2/91
c              More general fix for data types when input RINEX NOBTYP=5
c                (Replaces Scripps coding in SETTYP  rwk 91/4/9
c              Add Rogue swver 5.61 and 1.11 in SETTIM   rwk 91/04/09
c              Increase MAXEPH to 1300 in GETEPH - rwk 91/04/09
c        8.8b  SETTYP: count types of data
c
c version 8.10 Round TI ROM time tags to even second if produced by
c                translator TISTRX 10/90; Call this ROM 1.11.
c                Changes to MAKEX, RHEAD, SETTIM - rwk 91/5/9
c              Remove entry point nextsession to make compatible with fic2rx
c                Changes to MAKEX, MAKEK,OPENF, RSITED, RSCENO,RCOORD,
c                WXHEAD, FIC2RX rwk 91/5/10
c              Correct test on receiver type for K-file interval.
c                MAKEX   kurt 91/05/10
c              MAKEK: Add 'opnerr' (='nexts') to call of RCOORD -- rwk 910514
c              Change RRXHED in library to allow lower-case keywords
c                king 91/05/24
c version 8.12 Add COR version 4.13 to allow for both 59.00 and 59.08 at
c                Onsala (i.e. set slop =0.1)  king 91/05/31
c version 8.13 SETTIM: Add new Rogue software versions to message. yb/rwk
c              MAKEX:  Make 120s K-files for Ashtech as well as Trimble  yb
c version 8.14 RSITED: change comments to reflect F8.4 format for heights
c              WXHEAD: add extra decimal point for antenna height
c                bock/stark 91/07/10
c              Relinked for change in lib routine RRXHED
c                yb 91/07/13
c version 8.15 RCOORD,WXHEAD: Add extra decimal place for radius (0.1 mm)
c                yb 91/08/25
c              MAKEX: Allow low amplitude (1 or 2) for Ashtech data,
c                     zero amplitude is rejected for all data
c                yb 91/08/30
c              SETTIM: Add Ashtech v. 2.00 (P-code receiver)
c                yb 91/09/03
c              MAKEX: Patch for JPL RINEX screw up (no SNR for phases in
c                RINEX files)
c              MAKEX,ISGOOD check for pseudorange (omit for codeless receivers where clock
c                           drift can cause pseudoranges to appear unreasonable
c                           (i.e. Trimbles and Ashtechs)
c                yb 91/09/07
c              ISGOOD: check for zero pseudoranges
c                yb 91/10/07

c              Kinematic option  (Scripps version 9.1)
c              Added modkin.fti & modified MAKEX to generalize X-file epoch format
c                yb & kfg 91/10/07
c              For backwards compatibility modify X-file header, WXHEAD
c                yb 91/10/18
c              MAKEX: change kses to character*1 for station.info file
c                     (compatible with sceno file)  yb 91/11/20
c              MAKEK : reads station positions from lfile at each epoch
c              GETPR : changed to read expanded x-file format  jfg -- 911121
c              MAKEK : replace ucoord by ukinl
c              MAKEK : do nout output to k-file if kflag=2 (roving receiver)
c                     yb -- 911123
c
c              commented out logic to check on time span for kinematic stations.
c              needs to be removed or fixed jfg - 1911208

c              MAKEJ:  Write prns file to allow checing against scenario file.
c                      SIO Dec 91

c              included multiple antenna type/height option for static sites
c              in kinematic version (change in makex.ftn) jg 1/20/92
c
c              Fixed switched FFLAGS and FDFLAG in Sun Makefile  donnellan/king 1/22/92

c              Add receiver sampling time to X-file if it is greater than the
c                 requested interval.  Changes to MAKEX, WXHEAD, RHEAD, BLK101
c                 BLK1001, BLK1101, BLK1201, BLK1301.   King  92/1/23.

c              SETTIM: Add software versons 4.4, 4.42, 4.53 for Trimble.  king 91/1/23
c              Allow command-line argument for the batch file: MAKEX.  Oral/King 92/1/23
c              Pass the RINEX loss-of-lock flag to the X-file and vica versa
c                  Changes to MAKEX, RRINEX(lib)     Oral/King  92/1/23
c              Allow command-line argument for COUNTX.  Oral/King 92/1/23
c              Set nepoch = maxepc if too many epochs, in RSCENO.  King  92/1/23

c              MAKEX: Change file opens from 'readonly' to 'old', and 'append' to
c                   new to match FTN77 standard for the HP.   King 92/1/29

c version 8.16 Changes to use station.info and session.info files   King 92/2/10
c              RBATCH: Comments and change from GESAR default.
c              RBATCH, MAKEX:  Remove unused nepochs from call.   King 92/2/17
c              Add expanded data format to X-file header; changes in MAKEX
c                    and WXHEAD.     King 92/2/13
c              SETTIM: Add Trimble 4.60 (hypothetical) and 4.62. Bock/King 92/2/17
c              Set experiment type from station.info; put kinematic code in
c                    subroutine.  Changes to MAKEX, RBATCH, RSCENO; new routines
c                    KINCRD and library RSTNFO.   King 92/2/17
c              RSCENO: Read session as integer (but backward compatible).  king 92/2/23
c              SETTIM and MAKEX:  Allow for Mini-Mac sampling at 59.001 in German
c                    RINEX files, either tagged correctly (swver 1.59) or incorrectly
c                    (swver 1.89).    King  92/2/24
c              GETEPH: Add mprns, iwknes to SAVE list.   Donnellan/King 92/2/24
c              RXSCAN: Change to make a bar-graph.  Feigl 92/2/24
c              WXHEAD, MAKEX: Change antenna offsets to R*8.   King 92/2/25
c              MAKEX: Fix bug in opening second RINEX file.  Oral 92/2/25
c              Always use LREAD for site coordinates.  Remove RCOORD; changes
c                    to MAKEX, MAKEK, Makefile.    King 92/2/26

c version 8.17 MAKEX: Update Bob's kluge to get day of year corrrectly
c               when less than three characters.
c              Fix RSTNFO (in LIB) and recompile Bock 92/03/01
c              RHEAD: Update for Ashtech XII and P-12 receiver
c              MAKEX: Fix fend problem (would skip to new file even if did not exist
c              Makefile: Remove unused subroutines   Bock 92/03/02
c              MAKEK:  Add skd variable to call to XHDRED.  Bock 920306
c              SETTIM: Rogue version 6.11            Bock 92/03/09
c              RHEAD: (??Found no changes in YB version)
c                     Fix missing commas for FICA reads - King 92/3/21
c              MAKEX: Fix bug (line 790 too long) - King 92/3/21
c              RBATCH: If session number = 0, set = 1.  King 92/3/24
c
c version 8.16 Combine MVERSN, KVERSN, JVERSN. Change Makefile, MAKEJ, MAKEK. King 92/4/27
c              MAKEX, WXHEAD:  Add session number to X-file header; write receiver
c                              interval unconditionally.    King 92/4/27
c              Makefile: Remove references to FAKEK.  King 92/4/27
c              Remove unused variables:  BLK1080, BLK401, DOFICA, MAKEX, KINCRD,
c                   RBATCH, RSCENO, RXSCAN, STNCLK.  King 92/5/2
c version 9.1  Compiled at SIO by Yehuda 92/05/03
c version 9.11 RSCENO : if SN=0 => SN=1 Bock 92/05/06
c              MAKEK : Add variables to xhdred call
c                      Remove rcoord variable
c                      Bock 92/05/07

c version 9.12 RXSCAN : Fix call to RRINEX    King 92/6/4
c              MAKEK  : Remove unused variables  King 92/6/4
c              MAKEJ, JREADC  : Allow reading of multiple C-files for modeling SA
c                 New routines PICKFNS, LIFT_ARG, COUNT_ARG (put eventually in /lib)
c                 Change Makefile.   Dong/King 92/6/4
c              RSCENO : Match session numbers 0 and 1 correctly.  King 92/6/4
c              SETTIM : Add Rogue ver 7.00, Trimble 4.64, Ashtech 6.00
c                        Bock 92/06/06
c              SETTIM : Add COR ver 5.7, Trimble 3.24.  King 92/06/23
c              MAKEJ, MAKEX, GETEPH  : Add calling argument to CHECKE in library
c                      to make RINEX reads more robust.   King 92/06/25
c version 9.13 RHEAD  : Blank out header to avoid writing nulls.   King 92/07/06

c version 9.14 New MAKEJ to use multiple C-file input to handle SA modeling.
c              Changes to Makefile and MAKEJ; new routines J_FROM_E (renaming of
c                  FROME), J_FROM_C, RATE_EST, AVGCLK, ALLANV (dummy for now).
c                  Dong/King 92/7/10.
c                  **commented added 93/1/14: JREADC also updated 92/7/10  King
c              Fix time-tag bug in J_FROM_C.   King 92/7/15
c         9.15 MAKEX: Mismatched calling statements for second RHEAD call.  King 92/7/17
c         9.16 MAKEX.FTI  Increase MAXLIN to 100.
c              BLK70, BLK80, SETTYP.  Explicitly declare variables.
c              MAKEX, SETTYP, FIC2RX.  Put warning in for inconsistent TRM L2FACT.
c                     King 92/9/21
c              GETEPH, J_FROM_E : Remove isat from call to reade.  King 92/9/24
c version 9.17 Remove FAKEK from SIO version.
c              MAKEJ : Remove isat from call to reade Bock 92/9/30
c              SETTYP: Correct bug in warning for TRM L2FACT.   King 92/9/30
c              SETTIM, SETTYP:  Addition of SERCEL receiver.  Feigl/King 92/10/01
c              SETTYP: Change warning on TRM L2FACT to check existence of P2.
c                      Bock/King 92/10/05
c              MAKEJ : Replace with MIT version. No subroutines at bottom.
c                      Add j_from_c, allanv.ftn, avgclk.ftn and rate_est.ftn
c                      Bock 92/10/08
c version 9.18 SETTIM: Add information about Rogue Receivers in commment form
c              SETTIM/SETTYP: Add Turbo-Rogue (code "TRB"ROG) versions 1.00 -  2.50
c                      Bock 92/10/29
c version 9.19 MAKEJ: Add lib version to screen display.
c                     ReMake with lib/READE fixed to avoid CR problem.  King 92/11/10
c              J_FROM_C: Add stop if input and output J-file names the same.  King 92/11/27
c version 9.20 SETTIM/RHEAD/SETTYP:
c              Add ver 5.52 for Trimble SSE (first example,Caltrans resurvey). Bock 92/12/16
c              SETTIM: Add TRM versions 5.51, 5.53.   Bock/King 92/12/22
c              RHEAD : Make Trimble versions > 5.50 generic Trimble 4000 (FICA only).
c              SETTYP: Allow determining ndat from receiver software version (FICA only) only
c                      for TRM <5.5 and ASH <2.0.   King 92/12/22
c version 9.21 ISGOOD: Do not delete data if P2 observation is of zero value. Bock 92/12/29
c              SETTIM: Add TRM version 4.80 Bock 92/12/29
c version 9.22 MAKEJ, JREADC, FIC2RX : Update to 24 satellites Bock 93/01/11
c              SETTIM: Add Rogue version 7.30 & Modify Turbo-Rogue Output  Bock 93/01/17
c version 9.23 SETTYP, Makefile: Correct typo in setting L2FACT for X-to-RINEX.  Move to
c                      library for use with XTORX.  Feigl/King 93/01/25
c version 9.24 MAKEX, WXHEAD: Add antenna code to X-file header.
c              MAKEX : Return RINEX height from HISUB (in /lib).
c                         King 93/2/11
c              MAKEK : Add antcod to XHDRED call. (Remove non-library
c                        XHDRED from MIT Apollo Makefile and from
c                        MIT Apollo and Sun directories.)   King 93/02/12
c version 9.25 Makefile: New HP/Sun compatible.  Vigny/Fang/King 93/03/01
c              GETPR : Remove unused variables.  King 93/03/01
c version 9.26 SETTIM: Add WM102 receiver (WM2).  King 93/03/10
c              GETPR:  Fix typo from 93/03/01 change.  King 93/03/10
c              WXHEAD: Replace nulls with blanks at end of header.  Herring/King 93/3/11
c              SETTIM: Change WM102 channels from 5 to 8.   King 93/03/12
c              MAKEX:  Put in and remove debug for stdrel segmentation-fault problem;
c                        fix by deleting and remaking gs/makex3.a
c                      Remove debug print for site radius.   King 93/03/15
c version 9.27 RBATCH, STNCRD: Fix debug printout, add station coordinates to
c                        debug.   King 93/03/18
c              SETTIM: Remove comment about WM-102 channels.  King 93/04/02
c              MAKEX : Remove % from Apollo; fix STNCRD call in Sun/stdrel.
c                        King 93/04/02
c              MAKEX and KINCRD: Fixed code up to handle kinematic files
c                        JFG 93/06/22
c              Consolidate Apollo and SUN versions at SIO:
c               copy MAKEX&KINCRD from SUN to Apollo
c               copy SETTIM to SUN
c              Bock 930706
c version 9.28 SETTIM: Add new version number for Ashtech version 7B
c                      (changes documented in code) Bock 930709
c              KINCRD: Declare variables explicitly.    King 930714
c version 9.29 SETTIM: Add new firmware versions for TurboRogues Bock 930725
c              KINCRD: Remove unused variables.   King 930826
c version 9.30 STNCLK: Don't set error flag (but still issue warning) PR
c                      out of range (TRM clocks sometimes ok).  Herring/King 930910
c              MAKEK:  Fix calling sequence for STNCRD.   King 930924
c              WXHEAD: Add year to calling argument for NSNPRN.   King 921018
c              Increase formats for 32 satellites:  MAKEX, MAKEJ, FIC2RX,
c                        JREADC, J_FROM_C.   King 931018
c              Remove unused variables:  WXHEAD, KINCRD.   King 931018
c version 9.31 MAKEK: Fix second call for STNCRD.   King 931213
c version 9.32 RSCENO: Fix bug in reading scenario files where maxsat > number satellites.
c              SETTIM: Add new firmware version for Rogues (7.4).
c                      Add new firmware for Ashtech Z-12 (8.00).  Bock 940104
c              MAKEX: RSCENO Replaced by call to library routine RSESFO.  Bock 940104
c              Makefile:  remove RSCENO.  King 940107
c              MAKEX, STNCLK: Write clock warnings only to infor file, with a single
c                   summary warning to screen.
c version 9.33 SETTIM:  Add setup for Leica.  Morgan/King 940111
c              SETTIM:  Add Trimble version 4.81.   Feng/King 940114
c version 9.34 WXHEAD: Change call to NSNPRN--now a subroutine and more arguments.
c                   King 940119
c              MAKEX: Clean up output.  Bock 940122
c              SETTIM: Add Ashtech Z-12 version 8.01 (actualy IC01).  Bock 940125
c              Makefile: Add all dependencies explicitly.  Fang 940215
c version 9.35 SETTIM: Make all Trimble versions between 3.12 and 5.99 use the same
c                   offset (0.0) and slop (0.3s).  King 940322
c              MAKEX, RHEAD: Clean up printout and add test for mismatched receiver
c                   software.   King 940322
c              SETTIM: Add swver 7.10 for Ashtech P-code 7A.  Ferhat/King 940412
c version 9.36 MAKEX: Call KINCRD to read station.info and call HISUB if a change of
c                   antenna offset is detected in the RINEX file.    King 940502
c              MAKEK: Fix bug in debug.   King 940503
c              Makefile: Change HP compiler switch from +E1 +e to +U77.  Fang/King 940503
c              GETEPH: User smarter algorithm from BCTOT to find ephemerides.  King 940505
c version 9.37 MAKEX, KINCRD: Add calling argument to for library RSTNFO.  King 940511
c              MAKEX: Update antenna heights for static observations.   King 940517
c              SETTIM: Enter Rogue SNR-8 version 7.6 - 7.8.  Bock 940602
c              SETTIM: Add Ashtech swver 8.01 (1C01).  Bock 940623
c              Makefile: Shift ranlib to execute only once. Herring/King 940624
c              SETTIM: Add Ashtech 8.10 (1D00 firmware).  Bock 940702
c version 9.38 MAKEK: Use UTC/GPS time flag for UTC-GPST conversion.  King 940722.
c              MAKEX, WXHEAD: Write X-file in GPST, not UTC.  King 940728
c              MVERSN: Blank version to avoid non-ascii characters.  King 940728
c              SETTIM: Increase Trimble 'slop' to 1 s to account for poorly tuned
c                   oscillators and long tracks.  Herring/King 940805
c version 9.39 RHEAD: Add Ashtech (TOPCON) and Rogue receiver names.  Bock 940810
c              MAKEX,SETTIM: Add TOPCON receivers (equivalent to Ashtech).  Bock 940810
c              MAKEX: increase format dimension for ireject variable.  Bock 940810
c              MAKEJ: Blank header to avoid writing nulls.  King 940815
c              SETTIM: Enter Turbo-Rogue version 3.0.  Bock 941012
c              SETTIM: Enter Trimble version 6.04.  Bock 941104
c              GETEPH: Fix declaration of satok and initialization of first_eph. King 941117
c              MAKEX: Change x to 1x  and '' to ' ' in formats to satisfy IBM compiler.
c                       King 941118/941205
c              SETTIM: Ashtech version 1E00 (named 8.20).  Bock 941215
c version 9.40 Changes to compile with XL Fortran on RISC/6000:
c                 DOFICA, MAKEX, OPENF, RBATCH, TAGCHK King 950103
c version 9.41 WXHEAD: Replace call to NSNPRN by call to SVNAV_READ.  King 950403
c              WXHEAD: Add hrs/min to SVNAV_READ call.  King 950405
c              WXHEAD: Correct month,day variables passed to svnav_read PT 950407
c              COUNT_ARG, GETEPH: Declare variables explicitly.  King 950511
c version 9.42 JREADC: New C-file format.   King 950520
c              SETTIM: Enter Turbo-Rogue 3.20.  Bock 950608
c version 9.43 JREADC: Add antenna reference point offsets to C-file.  King 950717
c version 9.44 MAKEX: Remove unexecuted statements at end.  Sanli/King 950719
c              J_FROM_E: Add SV to debug to calling argument (mis-match with call
c                  in MAKEJ.  Sanli/King 950719
c              JREADC:  Change variable names itide/iut1pol=>ietide/isptide
c                  to agree with rest of GAMIT.  King 950720
c              SETTIM: Enter Ashtech firmware 8.30 & 8.70.  Bock 081195
c version 9.45 ISGOOD: Don't use data when L1 and L2 are identical (TurboRogue malady).
c                  Check on pseudorange since phases may have a bias.  Herring/King 950818
c version 9.46 RHEAD, RXSCAN: Pass 20-character receiver and antenna serial numbers
c                  correctly from RINEX to X-file header.   King 951019
c              MAKEK: Add receiver and antenna ids to lib/xhdred call.   King 951019
c version 9.46 MAKEX: Changed call to hisub. McClusky 951019
c              KINCRD: Changed call to hisub. McClusky 951019
c              MAKEX.FTI: Added antmod.dat filename,unit number,status to commons. McClusky 951019
c version 9.47 SETTIM:  Add Trimble firmware 6.10.  Fang/King 951125
c              SETTIM:  Allow through 6.29 .  King 951128/951207
c              Makefile: Variable ranlib for compatibility with Solaris 2. Fang/King 951208
c version 9.48 MAKEJ: Change message with J-file choice.   King 951215
c              SETTIM: Allow all TurboRogue versions 1.0-3.20; set nchan=12 for 3.20. King 951222
c version 9.49 MAKEX: For rejected data, printout PRN, not current-epoch index.
c                 King 960105/960122.
c              J_FROM_E: Print correct stop time when j-file entries are not time-ordered. King 960122
c              MVERSN: Remove redundant machine type containing nulls.  King 960122
c              MAKEJ, WXHEAD: Expand 'header' to get full version info.  King 960122
c version 9.50 Modified all routine in the makex directory with the report_status routine:
c              avgclk, closem, dofica, geteph, getpr, j_from_c, j_from_e, jreadc, kincrd
c              makej, makek, makex, mversn, openf, pickfns, rate_est, rbatch, rhead,
c              rxscan, settim, stnclk, stncrd, tagchk, wclock and, wxhead
c              GETEPH, STNCLK: replaced call to get_prog_name with call to rcpar, inline
c              with the new changes to report_stat. McClusky 960227
c version 9.51 FIC2RX: Make implicit none and remove unused variable.  King 960517
c              J_FROM_C: Fix bug in getting starting epoch (assumed C-files were UTC).  King 960517 
c version 9.52 Changes for new kf/gamit structure and Makefiles from unimake.  King 960718
c                   All routines:  removed trailing blanks and replaced lib/*.fti 
c                   includes by ../includes/*.h.  
c              Remove obsolete FAKEK from directory and Makefile.generic.  King 960719  
c              SETTIM: Allow Trimble firmware through 7.15.  Bennett/King 960802 
c              J_FROM_C, J_FROM_E: remove unused variables.  King 960813    
c version 9.53 RCOORD, RSITED, WCLOCK: Replace functions 'lowera' and 'uppera' with calls to 
c                 subroutine 'lowers' and 'uppers' to avoid segmentation fault when calling 
c                 argument not character*256.  King 960829
c version 9.54 SETTIM: Set UTC time for Rogue v2.2.  McClusky/King 961002
c              SETTIM: Add Ashtech firmware 1G00 = v.8.40.  King 961003
c              SETTIM: Add comma to format statement.  Tregoning/King 961017
c version 9.55 GETEPH, J_FROM_E:  Change calling arguments for (lib) READE for BCCCHECK updates.  King 961107
c              SETTIM: Add Ashtech firmware 1E24 (8.22), 1E95 (8.25).  Bock 961122
c              MAKEX: Change 'no obs' message from status to warning; write also
c                  in the infor file.   King 961203
c              MAKEX, KINCRD:  Change 'rcvers' => 'swvers' to match rstnfo.f.   King 970111
c version 9.56 G77 compiler changes:  DOFICA, GETEPH, INSCEN, J_FROM_E, JREADC, RXSCAN, STNCRD
c                  Initialised previously uninitialised variables, removed unused variables
c                  and explicitly typed all unexplicitly typed variables. McClusky 970113 
c version 9.57 Makefile, /includes/makex.h, GET_RXFILES, MAKEX:  New routine to pre-scan a set of input 
c                  directories and find all files with data within the requested span.  King 970211/14
c                  GET_FICFILES also added but currently dummy. 
c              GETEPH, OPENF: Minor change to print statements.  King 970212/970213 
c              MAKEK, MAKEX: Add messages for failed lib/getdir.  King 970212  
c              SETTIM: Remove debug print; fix check on swver for Trimble.  King 970213     
c              MAKEX, RHEAD: Write screen output only thru MAKEX.status/warning/fatal except 
c                   when in debug mode.  King 970213
c              STNCLK: Increase tolerance on obs. range to 40,000 km.  King 970213
c              MAKEX: Add antenna offset update for new RINEX file; remove redundant
c                   offset and coordinate variables (use all from modkin.h). King 970218 
c              OPENF: Fix calls to report_stat.  King 970218    
c              RBATCH: Change iflag 1 call to entry point for argument clarity.  King 970219
c              RHEAD: Don't print header info to screen unless debug mode.  King 970219 
c              GET_RXFILES: Initialize file names to blank.   King 970220
c              STNCLK: Loosen observed range tolerance to 10,000-40,000 km.  King 970225
c              GET_RXFILES: Fix trap of ../   King 970303   
c              MAKEX: Fix reading of multiple antenna heights.  King 970303
c              GET_RXFILES: Fix multiple declaration problem caught by g77. McClusky 970303
c version 9.58 MAKEX: Initialised uninitialised variables newfiles and ant_event. McClusky 970305
c              MAKEK: Initialised uninitialised variable nfiles2. McClusky 970305 
c              GET_RXFILES: Remove commented initialization and unsed variable.  King 970306/7
c              GET_FICFILES: Change from dummy to functional routine. King 970306
c              DOFICA:  Fix bug in calling BLK101.  King 970307
c              Makefile, FIC2RX, RFICAHD:  New routines to translate FICA to RINEX.  King 970311
c              DOFICA, FIC2RX, GET_FICFILES :  Add data-found logical to DOFICA calling sequence;
c                  change error handling to use report_stat.   King 970311       
c              GET_RXFILES: Fixed dimensioning of variables illi issi. McClusky 970311 
c              BLK101:  Print to screen and log file on on debug.  King 970317  
c              DOFICA:  Initialize svid array to zero.  King 970318  
c              GET_FICFILES: Add close of scanned files. 
c              Makefile, XTORX.  Move from /utils; subtract initial integer.  King 970321
c              GET_RXFILES: Remove duplicate declaration.  McClusky 970325 
c              MAKEX, SETTIM: Add parameters for Leica; echo slop.  King 970326  
c              ISGOOD:  Add Leica to receivers for which to skip pseudorange tolerance 
c                  test (clocks drift too much).  King 970326
c              XTORX: Fix double declatation if uantmod. McClusky 970327
c              FIC2RX: Fix uninitialised variable nchn. McClusky 970327
c              RXFICAHD: Fix multible declarations of several variables, and illegal 
c                  initialisation of logical nerr. McClusky 970327
c              STNCLK: Print debug for orbit info only the first time through unless
c                  a statement uncommented.  King 970430  
c              MAKEX: Fix debug print of kjd_stop.  King 970430  
c version 9.59 GET_RXFILES:  Allow searching for one week prior to start.  King 970502
c              STNCLK:  Fix report_stat call.  Murray/King 970630   
c              SETTIM:  Update for Trimble firmware 7.22.  King 970702     
c              SETTIM:  Update for Leica firmware 3.6.  King 970717 
c              RHEAD: Fix format in report_stat call for non-interger ircint.  King 970717
c              SETTIM: Shorten warning for unknown firmware version.  King 970718 
c version 9.60 SETTIM: Add Ashtech firmware 1F50 (8.35), 1F60 (8.36). Bock/Genrich 970728
c              GET_RXFILES: Correct format in report_stat message.  King 970809/970815
c              MAKEX:  Call lib/rrinex with non-dummy epoch.   King 970817  
c              XTORX: Fix bug in start time for rstnfo call.  King 970817 
c              RFICAHD:  Stop looking for header if data block encountered.  King 970829 
c              FIC2RX, RFICAHD: Set receiver type from BLK 0 if header block missing.  King 970831
c              BLK401: Write better diagnostics to .status and .infor files.   King 970831 
c              FIC2RX: Code cleanly for command-line or interactive input.  King 970902
c              SETTIM:  Add documentation for Trimble firmware versions.  King 970902 
c              FIC2RX: Correct bug in passing epoch of clock to lib/wrxnav.  King 970903
c              BLK9:  Minor correction to documentation.   King 970903    
c              GET_FICFILES, MAKEX:  Fix naming of FICA files, and error message.  King 970904 
c              FIC2RX: Set wavelength factors for TI = 1 1 if no BLK 101.  King 970905  
c              FIC2RX: Trap bad anth in internal read.  King 970917  
c version 9.61 GET_FICFILES:  Add checks and fatal message if no FICA files of the right
c                  name form exist.  King 971024/971029        
c              FIC2RX:  Remove redundant commas in format.  King 971030
c              MAKEX:  Don't open X- or K-files unless data exist (avoid writing null files). King 971215
c              SETTIM: Add Ashtech v.8.23 and consolidate Z-12 'if's.   Bock/King 971229 
c              SETTIM: Add Leica v.3.7.  King 980102  
c              SETTIM: Add Leica v. 3.3, 3.4, 3.6, 4.1, 5.2; add Ashtech 8.24, 8.32.  King 980109  
c              MAKEX: Add total number observations to .infor summary.   King 980114  
c              RFICAHD:  Remove redundant comma.  King 980115 
c version 9.62 DOFICA, GETEPH, MAKEK, MAKEX, RXSCAN:  Write times in GPST not UTC, using modified
c                  lib/wtime.f.  McClusky/King 980123 
c version 9.63 MAKEX:  Don't write data more frequently than ircint indicates it should be there.
c                  (To avoid problem with autcln putting bias flags in that solve/cview ignore in
c                  cases where we've set the RINEX sampling to, e.g. 120s but are writing 30s X-files.
c                  This occurs only when we have a RINEX file with mixed sampling.)  King 980129
c              SETTIM:  Fix bug for detecting Ashtech firmware > 8.0.  Morgan/King 980204
c              GETPR, MAKEK:  Fix call to report_stat.  King 980227
c              FIC2RX:  Remove extraneous comma in call to readj.   King 980227
c version 9.64 MAKEX, RBATCH: Remove batch-file unit number from calling sequence since
c                  makex.h is included; having the unit both places (but with different names)
c                  is a relic that seems to cause problems on some compilers.   King 980309
c              J_FROM_C:  Comment out redundant open of C-file (perhaps eventually keep
c                  this one and remove the one in jreadc).   King 980316 
c              SETTIM:  Add Trimble firmware 7.26.   King 980413
c              GET_RXFILES: Fix array-bound problem when using maximum number of files.  King 980520
c              MAKEX:  Add echo of # epochs, inverval, session length to STATUS report.  King 980520
c              RHEAD: Initialize to blank rcvr name from batch input; add Leica.  King 980521
c version 9.65 WXHEAD: Add SV antennas offsets to svanv_read call.  King 980615 
c version 9.66 MAKEX, FIC2RX, J_FROM_C, J_FROM_E, MAKEK, MAKEX:  Replace calls to gamit/lib sb pickfn.f 
c                with libraries/comlib function pickfn.c.  King 980720/22 
c              XTORX: Add code to use old-style X-files without time line; fix bug in UTC->GPST 
c                conversion; don't warning about non-std PRs unless a MiniMac.  King 980724/5
c              XTORX: Fix warning message about mismatched rcvr firmware.  King 980806/980807
c              RFICHD: Remove extraneous commas in formats (SGI complained).  Morgan/King 980901
c              SETTIM: Add Z-12 firmware version CB00->9.0.  Bock/King 980902 
c              FIC2RX, J_FROM_C, J_FROM_E, MAKEK, MAKEX, XTORX: Fix problem with pickfn calls
c                 using new C routine.  Herring/King 980905
c              JREADC:  Add 'norb' to calling arguments for lib/readc2.  King 980911
c              OPENF: Add L-file to list of fatal opens.  King 980918
c              JREADC: Add 'svclock' to calling arguments for lib/readc5. King 981002
c              RFICAHD: Fix call to report_stat.  King 990121  
c version 9.67 MAKEK: Allow command-line arguments.  King 990122
c              J_FROM_E: Allow RINEX as well as E-file name in wildcard.  Fang/King 990208
c              J_FROM_E: Blank out E-file name before setting.  King 990217
c              MAKEJ, J_FROM_E: Allow command-line arguments.  King 990219
c              J_FROM_E:  Add space to wildcard argument.  Fang/King 990305
c              MAKEJ: Correct format for 4-digit GPS week.  Fang/King 990308  
c              SCANX: Do not use comma as continuation character (IBM and f90 restriction).  Fang/King 990324
c version 9.68 SETTIM: Add single and dual-freq Geotracer (GEO) and UNAVCO CMC Allstar12
c                receiver (CMC).  King 990410/990415/990505     
c              SETTIM: Add Z-12 firmware version CC00->9.1.  Bock 990418
c              GET_RXFILES: Fix bug in search and add search directory to report_stat warning message.  King 990413
c              SETTIM:  Add Trimble firmware 7.27.   King 990422    
c              XTORX: Use double-read kluge to avoid a problem reading the X-file 
c                character from stdin after a call to pickfn.   Fang/King 990525
c              MAKEX,ISGOOD: Modify pseudorange check for modern Ashtech & Trimble receivers. Bock 990419
c              XTORX, SCANX:  Fix assignment of valid obs in RINEX header.  M King/ R King  990526
c              XTORX: Assign D-file unit (fix dangerous but passable bug).  King 990527
c              SETTIM: Add Trimble firmware 7.28.   King 990602
c              SETTIM: Add artificial firmware 1.96 for MiniMac sampled at GPST+5. Shen/King 990603
c              ISGOOD: Remove unused 'swver'.  King 990604
c version 9.69 WXHEAD: Add warning about no sampling interval on RINEX file.  MKing/RKing 990610
c              XTORX: Fix format to handle > 12 satellites.  van Domselaar 990708 
c version 9.70 Remove unused RSCENO, RSITED; change Makefile.  King 990727
c              KINCRD; Fix uninitialized epoch variables.  King 990727      
c              RBATCH: Change batch-file format.  King 990727
c              Changes for 4-digit years:  BLK1080, BLK670, J_FROM_C, J_FROM_E, JREADC,
c                MAKEJ, MAKEX, RBATCH, RDFICA, SCANX, WRHEAD, WRINAV, XTORX.   King 990728
c              MAKEX:  Add argument to lib/rsesfo call (for MAKEXP).  King 990817  
c              SETTIM: Add TurboRogue 3.30 firmware and minor changes.  Bock 990810
c              SETTIM: Add Ashtech CD00 firmware (9.20).  Bock 990818
c              SETTIM: Completely revise LEI section; Add Ashtech firmware UB00.  Bock 990821
c              SETTIM: Fix Y2K bug in code for Rogue version 2.30; fix selection of 
c                 firmware for unknown Rogue.  King 990830 
c              MAKEJ: Fix assignment of source nav file in header.  King 990917
c              RBATCH: Fix Y2K bug with old-style batch files.  King 990924
c              SETTIM: Add a mess of Ashtech firmware codes for SOPAC DB.  Bock 991026
c              MAKEX: fix for kinematic mode processing.  Bock 991023
c              MAKEX: merge kinematic and dynamic options for single-epoch processing.  Bock 991120
c              SETTIM: Remove / from 'message' format for Ashtech/Rogue report_stat calls.  King 991203
c version 9.71 GET_RXFILES:  Fix Y2K bug.  Fang 000101
c              GET_RXFILES: Allow ././ssssddds.yyo form.  Fang 000215
c              FIC2RX: Allow two possible FICA filename styles; fix y2k bug.  King 000301
c              BLK101: Fix declaration of 'auser'.  King 000301   
c              BLK1101, FIC2RX, RHEAD, RFICAHD: Fix passing of MiniMac values for FICA.  King 000302
c              FIC2RX: Fix observable types for MiniMac.  King 000302
c              SCANX, XTORX: Fix y2k problem with old x-files.  King 000307/000309
c              MAKEX: Fix bug in station.info read when RINEX file skips epoch and goes past
c                 the end of the scenario.  King 000318
c version 9.72 RHEAD: Fix bug causing array overrun with more than maxlin -7 comments lines 
c                 in the RINEX file.  McClusky/King 000608
c              SETTIM:  Add firmware 1.1-1.9 for Trimble 4700 total station.  King 000627
c              MAKEK:  Fix bug passing string to getdir causing segmentation fault on Linux.  Herring 000719  
c              STNCLK: Fix evaluation of broadcast ephemeris for clock terms (reduces error
c                  from about 0.5 microseconds to less than 0.1 microseconds).  Herring 00071 
c              MAKEJ: Change messages from 'e-file' to 'navigation file'.  King 000725
c version 9.73 SETTIM: Add firmware 2.x for Trimble 7400MSI.  Mattvd/King 001027    
c              SETTIM: Add firmware 4.05 for Geotracer 3220.  Wangqh/King 001106
c              BLK1101: Remove unused variable.  King 001204
c              REHEAD: Remove extra comma in BLK1101 call.  King 001208
c version 9.74 SETTIM: Add Javad JPS with open firmware.  King 001228
c              J_FROM_E, PICKFN: Fix or remove uninitialized variables to satisfy HP compiler.  King 010103
c              SETTIM: Add Ashtech Z-18 with firmware V6400 assigned as 9.8.  King 010109
c              SETTIM: Assign # channels when Ashtech firmware not matched.  King 010131
c version 9.75 FIC2RX, GET_RXFILES, RFICAHD, RHEAD, RXSCAN, XTORX: Changes for RINEX version 2.10:  version and interval 
c                      now real, extra decimal place for time of first observation.  King 010301
c              ISGOOD:  Remove test on 15,000-30,000 km pseudorange limit.   King 010302
c              GET_RXFILES, MAKEX, RXSCAN: Add SV-type id to RRINEX call; discard non-GPS SVs;
c                 change dimensions of issi/illi to (,2).  King 010313/010322
c              J_FROM_E:  Removed initialization start/end saved values at top (arguments
c                 no intialized themselves).  King 010314  
c              MAKEK:  Fix writing of site name to status file. King 010622  
c              MAKEJ:  Remove blank lines in screen print near end.  King 010622  
c              RFICAHD: Add code to detect ASDAP (TI-ROM) with missing header.  King 010705  
c              FIC2RX, BLK70, BLK80: Clean up print output. King 010706  
c              MAKEJ, MAKEK, GETEPH, STNCLK: Fix missing site name in status echo; don't write old screen messages.  King 010717
c              FIC2RX: Use rcpar rather than getarg to avoid problem with HP UX-11.  Herring/King 010718   
c              J_FROM_E: Comment out initial assignment of time_end_save from time_end.  King 010719 
c version 9.76 GET_RXFILES: Increase format for allowed number of searched files from 10 to 20. Vigny/King 010727  
c              STNCLK: Trap and omit clock offset values gt 1 s, to avoid k-file problems with TI 4100.  King 010727
c              MAKEK:  Add k-file name to message about bad clocks.  King 010727  
c              MAKEX: Remove extraneous parenthesis.  Morgan/King 010731    
c              RHEAD: Rephrase logical if's to avoid compiler complaint.  Morgan/King 010731
c              J_FROM_E: Fix names of undefined variables in debug print.  Morgan/King 010731  
c              FIC2RX: Initialize 'idflag'.  Morgan/King 010731  
c              RFICAHD: Add code to detect UNPACK (GESAR) with missing header.  King 010731     
c version 9.77 SETTIM:  Fix setting of 1 ms offset in search with MIN 1.49.   King 010815   
c              MAKEK:   Add check for pseudorange = 0 (already in MAKEX) to catch bad TI data.  King 010829
c              FIC2RX:  Write more warnings into the print file.  King 011005 
c version 9.78 MAKEX, KINCRD, XTORX,  Allow reading new-style station.info files
c                  Also change ../includes/modkin.h sname to character*16 to match other code.  King 020308 
c              MAKEX, KINKCRD, XTORX: Change 'swvers' from R*8 to R*4 to match rest of GAMIT.  King 020312 
c              SETTIM:  Add comments re Trimble BD-750 (MSAG) 'board' receiver; no change in code.  King 020320
c              GET_RXFILES: Do RINEX file sort by time, not name, to accommodate hourly file.  King 020321
c version 9.79 GET_RXFILES: Fix bugs in rwk changes of 020321.  Murray 020417   
c              SETTIM: Add Trimble 5700 (assuming firmware 1.03-1.04).  Jamason/King 020503.
c              SETTIM: Set nchan = 12 for firmware 1.0-2.99.  King 020531  
c version 9.80 MAKEX, MAKEK, FIC2RX: Allow reading of either l-file or apr file through modified 
c                     /lib/lread.f and  /includes/modkin.h.  King 020807
c version 9.81 SETTIM:  Add firmware code 9.81 for Ashtech Z-XII3 and 9.90 for Z-XII3T receivers.  King 020918
c              FIC2RX, XTORX: Change name of library routine geoxyz to geohms_xyz to avoid conflict with
c                  new, more general geoxyz.  King 021003
c              SETTIM: Add comments to better explain the Ashtech firmware versions.  King 021008 
c version 9.82 SETTIM: Add swver 9.90 for Ashtech Z-FX firmware ZC00.  King 021204  
c              XTORX: Fix bug in opening station.info.  King 030108 
c              MAKEX:  Fix bug in antenna update with old-style station.info.  King 030113 
c              SETTIM: Add comments for the Trimble MS750.  King 030127   
c              SETTIM: Separate Ashtech-derivative Topcons and add the Topcon Hiper firmware. King 030226 
c version 9.83 SETTIM: Make all Topcon except the Hiper derivative of other manufacturers (changes also to MAKEXP). King 030318 
c              XTORX: Use day-of-yr from x-file name rather than start time in naming RINEX file; issue warning. King 030320
c              SETTIM: Add Novatel Millenium (NOV 4.45-4.52) receiver card.  King 030331 
c version 9.84 MAKEX, KINCRD, XTORX: Remove obsolete 'icall' from call to hisub.  King 030418     
c              MAKEX: Remove debug for antmod.  King 030430.
c version 9.85 Makefile, FIC2RX, KINCRD, MAKEK, MAKEX, STNCLK, WXHEAD: Changes for reading velocity-dependent
c                coordinates from L-file.  Remove STNCRD.   King 030510 
c              Remove unused WRHEAD.  King 030510
c              MAKEX, MAKEK, WCLOCK:  Write to the k-file the observed PR rather than the theoretical
c                range (now match long-standing documentation). King 030520
c version 9.87 GET_RXFILES: Remove option to use header start,stop times since some RINEX files have
c                these wrong.  Herring/King 030820
c              SETTIM: Fix check for codeless Trimble firmware and P2 present (4700 has firmware 1.xx);
c                add comments for changed and new Ashtech firmware. King 030916
c              J_FROM_C: Fix problem with quoted argument for pickfn.  Champollion (Montpelier)/King 030925
c              MVERSN, FIC2RX, MAKEJ, MAKEX, WXHEAD, XTORX  : Increase length of version variable (ifc compiler).  King 031026
c              SETTIM: Fix line > 72 characters.  King 031026  
c              FIC2RX, XTORX: Fix mismatched character lengths. King 031026  
c              MAKEX:  Change 'back-door' debug option to use 'debug' rather than 'test' for the makex batch file name. King 031114
c              MAKEX:  Add NavCom receivers (assume swver 'NAV' and '1.8'.  King 031218
c              MAKEX, GET_RXFILES: Add command-line arguments to control number of days backwards and forwards
c                 searched for RINEX data.  McClusky/King 031230 
c              GET_RXFILES:  Fix bug in 031230 mod.  McClusky/King 040107
c              SETTIM: Add Z-12 firmware 1SWE1D0 as 8.86 (documentation only). King 040116
c version 9.88 SETTIM: Change default for search window ('slop') from 1 ms to 100 ms.  King 040214
c              SETTIM: Add Sokkia receivers.  King 040319
c              SETTIM: Add 4000MSGR firmware 7.58.  King 040412
c              J_FROM_E: Fix undefined start time in summary printout.  King 040614
c version 9.89 MAKEX: Correct nominal time tag by PR for checking within epoch window.  King 040624
c version 9.90 MAKEX: Fix bug in correcting nominal time tag (k-file wasn't written).  King 040719
c              SETTIM: Add 0.70 as firmware for the Trimble NETRS receiver.  King 040803
c              SETTIM: Add Leica GRX1200 (firmware CC00) as swver 7.0  King 040823
c version 9.91 MAKEX, RHEAD, WXHEAD, new SET_DCB_FLAG, Makefile: Deterine if P2 and C1 corrections needed 
c                 for cross-correlating receivers and write a new flag on the X-file.   King 041012
c              SETTIM: Change 'slop' for unknown Leica firmware from .01 to .1.  King 041031
c              SETTIM: Comments for NOVEU4 and LC1230.  King 041109
c              MAKEX: Fix bug in calling argument to rhead with second RINEX file.  King 041113
c              SET_DCB_FLAG: Add receiver type to report_stat message.  King 041203
c              WXHEAD: Trap case of RINEX files already altered by CC2NONCC.  King 041203
c              WXHEAD: MAKEX, SET_DCB_FLAG: Override DCB flag if CC2NONCC previously run. King 041208
c              WXHEAD: Add receiver code to X-file header.  King 041208
c              WXHEAD: Fix bug in testing for CC2NONCC (may not matter).  King 041209
c              SET_DCB_FLAG: Fix bug in calling arguments introduced 041209.  King 041216
c version 9.92 JREADC:  Add model parameters to calling arguments for lib/readc2. King 050129
c              MAKEK, XTORX: Remove 'extra' from call to lib/xhdred.  King 050129 
c              GET_RXFILES: Remove nulls from file names returns by getdir.  Herring/King 050131
c              SETTIM: Allow 'JNS' as well as 'JPS' for Javad receivers.  King 050303
c              SETTIM: Add comment on Topcon GB-1000 (TPSGB1) and allow firmware 2.4.  King 050311
c              SETTIM: Add Septentrio (SEP) receivers.  King 050520
c version 9.93 MAKEX, KINCRD, XTORX:  Add radome to read_rcvant and hisub calls.  King 050618 
c         9.94 MAKEX, KINCRD, XtORX:  Add warnings variable to hisub calls.  King 050719
c              KINCRD: Minor compile error caught by HP.  Herring 050826    
c version 9.95 MAKEX, WXHEAD, XTORX , Makefile: Changes to use lib/hisub.f w/ hi.dat (retain hisub.f option
c                for now); change antenna info written to X-file from phase centers to ARP.
c                Remove kincrd.f from Makefile and directory (unused since May 2003).  King 050927   
c version 9.96 MAKEX: If more channels on the RINEX file than indicated by rcvr/swver in settim, 
c                print a warning and increase 'nchan' to avoid losing data.  King 060307
c              SETTIM: Trap unknown Leica firmware 2.01-2.49 and set like 6.0.  King 060307
c              GET_RXFILES: Fix bug causing incorrect sorting and skipping a file in the case
c                  where the RINEX file has a blank between the header and data records.  King 060511
c version 9.97 MAKEX: Move screen messages about command-lines so that they appear only if the
c                  command-lines are missing. King 060706
c              MAKEX: Change arguments for rsesfo (see lib/lversn.f).  King 060815  
c              MAKEX, MAKEK, FIC2RX, Makefile: Change modkin.h to model.h.   King 060823   
c              MAKEX: variables dlat and dlon were real, should be int. McClusky 060829
c              Makefile, COUNT_ARG, LIFT_ARG:  Remove count_arg.f & lift_arg  (in /lib). King 061108
c version 9.98 MAKEX: Change arguments for rsesfo (see /lib/lversn.f). King 061207
c              MAKEX, CLOSEM: Correct typo in open for antmod.dat, and add close.  King 070103.  
c              FIC2RX, MAKEX, MAKEK: Change name of model.h variable from ukinl to iul.  King 070118
c              SETTIM: Add Thales Z-Max receiver with firmware 9.96 (=> Mxxx).  King 070118
c              SETTIM: Add TPS ODYSSEY_E receiver with firmware 2.5.  King 070129
c              SETTIM: Add Leica GMX902 (LCG902) with firmware 8.03.  King 070307
c              SETTIM: Add comment that Trimble R7 uses 2.x.  King 070315
c              MVERSN: Remove screen unit number from 'lversn' call.  King 070416
c              MAKEX, GET_RXFILES, RXSCAN: Remove unneeded variables from call to /lib/rrinex. King 070416
c              MAKEX, RXSCAN: remove unneed variable from call to /lib/gps_only.  King 070416
c              GET_RXFILES: Remove debug.  King 070416
c              XTORX: Fix bugs in using new station.info format.  King 070605
c              SETTIM: Add TPS GD3 72-channel receiver, firmware 3.1.   King 070607
c              SETTIM: Add artificial swver (2.99) for L1-only TR4600.  King 070607
c              SETTIM: Fix typo for early NETRS firmware (0.07 --> 0.70) (no effect). King 070705
c              ISGOOD: Restore check on PR reasonableness but only for gross cases.  Herring/King 070710
c              SETTIM: Allow firmware 2.6 for TPS receivers.  King 070817
c              FIC2RX, J_FROM_C, MAKEX: Remove unused 'iscrn' from call to /lib/readj. King 070910
c              J_FROM_E, JREADC, MAKEK, MAKEX: Remove unused statement labels.  King 070910 
c              BLK1001:  Add dummy statement for unused calling arguments.  King 070910   
c              ISGOOD, MAKEK, MAKEX, STNCLK: Remove unused 'rcvrsw' and 'swver' from calling argument. King 070910  
c              MAKEX, OPENF: Remove unused 'qbatch' from calling argument.  King 070910
c              MAKEX, MAKEK, MAKEJ, XTORX, MVERSN: Remove unused 'iscrn' from calling argument.  King 070910  
c              SCANX, XTORX: Remove unused 'ischan' from calling argument. King 070910 
c              GETPR, MAKEX, TAGCHK: Remove unused calling arguments.  King 070910
c              RXSCAN: Fix format type mismatch, identified by gfortran.  King 071112
c              SETTIM: Allow 14 SVs for Novatel receivers.  King 071126    
c              BLK101: Fix undeclared variable. Herring 080207
c              SETTIM: Add comment re JPS EGGDT firmeware. 080318
c              SET_DCB_FLAG: Fix logic for AS-off case and non-Trimble receiver.  King 080423
c              MAKEX:  Fix formats in debug print for kstart,kstop.  King 080501
c         9.99 MAKEX, CLOSUM, FIX2RX, : Remove references to old-sytle hisub and station.info.  Herring/King/McClusky 080509  
c              MAKEK, XTORX: Change antenna offset arguments for /lib/xhdred.  King 080509
c              RREAD: If wavenlength factor missing from RINEX header, set = 1.  Herring/McClusky/King 080522
c              MAKEX: Fix bugs in defining hi.dat and wrong case for site name.  King 080522
c              MAKEX: Fix longstanding bug initializing with maxchn instead of maxsat.  King 080523
c              SETTIM: Update #channels for Trimble SSE and SSi to 12.  King 080523
c              MAKEX: Fix bug causing array overrun; add trap for nepoch > maxepc.  King 080924    
c              SETTIM: Add trap for nchan > maxchn.  King 081006
c              SETTIM: Add comments on two Septentrio receivers and firmware.  King 081204  
c        10.00 XTORX; Change 'readd' to 'readdf' to match change in /lib.  King 090416 
c              MAKEX, SET_DCB_FLAG: Add site code to report_stat message.  McClusky/King 090508
c              MAKEX, GET_RXFILES; Add 'debug' to calling arguments.  King 090514   
c        10.01 SETTIM: Add firmware for more Novatel receivers. King 090630
c              SETTIM: Add comments on firmware assignment for Leica GRX1200PRO and GRX1200GGPRO. King 090703
c              MAKEX, MAKEK, XTORX: Remove 'trkcod' from calls to /lib/rsntfo; Remove call to check_old_stnfo 
c                (no longer needed); fix variables for new model.h  King 100209 
c              SETTIM: Allow 14 channels for the Leica GRX1200GGPRO firmware > 7.20.  King 100330
c        10.02 JREADC; Add c-file variables for both the Feb 2005 and Aug 2010 format changes. 100827 
c              MAKEX, FIC2RX, XTORX: Add rcvers, rcvrsn, antsn to rstnfo call sequence; remove
c                declarations conflicting with new these new variables in model.h. King 100908
c        10.03 SETTIM: Change mapping of Trimble firmware (treat all > 3.72 the same). King 101129
c              SETTIM: Increase Sokia channels to 24 (may still be low) to suppress warnings for the
c                 SOK GSR2700 RSX.   King 110930 
c              MAKEX; Allow 1000 rather than 100 rejected data points for message; clean up some debug. King 111019
c              GET_RXFILES; Fix bug in GAMIT status message when too many RINEX files found.  King 120112
c        10.04 MAKEX, MAKEK, FIC2RX: Change some variable names to match the new model.h.  Moore/King 120731
c              MAKEX: Another minor change for new model.h.   King 120802
c        10.05 JREADC: Add two variables to call for readc2.  King 130116
c        10.06 SETTIM: Add Ashtech firmware versions 9.97-10.01. King 130605
c              SETTIM: Add comment for SEPT POLARX4 vers 2.3.4.  King 130703
c              MAKEX: Add trap for # channels > maxchn.  King 140409
c              SETTIM: Add TRM4-3 to list of firmware 4.xx (comment only).  King 140612
c              SETTIM: Add Hemisphere (HEM) receivers. King 140811
c        10.10 MAKEX, FIC2RX, GETEPH JREADC, MAKEJ, MAKEK, OPENF, RBATCH, RXSCAN, WXHEAD, XTORX: Changes 
c                for non-GPS GNSS. King 141004
c              SETTIM: Add CHC receiverss. King 150422
c        10.11 WXHEAD: Add start/stop times to svnav_read call. King 150520
c              MAKEX:  Change itimdif to integer*8 kING 150528 
c              MAKEX:  Fix bug causing fatal when no data found in the RINEX file.  Wang/King 150608
c        10.12 MAKEJ, J_FROM_SP3 (temporarily dummy), READJ, Makefile.:  Allow getting SV clock terms f
c               from a fit to the sp3 clocks (not yet implemented).  King 150520
c              MAKEX, MAKEK, MAKEJ, CLOSEM, GETEPH,J_FROM_C, J_FROM_E, J_FROM_SP3. RBATCH: Change name of
c               broadcast file in makex.h from 'rdorb' to 'nav' and add 'sp3;; implement j-file from SP3, 
c               with nav-file as backup if SV missing or bad fit.  King 150730
c        10.13 MAKEJ, CLOSEM, CLKFIT(new), GETEPH, J_FROM_C, J_FROM_E, J_FROM_SP3, MAKEK, 
c              MAKEX, RBATCH, Makefile includes/makex.h: Implement fitting of SP3 file for J-file 
c               offset and rate (untested with breaks).  King 150804 
c        10.14: SETTIM: Add Leica firmware (artificial 9.0 + ) King 151009
c               MAKEX, GET_RXFILES, RHEAD, RXSCAN, SETTYP(move from /lib and edit), XTORX, Makefile:  
c                 Remove all FICA  code; implement reading of RINEX 3. 151014
c               MAKEX, SETTIM: No longer use receiver-specific 'nchan', just check for 'nprn' > MAXCHN.  King 151015
c        10.15  CLKFIT, J_FROM_SP3: Fix bug for all SVs following one that has no SP3 clocks. King 151019 
c        10.16  Makefile, remove all FICA routines to /fica and all references to FICA in /makex
c               (CLOSEM, MAKEJ, MAKEK, ISGOOD, OPENF, WSVCLK). King 151023
c               GET_RXFILES, MAKEX, RXSCAN, SETTYP: Make dimension of iobtypx 6 instead of 4 to allow backup
c                 L1 and L2 pseudorange observations if the primary missing (see /lib/lversn.f) King 151028
c        10.17  MAKEX: Fix error in calling arguments to 'rrhead' when more than one RINEX file
c                 for a site in the session.  King 151110                                      
c        10.18  MAKEX, MAKEK, GETEPH, STNCLK: Option for all GNSS navigation messages.  King 151120
c               RXSCAN: Correct calling arguments for /lib/sel_obtyp.  King 151130 
c               GETEPH: Change to using the last nav-file epoch rather than the closest one. King 151209
c               Makefile, MAKEX, MAKEK, GETSP3 (new), GETNAV (renamed from GETEPH), ../includes/makex.h.: 
c                 Allow using SP3 files for orbits to get the station clock.  King 151212
c               Makefile, MAKEJ, J_FROM_NAV (renamed from j_from_e), J_FROM_SP3, GET_NAVCLK: Code GNSS for 
c                 sp3 j-files.  King 151215
c               WXHEAD: Fix format in writing SV name to x-file header. King 151228
c               J_FROM_SP3: Fix logic for selecting SVs. King 160107
c               MAKEX: Write observation types to the infor file.  King 160109
c               GET_NAVCLK, J_FROM_SP3: Don't add the leap second to the Glonass offsets.  King 160113
c               MAKEX, SETTYP, XTORX, WXHEAD: Revert to a (1 2 3 4) scheme for data types on the x-file, 
c                 but add the C*3 observable codes (no change yet to xtorx to use these). King 160113  
c               GET_RXFILES: Add debug print to identify the rrxhed, rrinex debug. King 160113
c               GETNAV, J_FROM_NAV: Test for bad record from nav-file.  King 160218
c               SETTYP: Fix firmware warning for Trimble codeless (4.0 only, not 4.xx). King 160219
c               RBATCH, includes/makex.h: Increase file names to 120 characters. King 160223 
c               MAKEX: Change message, now a warning, on a RINEX read to say explictly 'RINEX',
c                   not 'input data file'.  King 160302
c               RXSCAN: Add gnss to calling argument. King 160303
c               XTORX: Add gnss arguments for lib/wrxhed.  King 160323 
c               WXHEAD: Change 'HEIGHT' to "RADIUS' and add "SPHERICAL" to coordinate descriptor for 
c                  clarification; write the rcvr code and firmware.  King 160420
c               GETSP3, J_FROM_SP3: Remove sp3gnss array from calling arguement for lib/rsp3hd. King 160602
c               SETTYP: Allow L9 and C9 (IRNSS) observables.  King 160808
c               GETNAV: Fix bug in error message; allow IRNSS.  King 160808
c               STNCLK: Allow geostationary altitudes in range check. King 160817
c               STNCLK: Add site coords to clock debug.  King 160822
c        10.19  GETNAV, MAKEX, MAKEK, WCLOCK: Fix bug in extracting nav-file orbits and clocks;
c                  change the k-file format to include a header.  King 160831
c               STNCLK: Flag to not use clock estimate if the pseudorange is near zero. Herring/King 160831  
c        10.20  GETNAV, MAKEX: Fix epoch used for clock and orbit values. King 160912 
c        10.21  MAKEX: Add missing argument to initial call of getsp3.  King 170411
c        10.22  MAKEK, XTORX: Add observable type for call to lib/xhdred. King 170525
c               MAKEX: Add writing of the GNSS to the screen and infor file. King 170531 
c               GETNAV: Increase number of ephemeris records allowed to 10000 to accommodate Galileo;
c                 add PRN to Kepler eqn non-convergence message. King 170601
c               WXHEAD: Change SATELLITE line to clarify PRN. King 170601
c               MAKEJ: Remove temporary force of gnss to 'G'.  King 170601
c               RXSCAN: Fix bug causing nsat > maxsat.  King 171212
c        10.23  GETSP3: Fix ordering of indices for clock and position variables. Herring 180309
c               RBATCH: Add trap for numsit > maxsit, King 180314 
c               MAKEK, MAKEX: Fix calling arguments for getnav and stnclk.  King 180319   
c               MAKEK, MAKEX: Add new common includes/units.h to provide unit numbers
c                    removed from includes/model.h. King 180320
c               MAKEK, XTORX: Change 'skd' to 'gnss' in call to lib/xhdred. King 180322 
c               MAKEX, XTORX: Remove obsolete 'kinflg' from call to lib/rstnfo. King 180322
c               GETNAV: Change Glonass error message.  King 180326 
c               RHEAD: If RINEX 2 multi-gnss, need to change 'P2' to 'x2' to get 'C2 for Glonass.
c                    Hercring/King 180326
c               STNCLK: Add missing comma in format statement. Pickle/King 180413 
c               RHEAD, SETTYP: For GLonass, change 'P2' to 'x2' only if C2 is present.  King 180521 
c               MAKEX, SETTYP, WXHEAD: Avoid unprintable characters in writing the observables 
c                 on the x-file header; add some debug. King 180531 
c        10.24  GET_NAVCLK:  Add 'I' to allowable gnss codes. King 180716 
c               OPENF: Change 'warning' to 'fatal' for missing j-file is specified as needed. King 181009
c               RBATCH: Trap 'numsit' > 'maxsit' for the new-format case.  King 181127
c               RBATCH: Correct trap of 'numsit' for both old- and new formats.  King 181129
c               J_FROM_SP3: Allow sp3 version 'd'.  King 190219/Floyd 190226
c               GET_NAV: Fix code to include the program name in error messages. King 180226
c        10.25  MAKEX, XTORX: Change obsolete 'maxcfl' to 'maxsit'.  King 190524
c               WXHEAD: Added antpwr to svnav_read call.  Herring 190702
c               GETSP3: Added test for sp3d format.  Herring 190703
c        10.26  JREADC: Updated for L1/L2 Satellite PCO values in c-file svantdx(2,3,nsat) Herring 200126.
c               MAKEX, XTORX: Added antdaz ground antenna azimut Herring 200205.
c               RBATCH: Added check for blank station name to stop input of sites.  Old behavior restoted Herring 200302.
c        10.27  MAKEX, MAKEK: Fix misnamed variable causing MAKEK to fail; fix 
c                 inconsistent dimensions in both routines, with no effect; fix printing
c                 of the observables to the infor file, also with no effect. King 200504 
c        10.28  MAKEX: Option added to pass lower frequency for G C E GNSS systems (-lfreq option.
c                 in sh_gamit/sh_makex/sh_preproc).  TAH 200511.
c                 Fixed format reporting frequencies. TAH 20526.
c               SETTYP: Allowed L8 for Galileo E8 measurements TAH 200511.
c               RBATCH: Added back feature to have blank line act as EOF in batch file. TAH 200526.
c               RHEAD: For Glonass processing: Removed code that set P2 to x2 when C2 available.  
c                 Made P2 be secondary choice to C2 so that files with C2 in header but no C2 
c                 data, will include the P2 data in the x-file. TAH 200526.
c               MAKEX: Updated 32I to 50I to allow for 35 Beidou satellites TAH 200618
c               J_FROM_SP3: Made no SP3 clock a warning and wrote zero for clock TAH 201003.
c        10.29  WXHEAD: Remove comma before output string on line 141. Floyd 20201201
c               MAKEX.H: Increased maxchn to 99 (RX3 files with 82 channels found) TAH 201215.

      RETURN
      END


