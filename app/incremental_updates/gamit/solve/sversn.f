      Subroutine SVERSN(VERS)

      implicit none   
 
      integer*4 nblen              

      CHARACTER*40 VERS
      character*45 libver
      CHARACTER*10 GETMAC,MACHIN
      character*256 message

      MACHIN = GETMAC(1)
      WRITE (VERS,5) MACHIN(1:nblen(machin))
    5 format ('10.57 2021/06/19 13:45 UTC (',a,')')
c   5 format ('%I% of %E% %U% (',a,')')
c     The above line contains magic SCCS keywords which will be
c     automatically updated by delta and get.
c     Please do not change them.

c     get library version
      CALL LVERSN(libver)
      WRITE (message,10) vers,libver
10    format('Started SOLVE ver. ',a40,' Library ver. ',a45)
      call report_stat('STATUS','SOLVE','sversn',' ',message,0)  

CVERSION V01.0001   SENT TO AERO SERVICE - APRIL 2, 1986
CVERSION V01.0002   INCLUDES NORMAL EQUATION NORMALIZATION - 4/8/86
CVERSION V01.0003   IONOSPHERIC CONSTRAINTS BY PROPORTIONAL ACCURACY -
C    (4/11/86)
CVERSION V01.0004   FIXED BUG IN COVST1 - PROPAGATION TO KEPLERIAN
C   COVARIANCES FROM STATE VECTORS - 4/13/86
CVERSION V01.0005   FIXED "SIGMA" BUG IN BIAS SEARCH - WAS NOT SORTING
C   BIASES CORRECTLY FOR MORE THAN TWO STATIONS -
C   4/18/86
CVERSION V01.0006   IN UPSIT - REMOVED CHI-SQUARED LIMIT FOR UPDATE
CVERSION V01.0007  IN DOUBLE- FIXED IUSE ARRAY ERROR, DOES NOT HAVE AN
C   EFFECT ON PREVIOUS WORK BUT COULD AFFECT IBM AT
C   VERSION (BECAUSE OF OVERLAY QUESTION) - 5/29/86
C
C   IN DOUBLE - FIXED IUSET ARRAY, ALSO ADDED THE
C   SATELLITE POINTERS (ARRAY ISATXX) - 6/10/86
C
C   SENT TO AERO SERVICE - 6/10/86
CVERSION V01.0008   FIXED BUG IN COVBLN - WRONG PROPAGATION INTO XYZ
C   COVARIANCES WHEN NUMBER OF STATIONS GREATER THAN 2
C   ADDED SUBROUTINE LIVE - 7/7/86
CVERSION V01.0009   ABILITY TO APPLY WEIGHTS TO SITE COORDINATES
C    CHANGED STATWT,LSQUAR & WSTAT - 7/8/86
CVERSION V01.0010   EXPANDED FIRST DIMENSION OF TPART ARRAY TO 20
C    SUBROUTINES LSQX,LSQUAR,NORMD,PCLOCK - 7/9/86
CVERSION V01.0011   DOUBLE-DIFFERENCE OPERATOR ALGORITHMS BYPASSED
C    IF OBSERVATION PATTERN IS THE SAME -
C    SUBROUTINES LSQX,LSQUAR,NORMD,OPERA,BSTAT 7/16/86
CVERSION V01.0012   ADDED ROUNDING BIAS SEARCH 7/17/86
C    SUBROUTINES LSQX AND NBIASR(NEW)
CVERSION V01.0013   DOUBLE-DIFFERENCE OPERATOR BYPASS UPDATE - 7/18/86
C    SUBROUTINES DOUBLE & NORMD (ARRAY IUSEBK)
C
CVERSION V01.0014   DOUBLE-DIFFERENCE OPERATOR BYPASS UPDATE - 7/26/86
C    ADDED VIRTUAL ARRAY W3 (SINCE W2 IS USED FOR
C    CONSTRAINT), LSQX,LSQUAR,NORMD
C   BYPASS IS COMMENTED OUT TEMPORARILY UNTIL BUG
C    CAN BE FOUND - 7/28/86
CVERSION V02.0014   COPIED TO AT - 7/29/86
CVERSION V02.0015   BISET BUG FIXED ON LSI - 8/15/86
CVERSION V02.0016   FIXED BUG IN DOUBLE - 9/12/86
CVERSION V02.0017   LOOP THRU BIAS ROUNDING AND SEARCHING IN LSQX-
C                           9/15/86
CVERSION V02.0018   MODIFIED LSQINT AND LSQIO - 10/7/86
CVERSION V02.0019   MODIFIED LSQINT - 11/3/86
CVERSION V3.10      MANUSCRIPTA GEODAETICA ALGORITHMS
C                    (WEIGHTING ALGORITHM) - 11/10/86
C                   SCHAFFRIN & BOCK ALGORITHM - 11/26/86
C                   INTERACTIVE FEATURES RESTORED - 1/7/87
C                   FOUND BUG IN CALL TO FXBIAS FROM NORMD
C                    (NPART INSTEAD OF NTPART) - 1/15/87
C                   MINOR MODIFICATIONS IN LSQERR, LSQDO1, UPDATE,
C                    UPSIT,LSQINT,NBIAS - 1/17/87
C                   NEW BSTAT SUBROUTINE - 1/18/87
C                   MODIFICATIONS IN QUICK ALGORITHM, BEFORE A NEW
C                    BIAS PARAMETER WAS BEING INTRODUCED REGARDLESS
C                    OF WHETHER THE SATELLITE APPEARED AGAIN.  ALSO,
C                    ALL BIASES WERE BEING REMOVED IMPLICITELY AT THE
C                    END OF A SERIES. THIS SHOULDN'T HAVE WORKED.
C                    NOW, A BIAS PARAMETER IS INTRODUCED ONLY IF THE
C                    SATELLITE RE-APPEARS, AND ONLY FOR AN INDEPENDENT
C                    SET OF BIASES AT THE END OF THE SERIES. ROUTINE
C                    BIAYES WAS MOVED INTO NORMD. BIAYES INTRODUCES
C                    BIASES FOR A SINGLE-DIFFERENCE AND QUICK SOLUTION
C                    AND THEN REMOVES THEM IMPLICITELY USING FXBIAS.
C                    1/24/87
C                   ADDED SUBROUTINE USELSS TO REMOVE PHASES THAT
C                    DO NOT CONTRIBUTE TO A DOUBLE DIFFERENCE - THIS
C                    FUNCTION HAS BEEN REMOVED FROM SUBROUTINE DOUBLE
C                    AND SHOULD SOLVE THE INFINITE LOOP EXPERIENCED
C                    PERIODICALLY IN FORMING THE DOUBLE DIFFERENCE
C                    OPERATOR - 1/27/87
C                   WENT BACK TO ALL BIASES REMOVED IMPLICITELY FOR
C                    QUICK AND SINGLE DIFFERENCE SOLUTIONS - 1/27/87
C                   FOUND BUG IN NORMD IN SETTING UP BIAS SLOTS
C                    (IPNTB'S) - 1/27/87
C                   FIXED BUG IN STATWT - 1/28/87
C                   CHANGED "SITCOR" TO "SESSION" FILE - 2/4/87
C                   IEDIT CHANGES IN ATWGHT,STATWT,SATWT,LSQIO,LSQINT -
C                   NOTE A BLANK AT END OF MASK FOR IEDIT - 2/5/87
C                   FIXED BUG IN ATWGHT - 2/12/87
C                    MINOR DEBUG MODS IN LSQMEN,LSQIO,LSQINT
C                   OUTPUT MODS IN LSQX,LSQINT,LSQDO1,LSQERR,NBIAS-2/12/87
C                   ADDED JSNGLE VARIABLE IN LSQX,LSQINT,LSQDO1-INDICATES
C                    THAT A SINGLE-DIFFERENCE SOLUTION HAS BEEN RUN
C                    PREVIOUSLY, IN COMMON BLOCK SFLAG - 2/14/87
C                   MINOR MOD IN LSQERR - 2/15/87
C                   Q-FILE OUTPUT MOD IN LSQDO1,ADDED ROUTINE CR - 2/16/87
C                   CHANGED NUMBER OF BIASES TO 72 IN BIAS-SEARCH
C                    ROUTINES (36 FOR L1 AND L2), ALSO EXAMINE THE
C                    10 LOWEST CHI-SQUARES (USED TO BE 5) BUT ONLY OUTPUT
C                    5 LOWEST, DID THIS SINCE SOMETIMES ALL 5 COLUMNS
C                    HAD THE SAME BIAS BUT 5TH CONTRAST WAS JUST SHORT
C                    OF CUTOFF - 2/17/87
C                   CHANGED CONTRAST CUTOFF TO 3 (INSTEAD OF 5) IN
C                    AUTOMATIC BIAS SEARCH MODE - 2/17/87
C                   FIX BIASES OF STATION THAT IS CLOSEST TO CENTROID OF
C                    STATIONS, MOD IN BIAYES AND NORMD, ADDED ROUTINE
C                    BCENTR - 2/18/87
C VERSION 3.11      SENT TO AERO - 2/18/87
C                   FIXED ROUTINE DOUBLE (AFTER IT HUNG UP, AGAIN), ARRAY
C                    IUSE REFERS TO THE ACTUAL PHASE MATRIX NOW - 2/19/87
C                   ADDED BIASES OPTION TO THE PARAMETER MENU IN LSQMEN,
C                    MODIFIED LSQHLP AND LSQX (LSQMEN NOW NEEDS ISLOT1)
C                    2/22/87
C                   RENAMED LSQX, LSQM - 2/28/87
C                   FORMAT CHANGES IN QHELP, LSQUAR, LSQINT - 3/3/87
C                   RENAME LSQM, LSQX FOR AERO SERVICE - 3/3/87
C                   ANOTHER BUG IN DOUBLE DIFFERENCE OPERATOR
C                    CONSTRUCTION, ADDED MATRIX IPHI - TEST IF ALL
C                    PHASES HAVE BEEN USED IN A DOUBLE DIFFERENCES,
C                    REINFORCES VECTOR CHECK - 3/6/87
C                   TO CHANGE GEODETIC DATUM TO WGS72 ELLIPSOID
C                    LINK ROUTINE COVBLN.WGS - 3/6/87
C VERSION 3.12      SENT TO AERO - 3/6/87
C                   FIXED ERROR IN PCLOCK - 3/9/87
C                   FIXED CALL TO GEOXYZ IN LSQIO - 3/18/87
C                   FIXED TWO BUGS IN BSTAT (DELETED JSTATX) - 4/9/87
C VERSION 3.13      TABS : COVST2
C                   IEDIT : ATWGHT,LSQIO,LSQINT,LSQUAR,SATWT,STATWT
C                   MENEDT : LSQMEN
C                   REMOVED IER FROM LSQUAR - 4/9/87
C                   ADD INQUIRE CALL IN LSQX, CALL TO LSQINT
C                   ADD DUMMY ARGUMENT IN IEDIT FOR APOLLO COMPATIBILITY
C                   ADD DUMMY ARGUMENT IN MENEDT FOR APOLLO COMPATIBILITY
C                   MODIFIED UPDATE FOR APOLLO COMPATIBILITY - 4/10/87
C VERSION 3.14      ADDED ATMOSPHERIC CONSTRAINTS
C                    ROUTINES APPLY,ATMCON,ATMDIF,ATMSTR,BASLIN,JEL
C                    MODIFIED NORMD
C                   CORRECT UPDATE FOR APOLLO COMPATIBILITY - 4/16/87
C                  FORM FEED (FF) UPDATE FOR APOLLO COMPATIBILITY
C                  ROUTINES LSQINT,LSQDO1,LSQERR,LSQX,NBIAS - 4/20/87
C                  ADDED GEODETIC DATUM OPTION
C                   LSQX,LSQINT,LSQIO,LSQDO1,COVBLN
C                   NEW : GDATUM - 4/2/87
C VERSION 3.20     SENT TO AERO - 4/22/87
C VERSION 3.21     ADDED INCLUDES FROM 'dimpar.fti' TO 39 SUBROUTINES
C                   D. DONG AND R. KING - 4/29/87
C VERSION 3.22     USE DONG VERSIONS OF DBCNT,NORMD(?),DOUBLE,NRMSCL
C                   NBIAS,BISET; CALCULATE MAXPRM FROM MAXSIT,MAXORB,
C                   AND MAXSAT (INCREASES MAXPRM FROM 160 TO 180);
C                   RECOMPILE ENTIRE PROGRAM - R. KING 6/23/87
C                   CONVERT TO USE WITH NEW C-FILE FORMAT;  ADD DONG
C                   CHANGES TO MENU AND PRINTOUT OF NRMS;  CHANGES TO
C                   15 SUBROUTINES - R. KING 6/27/87.
C VERSION 3.23     TRANSFER TO APOLLO, HARDWIRING GDETIC.DAT IN GDATUM
C                   FIX: ATMDIF,DELCON NON-INDEPENDENT DIFFERENCE
C                        OPERATOR FOR TROPO CONSTRAINT--MHM,JLD 870729
C                   FIX: NORMD INDICE BUG FOR TPART
C                   FIX: LSQDO1 IWFLAG AND IOFLAG SWITCHES
C                   FIX: LSQX TO ALLOW OPTION 3 (REPEAT WITH NEW DATA)
C                   MOD: NORMD TO SUPPRESS OUTPUT OF RESIDUALS IN BATCH
C                   MHM AND K. FEIGL 870730
C                   TROPOSPHERIC CONSTRAINT QUESTION ADDED, IATCON ADDED TO
C                    COMMON STWGHT:  LSQX,LSQUAR,NORMD -870731 MHM
C                   Fixed GM value in KEPXYZ - 7/5/87
C                   MOD: COVBLN, LSQDO1 Output NEU local on screen and Qfile
C                   MOD: BIAYES Correct Bias setting with missing data
C                        --YBock 870803
C VERSION 3.24     CHANGE: dimpar.fti TO FIT MAXSIT=20.
C                   ADD: INTEGER*4 IN JEL,ATMDIF,DELCON,NRMSCL,SOLVE1,SOLVE2,
C                        SOLVE3,LSQUAR,ATMCON,APPLY TO AVOID OVERFLOW.
C                   FIX: STATWT,SATWT,LSQMEN TO ALLOW 20-STATION MODE.
C                   ADD: COMMON BLOCK /NMAT/ INTO APPLY,DELCON,ATMCON TO FIT
C                        OPTION 2.
C                   CHANGE: IFORCE IN LSQX,LSQUAR TO SEPARATE OPTION 2 AND 3.
C                   ADD: CRITERIA OPTION IN BIAYES,NBIAS,LSQX.
C                   ADD: TROPOSPHERIC CONSTRAIN INFORMATION IN LSQUAR.
C                   ADD: COMMON BLOCK /PROTEC/ IN LSQUAR TO SAVE SOME VARIABLES.
C                   CORRECT: ERROR IN NORMD AND DOUBLE.DONG.
C                   ADD: COMMON BLOCK /PARA/ IN LSQINT.
C                   MOD: NORMD TO MAKE IMPLICIT BIASES MODE WORK AND TO TAKE
C                        OFF LARGE INTEGER IN FIRST DOUBLE DIFFERENCE EPOCH.
C                   ADD: SUBROUTINE GETUSR.
C                   MOD: UPSIT,LSQPRT.
C                       --- DND --- 870913
C VERSION 3.25     CHANGE: ADD CONSTANT TERM TO IONOSPHERIC CONSTRAINT,
C                          CHANGE DISTANCE-DEPENDENT VARIABLE PROK TO BKAPPA,
C                          AND NAME CONSTANT TERM AKAPPA
C                   MOD: IN LSQX MODIFY COMMON BLOCK PICKWT
C                   MOD: IN LSQINT MODIFY COMMON BLOCK PICKWT
C                   MOD: IN OPERA MODIFY COMMON BLOCK PICKWT
C                   MOD: IN OPERA2 MODIFY COMMON BLOCK PICKWT AND CALL TO CPMAT2
C                   MOD: REMOVE COMMON BLOCK PICKWT FROM NORMD, LSQUAR
C                   ADD: IN ATWGHT ENTER CONSTANT TERM AKAPPA
C                   NOTE: SEARCH FOR OTHER OCCURRENCES OF COMMON BLOCK PICKWT
C                       --- YB --- 870930
C VERSION 3.26      MOD: In NMU11, NM2233, NM2131, NU23, NM32, JEL,
C                         SOLVE1, SOLVE2, SOLVE3, NORMD, NRMSCL, LSQUAR
C                          - INT4 Problem
C                       --- YB --- 871019
c version 4.1       first version in SCCS, also mod o-file format
c                       ---kf 871021
c version 4.2       minor mods in LSQDO1 to fix mod o-file format
c                       ---kf 871022
c version 4.3       mod ATPA to initialize W array
c                       ---kf + dnd 871111
c version 4.4       mod KEPXYZ to avoid taking square root of negatvie QUAN1
c                       and crashing.  Now, stop, writing error message
c                       to screen and Q-file.  We think that QUAN1 only goes
c                       negative when ENORMOUS cycle slips are still present.
c                       Monday, November 16, 1987   6:50:48 pm (EST)
c                       dnd & kurt
c version 4.5       In NORMD call to LRGINT was skipped for single-differences, this
C                    has now been fixed - yb 11/23/87
C                   In LSQDO1 units for baseline vector (meters) are now correctly
C                    output to Q- and O-files - yb 11/23/87
c                   For session mode, moved call to BIAYES (from LSQUAR) back into
c                    the 1000 loop. The call to BIAYES is done session by session.
c                    BIAYES was modified to work session by session correctly.  The
c                    above changes necessitated a modification for the "2" option.
c                    In the "2" option, the normal equations are not refilled, only
c                    the parameters to be adjusted are varied.  Therefore, if the
c                    user uses the "biases" option in LSQMEN, a call to BIAYES must
C                    be made in LSQUAR outside of the 1000 loop which is skipped for
c                    the "2" option.  Therefore, have added an array FREESV that
c                    saves the previous set of independent biases set by BIAYES.
c                    BIAYES now has two modes.  The logical variable BIASES has been
c                    added to LSQMEN and BIAYES to distinguish between the biases
c                    being set (in LSQMEN) automatically with the "biases" option or
c                    individually by the user - yb 11/23/87
C                   In CHECK, output to unit 6, not seven.  In LSQX, dimension
C                     array SHFT(3) to be compatible with GDATUM (not used, though)
C                     YB 12/5/87
C                   In ATWGHT, BKAPPA is now correct (was PKAPPA before)
C                     YB 12/15/87
C                   Testing new CPHI and CKAPPA models, have new routines for
C                     NORMD,CPMAT2,OPERA and OPERA2 (.NEW), not implemented yet
C                     YB 12/15/87
c
C version 4.6       Implemented new CPHI and CKAPPA models, required modifications
C                     of NORMD,OPERA,OPERA2,CPMAT2,LSQINT, and ATWGHT - YB 12/17/87


c version 4.6 (according to SCCS)
c Thursday, February 11, 1988
c MAJOR CHANGES TO IMPLEMENT NEW WEIGHTING AND OPTIMAL D-OPERATOR ALGORITHMS:
C DONG AND YEHUDA GO WILD AND CHANGE ALL THESE ROUTINES:
C
C
C    atwght
C    biset
C    bsort
C    bstat
C    covbln
C    covst1
C    cpmat2
C    fxbias
C    lsqerr
C    lsqint
C    lsquar
C    lsqx
C    nbias
C    nbias1
C    nbiasr
C    nm2233
C    nm32
C    normd
C    nrmscl
C    opera
C    opera2
C    satwt
C    sigsrt
C    solve1
C    statwt
C    vinv2
C    wsat
C
C    these are new subroutines
C
C    *bcheck
C    *bisopt
C    *dbar
C    *dbias
C    *dopt
C    *nrmsc2
C    *remenu
C
C  VERSION 4.7          NEW NAME:  SOLVE  -- MHM 880409
C           Add Dong's fast bias-fixing algorithms and multi-session logic
C           Modify BCHECK, BISET, BISOPT, BSORT, COVBLN, DBIAS, LSQERR, LSQUAR,
C            NBIAS, NBIAS1, NBIASR, NORMD, REMENU and SOLVE
C           Add routines BNEW, BNEW1, BDECI, BEXPEC - yb and dnd 4/9/88
C           Hard-wire automatic bias search in BISOPT - yb 4/10/88
C           Format changes in LSQINT - yb 4/11/88
C           Add function NSNPRN to UPORB for g-file update - yb 4/11/88
C           Remove extraneous arrays from UPORB
C           In LSQMEN and DOPT, array OBFIL should be dimensioned MAXCFL not
C            MAXNET, however this will cause format (non-fatal) problem in LSQMEN in
C            multi-session simultaneous adjustment and should be fixed - yb 4/15/88
C           Modify MAXCFL = 40 and recompile all routines - yb 4/15/88
c
c  version 6.1   NEW RELEASE --  MHM880417
c  version 6.2   fix minor spelling errors is SVERSN
c  VERSION 6.3   Replace bad BCHECK: reimplement new fast bias-fixing algorithms
c                and multi-session logic -- MHM 880419
c  version 6.4   Accommodate NCHI=0 in NBIAS -- MHM 880420
c  version 6.5   ATWGHT:  Correct format for AKAPPA, BKAPPA read -- MHM 880424
c  version 6.6   ATWGHT:  Really correct format for AKAPPA, BKAPPA -- MHM 880426
c  version 6.7   NBIAS, NBIAS1: Put debug on certain output lines -- MHM 880502
c  version 6.8   NORMD: In quick solution, use 98 flagged points as good data and
c                 add a new bias parameter at the epoch of the flagged point - also
c                 remove amplitude test on IERR (done in MODEL) --- YB 880504
c  version 6.9   DOPT:  Correct bugs in creating the optimal
c                 bias parameters in the presence of missing data.
C                BDECI,NBIASR:  Institute a new criterion for fixing bias parameters
c                 DND, MHM 880509
c  version 6.10  Change algorithm of chi-square calculation to save CPU time.
c                 The modified subroutines are SOLVE2, SOLVE3.  --DND  880721
c  version 6.11  To fit multi-session mode.   --DND 880727
c  version 6.12  Allow PRN numbers output, and add END after last sat in UPORB;
c                 Remove debug from BCHECK -- rwk 880817
c  version 6.13  Fix bug in NORMD, change O-file format yb, kf 880911
c                 mod NORMD,SOLVE,LSQDO1
c  version 6.14  Fix L1 only bug in NORMD - yb 880915
c  version 6.15  Implement keyword output to Q and O files -- MHM 881111
c                 (batch and interactive).  See new routine KEYWRD
C                 for documentation.  Changes to LSQUAR,LSQERR,LSQDO1
c  version 6.16  BLOCK2 common block had too many variables in COVST1 & COVST2,
c                 moved variables to BLOCK8 common block
c                Add common blocks PROTEC, BICNT, AJUNK to SOLVE,
c                  update common block comments to reflect more accurately which common blocks
c                  are where
c                Modify variable NBIAS in common block BICNT to ibias (NBIAS is the name
c                  of a subroutine), this affects SOLVE,LSQUAR&NORMD, in NORMD change
c                  old variable IBIAS to JBIAS, then NBIAS to IBIAS
c                -- YB 881120
c  version 6.17  BDECI: Change rediculous [!] spelling mistakes.
C                LSQUAR: Add function to check dependent biases,Restore WSTAT.
C                INVERS: Add IJOB=4 for decomposition only.  -- DND 881122
c  version 6.18  USELSS: Avoid dangerous comparison with 0.0d0.
c                LSQDO1: Fix station name output bug in O-file with multi-session.
c                                       -DND  881202
c  version 6.19  INVERS: provide information when matrix ill-conditioned.
c                VINV2: provide information when matrix ill-conditioned.
c                WSAT : Fix integer*4 bug for NOBS
c                STATWT: Fit orbit-only mode.
c  version 6.20  Remove epoch loop in Lc after L1,L2 mode.
c                 Change :SOLVE, NORMD, LSQUAR, DBIAS, NMU11, NM2131,
C                         NM2233, NU23,APPLY,DOPT,LSQINT
C                 Add : LCNORM
C                DBAR : fix dimension bug in DR.
c                NBIAS : add LPART in BLOCK8
C                LSQIO : pass SEMI and FINV
c  version 6.21  NORMD : fix epoch option bug (nepoch-->iend)
c
c  version 6.22  should be the same as 6.21, except putting the
c                SCCS keywords back in this file
c                Kurt 890301
c  version 6.23  Ignore unweighted points (ierr = -1) in NORMD -- MHM 890319
c  version 6.24  Accept millimeter-level station constraints in
c                STATWT  -- MHM 890418
c  version 6.25  Correct integer*4 bug (NOBS) in WSTAT -- MHM 890421
c  version 6.26  5-decimal precision correlation output in LSQERR -- MHM 890505
c  version 6.27  Fix initialization bugs in NORMDI,LSQERR,KEPXYZ,NM2233 -- MHM 890611
c  version 7.1   Integer*4, Case Insensitive release -- MHM 890820
c                 Integer*4 changes to APPLY,ATMCON,ATMDIF,BASLIN,BCHECK,BISET,BISOPT,
c                  BNEW,BNEW1,COVBLN,COVST1,COVST2,DBIAS,DELCON,DOPT,FXBIAS,JEL,LSQDO1,
c                  LSQERR,LSQINT,LSQMEN,LSQPRT,LSQUAR,NBIAS,NBIAS1,NBIASR,NM2131,
c                  NM2233,NM32,NMU11,NORMD,NORMDI,NRMSCL,PCLOCK,REMENU,SOLVE1,SOLVE2,
c                  SOLVE2,SOLVE3,STATWT,VINV2,WSAT,WSTAT,SOLVE
c                 Case insensitive in BADSTA,GDATUM,KEYWRD,LCNORM,LSQDO1,LSQERR,LSQINT,
c                  LSQIO,LSQMEN,LSQPRT,LSQUAR,NORMD,NORMDI,SATWT,STATWT,UPORB,UPSIT,WSTAT,SOLVE
c                 Reverse ISLOT2 dimensions for column reads in APPLY,BCHECK,LSQUAR,NORMD
c                  NORMDI
C                 Reverse ITOR,TOR dimensions in EPTIME,LSQDO1,LSQUAR,PCLOCK
c                 Machine precision added to INVERS
c                 M- and C-file primitives in LSQIO,LSQUAR,NORMD,UPDATE,SOLVE
c                 Clean up fortran structure in NORMD,UPDATE
c  version 7.2    Correct logic for multi-session station and satellite constraints
c                  in STATWT, SATWT, and LSQUAR
c                 Correct bug in c-file read loop; change gap detection logic in NORMD
c  version 7.3    Generalized maxsit formats in STATWT -- MHM 890828
c
c  version 7.4    1. For minimac receivers, L2-L1 biases could be half-integer.
c                    Change NBIAS,NBIAS1,NBIASR,BDECI,REMENU to fit half-integer
c                    bias-fixing mode for minimac receivers.
c                 2. Automatic checking if there is any station or satellite which
c                    has no effective observation.  Then automaticly fix the
c                    corresponding parameters of orbits,coordinates,troposphere if
c                    they exist. (LSQUAR,NORMD,STATWT)   ---- DND 890909
c  version 7.5    1. put station and satellite constraints before epoch loop
c                    (LSQUAR,WSAT,WSTAT), so that the quick algorithm can also put
c                    station and satellite constraints.     --- DND 890924
c  version 7.6    1. solve additional bias parameters automatically
c                    (FXBIAS,NORMD)
c  version 7.7    1. fix bug in I > LPART + LBIAS1 (FXBIAS)
c  version 7.8    1. Add new approach: using LC mode solution to fix coordinates,
c                    atmos., orbits, then solve and fix L1,L2 separate mode WL
c                    biases. (BNEW,DOPT,GETWL,LSQINT,LSQUAR,NORMD,SOLVE)
c                 2. Output pseudo-range WL estimate (NORMD,PSEUWL) -DND 891026
c                 3. add common blck /biasod/ in REMENU to fix bug in HALF.
c                 4. change output title in LC first mode(GETWL).  -DND- 891206
c  version 7.9    1. move label 204 after getwl(LSQUAR)
c                 2. include l2flag=3 in KEYWRD           -DND- 891206
c  version 7.10   1. add baseline length constraint(all continental biases are
c                    always free) (NBIAS,NBIASR,DOPT)  (temporary ?!)
c                 2. fix dependent bias index bug in multi-session mode.(LSQUAR)
c                 3. fix index bug in multi-session mode.(FXBIAS,NORMD)
c                 4. change rounding area (1.0 --> 1.2) (BDECI)   -DND- 891219
c                 5. MAXSAT = 12  MHM 900104
c  version 7.11   1. fix site and sat weighting bug for LC mode. (LSQUAR,WSAT,
c                    WSTAT)
c                 2. 4 characters for site code. (LSQUAR,LSQIO,UPSIT,LSQDO1)
c                 3. fix keywird bug. (KEYWRD)
c  version 7.12   1. fix s-file updating bug. (UPSIT)
c                 2. Add nrms and factor of 3 comment. (LSQERR)
c                 3. Remove one debug show. (BNEW)
c  version 7.13   1. index can't use function upperc as dummy valiable.
c                    Fix s-file updating bug. (UPSIT)  -DND-  900207
c  version 7.14   Back port from sun: replace dimension a(1) with dimension a(*) in
c                 atmdif,covbln,covst1,covst2,lsqdo1,normd,wsat
c                 Kurt Wednesday, February 28, 1990   9:37:02 am (EST)
c  version 7.15   Delete routines proper and cross because they are
c                 in the library.
c                 Wednesday, February 28, 1990   10:57:03 am (EST)
c  version 7.16   Compile for dynamic storage and 68020 architecture.
c                 NBIAS,DBIAS,UPDATE,LSQUAR need static storage
c                 to fit in stack.
c                 Kurt Tuesday, March 6, 1990   10:30:57 am (EST)
c  version 7.17   Yehuda Bock April 29, 1990
c                 Modify UPSIT and LSQDO1 for new station coordinate
c                 conventions, added only comments to COVBLN
c                 Kurt Feigl May 6, 1990
c                 Also change format in LSQDO1 to allow long (global) baselines.
c                 NORMD Write sat numbers and C-file names in
c                 Q-file list of the number of double difference observations.
c                 Also, do not skip a line between epoch prints.
c
c  version 7.18   Convert baseline length of fixed station to meters
c                 Modification to LSQDO1 P. Morgan  900507
c
c  version 7.19   ?
c
c  version 7.20   Change format in LSQDO1.  Kurt Feigl 900507
c  version 7.21   Change format in LSQDO1.  Kurt Feigl 900508
c  version 7.22   Try to get hemispheres correct. Kurt Feigl 900509
c                 LSQERR: correct format of covariance printout.
c                 LSQDO1: do not write titles every 60 lines to O-file.
c                 UPSIT:  major clean up.
c  version 7.23   LSQIO: call SUICID if M-file cannot be found.
c  version 7.24   FXBIAS: description for diagonal term=0.
c                 INVERS: EPS=2.5D-10 to avoid the effect of cumulative error.
c                         DND 900618
c  version 7.25   change dimension of TOR and ITOR of
c                 EPTIME,LSQDO1,LSQUAR,PCLOCK,SOLVE    DND 900702
c  version 7.26   Collected small changes July 4, 1990 kurt
c         covst2  Print only sigmas, not correlations of orbital parameters.
c         invers  Change EPS from 2.5d-12 to 2.5d-10.
c         lsqdo1  Avoid Peter Morgan's problem of station mismatch.
c         satwt   Change format to match DRIVER.
c  version 7.27   Kurt 900714
c         satwt   IEDIT strings must be character*80.
c  version 7.28   Kurt 900815
c                 More sun changes.
c  version 7.29   Kurt 900815
c                 repair the inadvertant damage to ATPA
c  version 7.30   Found some more dimension (1) in ATPAL and MATPRT.
c  version 7.31   1. remove factor 3 from uncertainty calculation (CONBLN,COVST2,
c                    LSQDO1,LSQERR,
c                 2. global file output(GDATUM,GLOHED,LSQERR,LSQPRT,LSQUAR,SOLVE)
c                 3. remove debug output(DBIAS)
c                 4. change comments output(GETWL)
c                 5. using pseudo-range wide-lane(NBIAS,NBIASP,NBIASR,NORMD,PSEUWL,
c                    REMENU)
c  version 7.32   1. modify output format(statwt,atwght,getwl,lcnorm,lsquar,remenu,
c                           glohed,solve)
c                 2. modify criterion (getwl,nbiasp,keywrd)
c                 3. modify tor,itor for multi-session(pclock,normd,lsquar)
c                         DND 900901
c  version 7.33   Change calls to READC2, READC4, and READC5 for new C-file
c                    format, in LSQUAR, NORMD.  --rwk 901010.
c  version 8.1    New C-file format and Sun Backport. November 12, 1990 Kurt
c
c  version 8.2    1. add new option: tight constraint first, then with loose
c                    constraint (SOLVE,LSQUAR,GETWL,LCLOOS,LWSAT,LWSTAT,SATWT)
c                 2. reduce array's size to save space (BNEW,BNEW1,NBIAS,NBIAS1,
c                    NBIASR,NBIASP)
c  version 8.3    1. fix array inconsistancy (WSAT,LWSAT,NRMSC2)
c                 2. option 7: fix wide-lane biases with pseudo-range priority.
c                    (BISOPT,NBIAS,NBIASR,NBIASP,SOLVE,LSQINT,LSQUAR,KEYWRD,GETWL)
c                 3. Bob's pushing approach for bias flag. (NORMD)
c  version 8.4    1. fix bug which leads missing effective data occasionally(NORMD).
c  version 8.5    1. update LC normal mtrix in LC only mode (FXBIAS)
c                 2. fix a bug in DOPT              DND 910315
c  version 8.6    1. big cleaner work! reconstructure SOLVE
c                 2. change logic of double.ftn (previous logic could miss very
c                    few independent dd-combinations)
c                 3. add solve.fti
c                 4. input format of batch-file is changed(using new driver)
c                                                       DND 910321
c                  new subroutines :
c                      abnew        filobs       oline        stack
c                      add1         filomc       qhead1       tapdrv
c                      add2         filpar       qhead2       updat3
c                      add2i        fnddbi       qhead3       wlmode
c                      bakcfl       formn1       qhead4       zero1d
c                      bfwork       formn2       reodrd       zero1i
c                      copy1d       getcda       rmjunk       zero2d
c                      copy1i       gethed       sort1i       zero2i
c                      dbadst       getslt       sort4i
c                      filld        modelc       sortbl
c
c                  legends of current status :
c                      c = clean; h = half clean; p = part(30%) clean; d = dirty
c
c                  name        status      name        status   name        status
c                  abnew.ftn     p         filobs.ftn     c     nrmscl.ftn    p
c                  add1.ftn      c         filomc.ftn     c     nu23.ftn      d
c                  add2.ftn      c         filpar.ftn     c     oline.ftn     c
c                  add2i.ftn     c         fnddbi.ftn     c     opera.ftn     d
c                  apply.ftn     c         formn1.ftn     c     opera2.ftn    d
c                  atmcon.ftn    c         formn2.ftn     c     opera3.ftn    d
c                  atmdif.ftn    c         fxbias.ftn     h     pclock.ftn    d
c                  atmstr.ftn    c         gdatum.ftn     p     pseuwl.ftn    c
c                  atpa.ftn      c         geoxyz.ftn     h     qhead1.ftn    c
c                  atpdal.ftn    c         getcda.ftn     h     qhead2.ftn    c
c                  atwght.ftn    c         gethed.ftn     h     qhead3.ftn    h
c                  avgvar.ftn    c         getslt.ftn     c     qhead4.ftn    h
c                  badsta.ftn    h         getwl.ftn      h     qhelp.ftn     c
c                  bakcfl.ftn    c         glohed.ftn     h     raddeg.ftn    h
c                  baslin.ftn    c         gpgt.ftn       d     remenu.ftn    h
c                  bcheck.ftn    c         gtpgx.ftn      d     reodrd.ftn    c
c                  bdeci.ftn     c         invers.ftn     c     rmjunk.ftn    c
c                  bexpec.ftn    c         jel.ftn        c     satwt.ftn     c
c                  bfwork.ftn    c         kepxyz.ftn     d     solve.ftn     c
c                  biset.ftn     c         keywrd.ftn     p     solve1.ftn    d
c                  bisopt.ftn    c         lcloos.ftn     c     solve2.ftn    h
c                  bnew.ftn      c         lcnorm.ftn     h     solve3.ftn    h
c                  bnew1.ftn     c         live.ftn       h     solvlc.ftn    h
c                  bsort.ftn     c         lrgint.ftn     h     sort1i.ftn     c
c                  bstat.ftn     d         lsqdo1.ftn     h     sort4i.ftn     c
c                  copy1d.ftn     c        lsqerr.ftn     p     sortbl.ftn     c
c                  copy1i.ftn     c        lsqhlp.ftn    c      sphxyz.ftn     c
c                  covbln.ftn     c        lsqint.ftn    d      stack.ftn      c
c                  covst1.ftn     c        lsqio.ftn     p      statwt.ftn     c
c                  covst2.ftn     c        lsqmen.ftn    d      sversn.ftn     c
c                  cpmat2.ftn     d        lsqprt.ftn    p      tapdrv.ftn     c
c                  cr.ftn         c        lsquar.ftn    c      tformat.ftn    h
c                  csort.ftn      h        lwsat.ftn     d      trans.ftn      c
c                  dbadst.ftn     h        lwstat.ftn    d      transx.ftn     c
c                  dbar.ftn       h        matprt.ftn    h      updat3.ftn     c
c                  dbcnt.ftn      c        mmply.ftn     h      update.ftn     d
c                  dbias.ftn      h        modelc.ftn    c      uporb.ftn      d
c                  dcheck.ftn     c        nbias.ftn     p      upsit.ftn      d
c                  degfrm.ftn     d        nbias1.ftn    p      uselss.ftn     p
c                  degred.ftn     d        nbiasp.ftn    c      vinv2.ftn      c
c                  delcon.ftn     c        nbiasr.ftn    h      wlmode.ftn     c
c                  distnc.ftn     c        newfil.ftn    h      wsat.ftn       d
c                  dopt.ftn       h        nm2131.ftn    d      wstat.ftn      d
c                  dot.ftn        c        nm2233.ftn    d      xyzkep.ftn     p
c                  double.ftn     c        nm32.ftn      d      zero1d.ftn     c
c                  eptime.ftn     c        nmu11.ftn     d      zero1i.ftn     c
c                  filld.ftn      c        normd.ftn     c      zero2d.ftn     c
c                  fillj.ftn      c        nrmsc2.ftn    d      zero2i.ftn     c
c  version 8.7     1. clean NM2131,NM2233,NMU11,NM32,NU23
c                  2. add l2flag for call above subroutines(FORMN1,FORMN2,NORMD)
c                  3. reduce cumulated error(OPERA,OPERA2)
c                                                           Dong 910323
c  version 8.8     1. clean DBIAS, BNEW
c                  2. copy LCLOOS,LWSTAT,LWSAT,AVGVAR from test directory
c                     (missed copying when update version 8.6)
c                                                           Dong 910324
c  version 8.9     1. modify BNEW,NMU11,NU23,NM32,DBAR,OPERA,OPERA2 to save
c                     cpu time
c                  2. small modification (GPGT,GPGTX)       Dong 910326
c
c  version 8.10    SOLVE:  Add 'Normal stop in SOLVE'       King 910411
c
c  version 8.11    1. take out large integer after every bias flag
c                     (LRGINT,NORMD)                        Dong 910501
c  version 8.12    1. initialize wl1 in NORMD
c                  2. compare multi t-file names (GLOHED)
c  version 8.13    1. give information of non-exist G-file (UPORB)
c                  2. in the case of conflict wl estimates, take the better one
c                       (NBIASR)
c  version 8.14    UPORB:  Don't crash with non-existent G-file.  Dong/Kurt 910703
c                  SVERSB: Input version number in format statement YB 910707
c                  LSQINT: Cosmetic Q-file output change
c                  LSQERR: Cosmetic Q-file output change
c                  UPSIT:  lowercase check for station code
c                      YB 910711
c  version 8.15    LCNORM: real*8 scale
c
c  version 8.16    Major changes from Scripps for kinematic processsing:
c                  SATWT: Change format of radiation constraints to f7.2
c                      YB 91/07/24
c                  LSQPRT: Status=unknown for Q-file open
c                      YB 91/08/22
c                  UPSIT: Add significant place for radius (0.1 mm)
c                      YB 91/08/25
c                  Generalize multi-session option for different number of
c                   satellites in each session
c                  SOLVE,DOPT,GETCDA,LSQUAR,NORMD,QHEAD2,REMENU,GETHED,SATWT,BCHECK
c                      YB 91/09/30
c                  Update call to satwt in LCLOOS
c                      YB 10/2/91
c                  Remove erroneous comment in LSQMEN
c                  Do not redefine NTPART for L1 only solutions (in LSQINT)
c                      YB 10/27/91
c                  Fix L1 and add L2 only option
c                  L2FLAG=-1 (L1 only), L2FLAG=-2 (L2 only)
c                  L2FLAG=0 (L1 receiver)
c                  BCHECK,LSQINT,KEYWRD,NORMD,QHEAD2,FILOMC
c                  FORMN1,FORMN2,FILOBS,NMU11
c                      YB 11/1/91
c                  Add common block for jusit in LSQUAR
c                  New LRGINT,GLOHED,GTPGX (in SCCS version)
c                  Spelling in PSEUWL
c                      YB 11/3/91
c                  Initialize nobs in LSQUAR
c                      YB 11/10/91
c
c                  Dong's changes for bias flags and H-file format:
c                  1. avoid bias-flag confusion (LRGINT)
c                  2. modify h-file format (LSQUAR,GLOHED)     Dong 911102
c                  1. temporarily resolve bias-flag confusion (LRGINT,NORMD) 911104
c                  1. add ut1, x,y pole and nutation terms to h-file
c                  2. add antenna offset to h-file (SOLVE,LSQERR,GETHED,GLOHED)
c                                                   Dong 911217
c
c  version 8.17    merge Yehuda and Dong's later modifications  Dong 920122
c                  1. using modes instead of values to determine the bias status
c                    (SOLVE,LSQERR,KEYWRD)
c                  2. avoid singularity when no tight satellite constraint but has
c                     loose satellite constraint (LWSAT)
c                  3. add digital to ut1,x,y pole and antenna offset (GLOHED)
c  version 8.18    extend loose constraint solution to all modes
c                    (solve,lwsat,lwstat, +(new)loos12)   Dong 920130
c
c  version 8.19    Update UPDAT3 for l-file and i-file
c                  Add new subroutine UPL to update l-file
c                  Add new subroutine UPI to update i-file (incomplete for now)
c                  Yehuda Bock 2/2/92
c
c  version 8.20    Add ambiguity-resolution controls to input.  Changes to
c                    ATWGHT, BISOPT, BDECI, NBIASR, SOLVE.FTI.  Dong/King 92/2/26
c                  Remove junk variables from LSQUAR, LWSAT, WSAT, UPSIT, VINV2
c                    GLOHED, SOLVE.FTI.   Dong 92/2/26
c                  KEPXYZ: Return non-zero exit status.  Donnellan/King 92/2/26
c                  NORMD: Print only every 50th epoch or when action.  King 92/2/27
c                  New Q-file formats:  Changes to LSQDO1, QHEAD3, COVST2, LSQINT,
c                     QHEAD2,  King 92/2/27
c                  Add a space to the U coordinate (for grep) in Q- and O-files
c                     Changes to LSQDO1 and OLINE.   Oral/King 92/2/27
c
c  version 8.21    Fix bug in QHEAD3.   King 92/3/26
c                  SOLVE:   Get loose solution for case of LC_HELP.  Dong 92/3/30
c                  ATWGHT:  Fix bug in reading ambiguity-resolution parameters.  Dong 92/3/31
c         		   Correct bugs for loose constraints:  LWSAT, LWSTAT.  Dong 92/4/2
c                  NORMD:  Restore reml1, reml2 to LRGINT calls (lost in MIT Apollo
c                          version)  Dong/King 92/4/8
c                  New C-file format for skd and ircint:  Changes to GETHED. King 92/4/23
c                  DBADST: Fix undefined variable; comments added to LSQUAR.
c                  Eliminate unused variables in DELCON, DOPT, LSQERR, LSQMEN, REMENU
c                      SATWT, SOLVE2, SOLVE3, UPI, XYZKEP.   King 92/4/23
c                  LCLOOS:  Correct spelling error.
c                  SOLVE, LWSAT, LSTAT: Check l2flag=4 for LC_RANGE option.
c                      Bock/King  92/4/24
c                  GETHED:  Call to READC2 to include session number.  King 920429
c
c  Version 9.1     Changes copied to SIO by Bock  920505.  End of sccs.  King 920507.
c
c          9.11    GETHED:  Fix bug in isessn variable, affects multisession.  King 920507.
c                  NORMD :  Remove debug printout.   King 920507.
c                  NBIASR:  Ignore pseudorange value if very large.  Dong  920507
c          9.12    SOLVE, LSQDO1:  Read inputs correctly for LC option.  Dong 920519
c                  LWSAT:  Fix criterion for determining whether tight constraints
c                          were applied. Corrects erroneous loose-constraint solutions
c                          when 'tight' constraints applied first.  Dong 920604
c                  STATWT: Read station constraints for baseline option.  Dong 920604
c          9.13    HTOH, Makefile : Completely new HTOH to translate station names, rather
c                          than apply ties.  Murray 911021
c                         (Never put on Apollo; lost on sun in recent updates--King 920611)
c                  GLOHED: Change Earth-rotation-time format for clarity.  Herring/King 920611
c          9.14    SATWT : Make input format less restrictive.  Kurt/King 920706
c                  OLINE : Minor bug in printout.  Kurt/King 920708.
c          9.15    NBIASP: Bug in pseudorange bias check for option 7. Dong/Bock/King 920903
c                  BISOPT: Update MIT and SIO Suns to Apollo version.
c                  FORMN1: Aesthetic change.
c                  GETHED: Update MIT and SIO Suns to Apollo version.
c                  LRGINT: Replace MIT+SIO Apollo versions with MIT+SIO Sun versions. Bock/King
c                  FILOMC, GETHED, NORMD: Handle unmatched sampling intervals in
c                            quick solution.  Dong/King 920925
c          9.16    BDECI:  Change pause to SUICID.   Feigl/King  921001.
c                  QHEAD3: Add comments and generalize; convert zenith delays to meters.
c                            King  921005
c                  FILOMC: Correct bug in logic for unmatching sampling.  Dong/Bock/King 921006
c          9.17    SATWT : Upgrade menu to 24 satellites Fang/Bock 921008
c          9.18    STATWT: Reinstate yessat check (necessary for baseline solutions)
c                          Bock 921108
c          9.19    GLOHED: Update the output covariance to 500 lines
c                          (instead of 50), Bock 921128
c          9.20    GETHED: Fix multisession bias-fixing bug for half-cycle widelanes
c                          (copy lambda)  Dong/King 921120
c          9.21    GLOHED: Fix H-file label for pre/post nrms.  Herring/King 921203
c                  SOLVE:  default mode and original mode (with magic number > 1) for
c                            H-file output      Dong 921230
c                  GLOHED, KEYWRD, LSQERR: add solution mode for H-file name
c                            (GCR, GCX etc)   Dong 921230
c                  GLOHED, SATWT:  Update satellite selection output to 24 satellites. YB 930107
c                  SOLVE:  Add 'job' initialization.  Bock/Dong/King 930118
c                  BISOPT, SOLVE.FTI: Add control for bias-search criterion. Dong 930118
c                  BISET : Correct bug in chi-square search.  Dong 930118
c           9.23   NBIAS : Printout wide-lane ambiguity searches.
c                  ATWGHT: Change needed for bias-search control (lost by rwk).
c                  BISOPT: Printout ambiguity-resolution parameters.
c                            Dong/King 930204/930205
c                  ATWGHT: Fix bug in reading WL bias-search parameters.
c                  GETWL : Fix bug in using WL bias-search parameters.  Dong/King 930209
c                  ATWGHT, BISOPT : Change limits for detecting an input bias-search
c                          ratio.   King 930210
c           9.24   LSQPRT: Initialize h-file name
c                  NORMD,PSEUWL,REMENU: change 20 => maxsat
c                  QHEAD1: Add one more space for number of observations statistics
c                          debug in COVST1,COVST2
c                          Bock 930211
c                  GETWL: Change default search ratio from 3 to 5.     King 930215
c           9.25   Changes to common and calling arguments to pass separate ambiguity-
c                       resolution criteria for wide-lane and narrow lane.  Removed
c                       unused variable 'factor'.  Set stage for changing defaults:
c                       SOLVE.FTI, ATWGHT, BDECI, BISOPT, GETWL, NBIAS, NBIASR
c                  ATWGHT: Removed unused variables.
c                             King 930219
c                  Explicitly declare variables, changing all to R*8 if no reason
c                       otherwise:  BDECI, NBIAS, NBIASR, NBIAS1    King 930223
c                  LSQPRT:  Fix (Apollo) problem with bad characters in H-file name.
c                  GETWL, BISOPT: Set search ratio back to 3 (from 5) to match FIXDRV
c                        defaults (temporarily?)   King  930225
c           9.26   BISOPT: Fix format statement to read ambiguity parameters Bock 930306
c           9.27   ATWGHT, BDECI, BISOPT, GETWL, NBIASR:  Make bias-fixing defaults more conservative
c                        (old: 0.4, 0.4, 1000., 3. ; new: 0.15, 0.15, 1000. 10.)
c                        Changed also in FIXDRV.   King 930319.
c                  ATWGHT, GETWL : Move the printout for the wide-lane bias-fixing criteria from
c                        ATWGHT to GETWL to make sure display is accurate.  King 930322
c                  NBIASR, BDECI : Remove arbitrary factor of 3., scale wide-lane but not
c                           narrow-lane sigmas.   King 930322
c           9.28   NBIAS, NBIASR, QHEAD4, GETWL : Clean up internal and echoed documentation
c                       of ambiguity-resolution algorithms.   King 930331
c                  QHEAD4 : One more change to echoed ambiguity messages.  King 930405
c           9.29   SATWT : Incorporate Peng's change to generalize # of satellites YB 930531
c
c           MIT changes 930707 - 930728
c           9.30   SOLVE, LSQPRT, LSQIO, GDATUM, LSQINT, LSQMEN, STATWT, LCLOOS, LOOS12,
c                  LSQUAR, SATWT, WSTAT, WSAT, LWSTAT, LWSAT, UPDAT3, UPDATE, UPI, UPL, UPORB,
c                  LSQERR :   modified to fit new i/o system
c                  READ_BFILE : read all batch file options.    Dong 930707
c                  SOLVE, LSQUAR, STATWT, GDATUM, BISOPT, LCLOOS, LOOS12,
c                  SATWT, WSTAT, WSAT, LWSTAT, LWSAT:
c                       shift i/o system to FONDA style
c                  READ_BFL : read all batch file options with FONDA style.    Dong 930709
c                  SET_PARA,GET_SAT_APR,GET_SIT_APR : read repeated command.    Dong 930710
c                  Add FONDA library routines: GETCMD BLANK COUNT_ARG LIFT_ARG MATCH_NAME
c                       STRJST.   King 930713
c                  SOLVE:  Key file style on 'BATCH'   King 930713
c                  LSQUAR, new routine READ_BFL_SESSN (later READ_BFL_SESS) , Makefile: read
c                       session options of batch file   Dong 930716
c                  SET_PARA : Fix bug in selection of bias slots.  Dong/King  930722
c                  Makefile: Add new routines READ_BFILE, READ_BFL, READ_BFL_SESS,
c                       GETCMD, BLANK, COUNT_ARG, LIFT_ARG, MATCH_NAME, STRJST, SET_PARA.
c                  READ_BFL : Change keyword selection for quick solution.  King 930726
c                  READ_BFL : Trap missing keywords, use defaults for output file names.
c                               King 930727
c                  READ_BFILE : Skip reading G-file question if no orbit partials.   King 930728
c                  Debug in FNDDBI, READ_BFILE, BCHECK.
c
c           SIO changes 930728 - 930805
c           9.30   Multiple zenith delay parameter changes
c                    Add MAXATM=24 to dimpar.fti (i.e, a zenith delay
c                     per hour for 24 hour session).
c                  APPLY,KEYWRD,FXBIAS,FILPAR: Change islot1 check
c                  ZENWT: New routine to read zenith delay a prioris
c                  LSQUAR,RMJUNK,STAWT,ZENWT: add variable yesatm to
c                   SITSAT common block
c                  SOLVE: Add lowerc check for batch or interactive
c                  LSQUAR,NORMD,FILPAR,RMJUNK: Pass zenith delay info,
c                   including RLABEL
c                  FILPAR: Add common block sitsat
c                  QHEAD1: Expand space for number of observations
c                  WZEN: New subroutine to apply zenith delay constraints
c                        to normal equations
c                  LSQUAR: New call to WZEN
c                  LWSTAT,WSTAT: Correct dimension for COVPLR
c                  ATWGHT: Informative comment added
c                  BCHECK: Change islot1 check for multi zenith delays,
c                           remove some misleading comments
c                  DOPT,DBADST: correct spelling error
c                  GETHED: correct debug statement (remove nflag)
c                  LSQUAR: Initialized ALC array by maximum dimension
c                  YB 930728
c                  LSQUAR: Do not call wzen if not solving for zenith
c                          delay parameters YB 930804
c                  FILPAR,RMJUNK,NORMD: check zenith delay flag (izenw)
c                  ZENWT: read in zenith delay info, even if no zenith
c                         delay information.   YB 930804
c             9.31 LSQINT: Fix L1 and add L2 only option
c                     L2FLAG=-1 (L1 only), L2FLAG=-2 (L2 only)
c                     L2FLAG= 0 (L1 receiver)
c                  LSQUAR: declare mbias0 variable
c                      YB 930804
c                  SOLVE: Correct spelling
c                  LWSTAT: Case when no station constraints in constrained
c                          solutions
c                  LCLOOS: commented out temp1&temp2 test
c                      YB 930805
c                  Many routines: remove % for include statements YB 930807
c
c             9.32 Merge MIT and SIO changes.   King 930818
c                  LSQINT: SIO changes merged with MIT LSQINT and new READ_BFILE.
c                  Merge changes: LSQINT, LSQUAR, LWSTAT, SOLVE, STATWT, WSTAT,
c                  Accept MIT changes:  BISOPT, BLANK, COUNT_ARG, FNDDBI, GDATUM,
c                      GET_CMD, GET_SAT_APR, GET_SIT_APR, LIFT_ARG, LOOS12, LWSAT,
c                      LSQERR, LSQIO, LSQMEN, LSQPRT, MATCH_NAME, READ_BFL,
c                      READ_BFL_SESS, SET_PARA, STRJST, WSAT, UPORB, UPL, UPI,
c                      UPDAT3, UPDATE, LCLOOS (ignore SIO check for now)
c                  Accept SIO changes: APPLY, ATWGHT, BCHECK, DBADST, DOPT, FILPAR,
c                      FXBIAS, GETHED, KEYWRD, NORMD, QHEAD1, RMJUNK.  King 930820
c                  Incorporate Dong version of common /constr/ and add zencv (YB czen)
c                      and zencv2.  LCLOOS, LOOS12, LSQUAR, LWSAT, LWSTAT, READ_BFILE,
c                      READ_BFL.
c                  Change common /sitsast/ variable 'yesatm' to 'yeszen' in FILPAR,
c                      LSQUAR, RMJUNK, and STATWT, and add to READ_BFL and READ_BFILE.
c                  Add izenw to common /csat/ in LOOS12, LCLOOS, READ_BFILE, READ_BFL, SOLVE.
c                  New routine GET_ZEN_APR.  King 930823
c                  Makefile:  Add SIO WZEN, ZENWT
c                             Add MIT READ_BFILE, READ_BFL, READ_BFL_SESS, BLANK
c                                     COUNT_ARG, GETCMD, GET_SAT_APR, GET_SIT_APR,
c                                     GET_ZEN_APR, LIFT_ARG, MATCH_NAME, STRJST, SET_PARA
c                  Include zenith-delay weighting for loose solution.
c                       LCLOOS, LOOS12,
c                  Don't read zenith-delay constraints if not turned on. ZENWT
c                  Don't use zenith-delay partials if not turned on.  FILPAR
c                        King 930825
c                  ZENWT: Correct bug in scaling input zenith delays constraints.  King 930827
c                  LCLOOS, LOOS12: Fix dimension bug (izen).   King 930828
c                  Debug in LSQUAR, READ_BFL, REMENU, FNDDBI, LCLOOS.
c                  LSQUAR:  Call ATWGHT only if old-style batch file.  King 930908
c                  READ_BFL_SESS: Read number of zenith-delay parameters; require
c                     to be the same for all sites, for now.   King 930908
c                  GET_ZEN_APR: Printout zenith info correctly.  King 930908
c                  LSQUAR,   : Use common zencv variable.  King 930909
c                  For new-style batch, move zenith-delay number to a priori read and
c                     fix a priori writes to screen and q-file:  READ_BFL, READ_BFILE,
c                     READ_BFL_SESS, GET_SIT_APR, GET_SAT_APR, GET_ZEN_APR.   King 930909
c                  Add izen to common/constr/ : LCLOOS LOOS12 LSQUAR LWSAT LWSTAT
c                     READ_BFILE  READ_BFL.    King 930909
c                  Move Q-file printout  (LSQPRT, LSQINT) and GDATUM inside READ_BFL
c                     and READ_BFILE; changes to these latter two, plus SOLVE. King 930910
c                  Change zenith-delay internal units from centimeters to meters
c                     (changed in MODEL):  GET_ZEN_APR, QHEAD3, ZENWT.  King  930914
c                  Allow  stepwise or piecewise-linear model for multiple-zenith delays.
c                     New routine PZENTH (change Makefile); new common /zendly/; changes to
c                     LSQUAR, LCLOOS, LOOS12, NORMD, READ_BFL, GET_ZEN_APR.  Hardwire model
c                     (currently = constant) with old batch files, allow control with new.
c                        King 930915
c                  Remove old-style batchfile and single-difference options.
c                        SOLVE, BISOPT, LCLOOS, LOOS12, LSQUAR, READ_BFL Makefile
c                        (ATWGHT, LSQMEN, READ_BFILE, SATWT, STATWT, ZENWT removed). King 930918
c                  Restructure main and (non-solve.fti) commons.  Changes mainly to SOLVE,
c                        READ_BFL, and new routine WRITE_SOLN, but minor changes to
c                        48 subroutines.  King 930831
c                        SOLVE.FTI :  remove iounit, move iwflag to local common /flags/
c                  BIAYES: Removed from MIT Apollo directory (not used, not in /stdrel).
c                  Rationalize units of apriori station, satellite, and zenith-delay  errors
c                        (now universally sigmas, not variances:  GET_SAT_APR, GET_SIT_APR,
c                        GET_ZEN_APR, LCLOOS, LOOS12, LWSAT, LWSTAT, READ_BFL, WSAT, WSTAT, WZEN,
c                        SOLVE ; new routine LWZEN for loose solutions; Makefile.   King 930924
c                  Remove excess printout and form-feeds from Q-file:   LSQERR, COVST2,
c                  SET_PARA:  Add multiple zenith delays.   King 930925
c                  BISOPT: Fix printout of 'narrow'/'wide' criteria for L1L2_IND case.  King 930929
c                  READ_BFL: Fix test on zenith-delay parameter. King 930929
c                  LSQUAR: write 'free', 'r2sum', and 'ntpart' to scratch unit 27 for
c                         use in L1L2_INDEP-mode loose solution.  Add warning to
c                         LOOS12.    King 930929
c                  Remove % before 'include'  and change 'D' => 'CD' in all subroutines.
c                         King 931001.
c                  Correct mistyped and remove undeclared variables.  Shorten
c                         READ_BFL_SESSN to READ_BFL_SESS to accomodate Sun. King 931002.
c                  Makefile:  remove reference to makex.fti (needed for kinematics?) King 931005
c                  NRMSCL:  Correct spelling.  King 931005
c                  FNDDBI:  Modify debug to catch problems.
c                  GET_SAT_APR, FNDDBI, LSQERR:  Fix printout for many stns, params.  King 931007
c                  FILOMC:  Don't estimate bias at gap for QUICK (assume AUTCLN
c                         flags all non-fixed gaps).    Herring/King 931008
c                  NORMD, PZENTH: Fix logic for choosing tabular points and forming partials
c                         using piecewise-linear model.   Herring/King 931011
c                  GET_ZEN_APR: Force constant model if nzen = 1.   Herring/King 931011
c                  LSQERR: Don't write version number twice to o-file.  King 931011
c                  GLOHED: Skip zenith-delay parameters in writing H-file.  King 931012
c                  LCLOOS, LSQUAR, NORMD, PZENTH, READ_BFL : fix dimensions of idtzen.   King 931012
c                  GET_SAT_APR: Fix reading of tight constraints.  King 931013
c                  LSQDO1, ZENOUT (new), Makefile: Write a dated summary of zenith-delay
c                        estimates to the O-file.   Herring/King 931014
c                  WSAT :  Comment (debug) printout of covariance matrices.  King 931014
c                  LSQERR: Add session times to screen,q,o output.  King 931015
c                  UPORB:  Add year to calling arguments for NSNPRN.  King 931018
c                  Increase formats to 32 satellites: GLOHED, LSQERR, NORMDI, QHEAD2,
c                        READ_BFL_SESS.  King 931018
c                  ZENOUT: Fix declaration of zenmod, and use midpt for constant model.
c                        King 931018
c                  PSEUWL: Dimension idwl in /wl/ correctly with maxsat.  King 931020
c             9.33 GET_ZEN_APR: Set defaults for numzen and zenmod if no input; print
c                        model in qfile. King 931214/931217
c                  Makefile:  Remove ZENWT.  King 931214
c                  LSQERR: Don't write loose-solution correlations in log file. King 931214
c                  FILOMC:  Restore for now flagging gaps in quick. King 931216
c                  FILOMC, FXBIAS, KEYWRD, LSQINT, LSQUAR, NORMD, READ_BFL, SOLVE:
c                        Make permanent the option of quick but without flagging gaps.
c                        Set by (new) 'Biases = implicit' in batch file and implemented
c                        by lquick = 2 throughout SOLVE.  All but one test in SOLVE changed
c                        from lquick.eq.1 to lquick.ge.1.  King 931220/931221
c                  ZENOUT: Skip output of zenith-delay info if single parameter per site.
c                        King 931220
c                  GDATUM: Comment out echoing of datums.   King 931222
c                  FNDDBI: Remove debug.   King 940111
c                  GETHED: Stop if sites or sats exceed dimensions.  King 940111
c                  SVERSN: Add library version to header.   King 940111
c                  FXBIAS: Delete extraneous printout.   King 940113
c             9.34 UPORB: Change call to NSNPRN (now a subroutine).  King 940120
c                  GPGT:  Tested (ok) and cleaned up a bit.
c                  MMPLY: Tested (ok), no change.
c                  NRMSC2: Modified to look more like NRMSCL.
c                  WSAT:   Added some comments.
c                  INVER2, Makefile: New routine, modification of VINV2
c                  COVST1: Call INVER2
c                     Bock 940128
c                  READ_BFL: Fix narrow lane ambiguity read.  Bock 940201
c                  DBIAS,DELCON,LWSAT,OPERA,OPERA2,OPERA3,WSAT: replace INVERS by INVER2
c                  INVER2: Write output to unit 6.  Bock 940201
c                  KEPXYZ: Tested (ok), some new comments.
c                  COVST1: Use GEODYN routine ELEM for Keplerian elements propagation.
c                  ELEM, Makefile: Modified GEODYN routine for Keplerian propagation.
c                     Bock 940202
c                  READ_BFL_SESS: Fix excluding sites and satellites.  McClusky/King 940204
c                  BDECI: Remove 'sngl'. In ERF & ERFC new variable 'half'.  Bock/Fang 940206
c                  COVST1: Comment dimensioning of unused variables.  King 940207
c          9.35 LSQUAR, RMJUNK: Fix bug in adjustments when no observations on a site. King 940321
c          9.36 SOLVE,READ_BFL,KEYWRD: Restore L1_L2_INDEPENDENT option.
c                 Bock/Genrich 940406
c               SOLVE, BISOPT: Restore mopt=3 option
c                 Bock/Genrich 940607
c               READ_BFL: Allow 'LC    ' as well as 'LC_ONLY' for LC mode.  King 940419
c               LSQDO1, QHEAD3: Fix bug in converting longitude adjust and sigma to meters.
c                 King 940419
c               LCLOOS: Fix SUICID print to identify subroutine correctly.   King 940419
c               REMENU: Set wavelength factor correctly for L2_ONLY biases.  King 940420
c               NRMSCL, NRMSC2: Set diagonal to 1 if zero to provide more graceful
c                    failure (probably not necessary).  Herring/King 940422
c          9.37 Makefile: Change HP compiler switches from +e +E1 to +U77.  Fang/King 940506
c               LSQINT: Record L1,L2 INDEPENDENT observations.  Bock 940508
c               SOLVE: Clear IEEE exceptions (commented for now).  King 9405010
c               Changes to incorporate earth-rotation parameters.  McClusky/King  11 May 94
c                 FILPAR: added code to add earth orientation partials to normal equations
c                         ["c" matrix.]. Added yeseop to sitsat common
c                 GET_EOP_APR: new subroutine to get a priori constraints for pole and ut1
c                              (earth orientation) parameters from new solve batch file.
c                  LOOS12: added call to lweop to apply loose earth orientation constrainsts
c                          to normal equations. Added eop_apr(6) and eop_apr2(6) to constr common,
c                          ieopw(2) to csat common and yeseop to sitsat common.
c                  LCLOOS: added call to lweop to apply loose earth orientation constrainsts
c                          to, normal equations. Added eop_apr(6) and eop_apr2(6) to constr common,
c                          ieopw(2) to csat common and yeseop to sitsat common.
c                  LWEOP:  new subroutine to calculate and apply loose constraints to
c                          normal equations.
c                  LSQUAR: call to new routine weop to apply earth orientation parameter
c                          constraints to normal equations. Added iweop(2) to csat common,
c                         yeseop to sitsat common, and eop_apr(6) and eop_apr2(6) to constr common.
c                  NORMD:  add iweop(2) to csat common statement.
c                  READ_BFL: added code to convert free matrix from (-3,-2,-1,0,1,2,3) to 0 or 1
c                            (see changes to set_para.f). Also added call to new routine get_eop_apr.f
c                            Added eop_apr(6) and eop_apr2(6) to constr common, ieopw(2) to csat common
c                            and yeseop to sitsat common.
c                  SET_PARA: changed to ensure all_parameters statements do not overwrite
c                            either group set flags (2/-2) or station/satellite set flags (3/-3).
c                            Also added code to read global settings.
c                  SOLVE:  added iweop(2) to csat common, yeseop to sitsat common, and eop_apr(6) and
c                          eop_apr2(6) to constr common.
c                  WEOP:   new subroutine to apply apriori tight constraints to earth orientation
c                          parameters.
c                  LWSAT, LWSTAT, LWZEN, WSAT, WSTAT, WZEN: added eop_apr(6) and eop_apr2(6) to
c                          constr common
c                  RMJUNK, ZENWT: added yeseop to sitsat common
c version 9.38  RMJUNK: Fix bug when two sties are excluded.  King 940607
c               SOLVE, LOOS12:  Fix L1/L2 loose solution bug and pass phase-obs keyword correctly
c                          for L1/L2 loose solution; get SOLVE to call LOOS12 to do L1 or
c                          L2-only solutions.  King 940607
c               GET_SAT_APR: Increase field size for printing a priori constraints.  King 940608
c version 9.39  Modifications to use Gauss Markov weighting of atmospheric parameters. McClusky 940616
c                  GET_ZEN_APR, READ_BFL : Get zen_mar, zen_mar2, tau, and tau2 from batch file.
c                  WZEN, LWZEN: Apply Gauss-Markov process to weight matrix.
c                  Add zen_mar, zen_mar2, tau, and tau2 to common /constr/.  Changes to LCLOOS, LOOS12
c                       LSQUAR, LWEOP, LWSTAT, LWSAT, LWZEN, READ_BFL, SOLVE, WEOP, WSAT, WSTAT, WZEN.
c               ZENOUT: Print both total and adjustments for zenith delays in O-file.   King 940620
c               READ_BFL: Check agreement between M-file and SOLVE batchfile numzen.  King 940624
c               Makefile: Shift ranlib to be executed only once.  Herring/King 940624
c               LWZEN: Remove debug for zenith weight matrices.  King 940624
c               SOLVE: Clean up by moving L1L2 Independent option to its own loop.  Bock 940625
c               BISOPT: Declare explicitly variable 'lane'.
c               NBIAS: Add a comment.
c               NBIASR: Remove option to replace phase estimate by pseudorange estimate in
c                widelane search.
c                 Bock 940702
c               NBIASR: Change goto to comments to avoid compiler warning.  King 940705
c               Remove unused routines from directories and Makefile:  LSQHLP, UPSIT.   King 940707
c               Remove unused routines from directories: DAYNUM, GETUSR, NSNPRN, UPPERC.   King 940707
c               READ_BFL: No longer read G-file line from input stream; create
c                   in SOLVE to reflect actual run conditions.  King 940722
c               UPORB: If time type blank, write 'UTC' for documentation.  King 940803
c               SOLVE: Add format to stop statement.  Bock 940810
c               UPORB: Change message for missing END statement.  King 940818
c version 9.40  READ_BFL: Declare argument for /lib/getusr.   King 940927
c               NBIAS: Remove call to pseudorange search (NBIASP).  Bock 941003
c               UPORB: Take q-file comment lines out of loop.  Bock 941008
c version 9.41  Fix dimensioning bug for biases.   Use maxbis to dimension bdeci and nfix in
c                 BISET, BNEW, BNEW1, BSORT, GETWL, LCLOOS, LOOS12,NBIAS, NBIAS1, NBIASP, NBIASR.
c                 Note:  A cleaner fix would be to put these in common, but since this would
c                 change a lot of calling arguments, I resisted the temptation.   Temporarily,
c                 at least, restore call to NBIASP in NBIAS.  King 941115
c               Remove unused ATPAL, ATPALX.   King 941115
c               Called INVER2 rather than INVERS in BNEW.  King 941115
c version 9.42  Allow input control for pseudorange widelane ambiguity resolution.
c                 BISOPT, NBIASP, NBIASR, READ_BFL_SESS, QHEAD4, solve.fti.  King 941118
c               Fix IEEE underflow in BDECI.   King 941118
c               GETHED, GLOHED, GDATUM: Add comment information to H-file.   King 941118
c               NORMDI:  Obsolete, remove from directories.  King 941123
c               Makefile:  Remove unused ZENWT.   King 950102
c               GET_EOP_APR, GET_SAT_APR, KEYWRD,  OLINE: Minor edits to make ANSI standard
c                  Fortran for XL/RISC compiler.  King 950103
c version 9.43  UPORB: Replace call to NSNPRN by call to SVNAV_READ.  King 950403
c               UPORB: Add hrs/min to SVNAV_READ call. King 950405
c               GETWL: Correct bug with printing of * symbol.  King 950420
c version 9.44  GETHED, GLOHED: Add reference frame and precession constant to h-file
c                 header; get temporarily from T-file.    King 950512
c version 9.45  Move EOP slots from 1900-1906 to 14001-14006 to make way for once-per-rev
c                 non-gravitational parameters in 1401-2000.  Add C-file partial number to
c                 common /parts/ in solve.fti.  Changes to COVST1, COVST2, FILPAR, GET_EOP,
c                 GET_SAT_APR, GLOHED, KEYWRD, LWEOP, LWSAT, QHEAD3, READ_BFL, SET_PARA,
c                 WEOP, WSAT.  Change common /constr/ also in SOLVE, LCLOOS, LOOS12, LSQUAR,
c                 LWSTAT,LWZEN UPDAT3, UPORB, WSTAT, WZEN.  King 950516/19/25
c               READ_BFL:  Don't allow a mismatch between M-file and batch file (temporary?)  King 950518
c     include '../includes/solve.h' to common/block2/
c version 9.46  New C- and M-file formats: GETHED, GLOHED, GETCDA, READ_BFL, UPDATE
c                  King  950523
c               GETHED, NORMD:  Remove calculation of zenith-delay tabular points; now done in CFMRG. King 950524
c               Change dimensions of izen in common/zendly/ from maxsit to maxcfl--this is academic
c                  at present since we force numzen to be the same for all sites, but it may cause
c                  confusion later.  Changes to GETHED,GET_ZEN_APR,  GLOHED, LCLOOS, LOOS12, LSQUAR,
c                  NORMD, PZENTH, READ_BFL, RMJUNK, SOLVE, WZEN, ZENOUT.  King 950524
c               QHEAD1:  Replace exlicit commons with included solve.fti.  Fixes mismatch of variables
c                  in common /block2/ leading to erroneous screen print for # of live parameters.  King 950525
c               SET_PARA:  Fix bug in setting free flag for EOPs.  King 950529
c               WSAT:  Allow Cartesian constraints for orbits.  Temporarily hardwired to
c                   Keplerian.   King 950530
c               GET_SAT_APR, WSAT, LWSAT:  Add control for Keplerian or Cartesian weighting of satellite
c                   initial conditions.  Change common /constr/also in LCLOOS, LOOS12, LSQUAR, LWEOP
c                   LWSTAT, LWZEN, READ_BFL, SOLVE, WEOP, WSTAT, WZEN
c               GETHED, GLOHED:  Change dimensions in /modcom/ from maxsit to maxcfl.  Set
c                   correctly elevation cutoff, antenna model, and number of zenith delays
c                   for H-file and write them with session rather than site parameters since
c                   they are stored on the C-file rather than the M-file.  King 950525
c               LSQUAR, READ_BFL_SESS: Issue warning and turn off tropsopheric spatial
c                   constraints if using multiple zenith delays.   King 950531
c               GETHED, GLOHED, UPDAT3, UPORB: Write radiation pressure model on H-file and
c                   correctly update parameter values on G-file.   King 950613
c               READ_BFL: Fix declaration of 'dattim', used for G-file header.  King 950613
c version 9.47  GETHED, UPDATE:  Add decimation factor to M-file record 2.  King 950616
c               UPORB:  Fix format for output G-file.  King 950619
c               DOPT,READ_BFL, BISOPT,NBIASR (SOLVE.FTI):  Allow input of maximum baseline
c                     length for bias fixing; increase default from 500 to 1000 km.   Bock 950604
c               READ_BFL, READ_BFL_SESS: Restore 500 km default for maximum baseline
c                     length for bias fixing.  King/Bock 950620
c               GET_SAT_APR: Removed bug in format statement   Tregoning 950630
c version 9.48  GETHED, GLOHED, SOLVE:  Change C-file record 2 to include Add antenna reference
c                   point offset; change dimensions of common/ante/anoff.  Add new common
c                   /rxant/ for rcvr/antenna type/sn.   Remove /ante/ from LSQERR (not used). King 950718
c version 9.49  READ_BFL: Remove duplicate declarations.   King 950719
c               GET_SAT_APR: Remove extraneous comma in format.   Herring/King 950719
c               GET_EOP_APR: Fix logic with missing constraints when parameters estimated.
c                    Sanli/King 950719
c               LSQERR: Open statement requires status for DEC.  Sanli/King 950719
c               LSQUAR: Fix typo in check for bad inversion.  Sanli/King 950719
c               READ_BFL, READ_BFL_SESS: Set default maximum for bias fixing if missing
c                     in batch file.   King 950721
c               GET_ZEN_APR: Stop if number zenith delays exceeds dimensions.  King 950721
c               GET_EOP_APR: Make variables read from batch file double precision.  King 950721
c version 9.50  GETCDA: Fix mismatched calling arguments (isnr/amp) to READC5.  King 950725
c     include '../includes/solve.h'pseudorange input variables to common/bcrit/. Sanli/King 950728
c               UPORB: Fix format for output g-file header.  King 950803
c version 9.51  READ_BFL: Set lrecl=250 (required for DEC to handle long records.  Sanli/King 950807
c version 9.52  GLOHED: Add time type to H-file.  King  9508011
c version 9.53  NORMD: remove screen-write of 'epoch < start epoch'.  King 950812
c version 9.54  Add satellite cartesian covariances to h-file.  Principal changes in
c                  GLOHED, LWSAT, WSAT.  Common changed in LCLOOS, LOOS12, LSQUAR,
c                  LWEOP, LWSTAT, LWZEN, READ_BFL, SOLVE, WEOP, WSTAT, WZEN.  King 050819
c               GET_SAT_APR, READ_BFL: Printout tight and loose satellite constraints consistently
c                  in km, km/s and dimensionless units (not ppm or %).  King 950822/23
c version 9.55  GET_ERR_APR: New routine to read station and satellite apriori data errors. McClusky 950823
c               READ_BFL: Added subroutine call to new routine GET_ERR_APR, plus added
c                         new common block called /data_noise/. McClusky 950823
c               SOLVE: New common block called /data_noise/ added. McClusky 950823
c               GLOHED: Remove blank line before a priori covariance.  Herring/King 950823
c               LWSAT: Remove multiple declarations, and 1 format bug  Tregoning/Fang 950824
c               GETCDA: Modified to read satellite elevation angles into an array. McClusky 950829
c               NORMD: Modified calling argument lists of GETCDA, FILOMC, and OPERA
c                      to pass satellite elevation angle arrays. McClusky 950829
c               OPERA: Modified calling argument list to pass satellite elevation angles
c                      to the routine CPMAT2 where they are used. McClusky 950829
c               OPERA2: Modified calling argument list to pass satellite elevation angles
c                      to the routine CPMAT2 where they are used. McClusky 950829
c               OPERA3: Modified calling argument list to pass satellite elevation angles
c                      to the routine CPMAT2 where they are used (I think this routine is obsolete?). McClusky 950829
c               FILOMC: Modified to fill epoch site/sat elevation angle array. McClusky 950829
c               CPMAT2: Modified to allow individual sites and satellites to have
c                       different data weights, and data weighting models. McClusky 950829
c               SOLVE.FTI: Added variables decim, and elvcut_solve to the /cutoff/ common
c                          Added a new common called /constants/.   McClusky 950829
c               READ_BFL: Modified to allow reading of a decimation interval, and solve
c                         cutoff angle from batch file. McClusky 950830
c               FILOMC: Modified to allow data decimation and new solve cutoff angle. McClusky 950830
c               NORMD: Modified to allow skipping of data when decimation is applied McClusky 950830
c version 9.56  GET_SAT_APR: Fix HP compile bug   Fang/Tregoning   950823
c               SOLVE: defined constants pi and convd. McClusky 950831
c               READ_BFL,LSQPRT,GDATUM: Fix hfile output for quick solutions.  Bock 950824
c               LSQPRT, GDATUM.  Continue hfile fix for quick solutions.
c               READ_BFL, LSQPRT, LSQERR.  Pass Owner name to h-file.  Bock 950825/950826
c               READ_BFL: Remove 'recl' in h-file open, added for DEC but prohibited
c                   by Sun.  Need to look for another workaround.  King 950830
c               READ_BFL, FILOMC,: minor cosmetic changes. McClusky 950831
c               NORMD: Fixed bug in decimation code. McClusky 950831
c               GLOHED: Modified to write out data decimation interval used, and
c                       the new measurement error model information. McClusky 950831
c               LSQERR: Modified to write decimation interval with epoch range. McClusky 950831
c               GETHED: Modified so that the elevation angle cutoff used at each
c                     station in the solution is stored in array elevcut. McClusky 950831
c               QHEAD2: Modified so that elevation angles used at each station
c                      are written to the qfile. McClusky 950831
c               GET_SOLVE_CUT: New routine to read cutoff angles from solve batch file. McClusky 950901
c               SOLVE.FTI: Modified the /cutoff/ common parameter list (solve_cutoff -> elvcut_solve
c                          this is in keeping with previous convention. McClusky 950901
c               GET_ERR_APR: cosmetic changes. McClusky 950901
c               READ_BFL: Added code to read cutoff elevation angle from solve batch file. McClusky 950901
c               LSQINT: Modified format of error model parameters written on qfile. McClusky 950901
c               CPMAT2: Modified to use a constant and an elevation angle dependant part in
c                       the computation of the a priori data VCV matrix. McClusky 950903
c               CPMAT2: Modified to pass in a variable describing whether data_error or
c                       ionospheric VCV is being formed. McClusky 950905
c               OPERA: Added icall to the calling argument list of cpmat2. McClusky 950905
c               OPERA2: Added icall to the calling argument list of cpmat2. McClusky 950905
c               READ_BFL: Initialised error1 and error2 variables. McClusky 950905
c               LSQERR, LSQINT: Remove printout of error terms (now in READ_BFL).  McClusky 950905
c               QHEAD2: Fixed bug in qfile output format. McClusky 950906
c               READ_BFL: Modified to allow baseline mode in new batch file format.  McClusky 950907
c               GET_ERR_APR: Modified to allow baseline mode in new batch file format.  McClusky 950907
c               UPDATE: Modified to write elevation angle used in solve solution to mfile.  McClusky 950907
c               READ_BFL,LSQPRT,LSQERR. Fix Quick solution h-file problem.  Bock 950911
c               GETHED: Added kludge to get around erroneous radians value for cutoff
c               in cfile header written by autcln. McClusky 950907
c version 9.57  GLOHED:  Correct conversion of coordinate variances from km to radians in case where a
c                    station is skipped.  Tregoning/King  950921
c               GETHED: Removed kludge to get around erroneous radians value for cutoff
c                     in cfile header written by autcln. McClusky 951002
c               GLOHED: Slight change in headers for rcvr/ant info.  King 951024
c version 9.57  GETHED: Modified code to stop null characters for receiver and antenna
c                     information getting written into hfile. McClusky 951109
c               RMJUNK: Fix bug in removing satellites with no observations in case
c                     when the number of satellite parameters > 9.  King 951110
c version 9.58  Makefile: Variable ranlib for compatibility with Solaris 2.  Fang/King 951208
c               GET_ERR_APR: Rewrite free format character read (for DEC).  Tregoning 960102
c               GET_ERR_APR: Fix bug in call to read_line.  Herring/Tregoning 960112
c               GET_ERR_APR: Change pointer used in decoding character string.  Tregoning 960129
c version 9.59  Changes to 57 routines to add REPORT_STAT.  King 960226-28
c               NBIAS1, Intercept interactive calls
c               Remove unused CHECK and ZENWT from directories, and CR, QHELP  from directories
c                   and Makefile.  King 960226
c version 9.60  LSQERR, QHEAD4:  Add status reports.  King 960318
c               DOPT: Correct subroutine name in call to stat_report.  King 960319
c               DBIAS, FNDDBI, LSQUAR: Add and rearrange status checks.  King 960319
c               GET_SAT_APR:  Trap case of 9 parameters with norb (from m-file) = 15.  King 960320
c     include '../includes/solve.h'BIAS1, NBIASR: Declare all variables explicitly to begin
c                  process of using none in code.   King 960320/960401
c               NBIAS, NBIASR: Remove redundant declarations.  King 960509
c               GET_EOP_APR: replace '' with ' ' in report_stat call (DEC bug) Tregoning 960527
c version 9.61  READ_BFL: Use 2-digit year in h-file name; temporarily get it by reading
c                 the second record of the m-file; make implicit none.  King 960607/960612
c               VINV2: Correct message format for report_stat.   King 960619
c version 9.62  NBIAS1: Comment out iflag variable (DEC complained).  Bock 960609
c               GET_ERR_APR, GET_SAT_APR, GET_SIT_APR, GET_SOLVE_CUT, GET_ZEN_APR,
c                   LWSAT, LWSTAT, LWZEN.  Add station and satellite names to
c                   printout of constraints.  King 960709  
c version 9.63  Changes for new kf/gamit structure and Makefiles from unimake.  King 960718
c                     All routines:  removed trailing blanks and replaced lib/*.fti 
c                     includes by ../includes/*.h. 
c               FILOMC: Temporarily put in code to allow high elevation angle cutoff. McClusky 960808 
c               BAKCFL: Correct number of arguments in calling report_stat; change missing
c                     menu-item message to fatal and add comments about the likelihood that
c                     'Choice of experiment = ORBIT' no longer works.  King 960823  
c version 9.64  GET_ZEN_APR, SET_PARA:  Replace functions 'lowera' and 'uppera' with calls to 
c                     'lowers' and 'uppers' to avoid segmentation fault if calling argument not 
c                      character*256. King 960829
c               LSQUAR:  Remove debug print 'calling getslt'.  King 960904
c version 9.65  SOLVE: Stop if GAMIT.fatal exists.  King 960912
c version 9.66  Extra parameter slots added for the estimation of 1/session N/W and E/W 
c               atmospheric gradient parameters. McClusky 961001
c               Additionally the EOP global parameter slots have been moved from 14001-14006 to 
c               80001-80006 allowing room for an additional 60,000 parameter slots. McClusky 961001
c               Routines associated with the EOP slot changes were: FILPAR, GET_EOP_APR, KEYWRD,
c               LWEOP, READ_BFL, SET_PARA, WEOP. Mostly the changes involved modifying "if" statements 
c               checking slot numbers from 1400? to 8000?. Except READ_BFL which checks to ensure that 
c               the mfile contains the updated islot1 information. McClusky 961001
c               Added yesgrad to sitsat common block. Routines modified include:
c               FILPAR, LCLOOS, LOOS12 LSQUAR, READ_BFL, and RMJUNK. McClusky 961001
c               Added apr_grad and apr_grad2 to constr common block. Routines modified include:
c               GLOHED, LCLOOS, LOOS12 LSQUAR, LWEOP, LWGRAD, LWSAT, LWSTAT, LWZEN, READ_BFL, 
c               SOLVE, WEOP, WGRAD, WSAT, WSTAT,and WZEN. McClusky 961001
c               Added igradw to csat common block. Routines modified include:
c               LCLOOS, LOOS12 LSQUAR,NORMD,READ_BFL,SOLVE. McClusky 961001
c               New routines added are: LWGRAD, WGRAD, GET_GRAD_APR, ATM_GRAD_PART. McClusky 961001
c               FILPAR: Added atm grad parameters to solution sub matrix "c".  McClusky 961001
c               LSQUAR: Added call to new routine WGRAD. McClusky 961001
c               RMJUNK: Added atm grad parameters to index list.  McClusky 961001 
c               LOOS12: Added call to LWGRAD, applys loose constraints to atm grad parameters. McClusky 961001
c               LCLOOS: Added call to LWGRAD, applys loose constraints to atm grad parameters. McClusky 961001
c               KEYWRD: Added GRAD and EOP to soultion keyword list
c               GETCDA: Added call to new routine ATM_GRAD_PART to calculate atm grad partials. McClusky 961001
c                       Also create the gpart array which contains the atm grad partials for each epoch. McClusky 961001 
c               READ_BFL: Added call to new routine GET_GRAD_APR which extracts the atm grad a priori 
c                       constraints from the solve batch file. McClusky 961001
c               SET_PARA: Modified to allow the atm grad parameters to be fixed or freed in the solution. McClusky 9611001
c               GLOHED: Modified to stop atm grad elements getting written into the hfile. McClusky 9611001 
c               QHEAD3: Rescale atm grad parameter estimates from, cycles @ zenith to, meters @ 10deg. McClusky 9611001 
c               DIMPAR.h: Added maxeop and maxgrad parameters. Modified the maxdm and maxprm parameters. McClusky 96100
c               SOLVE.h: Added the variable gpart, and gparts common which contain atm grad partials. McClusky 96100 
c               GET_EOP_APR, GET_ERR_APR, GET_SAT_APR:  Remove carriage control from 'message' formats to avoid 
c                    crash in REPORT_STAT.  King  961008
c               GET_ZEN_APR, OLINE: Fix format statements for gcc compiler.  Tregoning/King 961017 
c               GET_GRAD_APR:  Fix assignment of loose constraints.  McClusky/King 961017
c               GET_EOP_APR: Fixed Solaris 2 problem with indexing of temp array (again). McClusky/King 961029 
c version 9.67  RMJUNK: Correct bug in diminishing zenith-delay parameter count when a station has no observations;
c                   add code for gradient parameters.  King 961102             
c               INVER2:  Remove carriage-control characters from message passed to report_stat.
c                   King 961107 
c               ATM_GRAD_PART:  Remove duplicate 'implicit none'.  King 961119
c version 9.68  NRMSCL:  Minor correction to warning format.   King 961123
c               RMJUNK: Fix bug in removing gradient parameters when station has no observations.   King 961123
c               BISOPT, BNEW, BNEW1, GETWL, NBIAS, NBIAS1:  Remove extraneous print, fix up formats.   King 961123
c               BAKCFL:  Fix declaration of 'message' leading to eof abort.  King 961212`
c version 9.69  Fix and streamline logic for indicating presence of partials and estimation of parameters:
c                 Replace yessit, etc in /csat/ and istatw etc. in /sitsat/ with new common /parmflags/
c                 using sitpar, sitest, sitwgt, etc.  Assign default constraints for all parameters (except
c                 still allow no constraints on coordinates and orbits) and tighten and simplify traps for
c                 improper input.  Changes to solve.h, FILPAR, GET_EOP_APR, GET_ERR_APR, GET_GRAD_APR, 
c                 GET_SAT_APR, GET_SIT_APR, GET_ZEN_APR, LCLOOS, LOOS12, NORMD, READ_BFL, RMJUNK, SOLVE.
c                   --King 961230 
c               WGRAD:  Do not write a message to the screen.  King 961231 
c               LWGRAD: Write to screen and q-file the loose constraints.  King 961231
c               LSQERR, NBIAS, QHEAD1:  Minor trims of q-file print.   King 961231
c               FNDDBI, LCLOOS, LOOS12, LWGRAD: Move n-s/e-w loop inside to facilitate printout.  King 961231
c               LCLOOS, LOOS12:  Fix spelling error, declare implicit none, remove unused
c                  commons.  King 970103 
c               READ_BFL:  Fix test for obsolete m-file.   King 970108     
c               GLOHED:  Fix error handling on failed t-file open.  King 970109 
c               GET_GRAD_APR, GET_ZEN_APR, READ_BFL:  Fix logic in detecting missing entries.  King 970109
c               READ_BFL_SESS:  Fix format in call to report_stat.  King 970109
c               LSQUAR, UPDATE:  Close the M-file before solving normal equations to avoid
c                  its mysterious clobbering on GSFC Solaris system.  King 970109 
c version 9.70  Changes to initialize and declare all variables for G77 compiler: 
c                   BDECI, BISOPT, BNEW, BNEW1, BSTAT, COVBLN, CPMAT2, DEGFRM, DOPT
c                   DOUBLE, FILOBS, FXBIAS, GLOHED, INVER2, INVERS, KEYWRD, LSQIO, LWGRAD
c                   LWSTAT, NBIAS1, NBIASR, NM2233, RMJUNK, VINV2, WGRAD, WSTAT, NORMD, GLOHED,
c                   FNDDBI, DOPT, DEGRED, QHEAD4, LCLOOS.  King/McClusky 970121
c               GET_GRAD_APR, LWGRAD, QHEAD3: Modified to remove non standard intrinsic function
c                   (sind) not allowed by the G77 compiler. McClusky 970121
c               READ_BFL:  Remove test on m-file version (now caught in lib/readm1).  King 970123
c version 9.71  LRGINT:  New simpler version which also handles correctly residuals flattened
c                    by autcln's estimates of phase clocks.  Herring 970213
c               NORMD: Print out only every 100th epoch, not 50th.  King 970429
c               READ_BFL:  Correct format for writing epoch interval in q-file.  King 970429  
c               BDECI, ATM_GRAD_PART: Add 'd0' to constants to satisfy LINUX (F77 standard). King 970508
c version 9.72  GET_ERR_APR, GLOHED, LOOS12:   Avoid splitting Hollerith strings, disliked by picky DEC 
c                     compiler.  King 970606    
c               GET_ERR_APR:  Fix bug when station coordinates fixed; replace print/stop by
c                   calls to report_stat.  King 970606 
c               LSQERR: Fix bug in status print for bias-fixed solution.  King 970620    
c               LSQERR: Fix omitting of correlation printout for loose solution.  King 970829
c version 9.73  GETSLT, LSQUAR: Remove debug.  King 971014
c version 9.74  ATPA:  Change code in integer computation to avoid HP compiler bug.  King 971016
c version 9.75  READ_BFL: Fixed bug in weighting data by satellite. (nsat not norb). McClusky 971021  
c               GLOHED:  Fix format for satellite errors.  McClusky 971024
c version 9.75H FNDDBI:  Replace integer k(k-1)/2 by real computation to avoid HP 10.2 compilere
c                 bug with optimization level 2 or 3.  Herring/King  971204
c version 9.76  LSQERR, READ_BFL: Add control of printout of correlations; add variable to common /flags/.
c                  To tidy things up, moved /flags/ to solve.h, requiring changes to **45** others routines.
c               GET_GRAD_APR: Remove extra commas in two report_stat calls.  King 971223 
c               LSQUAR: Fix saving of 'free' array for loos12.  King 971223  
c               SOLVE:  Open temporary files in working directory, not /tmp space.  Herring 971111/9712130
c version 9.77  Makefile, ATPA, DBIAS, FNDDBI:  Replace i(i-1)/2 calculations with a call to
c                 (new) /lib/jelf to cope with HP compiler bug.  Move JEL to /lib for 
c                  consistency.  Herring/King 980121  
c version 9.78  READ_BFL:  Fix 'owner' bug from moving gloprt to common at v.9.76.  King 980204
c               READ_BFL:  Change default end epoch from 225 to 2880; remove debug from 971223 and
c                   980204 changes.  King 980217   
c               BDECI:  Declare erfc external to avoid (Linux) conflict with intrinsic.  King 980227 
c version 9.79  SOLVE: Make the scratch file names unique so that you can run more than one job in 
c                   a directory simultaneously.  Herring/King  980309   
c               COVST2: Write satellite covariance to o-file, not screen/log-file.  King 980501
c version 9.80  READ_BFL:  Correct bug in assigning apr values for gradient parameters.  King 980605
c               LSQINT, GET_GRAD_APR, GET_ERR_APR, GET_SIT_APR, GET_SOLVE_CUT,GET_ZEN_APR, LWGRAD,
c                  LWSTAT, LWZEN:  Add 4-char site code to station names in list.  King 980605
c               LSQDO1:  Replace station name by 4-char code in baseline printout.  King 980605
c               UPORB: Change argument list for svnav_read to include antenna offsets. King 980615
c version 9.81  UPORB: Add BERN2 radiation pressure model.  King 980708  
c version 9.82  LWZEN: Fix bug in indexing (garbage results but concidentally survives if arrays 
c                  are initially zero).  McClusky/King  980716    
c               UPORB: Add BERN1 radiation pressure model.  King 980720
c version 9.83  SOLVE, GET_ERR_APR, READ_BFL, READ_BFL_SESS:  Data weights now read from the session
c                  part of the batch file; N-file is read (unit 13) for overriding weights; pre-1995
c                  style of batch-file weights no longer supported.  King 980729   
c               READ_BFL_SESS: Change 'exclud ' to exclude' in getcmd call.  King 980804  
c version 9.84  includes/solve.h:  Add svantest and svantwgt to common /parmflags/.  King 980910  
c               includes/dimpar.h:  Increase 'maxorb' to 18 to allow fof SV antenna offsets. King 980910
c               FILPAR, GET_SAT_APR, LWSAT, READ_BFL, SET_PARA, WSAT:  Include SV antenna offsets 
c                  ('svant') in satellite parameter list.  King 980909/980911                 
c               BCHECK, KEYWRD, QHEAD3:  Correct comments.  King 980911      
c               GETHED: Add 'norb' to calling arguments for lib/readc2, and make explicit none. King 980911
c               solve.h, GETHED,GET_ZEN_APR,GLOHED,LSQUAR,LWGRAD,LWZEN,NORMD,PZENTH,READ_BFL,RMJUNK
c                   ,SOLVE,WGRAD,WZEN,ZENOUT: Replace common/zendly/ with common /atmprms/, in solve.h.  King 980916 
c               Make implicit none:  GLOHED, LSQUAR, PZENTH, SOLVE, WZEN, WGRAD, ZENOUT  
c               Makefile, ADDWGT:  New routine to add weight matrix to normal eqns.  King 98916 
c               WSAT: Use addwgt.  King 980916     
c               FILPAR, UPDATE:  Fix to accomodate new labels from model/wrthed.   King 980917            
c               solve.h, GET_GRAD_APR, GET_SAT_APR, GET_SIT_APR, GET_EOP_APR, GET_ZEN_APR, GLOHED, LSQUAR, LCLOOS, 
c                  LOOS12, LWEOP, LWSAT, LWSTAT, LWZEN, READ_BFL, SOLVE, WEOP, WGRAD, WSAT, WSTAT, WZEN:  
c                  Move common /constr/ to solve.h and add variables for multiple gradient parameters.  King 980917
c               Add multiple gradient parameters:  FILPAR, GET_GRAD_APR, GETHED, GRADOUT(new), LSQUAR, 
c                  LSQDO1, LWGRAD, PGRAD (new), WGRAD, RMJUNK.   King 980917/980929
c               GLOHED: Fix label change for LAT and LONG using new labels.  King 980921      
c               WSTAT, WZEN, WGRAD, WEOP, LWSAT, LWSTAT, LWZEN, LWGRAD, LWEOP : Use new subroutine ADDWGT.  King 980929
c               FILPAR: Replace with more efficient version. de Jonge/King 981001  
c               GETCDA: Add 'svclock' to lib/readc5 call; make implicit none.  King 981003   
c               LCLOOS, LOOS12, WRITE_SOLN:  Write better divider between solutions in Q-file.  King 981003
c               GETHED: Check match of M-file and solve batch file for # zenith and gradient parameters.  King 981006
c               GETCDA, GLOHED: Remove duplicate declarations.   King 981007   
c               GET_ERR_APR: Put fatal stop if sat error in baseline mode (need to fix code).  King 981007
c               GLOHED: Fix omission of SV antenna offsets in h-file.  Herring/King 981021    
c               GET_SAT_APR: Fix bug in SV antenna constraints for loose solution.  Herring/King 981022
c               GLOHED: Fix bug in writing SV antenna offsets to h-file.  Herring/King 981026
c               GET_ERR_APR: Remove faulty stop when no orbit parameters estimated.  Tregoning/King 981027    
c               GDATUM, LSQERR, LSQPRT, READ_BFL, WRITE_SOLN, solve.h:  Remove common /hfilmode/ and
c                  put variable 'ihmode' into common /flags/ in solve.h.   King 981127
c               GDATUM, LSQERR, LSQPRT, READ_BFL: Don't write h-file for prefit solutions unless ihmode>0. 
c                  King 981112/981125/981201 
c               RMJUNK: Fix bug in clock-parameter count (problems when SV excluded).  Tregoning/King 981222
c version 9.85  BADSTA, GLOHED, GRADOUT, LSQDO1, LSQIO, LSQUAR, OLINE, REMENU:  Use C-file name, not X-file 
c                  name to get station id for multisession, in order to work with simulation mode.  King 981230
c               ABNEW,  : 
c                  Replace 'dimension' with explicit typing to satisfy SGI compiler.  Morgan/King 981230
c               BADSTA: This routine already dummy; comment most code to avoid compiler warnings.  King 981230
c               DBADST: Fix dimensions on /protec/iuse.  Is this routine used?  Morgan/King 981230
c version 9.86  Explicitly typed variables in 78 subroutines to satisfy SGI compiler.  Morgan/King 981231
c               Makefile: Add -u to FFLAGS for Sun.                
c               Remove to /old_gamit/active unused routines:  ATPAL, ATPALX, BCENTR, GTPALX, OPERA3.  King 981231
c               Make all routines implicit none.   King 990105     
c               READ_BFL:  Set number of zenith and gradient parameters (izen,igrad) for the
c                 case where they are on the m-file but not estimated; fixes problem in rmjunk.f.  King 990125
c               GLOHED:  Fix dimensioning of anoff in common  /ante/.  Fang/King 990126       
c               NBIAS:  Restore lost call to jelf replacing i*(i+1)/2 calculation (HP fix)  King 990126
c               NBIASP: Fix serious declaration bug.  King 990127  
c               Fix duplication dimensions and/or order of declaration/dimensioning for
c                  calling arguments (g77 compiler with implicit none) in 19 subroutines:
c                  ADD2I, BSORT, COPY1I, DOPT, DOUBLE, ELEM,  FILOMC, GETHED, GETWL, LIVE, REMENU
c                  RMJUNK, SORT1I, SORT4I, SORTBL, SPHZYZ, UPDAT3, ZERO1I, ZERO2I.  King 990130
c version 9.87  LSQDO1, GRADOUT ZENOUT: Fix indexing for sigmas zenith-delay and gradient sigmas
c                  for o-file.  King 990209   
c               GET_ZEN_APR:  Fix report_stat call.  King 990209    
c               RMJUNK:  Remove SV antenna offset parameters from live list when SV
c                   excluded.  King 990209     
c               GRADOUT: Fix units and separately identify N/S, E/W gradients in o-file summary.  King 990213  
c               GET_ERR_APR: Check both constant and elevation-dep term for zero before stopping.  King 990320
c               GLOHED: Do not use comma as continuation character (IBM and f90 restriction).  Fang/King 990324
c               READ_BFL, UPL, solve.h:  Add tolerance for updating L-file coordinates (default 0.3 m).  King 990330
c               QHEAD3:  Fix failure to multiply by cosine of latitude in getting longitude adjustment 
c                  in meters.  King 990331 
c version 9.88  GET_SIT_APR: Removed a bogus piece code that caused loose station constraints to 
c                            be corrupted when observed sites == maxsit. McClusky 990402
c               GLOHED: Fixed indexing problem when writing apr variances to h-file. McClusky 990402
c               GETWL, LCNORM, SOLVE1: Add iostat to reads of temp files to catch occasional
c                 bus error on HP.  Herring/King 990405                                       
c               UPL:  Fix bug in checking L-file update tolerance.  MKing, RKing  990416   
c version 9.89  UPL: Add comment line for history; remove extraneous blanks in line.  Herring 990624
c               FILPAR: Fix bug in setting tpart when parameter not estimated.  King 990705   
c version 9.90  Changes for 4-digit years;  GLOHED, GRADOUT, OLINE, READ_BFL, UPORB, ZENOUT.  King 990728/0818
c               LSQERR: Fix data statement to avoid nulls in o-file.  King 990805             
c               INVERS: Fix format in message to report_stat.  Murray/King 990830
c               BCHECK, READ_BFL:  Always assume that the M-file bias count includes both L1 and L2-L1 
c                  biases, even if the receivers are L1-only; trap L1_RECEIVER option.  King 990906
c               GLOHED: Write out a priori covariances for SV antenna offsets.  King 991029
c               LWSAT:  Change integrated-parameter loops from maxorb_cov to norb_cov to avoid
c                    subtraction of undefined variable.  McClusky/King 991101
c version 9.91  LSQDO1, UPDAT3, UPL: Don't updates a station's coordinates in the output L-file 
c                    if the sigma exceeds 10 m.  Herring/King 991124   
c               READ_BFL: Fix bug in setting logical for SV antenna estimates.  DeJonge/King 991201
c version 9.92  READ_BFL: Y2K bug.  Fang 000101                                               
c               READ_BFL: Too many arguments for get_sat_apr.  Fang 000215   
c version 9.93  NEU baselines now refer to WGS84 ellipsoid
c               COVBLN: Calls new subroutine DXYZ_TO_NEU
c               LSQDO1: Baseline output in o- and q-files to 5 significant places.  Bock 000301
c               LCNORM: Fix report_stat call.  DeJonge 000303                                 
c version 9.94  UPORB: Declare 'iyr2'.  King 000816      
c               solve.h: Change common /stime/ to /times/ to avoid conflict with IRIX system routine. King 000816
c               SOLVE, MODELC, GETWL, QHEAD4:  Use common /params/ instead of calling arguments for
c                  arrays, in order to circumvent (not yet understood) segmentation violation with Linux. King 000928
c               NBIAS, NBIAS1:  Initialize 'isig' in NBIAS1 and test to avoid array overrun from undefined loop
c                   limit in NBIAS1.  Avoids occasional segmentation violation with Linux.  King 001019  
c               BISOPT: Fix dimensions of 'nfix' and 'bdev' to match other routines.  King 001019
c version 9.95  UPL: Fix bug in logic for updating L-file when some station parameters not adjusted.  King 010529
c               GRADOUT: Define 'gradkey'.   Morgan/King 010731     
c               UPORB: Fix format statement.  Morgan/King010731        
c               SOLVE, GETWL: Add iostat check for opening and reading of temporary files.  King 010815  
c               SOLVE: Fix missing assignment of phase_obs for LC-only case (pre-fit q-file).  King 010820  
c version  9.96 GLOHED: Change t-file missing from fatal to warning; set frames blank.  King 020311 
c version  9.97 Rearrange and rename commons and variables for counting partials and observations to
c                 make clearer, and make the m-file counter an array to allow for mixing single- and
c                 dual-frequency c-files:  solve.h,, BCHECK, DBADST, DBIAS, BSTAT, DOPT, FILPAR,  FXBIAS, GETCDA,
c                 GETHED, LSQUAR, NORMD, QHEAD1, QHEAD2, READ_BFL, READ_BFL_SESS, REMENU, RMJUNK, SOLVE, UPDATE.  King 020531 
c version  9.98 LSQUAR:  Fix bug introducted in 9.97 changes, miscounting live parameters. Herring/King 020607 
c version  9.99 UPL: Allow for updating apr file instead of l-file.  King 020807/020923   
c version 10.01 Makefile, READ_BFL, GLOHED (comments only): Replace solve/gdatum.f with new lib/read_gdatum.f
c                 which reads new (as well as old) style gdetic.dat.  King 021002  
c               Makefile, COVBLN, DXYZ_TO_NEU, LSQIO: Remove GEOXYZ in favor of library version; change name of 
c                  N+h from 'height' to 'geodrad' for clarity.  King 021002                   
c version 10.02 UPL:  For apr-style file, add epoch calculation and relax requirement that site name
c                  begin in column 2.  King 030520                                            
c               READ_BFL, LSQPRT: Clean up code writing the datum line in the h-file.  King 030606   
c version 10.03 GETWL:  Open and close units 27 and 29 to avoid too-fast reading with some machines.
c                 Herring 030638                                                              
c version 10.04 REMENU: Correct declaration of 'idwl' in common /wl/ to agree with normd.f and pseuwl.f. Agnew/King 030902
c               SOLVE:  Correct dimensions of 'anoff' in common 'ante' to agree with glohed.f.  Agnew/King 030902  
c               ABNEW: Remove unreference label.  King 031027  
c               LSQERR: Remove blanks from too-long literal. King 031027
c               BCHECK, BISET, DBADST, DOPT. GET_EOP_APR, GLOHED, LSQERR, LSQINT: Remove unused statement (ifc warnings). King 031027
c               Makefile: Add ifc compiler flags.  McClusky 031029                            
c version 10.1  NBIASR, NBIAS: Add extra printout to trace bias fixing.  Herring 031205
c               NBIAS:  Comment out bias search code.  Herring 031205          
c               READ_BFL:  Replace 'LC_RANGE' with 'LC_AUTCLN" for WL bias-fixing.  King 031208    
c               LRGINT, NORMD: Do not remove initial integer if LC_AUTCLN option selected.  Herring/King 031209 
c                  Temporary for testing: Also do not remove initial integer if 'acbias.dat' present.  King 040202
c               Makefile, SOLVE, GET_WIDELANE (new), RESOLVE_WL (new): Replace MODELC and GETWL by GET_WIDELANE,
c                   RESOLVE_WL, and LC_SOLUTION to streamline and handle new autcln WL option.  King031208  
c               READ_BFL : Disallow any bias-search option except decision function.  King 031208   
c               LSQUAR, FILLD: Add comments.  King 031208                                     
c               SIGSRT: Remove unused routine.  King 031211    
c               BCHECK, DBADST, DBIAS, DOPT, FNDDBI, LSQUAR, REMENU: Change name of 'lastbi' variable to 
c                 'last_nonbias' to clarify.  King 031216  
c               NBIAS: Add one more sig fig in interation printout.  Herring 031218 
c               DOPT, READ_BIASES, Makefile:  New routine to read DD ambiguities from autcln and set
c                  pointers consistent with the code in DOPT; add comments and clean up some statement
c                  numbers in DOPT.  King 040102  
c               NBIAS: Add one more sig fig in interation printout.  Herring 031218
c               NBIASR:  Add temporary code to use WL biases fixed by AUTCLN (superseded by code
c                 in new routines below).   Herring 031203-040107. 
c               RESOLVE_WL, SOLVE, Makefile:  Combine code from BISOPT and NBIAS, removing obsolete approaches
c                  to ambiguity resolution.  King 040109.
c               ROUND_BIAS, Makefile:  Rename NBIASR (kept temporarily for TAH code), removing
c                  chi2-search code (previously commented out).  King 040109
c               READ_BFL_SESS, GET_ERR_APR, READ_BIASES, solve.h: Add n-file name to common.  King 040126 
c               QHEAD4, GET_WIDELANE, RESOLVE_WL: Remove LC_RANGE option, and reorganize bias-summary printout 
c                  to handle both LC_HELP and LC_AUTCLN.  King 040202    
c               WRITE_SOLN: Add an overall solution header for free and fixed output.  King 040204
c version 10.11 UPL: Fix format problem with velocities written on L-file. King 040303        
c               READ_BIASES: Increase length of read buffer to 128 and rearrange EOF logic.  Herring/King 040304
c               GET_WIDELANE, GET_NARROWLANE: Add and/or reword echo to q-file of # biases fixed.  King 040304
c               GLOHED, READ_BFL, RMJUNK, WZEN, LWZEN, WGRAD, LWGRAD  : Include in the solution the average zenith delay along 
c                 with the knots of the the piecewise linear function, and in the H-file the average only.  King 040304  
c               WZEN:  New formulation of FOGM to get the mean near zero.  Herring 040406  
c               LSQERR, ZENOUT: Write to the o-file the combination of avg and PWL zenith delays.  King 040423  
c               LWZEN: Add new formulatoin for loose solution.  King 040511                   
c               GLOHED: Fix bug in writing a priori variances to h-file.  King 040511    
c               DOUBLE: Remove '/' from report_stat message.  Tregoning. 040526.   
c               DOPT: If refererence SV = 0 from acbias.dat (LC_HELP mode), let SOLVE choose the reference.
c                  Herring (040305 in /active) King 040625.                                   
c               ZENOUT: Fix bug in dimensioning of two new arrays.  King 040628               
c version 10.12 LSQINT: Write different title for LC_AUTCLN at beginning of q-file.  King 040707 
c               RMJUNK: Fix bug induced by new zenith-delay code.   King 040707               
c               LSQUAR, DBIAS, FNDDBI: Rename last_nonbias (in these routines) to last_bias after
c                 its meaning changed (clarification only--no effect on logic).  King 040707  
c               DOPT:  Fix bug in setting stations list when station missing.  King 040708    
c               READ_BIASES, DOPT: Fix bugs in setting pointers when sites or sats missing.  King 040802
c               LC_SOLUTION: Remove common /acbiases/ (not used).  King 040802
c version 10.13 LCLOOS, ROUND_BIAS: Fix problem of erroneous loose bias-free solution.  King 040803
c version 10.14 ZENOUT: Restore 'q-file' key in o-file zenith delays, so that sh_gamit_atmos will work.  King 040913
c version 10.15 SOLVE, READ_BFL, GET_WIDELANE, LC_SOLUTION, solve.h: Add input option for using local scratch.  King 041018
c               GLOHED: Increase size of integer format for obs on h-file. Shimada 041105     
c version 10.16 SOLVE, GET_WIDELANE, LC_SOLUTION, solve.h: Use full q-file name, not just the 4-character experiment 
c                 in naming scratch files; put the scratch file names in common in order to close and reopen them
c                 to avoid timing problems with the AMD machines.  King 050105                
c version 10.17 GETHED, GLOHED, FILOMC, QHEAD2, SOLVE, UPDATE, UPORB: Block number, svant model, atm and ocean 
c                 models from c-file to common/modcom/ to be written on the h-file.  King 050129   
c               GLOHED, READ_BFL_SESS, REMENU:  Convert site code to uppercase after setting from c-file name. King 050211
c version 10.18 READ_BIASES: Add trap for 0 biases from N-file.   Herring/King 050222         
c               FILPAR, NORMD : Remove unused 'delta' variable in calling sequence.  King 050224
c               AVGLOAD (new), GLOHED, NORMD: Estimate average loading values for h-file. King 050228
c               FNDDBI:  Add check on large value as well as small value in determining whether 
c                  matrix is singular for dependent biases.  Herring 050307                   
c version 10.19 WRITE_SOLN: Fix syntax error in format statement.  Herring 050314             
c               AVGLOAD: Change program id for report_stat call in solving normal equations for average 
c                 atmospheric load from 'ORBFIT' to 'SOLVE'; comment out most messages.  Tregoning/King 050504    
c               READ_BIASES: Read widelane as real*8 to support codeless receivers. Herring 050713  
c version 10.20 OLINE, ZENOUT, GRADOUT: Convert site codes to uppercase for o-file after setting from 
c                 c-file name.  King 050720                                                   
c version 10.21 AVGLOAD: Fix bug in loading calculating h-file avg loading for #sites > 30. Shimada/King 050818
c               FILOMC: Skip obs flagged bad model (igmodl), used for PCV model cutoff angle. King 050906
c version 10.22 QHEAD3: Set PWL zenith delay values to zero since adjustments are wrt average value. King 051029 
c               FILPAR, GET_WIDELANE, LC_SOLUTION: Removed unused variables.  King 051029     
c version 10.23 GETHED, GLOHED: Set solid-E, ocean-tide, and atmospheric tides application bits for the
c                 H-file from the sestbl/c-file 'tides applied' variable.  King 051115        
c               LSQINT: Fix bug in printing observable/ambiguity type in q-file.  Agnew/King 060105 
c version 10.24 READ_BFL: Fix bug in counting zenith delays when numzen = 1.  King 060206     
c version 10.25 UPL: Expand format velocity in output apr-style l-file to allow for motions on ice.  Herring 060524 
c version 10.26 READ_BFL: Fix bug in assigning slots for a priori values of gradient parameters, introduced
c                  by 060206 mod.  King 060620
c               READ_BFL_SESS: Make format length C*16 instead of C*14 (possible fix for compile problem
c                 with gcc-2.95.3 on AMD Turion64 .  King 060620   (change kept but probably unnecessary--
c                 see below:)
c               READ_BFL_SESS: Fix dimension of suse from '50' to 'maxsit'.  Palamartchouk/King 060626
c version 10.27 GET_ZEN_APR: Fix bug in counting zenith delays, introduced by 060620 mode.  King 060718/060730
c               SOLVE: Trap 'quick' option, which no longer works.   King 060718ll
c version 10.28 BNEW1, SOLVE1, LOOS12, LCNORM, FNDDBI: All inversions done by gamit/lib routine INVER2. McClusky 060908
c               VINV2 and INVERS: now obsolete. McClusky 060908
c               DOPT: Trap missing n-file before calling read_biases.  King 060929
c version 10.29 READ_BFL: Replace 'quick' and implicit bias echo with a trap since these options no
c                 longer supported.  King 061007       
c               LSQUAR, LCNORM, READ_BFL, WBIAS(new), ADDWGT,solve.h. Makefile:  Apply an apriori constraint  to the 
c                   bias parameters.  King/Herring 061007/061012
c               FNDDBI, FXBIAS, GET_ERR_APR, GET_GRAD_APR, GET_SAT_APR, GET_SIT_APR, GET_ZEN_APR, GET_NARROWLANE,
c                 GET_WIDELANE, LSQDO1, LSQERR, LSQINT, LSQUAR, LWGRAD, LWSTAT, LWZEN, QHEAD1, QHEAD2, QHEAD3, 
c                 QHEAD4, READ_BFL, READ_BFL_SESS, ROUND_BIAS SOLVE1 SVERSN, SOLVE_CUT, WRITE_SOLN, UPDATE, UPL, 
c                 WRITE_SOLN, solve.h:  Add batch-file flag ('logprt') controlling print to screen.  King 061012   
c               READ_BFL_SESS: Fix bug in reading ion constraints.  Shimada/King 061102
c               COVST2, GET_EOP_APR, LWSAT, RESOLVE_WL, UPORB: Suppress remaining screen print unless
c                 'logprt' true. King 061112
c version 10.30 FNDDBI: Use condition number (ier 140) rather than 'info' number (ier 130) to detect
c                  ill-conditioned matrix and hence dependent biases. McClusky 061115
c               BCHECK, FILPAR, FXBIAS, READ_BFL, SET_PARA: Expand L1 and L2-L1 bias slots to allow
c                  45 SVs (remove implicit bias slots).  King 061127 
c version 10.31 CPMAT2, READ_BFL_SESS, RESOLVE_WL: Fix long-standing error in ionospheric constraint 
c                (LC_HELP) (sigma too large by factor of 50).  deJonge/King 061129
c               GET_ERR_APR: Fix spelling report_stat call. King 061130
c               LSQUAR: Enlarge format field for echoing bias constraints in status file.  King 061201
c               GET_SIT_APR: Add another decimal place to apr coordinate echo in q-file.  King 061201
c version 10.32 Makefile: Fix compiler flags.  Shimada/McClusky 061229     
c               DOPT: Correct bug (typo) using WL-baseline-max for limiting NL resolution. King 070112
c               Makefile, CHECK_ADJUST (new), GET_SIT_APR, GET_ERR_APR, LSQDO1,  READ_BFL, READ_BFL_SESS: 
c                 read the  N-file for overrides of tight-solution constraints.  King 070118
c version 10.33 RESOLVE_WL: Change summary of WL criteria, adding grep'able lines for sh_gamit.  King 070209
c version 10.34 READ_BIASES: Fix bug in setting baseline dimensions (was no harm). King 070221
c               Makefile, GET_BIAS_SCALE (new), GET_WIDELANE, GET_NARROWLANE, RESOLVE_WL, QHEAD4, solve.h:
c                  Scale bias uncertainties by nrms of deviations from integer for each baseline. King 070302 
c               ROUND_BIAS: Remove session loop.  King 070302
c               BDECI, GET_BIAS_SCALE, RESOLVE_WL, ROUND_BIAS:  Add debug variable.  King 070302 
c               GET_BIAS_SCALE, READ_BFL, solve.h: Hidden feature to turn off NL rescaling.  King 070309
c               BNEW, BNEW1, COVST1, DBIAS, DELCON, FNDDBI, LCNORM, LOOS12, LWGRAD, LWSAT, OPERA, OPERA2
c                   SOLVE1, WGRAD, WSAT, WZEN: Add 'rcond' to calling argument for lib/inver2.  McClusky/King 070319
c               FNDDBI: Added bias matrix condition number ratio checking to FNDDBI. McClusky 070320
c               AVGLOAD, BDEC,  BEXPEC, BSTAT, CPMAT2, DBIAS, DEREAD, FNDDBI, GETCDA, GET_EOP_APR,
c                  GET_ERR_APR, GETSLT, GLOHED, GRADOUT, LCLOOS, LOOS12, LRGIST, LSQDO1, LSQIO, LSQUAR, LWSTAT, 
c                  LWZEN, NBIAS1, NBIASP, NBIASR, NM2131, NMU11, NORMD, NRMSC2, OPERA, PCLOCK, PSEUWL, 
c                  READ_BFL, READ_BFL_SESS, SET_PARA,  SOLVE1,  SOLVE2, SOLVELC, UPDAT3, UPI, USELSS, ZENOUT; 
c                  Makefile (remove badsta).  Clean up code to avoid warning messages with gfortran. King 070327.
c               READ_BIASES: Add printout to list biases when number doesn't match autcln.  Herring 070402    
c               Makefile, ADDWGT, ADDWGTB (new), WBIAS : Move bias-constraint code entirely to new routine since
c                  matrices to be filled are different.  King 070403.
c               READ_BFL: Change default for bias constraints to 1000. rather than 0.  Herring/King/McClusky 070409
c version 10.35 BFWORK, GET_BIAS_SCALE: Fix bug, now excluding from NL nrms deviations when WL not fixed.  King 070412
c               BCHECK, GETHED, LSQERR, LSQUAR, READ_BFL, READ_BFL_SESS, WL_FIXED (new), Makefile (remove BAKCFL):  
c                  Move opening and  reading of m-file and o-file  earlier to get session time information for the 
c                  q- and o-file print;  remove multisession c-file names; remove other multisession code.  King 070413
c               GET_BIAS_SCALE: Change variable name from 'include' to avoid problem with unimake. King 070416
c               READ_BFL: Remove call to clrscr.  King 070416. 
c               LSQINT, SVERSN: Remove print unit number from 'lversn' call.  King 070416  
c               SVERSN, LSQINT: Remove call to 'clrscr'.  King 070418 
c               FNDDBI, READ_+BFL, solve.h: Add a flag for bias debug.  McClusky/Herring/King 070430
c version 10.36 READ_BFL: Set the default rcond high so that no independent biases get removed.  Herring/McClusky/King 070501
c               FNDDBI: Fix format statement for writing out rcond.  King 070502
c               FNDDBI: Write debug format to q-file and log file, not screen.  Herring/King 070507
c               READ_BFL_SESS: Fix format bug with maxsit=100.  King 070530
c               LOOS12: Remove warning about untested L1/L2 solution; change comment; add commented debug.  King 070724
c               GET_SIT_APR: Add warning about site constraints loosened by prefit adjust.  King 070724      
c               LSQINT, READ_BFL, SOLVE, SVERSN: Remove unused calling argument.  King 070910  
c               NRMSCL, NRMSC2, OPERA, PCLOCK, SOLVE: Remove unused statement label.  King 070910
c version 10.37 GLOHED: Add nutation and gravity model to h-file.  King 071221   
c               SOLVE:  Change comment about unit numbers for c-files.  King 080121
c               CHECK_ADJUST: Write 'NONE' in q-file if no adjustments gt tolerance. King 080214    
c               GET_BIAS_SCALE:  Add dummy branch to satisfy compiler.  King 080509
c               GET_WIDELANE: Use autcln WL even if not fixed.  Herring 080529
c version 10.38 GET_WIDELANE:  Reset WL to zero if dependent bias.  Herring 080609
c               LSQDO1, UPDAT3, UPL: Distinguish which 8-character l-file sites to update based on
c                 the agreement of their coordinates with the apriori used in the run.  Herring/King/McClusky 080611
c version 10.39 DOPT, LRGINT, NBIASR, RESOLVE_WL: Remove all references to file acbias.dat.  King 080724            
c               SOLVE, READ_BFL: Add option to skip the loose solutions.  King 080905     
c               UPL:  Use velocity to match sites to update (needed for fast-moving sites). King 081113
c               NORMD: Initialize atml and hyrological loading averages to zero (for case of no data). King 081212
c version 10.40 GET_BIAS_SCALE: Fix bug in computing scale factor for half-wavelength L2.  King 090107 
c               BDECI: Add cutoff threshold to debug print.  King 090107
c               GET_SIT_APR: Fix logic for warning about coordinate constraints loosened.  King 090128/090708 
c               Makefile: Change compiler flag for Intel from -xN to -Xt.  McClusky/King 090708 
c               BISOPT, GET_NARROWLANE, NBIAS, NBIASR: Add 'logprt' check to writes to screen.  King 091021  
c               AVGLOAD: Initialize 'adjload' (not clear that fixed problem).  King 091106
c version 10.41 SOLVE,  GETHED, GLOHED, LSQINT, READ_BFL, READ_BFL_SESS, solve.h: Move common /rxant/ to solve.h; remove 
c                 supeflous listing of stations and  SVs in LSQINT; add print of station rcvr/ant info 
c                 to q-file in GETHED; remove obsolete single-difference code from lsqint. King 100209
c               GET_BIAS_SCALE: Add a blank line after the rescaling summary. King 100824
c version 10.42 SOLVE, AVGLOAD, FILOMC, GETCDA, GETHED, GLOHED, LSQPRT, NORMD, QHEAD2, UPDATE, UPORB:  Change h-file format to 
c                 accommodate longer rcvr and SV antenna models, and add h-file version number (2.0). King 100827
c               GLOHED: Add another digit to the number of observations.  Herring/King 100830
c               GLOHED: Small change to column headers.  King 101007/101104     
c               FILPAR(comment only), UPORB : Add UCLR1 radiation-pressure model. Petrie/King 110124
c version 10.43 GETHED: Change mapping of ietide for h-files; mean-pole-tide removed is now bit 4 (8),
c                  rather than bit 5 (16). Herring/King 110524
c               QHEAD3; Fix parameter-name test so that no confusion between orbits and
c                  sites named 'RADx' or 'ORBx'. King 111227
c               GETCDA, AVGLOAD: Fix R*4/R8* confusion in reading atm load from c-file, and several
c                  initialization and units errors in generating the weights for the (parallel, 
c                  not-used) straight-average calculation. King 120111
c version 10.44 QHEAD3: Fix bug from 111227 change.  King 120418
c               AVGLOAD: Fix bug in computing average hydrological load.  King 120524
c               GLOHED: Fix bug in setting iaload for h-file when unfiltered ATML is used.  Herring/King 120727
c version 10.45 AVGLOAD, FILOMC, GETHED, GLOHED, NORMD, QHEAD2, SOLVE, UPDATE, UPORB:  Add dryzen and wetzen to
c                  model information passed to the h-file (temporary setting until c-file format revised). King 121219
c version 10.46 GETHED: Replace temporary setting of dryzen and wetzen by c-file values (new format). King 130115
c version 10.47 AVGLOAD, FILOMC, GETHED, GLOHED, NORMD, QHEAD2, SOLVE, UPDATE, UPORB: Add E-radiation and
c                 antenna-radiation to the c-file and h-file.  King 140327    
c               AVGLOAD, FILOMC, GETHED, GLOHED, NORMD, QHEAD@, SOLVE , UPDATE, UPORB: Add ion source and
c                 magnetic field to the c-file and h-file. King 140401 
c               GET_BIAS_SCALE: Change vairable name 'include_dev' to 'use_dev' to avoid an erroneous
c                  grep by unimake. King 140409
c version 10.48 AVGLOAD, GETCDA, GLORED, NORMD: Use arithmetic average of atmospheric and hydrologic loads
c                  when the site is within 5 degrees of a pole to avoid ill-conditioned longitude. King 140515
c version 10.49 Most routines: move solve.h to /solve from /includes. King 141004
c                 30 routines: move common /params/ to parameters.h.
c                 7 routines: move common /modcom/ to models.h 
c                  models.h: Replace I*4 iblk with C*20 svantbody. 
c               GLOHED, LSQPRT: Replace iblk on h-file with svantbody. . 
c               UPORB:   Remove extraneous code for old-style g-files. King 141004
c version 10.50 Changes to fix problems with 141004 mods and for new c-file format: AVGLOD, 
c               BISOPT, FILOMC GETCDA, GETHED, GET_NARROWLANE, NBIAS, NBIASR,  NORMD, 
c               RESOLVE_WL, ROUND_BIAS, UPORB. King 141206
c               UPORB:  Fix failure to recognize the new GNSS g-file format.  King 150423
c version 10.51 FILOBS, FORMN2, FORMN2, GETHED, GET_WIDELANE, NORMD, OPERA, OPERA2, SOLVLC, solve.h: 
c                 Set the frequency gear ratio from the c-file frequencies; write the GNSS code 
c                 and frequencies on the q-file.  King 160114 
c               GETHED: Fix SV frequency printout so only once for non-Glonass. King 160711 
c version 10.52 FILPAR: Change reading of ATMZEN label to accomodate nzen > 99. King 170324
c               LSQD01: Replace obsolete (and buggy) code for setting 'id' for multisession; remove
cc                 equivalence of 'id' and 'work1' arrays. King 170328
c               ZENOUT: Fix bug in initiailizing arrays (maxsit not maxsat).  King 170330
c version 10.53 SOLVE: Use h-file rather than q-file names as basis of the tmp files to avoid 
c                 ambiguity across years. Pickle/King 171218 
c version 10.54 GETHED: Added gear calculation for GLonass. Herring 18319
c version 10.55 BCHECK, COVST1, FILPAR, GET_ERR_APR, GLOHED, KEYWRD, LSQERR, LWEOP (comment only)
c                 READ_BFL, SET_PARA, UPDAT3, UPORB : Change islot checks for additional solar
c                 radiation-pressure parameters. 
c               BAKCFL, Makefile: Remove now usused routine bakcfl. King 190504
c               Move many calling arguments to commons solve.h and parameters.h, put extra commons 
c                  into solve.h, and remove all references to multi-session. 109 routines. King 190524
c               NORMD Re-instated isame initialization Herring 190612
c               REMENU Re-instated 1 for session number in bias flags (backwards comparatablity) Herring 190616
c               SET_PARA: Fixed .gt. to .ge. for test of 301. Herring 190615
c               RMJUNK: Used correct variable in computing ilast for mode=2 jzen not nzen.  Herring 190615
c               NORMD: Removed extra argument from call to formn2.  Herring/King 190827 
c version 10.56 GETHEAD, GLOHED:Updated for L1/L2 Satellite PCO values in c-file.  Headers contain
c                 L1/L2 values but apriori for PCO estimate is LC (from WRTHED in model).
c                 svantdx(2,3,nsat).  Herring 200126.
c               LSQPRT: Increased h-file version from 3.0 to 3.1 for svant and antdaz, Herring 200205. 
c               GLOHED, GETHED: Increased format on ietide and evaluation of ietide for IERS20 pole
c                 tide model. Added to bit mapping for E-tide. Herring 200218.
c               LSQERR, QHEAD2: Updated 32I to 50I to allow for 35 Beidou satellites TAH 200618
c version 10.57 BCHECK, FILPAR, FXBIAS, GLOHED, GRADOUT, KEYWRD, LSQUAR, READ_BFL, SET_PARA, ZENOUT:
c               Reassign islot1 bias, atm, and grad parameters to allow for 45 SVs.  King 210619 

      RETURN               
      END














