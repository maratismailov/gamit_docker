Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego, 1995. All rights reserved.
      SUBROUTINE AVERSN (vers)

      implicit none

      integer*4 nblen
      CHARACTER*40 VERS  
      character*45 libver
      CHARACTER*10 GETMAC,MACHIN
      character*256 message
C.....

      MACHIN = GETMAC(1)
      WRITE (VERS,5) MACHIN(1:nblen(machin))
    5 format ('9.98 of 2021/01/13 19:52 UTC  (',a,')')
             
c     get library version
      CALL LVERSN(libver)

      WRITE(*,'(a)')' '
      WRITE(message,12) VERS,libver
   12 FORMAT('Started ARC, Version ',A40,'Library ver. ',a45)
      call report_stat('STATUS','ARC','aversn',' ',message,0)


C VERSION V.002    LAST LSI VERSION
C VERSION 3.10     IBM PC AT VERSION - 1/26/87
C                  MODIFIED SUBROUTINE TIMRED - G-FILE IC EPOCH MUST
C                    BE AT MULTIPLE OF 1350 SECONDS (22.5 MINUTES),
C                    OTHERWISE STOP. ADDED NVERSN,PROPER,HILITE
C VERSION 3.11    SENT TO AERO - 2/18/87
C                   TEMPORARILY REMOVED ORBITAL PARTIAL OPTION FOR
C                    AERO SERVICE , COMMENTED OUT LINES IN ROUTINE
C                    INIT (DENOTED BY CTEMP) - 3/3/87
C VERSION V3.20     BOB'S CHANGES
C                    MODIFIED INIT, ADDED ROUTINES WGS72, WGS84, MERIT,
C                    ADD REFERENCE SYSTEM TO INPUT FILE
C                    ORBIT PARTIAL OPTION RESTORED - 3/5/87
C VERSION V3.21     Rewrote subroutine TIMRED to simplify logic and
C                     allow arbitrary (i.e., not 1/64th day) IC epoch.
C                     Changed COMMON /TIMXTR/ to use INTEGER*4 JDs,
C                     eliminate FJDE and NEPOCH, rename CTUTC to TDTUTC.
C                     /TIMXTR/ changed in subroutines EPHDRD, ROTURT,
C                     SBFN, and TIMRED, and removed from REDSAT and SBOUT.
C                     Changed arguments to TAIUTC and UT1RED to INTEGER*4.
C                     Removed subroutines CTIME (TDTUTC calculated
C                     directly in TIMRED), FJLDAY (replaced by JULDAY),
C                     and CONVRT (not used). Subroutines MDYJUL and TIMINC
C                     added.  Commons /TIMAUX/ and /HEDAUX/ eliminated.
C                     Dummy routine RXREAD added, eventually to read
C                     observation start and stop from an R- or X-file.
C                     RWK 3/6/87
C VERSON V3.22      Changes according to Mark - YEHUDA 3/12/87
C                    Removed extraneous info from ADAM, SBFN, CALCOF,
C                     ERTORB, EVAL, HOFRAD, NUTRED, SBOUT, NINT
C                    Removed integer*4 NZONE,NTESS in LEGND2
C                    Added CONTYP, moved satellite loop in NAVSTR
C                    I*4 JD & IJD in TIMINC
C                    Added CONTYP to INIT
C                    Commented IAU(1976) System - PRCES
C VERSION 3.23      Fix bug in RADCON header write: declare LTCON etc.
C                     REAL*8 in WRTHED - rwk 3/27/87
C                   Add sb SHWPRT in ERTORB to printout eclipsing info
C                     only at beginning and end of eclipse. Add LAMPRT
C                     to COMMON SOLRAD in ERTORB,INIT,WRTHED.
C                   Add FJD to calls of ERTORB in SBFN and SBFN1.
C                   Remove standards write to screen in MERIT, WGS72,
C                     WGS84 -rwk 4/1/87
C VERSION 3.24      Set up general UT1 reading routines to
C                     allow bih and iris formats in UT1RED,
C                     Replace TAIUTC with version created by KF.
c version 4.1       Thursday, February 11, 1988 Kurt
c                     First version in SCCS: same as 3.24
c version 4.2       Tuesday, February 16, 1988  Kurt
c                     Read general tables in local directory.
c                     Add error-checking CHKERR to AT_SUBS
c                     Look in local directory for tables (FILOPN)
C                     Read general UT1 table (UT1RED).
C version 4.3       Generalize for variable number of satellites
C                     routines NAVSTR,IPKSAT,NAVMRG,REDSAT,WRTHED - YB 3/31/88
C                   Initialize IOPEP in NAVSTR - YB 4/5/88
C
C VERSION 4.4       NEW NAME:  ARC  -- MHM 880412
c VERSION 6.1       NEW RELEASE     -- MHM 880417
c VERSION 6.2       Generalize interpolation arrays for variable number of
C                     orbital parameters.  Replace small ARC dimpar.fti with
C                     large version used for all other models.  RWK 5/4/88
C                   Add a third radiation pressure coefficient along the SV
C                     Z-axis.  Change sign of y-axis force to +y-axis.
C                     Changes to ARC, ARCMRG, INIT, ERTORB, REVREC
C                     and SBOUT.  RWK 7/1/88
C                   Change SNAME and SATNAM in /SATLMN/ to C*16 so that
C                     input names can be 'PRN nn'.  Changes in ARC, INIT,
C                     IPKSAT, REDSAT, and WRTHED. New routine CHGNSN.   RWK 7/1/88
C                   Update SCCS and COM directories with above changes.
C                     rwk, mhm 880726.
c VERSION 6.3       Re-MADE with -SAVE option in all routines to fix problem
c                     with zero third non-grav. partial.  RWK 8/10/88
c version 6.4       Port to the alliant.  Minor changes in arcmrg, wrthed, sbfn, eval
c                     and adam to allow compilation with of D-lines.
c version 6.5       In REDSAT, pad SNAME with blanks to compare with SATNAM. MHM 880825
c                     (Previously, SNAME read from the Gfile had nulls after the
c                     last character, whereas SATNAM read from unit 5 had blanks)
c version 6.6       Fix bug in ARCMRG:  restore dropped line that reads coordinates from
c                     buffer file in no partials case.  RWK 12/15/88
C version 6.7       Recompile entire module without I*2 option.  Explicitly
c                     declare 6 header variables I*4 in WRTHED & ARCMRG.  RWK 8/14/89
c version 7.1       Official Release 7 -- MHM 890825
c
c version 7.2       Port to suns:  Kurt 900301
c                   1) Use library routines for DACOS, DASIN, SUICID, TAIUTC, DOT, CROSS
C                   2) Reorder common blocks with integers last:
c                      SOLTAB,LUNTAB,TIMEXT
C                   3) Create routine OSCRAT.
c version 7.3       Kurt 900504
c                   Start numbering Scratch file units at 40, to allow more than 10 sats.
c                   Change to ARC only.
c version 7.31
c                   10.2 change in ADAM,ARC,ARCMRG,CALCOF,EVAL,FILOPN,
c                   INIT,IPKSAT,NINT,OPNCHK,REDSAT,REVREC,SBFN,SBOUT,TIMRED,
c                   WRTHED
c                   Yehuda Bock 14 September 1990
c version 8.2
c                   Initialize IENDF in SOLRED - Kurt Feigl 901102, per Burc's request.
c version 8.3
c                   Change GM in MERIT to match CSR value (change of .0084, 2E-8)
c                     King 910313
c version 8.4       REDSAT : Fixed a 10.2? problem
c                     Yehuda Bock 4/2/91
c                   ARC:  Add 'Normal stop in ARC'    king 91/04/11
c version 8.5
c                   TIMRED:  Allow input times to be yy ddd    king 91/7/1
c version 8.6       SHWPRT:  Nicer printout of times     oral/king 92/1/6

c Version 9.1       AVERSN:  Update version number.   King 920507

c Version 9.11      FILOPN: Call TOPENS, viz OPNCHK; explicit declarations.
c                   Makefile:  Remove OPNCHK.    Bock/King 920905
c                   FILOPN, ARC:  Fix bug in above.  King 920908
c                   REDSAT : Restore padding of blanks in SNAME.  King 920921
c                   Remove unused variables: EPHDRD, INIT, NUTRED, PRCES, ROTURT
c                      UT1RED, WRTHED.   King 920925
C                   Makefile: Change ranlib order in Sun/HP version.  Feigl/King 921015
c                   AVERSN: Reconcile Apollo and SUN versions. Bock 930706
c
c Version 9.12      Add ROCK4 and ROCK42 radiation pressure model; clean up
c                      sign conventions; add SV mass to svnav.dat table, modifying
c                      common / SATLMN/.   C.S. Ho / King 930901
c                   Changes to ARC, ERTORB, INIT, MERIT, REDSAT, SBFN1, WGS72, WGS84

c Version 9.13      Reparameterize ROCK4 radiation pressure model and make 15-min
c                   (900 sec) tabular t-file available to the user by switch.
c                   Secret back door can in fact allow the user to select
c                   any tabular interval and integration interval as long as
c                   integration interval is integer multiple of tab interval.
c                   Also fixed bug in eclipsing satellite printout code and
c                   printout.  Simon McClusky / King 20 930820- 940103
c                   ARC: Add code to read tabular, and integration interval
c                        2 common list changed /INCON/ & /ADAMS/.
c                        1 new common added /ECLIPSE/.
c                   ADAM: Modified to not use power of 2 intervals.
c                         Commons changed /INCON/ & /ADAMS/
c                   CALCOF: modified common /ADAMS/
c                   ERTORB: modified parametisation of rock 4 and 42 models.
c                   EVAL: modified common /INCON/.
c                   INIT: Modified to not use power of 2 intervals, use
c                         user selected ones.
c                         Modified commons /INCON/ & /ADAMS/
c                   IPKSAT:  Modified common /SATLMN/
c                   NINT: modified commons /INCON/ & /ADAMS/
c                   REDSAT: Changed printout format. modified common /SATLMN/
c                   REVREC: modified common /INCON/
c                   SBFN: modified common /INCON/
c                   SBOUT: modified common /INCON/
c                   SHWPRT: modified printout format, and subroutine pjdhms
c                           so that fjd is not modified to pjd and passed back
c                           to sbfn1. Modified common /SATLMN/; added common /ECLIPSE/.
c                   TIMRED: removed hard-wired 1350-s tabular interval and modified
c                           output file format.
c                   WRTHED: modified commons /INCON/ and /SATLMN/ and output file format.
c                   *****
c                   ARC      Cosmetic changes.   King 940114.
c                   Declare variables explicitly:  All routines above.  King 940114/940117
c                   REDSAT: Allow use of old-style svnav.dat table.   McClusky/King 940119
c         9.14      REDSAT: Use modified library routine NSNPRN to get block number
c                           and mass.   King 940121
c                   ARC: Change warning to make 75 sec the standard stepsize.  King 940207
c                   Makefile: Change HP compiler switches from +e E1 to +U77.  Fang/King 940506
c         9.15      Makefile, INIT, IGS92 (new):  Add option for IGS (IERS92) constants.
c                        Changes to common /HARCOF/ in ARC, INIT, MERIT, SBFN, SBFN1,
c                        SCLCOF, WGS72, WGS84.    King 940525
c                   Replace arc/ROTURT by lib/ROTSNP.  Remove now unused NUTRED, NUTTAB, PRCES,
c                        SIDMAT, SIDTIM, UT1RED.  Changes to Makefile, ARC, FILOPN, SBFN, and TIMRED.
c                   Remove MONDAY and TIMINC from directory and Makefile (in /lib).  King 940526
c                   REDSAT: Remove unused variables.   King 940526
c                   ARC, EPHDRD, FILOPN, REDSAT: Simplify file opening.   King 940528
c                   Makefile:  Shift ranlib to be executed only once.  Herring/King 940622
c                   ARC: Bug fix now rock4 solar radiation pressure model can be used. McClusky 940712
c                   ERTORB: Fix bug in spherical formulation of rock4 model, errors in
c                           both acceleration and partials. McClusky 940712
c         9.16      ARC, ARCMRG, TIMRED, WRTHED: Add flags for time time (UTC/GPS) and
c                           reference frame (INERTIAL/EARTH-FIXED).  King 940713
c                   ARC, EPHDRD, REDSAT, SBFN, TIMRED, WRTHED: Change name of variable
c                           in common /timxtr/ from tdtutc to tdtoff; make altered
c                           routines lowercase.  King 940713
c                   SBFN, TIMRED: Changes units of /lib PNROT calling argument TDTGPST. King 940729
c                   TIMRED: Correct bug in converting G-file time-type to integration time-type.
c                         King 940722
c         9.17      ERTORB, SBFN, SBFN1, SHADOW: Add solid-Earth tides to accelerations.
c                      **temporarily comment out--not right yet. King 940930
c                      Change variables in ERTORB and try again.  King 941027
c                   ARC: Remove duplicate declarations; change screen message.  King 941118
c                   ERTORB: solid-Earth tide appears to be right.  Keep.   King 941122
c                   ARC: Declare character variable.   King 941130
c         9.20      ERTORB: Add new Bernese 9 parameter solar radiation pressure model Ref M.G. 1994.
c                           Change /srprok/ common: irock4 --> modrad, add variables icsnam,u,srpmod. McClusky 941216
c                           Change /coraux/ common: sbcor(3) --> sbcor(6). McClusky 941216
c                           Change /solrad/ common: change radcon(3) --> radcon(9). McClusky 941216
c                           Add /stwopi/ common: to routine. McClusky 941216
c                   SBFN: Change /coraux/ common: sbcor(3) --> sbcor(6) now contains satellite pos and vel. McClusky 941216
c                         Change /paraux/ common: dsbcor(3,10) --> dsbcor(3,15) to fit new srp partials. McClusky 941216
c                   SBFN1: Change /coraux/ common: sbcor(3) --> sbcor(6). McClusky 941216
c                          Change /paraux/ common: dsbcor(3,10) --> dsbcor(3,15) to fit new srp partials. McClusky 941216
c                   ARC: Change /coraux/ common: sbcor(3) --> sbcor(6). Change srprok common: irock4 --> modrad .
c                        Change /srprok/ common: add variables icsnam,u,srpmod. McClusky 941216
c                        Change /solrad/ common: change radcon(3) --> radcon(9). McClusky 941216
c                        Add preliminary GM earth value to declarations here instead of in KEPLR. McClusky 941216
c                        Change /paraux/ common: dsbcor(3,10) --> dsbcor(3,15) to fit new srp partials. McClusky 941216
c                   SHADOW: Change /coraux/ common: sbcor(3) to sbcor(6). McClusky 941216
c                   REDSAT: Change /srprok/ common: irock4 --> modrad, add variables icsnam,u,srpmod. McClusky 941216
c                           Modified to detect when gfile and requested radiation pressure models differ. McClusky 941216
c                           Change /solrad/ common: change radcon(3) --> radcon(9). McClusky 941216
c                   KEPLR: Add /const/ common to routine so that GM is consistant with user selected model. McClusky 941216
c                   INIT:  Add /srprok/ common to routine so that the user selcted srp model can be checked against the
c                          srp coefficients found on the gfile....
c                          Change /solrad/ common: change radcon(3) --> radcon(9). McClusky 941216
c                          Allow neq to be either 60 ( 9 param integration), or 96 (15 param integration). McClusky 941216
c                   WRTHED: Change /solrad/ common: change radcon(3) --> radcon(9). McClusky 941216
c                           Add /srprok/ common to routine. McClusky 941216
c                           Add code to write orbit parameter names in the comment line comt(2). McClusky 941216
c                   SBOUT: Modify output of positions and partials to accomodate variable length srp models. McClusky 941216
c                          Add neq to routine argument list. McClusky 941216
c                   ARCMRG: Just some commented debug added. McClusky 941216
c                   REVREC: Modify reordering of positions and partials to accomodate variable length srp models. McClusky 941216
c                           Add neq to routine argument list. McClusky 941216
c                   ADAM: Modified argument lists for sbout, and revrec, to pass in neq (number of equations). McClusky 941216
c                   REDSAT: Stop if # parameters exceeds MAXORB.   King 950106
c         9.21      ARC, SHWPRT:  Generate a one-line eclipse summary for AUTCLN.  King 950214/950217
c         9.22      SHWPRT: Correct bug in summary when end-of-eclipse outside of arc.
c                        King 950328
c                   WRTHED:  Remove debug.  King 950328
c         9.23      ERTORB:  Fix order of include and declaration.  King 950330
c         9.24      SVNAV_READ: New routine to replace NSNPRN now reads new svnav.dat format containing yaw
c                          information. Temporary: read svnav.new vice svnav.dat.  McClusky 950330.
c                   REDSAT: Call to nsnprn replaced with call to svnav_read. Also added new common called
c                          common/yawinfo/ McClusky 950330
c                   ARC: Added new common called common/yawinfo/ McClusky 950330
c                   ECLOUT: Added new common called common/yawinfo/. McClusky 950330
c                         : Added / Revised code to print eclipse/noonturn info into yaw file and archive file. McClusky 950330
c                   FILOPN: Modified to open yaw file, added iarh and iyaw variables to common/units/  McClusky 950330
c                   EPHDRD,ARC,SBFN: common/units/ had, iarh and iyaw variables added.  McClusky 950330
c                   SVNAV_READ:  Change read file to svnav.dat and move routine to lib.  King 950403
c                   FILOPN: Message and stop if not newsytle svnav.dat file.  King 950403
c          9.25     SBFN: Restructure similar to NOAA version.  King 950404
c                   IGS92: Replace MIT constants with NOAA IERS definitions--effectively equivalent.
c                      Change aultlc=>aunit in common/const/ in ARC, INIT, KEPLR, MERIT, SBFN, SBFN1
c                      SHADOW, ERTORB, WGS72, WGS84.    King 950405
c                   SBFN: Add common /timfrm/ to identify UTC integrations and do proper
c                      conversions for calling ROTSNP.   Fixes 10s sidereal time error in UTC
c                      integrations with versions 9.16ff.  King 950418
c                   ARC, IGS92, INIT, MERIT, SBFN, SBFN1, SCLCOF, Makefile:
c                      Add effects of solid-Earth tides using Jim Ray's code from IERS standards;
c                      Clean up comments and structure in SBFN1.  King 950419
c                   SBFN1, SCLCOF: Rationalize signs for gravitational harmonic accelerations--
c                      Change sign of force addition for zonals plus tesserals in SBFN1 (- => +),
c                      balanced by change of czone sign (+ => -) and tesseral scaling in SCLCOF
c                      ( - => + ).  No changes needed for partials.  King 950419.
c                   ERTORB:  Correct sign in vector formulation of tides (to be used only to
c                      check the harmonic formulation in SBFN1.   King 950419
c                   ** Merge McClusky/Tregoning changes below from 950406-950419 with King changes for 9.25 **
c                   ECLOUT: Fixed bug occuring when satellite was eclipsing at the IC time. McClusky 950406
c                   ECLOUT: SVNAV.DAT read in eclout so that sub daily changes in yaw bias status are place
c                          in yawfile. McClusky 950406
c                   REDSAT: Modified to pass yr,mo,day,hr,min into svnav_read. McClusky 950406
c                   TIMRED: Extended integration length to 11 epochs past the input start and stop times
c                           to allow interpolation of the tfile up to 1 hr prior to session start time. McClusky 950407
c                   Makefile: Add dimpar.fti dependencies for arc.o, ertorb.o, init.o.  King 940407
c                   ECLOUT: Finally put in implicit none after another bug associated with undeclared variables
c                            bit me!!! McClusky 950407
c                   TIMRED: reads '1950.0' on 1st line gfile header to determine precession model used PT 950414
c                   WRTHED: buf45 changed to buf54 to read '1950.0' from batch file PT 950414
c                   WRTHED: pass contyp in, to write out on tfile header PT 950414
c                   REDSAT: correct gfile IC's if wrong 1950 precession used PT950414
c                   SBFN: Add .true. argument (ie correct prec. now used) to rotsnp call PT 950414
c                   ROT_GFILE: new routine to rotate IC's using IAU76 precession model PT 950418
c                   REDSAT: undo changes of 950414 PT 950418
c                   ARC: call ROT_GFILE if there is no '1950.0' on gfile header. PT 950418
c                   ROT_GFILE, ARC: Add comment line to new g-file if converted; write messages
c                       to archive file as well as screen.  King 950419
c                   REDSAT: Fix problem in calling SVNAV_READ.  King 950420
c                   IGS92: Fix Sun radius problem introduced from NOAA iers.f (never in /active
c                       or stdrel.  King 950420
c                   SHADOW, KEPLR:  Change aultsc to aunit in common /const/ (no effect).  King 950420
c                   ECLOUT: Fixed another bug in eclipse summary code, added COMMON/INCON. McClusky 950421
c                   ROT_GFILE: Add precession velocity term to rotated velocities.   King 950421
c                   SIDTIM: Change sidereal time computations to IERS formulas PT 950424
c                   SBFN: Correct bug in getting body-fixed longitudes for Moon and Sun for tides.
c                       Ray/King 950425
c          9.26     ARC, EPHDRD: Allow for J2000 reference frame.  King 950426
c                   ARC,SBFN: hardwire for IAU68 to allow tests of solid earth tides only PT950429
c                   ARC, SBFN: code restored to IAU76 precession PT950429
c                   ARC,REDSAT,SHWPRT/ECLOUT: Added yaw_entries counter to yawinfo common PT950429
c                   ARC: Added number of yaw entries to 2nd header line of yaw file PT950429
c                   FILOPN: changed unit number of yaw file to 10, opened temporary yfile PT950430
c                   ARC,FILOPN: added unit number of yaw.tmp to units common PT950430
c                   ARC: add yaw_entires to yfile header, copy to correct filename PT950430
c                   ECLOUT: fix another bug in eclipse summary code!!! PT950502
c                   ARC,REDSAT,SHWPRT/ECLOUT: fixed yawinfo common PT950502
c                   ARC,WRTHED: only write IAU76 on tfile header if new precession used in arc PT950508
c                   ARC,TIMRED,ROT_GFILE,SBFN: Allow options for inertial frame and precession
c                      model through arc input file  PT950509
c                   ARC, SBFN, WRTHED: Add hidden control for turning off tides to input stream
c                      and common /models/; fix bugs in invoking IAU68 and turning off tides.
c                          King 950510
c                   ARC, FILOPN, INIT, IGS92, MERIT, WGS72, WGS84:  Add model printout to output
c                       file; print gravity model only once for all SVs.  King 950510
c                   Declare variables explicitly ARCMRG, INIT, MERIT, REDSAT, REVREC, SBOUT.  King 950511
c                   ARC,ERTORB,INIT,REDSAT: Add SRDYB rad. press model (dir., Y,B-biases - no Rock4)
c                      Tregoning 950511 / King 950524
c                   TIMRED: Change reading of frame to match FIXDRV.  King 950512
c                   ERTORB:  Fix compile bug.   King 950512
c                   ARC:  add comment to yfile header. Tregoning 950512
c          9.27     ARC, SBFN, WRTHED: Add gravitational field common /models/ and both this and radiation-
c                      pressure model to T-file header.  King 950517
c                   BINOML,HMS,HOFRAD,KEPLR,LEGND2,LEGNDR,LUNRED,MATMPY,MATTRN,MDYJUL,SCLCOF
c                     SHADOW,SOLRED: declare integer variables Tregoning 950605
c                   ECLOUT: fixed bug in writing yfile. Tregoning 950605
c                   AVERSN:  Remove call to CLRSCR.   King 950606
c                   ARC, ROT_GFILE:  Remove extra commas in format.  King 950606
c                   SHWPRT:  Change idint(jd) to jd to satisfy HP compiler.  King 950606
c                   ARC: Change srdyx rad. press.  model to srdyb  Tregoning 950606
c                   REDSAT: Fixed bug in srdyb model  Tregoning 950606
c                   FILOPN: Change svnav.dat unit number to 69 to match svnav_read.f, closed
c                             svnav.dat, changed status of yaw.tmp to scratch  Tregoning 950606
c                   TIMRED: Changed length of char variable to match character read Tregoning 950606
c                   ROT_GFILE: Remove extra comma in format statement.  King 950609
c                   ARC, ERTOB, REDSAT: Fix naming and documentation for radiation-pressure models.
c                       Temporarily add suicid call for undefined constant in body-axis (SRXYZ, modrok=2)
c                       model.  King 950612
c                   Move definition of T-file limits (maxytp,maxyt2) to lib/dimpar.fti:  ADAM, ARC, ARCMRG,
c                      CALCOF, EVAL, INIT, REDSAT, REVREC, SBFN, SBOUT, TIMRED, WRTHED.  King 950612
c                   ADAM, ARC, INIT, NINT, TIMRED: change name of initial integration values
c                       from V0 to satics in common/stint/ to match other routines.  Change
c                       diminsions in ARC from fixed 66 to variable maxyt2.  King 950612
c                   ARC, INIT, REDSAT, GHDRED (new),WRTHED,  Makefile:  Move reading of 2d g-file
c                       record from REDSAT to TIMRED, renamed GHDRED.  King 950613
c          9.28     ARC, EPHDRD:  Allow and check J2000 lunar/solar ephemeris.  King 950616/950620
c                   GHDRED:  Fix format for G-file read.   King 950619
c          9.29     ERTORB: Simplify constants and other code for radiation pressure models.
c                        Remove commons /ertaux and /srprok, containing many now unused
c                        variables; consolidate remaining ones in common /solrad--changes to
c                         in ARC, ERTORB, and SBFN, and
c                        'ltcon', 'ltcon2', 'ltconz', and 'sbarea' from common /solrad in ARC,
c                        ERTORB, GHDRED, INIT, REDSAT, SBFN, and WRTHED.  King 950621.
c                   FILOPN: Removed naming of yaw scratch file (tmp.yaw) to make HP happy. McClusky 950701
c                   ROT_GFILE: Removed naming of gfile scratch file (gfile.tmp) to make HP happy. McClusky 950701
c                   ARC, ERTORB, LUNRED, ROT_GFILE, SOLRED, WRTHED:  Remove duplicate declarations.   King 950707
c                   Makefile, ARC:  Remove GHDRED, ROT_GFILE to library and Modify calls.  Tregoning 950725
c                   ARC: add arguments to calls to GHDRED,ROT_GFILE   Tregoning 950727
c                   LUNRED, SOLRED: Remove local misspelled alternate of common 'iendfm' and 'iendfs'.
c                        King 950728
c                   ARC, ADAM: Change dimensions of /stint/ and /nordcf/ variables from
c                        66 to 'maxyt2'.  King 950729
c                   ARC, DTWOPI, EVRTCF: Convert the latter two from block data to subroutines and assign
c                        the constants (to accommodate the DEC).  Sanli/King 950801
c                   ERTORB: Save additional variables; multiply each term by lambda instead
c                        of constant so that some terms can be kept on during eclipse. King 950803
c          9.30     ARC:  Fix minor bug in test for radiation model.  King 950807
c          9.31     ERTORB: Remove ROCK4 model from all but the SRDYZ and SRXYZ models, and remove
c                         shadowing effects from all terms except direct solar and the equivalent
c                         ROCK4 accelerations. Tregoning/King 950909
c          9.32     ERTORB: Minor correction to Blk I 'radforce' (4.545=>4.54.  McClusky/King 950920
c                   IGS92:  0.25 ppb change in GM to geocentric value (.4418=.4415).  Ray/King 950920
c                   REDSAT: Trap zero ICs from g-file.  King 951005
c          9.33     Makefile:  Variable ranlib to make compatible with Solaris 2.  Fang/King 951207
c          9.34     EPHDRD: Minor change to ensure the correct reference frame for the
c                         solar and lunar ephem tables is obtained. McClusky 960112
c          9.35     Modified many routines to allow better error handeling using report_stat:
c                   arc, aversn, arcmrg, ephdrd, filopn, init, ipksat, mdyjul, oscrat, redsat,
c                   shwprt and,  wrthed.f. McClusky 960223
c                   SBFN:  Add comment (only) about leap seconds.  King 960517
c                   REDSAT: Add PRN name to message for missing satellite.  King 960705 
c           9.36    Changes for new kf/gamit structure and Makefiles from unimake.  King 960718
c                        All routines:  removed trailing blanks and replaced lib/*.fti 
c                        includes by ../includes/*.h.  
c                   KEPLR, MATMPY, MATTRN: Declare variables explicitly.   King 960815
c           9.37    ARC: Stop if GAMIT.fatal exists.   King 960828
c           9.38    Rename NINT=> START_INT to avoid conflict with intrinsic; changes
c                      also in ADAM, ARC, SBFN.  Tregoning/King 961017   
c                   SBFN:  Add yaw info to common /units/.   Tregoning/King 961017
c                   SHWPRT: Fix bug in eclipsing array logic.   Tregoning  961210    
c           9.39    SBFN1: rewrite some logic to allow -O3 optimisation using g77.
c                       Tregoning 961022/King 970102
c                   ARC: Record the y-file name in the arcout file.  King 970103   
c                   ARCMRG: Clean up logic to satisfy Gnu compiler.  King 970109 
c           9.40    Changes for G77 compiler and make all routines implicit none:
c                      ARCMRG, BINOML, HMS, HOFRAD, IGS92, LUNRED, SCLCOF, SHADOW, 
c                      SOLRED,  WRTHED.  McClusky  970110.
c           9.41    Fix eclipse print problem with simpler code.  King 970228  
c           9.42    Modified ARC to allow the creation of a T-file containing velocities, and 
c                      partial derivatives of orbit parameters wrt velocity (use "V" in partials
c                      line of arc input file): Main routines modified were ARCMRG, INIT, SBOUT, 
c                      REVREC, and WRTHED. Because the incon comon block was modified the following 
c                      routines were also changed: ARC, ADAM, EVAL, SBFN, SHWPRT ,START_INT. 
c                      McClusky 960424/King 970522  
c                   EVAL: Write interval check message to GAMIT.fatal, not unit 6.  King 970707
c                   SHWPRT: Correct misformatted report_stat call.  King 970707
c           9.43    ERTORB:  Add Rock4 (T30) accelerations for Block IIR satellites and
c                     change T20 code to check on iblock=2 (Blk II) or 3 (Blk IIA).  King 970806
c                   REDSAT: Change svnav.dat failure to fatal and revise message.  King 970808
c           9.44    INIT:  Fix call to report_stat.  King 980227
c           9.45    FILOPN, REDSAT, SHWPRT: Add sv antenna offsets to call of lib/svnav_read.  King 980603
c           9.46    ARC, ERTORB, INIT, REDSAT, WRTHED:  Add CODE 9801 radiation pressure model (BERN2); 
c                      change variable names and values for units for all models.  King 980701/980716
c                   ARC: Convert rad model to uppercase.  King 980708/10
c                   REDSAT: Enhance display of svnav.dat info in arcout file.  King 980709
c                   ERTORB: Fix duplicate declaration.  King 980722
c           9.47    ERTORB: Change sign of B-axis back to original so that BERNE, BERN1, and BERN2
c                      are all the same.  King 980805
c           9.48    ARC, INIT, ERTORB, WRTHED: Switch modrad values (3,4) for BERNE and SRDYB for 
c                      convenience in coding the next change.  King 9810016
c                   ERTORB: Fix bug of omitted adjustable  terms for BERNE.  King 981016
c                   ERTORB: Fix bug in partials for BERNE.  King 981026
c                   Remove implicit statement to statisfy SGI compiler, and add implicit none;
c                     all routines compiled ok previously with explicit compiler switch:
c                     ADAM, EVAL, REVREC, START_INT, CALCOF, MERIT, SBOUT, SBFN1, WGS72, 
c                     WGS84:  Morgan/King 981231   
c                   WRTHED: Correct format statement number in writing arcout header modrad=3,4
c                     (SRDYB, BERNE), and correct description of SRDYB.  King 990111
c                   ERTORB: Do not double '+' in arithmetic expression (IBM and f90 restriction).  Fang/King 990324
c           9.49    SBFN, EPHDRD: Add 'iyawtmp' to common /units/ (alread in ARC, FILOPN.  McClusky 990813
c           9.50    SHWPRT, WRTHED, remove MDYJUL, Makefile:  Changes for 2-digit years.  King 990727  
c           9.51    SHWPRT: Trap case of start but no end of eclipse when cause is integration
c                     startup; fixes bogus entries in arcout and y-file.  King 000328
c           9.52    ERTORB: Remove extraneous + in modrad=5 terms.  Bosy/King 001208
c                   AVERSN, ARC: Don't print copyright to screen; use report_stat for integration message.  King 001220
c           9.53    ERTORB: Comment out extraneous calculation of u.  King 010103
c           9.54    WGS72  WGS84  IGS92  SCLCOF  MERIT  ERTORB  SBFN  SBFN1  INIT  EGM96: Increased gravity field to 
c                   10x10 field for the EGM96 (IERS2000) field.  Introduced time dependent garvity coefficients for
c                   second degree harmonic terms in EGM96. Added PEP MJD of IC epoch (jde) to call to init. Herring 030304.
c           9.55    EGM96: Shorten in-line comment to avoid ifc warning.  King 031027
c                   ARC: Fix typo in string assignment.  King 031027
c                   WRTHED, Makefile: Remove sb hms.f in favor of lib/ds2hms.f to fix round-off problem.  King 031111
c           9.56    ARC, EGM96, ERTORB, IGS92, INIT, KEPLR, MERIT, SBFN, SBFN1, SHADOW, WGS73, WGS84: Add
c                      lunar radius to common /const/.  King 040511
c                   SHADOW, SHWPRT: Check for lunar as well as terrestrial eclipses.  King 040511
c                   ARC, EPHDRD, FILOPN, INIT, REDSAT, SHWPRT, WRTHED  : Use variable name for units.   King  050103
c                   SHWPRT: Mark lunar eclipses 'M' rather than 'E' in y-file (already in print file). King 050103
c           9.57    ARC, ECLOUT, ERTORB SHADOW, SHWPRT: Save eclipse type for summary, fix bugs in summary times.  King 050223
c                   ARC, SHADOW: Add hidden flag to turn off lunar eclipses.  King 050223
c           9.58    REDSAT: Enlarge character string for block number and allow iblk=5 from svnav.dat.  King 050706
c                   ERTORB: Set radiation pressure scaling for Block IIR-B (iblk=5) to be same as IIR-A iblk=4)
c                     for most models; trap for missing BERN1 BERN2 terms.(iblk=5). King 050707
c                   ERTORB, REDSAT: Allow iblk=6 for BLOCK IIR-M.  King 051027
c                   ERTORB: Save 'ecltyp' version.  King 051107
c           9.59    ARC, ARCMRG, ECLOUT, EGM96, EPHDRD, ERTORB, GET_LAMBDA, FILOPN, INIT, REDSAT, SBFN, SHWPRT, 
c                      WRTRED:  Open a file for debug and add its unit number to common/units/.  King 051216.
c           9.60    ARC, FILOPN: Remove unused variables.  King 070130
c           9.61    ARC: Remove unneeded arguments from call to rot_gfile.  King 070416
c           9.62    ARC, OSCRAT: Assign unique names to scratch files for cluster use.  King 070709
c                   SBFN: Remove unused 'iscrn' from call to /lib/rotsnp.  King 070906
c                   SBOUT: Add debug code for unused calling argument. King 070906     
c                   ARC: Fix obsolete loop structure.  King 070906       
c                   ADAM, EPHDRD, ERTORB, EVAL, FILOPN, : Remove unused statement labels. King 070910 
c                   ARC, FILOPN: Remove unused calling arguments.  King 070910
c                   EPHDRD: Change message reporting error opening luntab and soltab.  King 071024
c           9.63    ARC, SBFN, SHADOW, WRTHED: Add nutation model to batch file, g-file, and t-file;
c                       gravity field to g-file.  King 071221   
c                   ARC, ECLOUT, ERTORB, FILOPN, SHWPRT: Write beta angle to print file and yaw file.  King 080129
c                   EGM96: Minor change in print format.  King 080130
c           9.64    ARC: Set 'nics' before ghdred call (undefined at this point).  King 090417    
c           9.65    ERTORB, REDSAT: Allow iblk=7 for BLOCK II-F.  King 100529
c           9.66    ERTORB: Change radiation pressure constants for BLOCK II-F so that the 
c                     direct coefficient is closer to 1.0 (acceleration is 2.13 times IIR). Herring/King 101014
c           9.67    Makefile, ARC, EGM96, EGM08(new), IGS92, INIT, MERIT, SBFN1, SBFN, SCLCOF, WGS72, WGS84: Increase 
c                     dimensions  of gravity field in common /harcof and add EGM 2008. King 101222
c           9.68    Put all commons in new includes/arc.h (changes to almost all routines; reconcile variable
c                      names).  Remove merit.f, wgs72.f, and wgs84.f from Makefile.  Change name of sb 
c                      ipksat to get_sats.  King 101230/110114
c                   ARC, EFLUX(new), EARTHRAD(new), INIT, SATPROP (new), WRTHED, Makefile: Add new  routines 
c                      and code to handle antenna thrust and Earth radiation.   Petrie 101005/110114
c                   COEFFS(new), SRPFMOD(new), Makefile:  New routines to use U College London (UCL) fourier
c                      coefficients for radiation-pressure model.  Petrie 101229                        
c                   Merge King and Petrie changes.  King/Petrie 110114
c                   ARC, WRTHED: Assign 'UCLR1' to modrad=7; temporarly write 'BERNE' on t-file to satisfy the 
c                      rest of GAMIT.  King/Petrie
c           9.69    ERTORB: Fix bug in calling arguments for SHWPRT (serious error in yaw).  King 110523   
c                   ARC: Change numerical unit number '9' to 'iyawtmp' (no effect). King 110523
C           9.70    ARC, ERTORB, SRPMOD, : Updated with latest experiments from the arckp version: ertorb.f, srpmod.f. Petrie 
c           9.71    Makefile, ARC, EARTHRAD, ERTORB, INIT, WRTHED, BILEN8(new), BINOML(new) : Added UCLR2 SRP model and preparing 
c                      for yaw as part of arc. Also calling ../model/satatt.f Petrie 120404
c           9.72    Makefile, ARC, ARCMRG, ERTORB, FILOPN, INIT, REDSAT, SBFN, WRTHED; new routines CHECK_GMODELS,
c                      GET_SAT_INFO, READ_ARC_BATCH, SET_MODEL_CONTROLS, WRITE_ARCOUTHD;  remove GETSAT, REDSAT; 
c                      Simplify logic and variable naming, more stuff in common, add Earth radiation and antenna 
c                      radiation model names, change hidden controls for turning off tides and lunar eclipses from 
c                      the model line of the batch file to the binary-coded debug bit on the print-file line. 
c                   ARC, CHECK_GMODELS: Fix the logic for setting the radiation-pressure parameter names when
c                      the g-file model doesn't match the requested model.  King 120531
c                   CHECK_GMODELS: Fix model name in warning. Petire 120606
c                   GRIDVAL: Fix check on bad lat/lon.  Petrie 120606  
c                   ERTORB: Fix check on antenna radiation and yaw comments. Petrie 120608
c                   EGM08: Use JD start instead of undefined JD epoch for harmonics and mean pole. King 120610
c                   ARC, SET_MODEL_CONTROLS: Move call to dtwopi earlier.  King 120610
c                   EGM96: Fix bug in undefined epoch for 2nd order gravity coefficients. King 120730
c                   SATPROP: Add information for block IIF (approximate only) Petrie 120905
c                   ERTORB: Add a block-specific antenna radiation model. Petrie 120906
c                   EARTHRAD, EFLUX, SBFN: Changes to comments and commented-out debug. Petrie 120906
c                   ERTORB, Makefile.generic, new routines EARTHRADTUM, ERPFBOXW, PROPBOXW, SURFBOXW. Petrie 120910
c                   ERTORB: Fix bug in implementing TUM ERad model.  King 130327   
c             9.73  RDSRPGRD: Remove extraneous commas. King 131011
c                   CHECK_GMODELS: Clarify wording of model differences.  King 140228
c             9.74  SBFN, SBFN1, SET_MODEL_CONTROLS, GENERALREL (new), arc.h, Makefile:  Add general relativity
c                     accelerations, allow gravity model 'EGR08' for testing these.  McClusky/King 140402
c             9.80  Major mods to change svnav.dat format and to support GNSS:  Makefile, EARTHRAD, ECLOUT, 
c                     ERTORB. FILOPN, GET_SAT_INFO, GPSBLOCK (new, temporary), READ_ARC_BATCH, SATPROP, SRPFGRID,
c                     SRPFMOD. King 150107 
c                   Move GPSBLOCK to /lib.  King 150316
c             9.81  ECLOUT, GET_SAT_INFO: Add start/time times to svnav_read call.  King 150520 
c                   ECLOUT: Fix missing 'gnss' argument in last call to svnav_read.  Herring 150529
c                   ERTORB: Fix string equate length and "BLoCK IIF" typo to uppercase O.  Floyd 150725
c             9.82  ARC: Add satellite name to report_stat status. King 150902
c                   GET_SAT_INFO: Read svnav.dat with the IC epoch rather than the start time.  King 150902
c                   ECLOUT: Skip call to svnav_read (not needed) to avoid PN change over day boundary. King 150902
c                   ERTORB: Add temporary scaling for Glonass.  King 150919
c             9.83  ERTORB,INIT, WRTHED, SET_MODEL_CONTROLS: Add ECOM2 (alias BERN2) radiation pressure model;
c                     remove all models except BERNE(ECOM1), NCLR1, and NCLR2.  King 151017
c                   ERTORB, EARTHRAD, EFLUX, SATPROP, SRPFMOD, SRPFGRID: Replace 'iblock' with 'antbody' throughout 
c                     code; replace check of PRN 25 by check of SV 62 for antenna radiation thrust, and add Glonass. King 151019
c                   ERTORB: Limit thrust-model warnings to one per SV.  King 160107
c             9.84  GENERALREL: Fix two errors J2 relativity terms.  King 160201
c             9.85  SATPROP: Change warning to fatal for unrecongized spacecraft. King 161005
c                   READ_ARC_BATCH: Allow reading of old-style satnam.  King 161005
c                   ERTORB: Add test on idbprt>0 for writing 'earth accel' to debug print file. King 161017
c             9.86  SBFN, SBFN1, SCLCOF1(new), SET_MODEL_CONTROLS, READ_OTIDES, DOODSON_ANGLE, Makefile: Add new 
c                     code for solid-Earth and ocean tidal accelerations.  King 170420
c             9.87  READ_ARC_BATCH, SBFN, SET_MODEL_CONTROLS, includes/arc.h: Allow selection of the harmonic degree 
c                     for static gravity, solid-Earth tides, and ocean tides via integers in the batch file. King 170912 
c                   FILOPN: Open the ocean tide table only if ocean tides are called for; add a message for a 
c                     successful open. King 170925
c                   READ_ARC_BATCH.f, Makefile: New routine to print gravity coefficients.  Herring 180229
c                   EGM08: Change C20 to SLR value for tide-free adjustment, and add time-variable terms
c                      for C30 and C40.  Herring 170228
c             9.88  ARC, FILOPN, SBFN, arc.h: Read a PEP or JPL planetary ephemeris instead of soltab. and luntab. King 180301 
c                   EPHDRD: Fixed assignment of unit numbers for luntab. and soltab.  King 180319
c                   ARC, ECLOUT, EGM96,EGM08, INIT, OSCRAT, READ_ARC_BATCH, WRTHED: Change start-time variable names
c                      here and in includes/arc.h from jd0,t0 to jdb,tb to be less ambiguous and consistent
c                      with other modules. King 180320 
c                   ARC, ARCMRG, CHECK_GMODELS, EATHRADTUM, FILOPN, SBFN: Add new common includes/units.h 
c                      to provide global unit numbers removed from includes/arc.h.  King 180320
c                   CHECK_GMODELS, EARTRADTUM, ECLOUT, GET_SAT_INFO, READ_ARC_BATCH, SBFN, WRITE_ARCOUTHD,
C                      WRTHED: Add new common includes/global.h with revised arc.h . King 180321
c                   Makefile: Move EVRTCF from /arc to /lib. King 180324
c                   CHECK_GMODELS, FILOPN: Open the nutation file only if it is needed, King 180328
c                   ARC, SBFN, SBFN1: Remove data/neweph/ and test on the existence of 'nbody'
c                      for whether to use soltab/luntab or nbody.  King 180328
c             9.89  ERTORB: Add radiation scaling for Glonass, Beidou, and Galileo and change the
c                      default to keep the area/mass similar to Block IIF rather than just the area. King 180405
c                   ERTORB: Correct typo in effective area value for Glonass K-1.  King 180428
c             9.90  ARC, READ_ARC_BATCH: Add flag for Jupiter and Venus to the batch file.  Herring/King 181205
c             9.91  ERTORB: Add GPS 'BLOCK IIIA' for radiation pressure scaling and antenna power,
c                     with preliminary values. King 190109 
c                   EARTHRADTUM, SATPROP: Temporarily set BLOCK IIIA body values to be the same as BLOCK IIF. King 190129
c             9.92  ERTORB, GET_SAT_INFO, INIT, SBFN, SBFN1, SET_MODEL_CONTROLS, WRTHED includes/dimpar.h, 
c                     includes/arc.h: Allow 13 adjustable SRP coefficients for ECOM2. . King 190430
c                   ERTORB Replace u argument with u-u0 (uu0) in the BERNE OPR terms and in the partials
c                     This change corrects the partials error in the ECOMC OPR terms and make BERNE and ECOMC
c                     the same for OPR terms.  Pre-9.92 g-files will not integrate correctly with this change.
c                     Added distfct scaling to 2U and 4u partials. Herring 190615.
c                   WRTHED: Write overflow IC names into the 3rd line of the t-file comment. King 190618
c                   READ_ARC_BATCH: Cleaned up reading lbody code and reporting Herring 190622
c                   GET_SAT_INFO, READ_ARC_BATCH: Added antpwr to svnav_read calls Herring 190702
c                   ERTORB, ANTPWR_BASE: New routine to return antpwr based on satellite type (removed special
c                     case for G062 after 2011/04/05; not in igs_metadata.snx.  Herring 190702.
c                   EARTHRADTUM,PROPBOXW.f updated for GLONASS, GLO-M and GLO-K satellites Herring 190702.
c             9.93  ARC.H, GET_SAT_INFO, ERPFBOXW: Added SVN (satsvn) to arc.h in get_sat_info for use
c                     in Galileo albedo model mass.
c                   ERPFBOXW, PROPBOXW, SURFBOXW: Added new arrays to handle Galileo albedo but only for
c                     EARTH RADIATION PRESSURE (ANALYTICAL)=ERM=1 model.  Herring 190722.
c             9.94  CHECK_GMODELS, ERTORB READ_ARC_BATCH: Added code to allow SOFA IAU inertial frame integration.
c                   GET_SAT_IBFO: Remove used declaration Herring 200205.   
c             9.95  EARTHRADTUM: Modified to warn if not Albedo model rather than fatal.  Allows same sestbl.
c                     for all multi-gnss constellations. Herring 200225.
c                   READ_ARC_BATCH: Modified reading of arc output file to be more flexible in handling
c                     (non-documented) debug and yaw flags.  Herring 200331.
c             9.96  READ_OTIDES: Fixed logic in reading files with more than 18 frequencies (FES2014b) and
c                      added output of model name to GAMIT.status, cleaned up debug output. Herring 200503
c                   SATPROP: Added test so that BlockIIIA not coded warning is only printed once. Herring 200504
c                   EARTHRADTUM, SBFN: Replaced arbitrary setting of iut1pol bit mapped diurnal/semidiurnal 
c                      models with value read from sestbl. TAH 200505
c             9.97  EARTHRADTUM: Updated to handle GLONASS M+, K1A and K1B satellites (change to 
c                      name_to_blk changes mapping to BlkNum used in routine. TAH 200603.
c             9.98  ERTORB: Added "BEIDOU-3I" SV body type to block defining
c                   d0 for "BEIDOU-2I" and "BEIDOU-2G". MAF 20210113
c

      RETURN                                                                                     
              
      END                                                                                                                      



