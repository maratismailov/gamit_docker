Copyright (c) Massachusetts Institute of Technology and the University of
California at San Diego, 1995.  All rights reserved.

      SUBROUTINE CVERSN(iprnt)

      implicit none

      CHARACTER*10 GETMAC,MACHIN
      CHARACTER*40 VERS
      character*45 libver
      character*256 message

      integer*4 iprnt,nblen

      MACHIN = GETMAC(1)
      WRITE (VERS,5) MACHIN(1:nblen(machin))
    5 format ('9.64 of 2021/06/19 13:25 UTC (',a,')')
c     The following is old for sccs
c   5 format ('%I% of %E% %U% (',a,')')
c     The above line contains magic SCCS keywords which will be
c     automatically updated by delta and get.

c     get library version
      CALL LVERSN(libver)

c     write (iscrn,'(/)')
      WRITE (message,10)vers,libver
10    format('Started CFMRG ver. ',a40,' Library ver. ',a45)
      call report_stat('STATUS','CFMRG','cversn',' ',message,0)
c     rwk 061112: no longer print this to screen (status msg enough)
c      WRITE (iprnt,*) 'CFMRG v. '//vers
c      WRITE (iprnt,*) 'Library v. '//libver
c      CALL PROPER(iscrn)
      CALL PROPER(iprnt)
C
C VERSION V01.0001 LSI VERSION
C VERSION V02.0001 IBM PC AT VERSION - BATCH
C VERSION V03.0001 IBM PC AT VERSION - INTERACTIVE , 1/23/87
C VERSION V3.01    APOLLO COMPATIBLE, dimpar.fti VARIABLE PARAMETERS
C                    OPNCHK REMOVED, RECL INQUIRED  5/5/87
C VERSION 3.22  Use new C-file format - rwk & mhm 6/26/87
c version 4.1   Put into SCCS and handle 20 stations
c version 4.2   Corrected explicit E reading
c version 6.2   Modified to fit multi-session mode. DND 880726
c version 7.1   Structured fortran version.  MHM 890817
c                Added APRIOR,FILLS1,FILLS2,GETCFS,PROPER,QUERYB,QUERYP,
c                and REORDR subroutines
c                Integer*4, Case Insensitive, READ/WRITM1-3,MOPENS Primitives
c         7.2   MAXSAT = 12
c version 7.3   4 characters for site code
c version 7.4   Delete PROPER; use libary version.
c               Update this routine to handle machine type.
c version 8.1   New C-file format. 901112 Kurt
c version 8.2   Get multi-session correct for varying no. of satellites
c                (CFMRG,APRIOR,FILLS1,FILLS2)
c               Note CFMRG input batch file has changed
c               Yehuda Bock 910929
c               New C-file format:   Call to READC2 in GETCFS.
c                  King 920423/920430
c version 9.1   Copy to SIO and update version.  King/Bock 920502/4.
c version 9.11  Correct MAXSIT => MAXNET in MIT and SIO Sun versions
c                 APRIOR, CFMRG, FILLS1, FILLS2.  Bock/King 920905.
c               Remove unused variables:  CFMRG, FILLS1, FILLS2.  King 920925
c               CVERSN: Reconcile Apollo & SUN versions Bock 930706
c version 9.12  FILLS1:  Add '0' to single-digit PRN number for label.  Dong/King 930708
c               FILLS1:  Remove % from include statement Bock 930717
c version 9.13  CFMRG,QUERYP,FILLS1,FILLS2:
c               Add multiple zenith delay parameter array NZEN
c               Add MAXATM=24 to dimpar.fti (i.e, a zenith delay
c                 per hour for 24 hour session).
c                 Bock  030719
c               All routines: remove % from include statements
c                 Bock 930807
c version 9.14  Remove interactive prompts and divert screen output
c                 to file 'cfmrg.out':  CFMRG, GETCFS, QUERYB, QUERYP
c                 King 940107/940112
c               Makefile: Change HP compiler switches from+e +E1 to +U77.  Fang/King 940506
c version 9.15  Add earth orientation parameters/partials
c               APRIORI: add apreop(6) to MERCOM common, and fill array with
c                        earth orientation parameter prevals. McCLUSKY 940429
c               CFMRG: Add apreop(6) to MERCOM common, add variable ieochk
c                      to queryp call argument list. Indicates if eop partials
c                      present on the cfiles. McCLUSKY 940429
c               FILLS1: Assign eop parameters islot1 values of 1901 --> 1906,
c                       bin numbers 24 --> 29 and add X POLE YPOLE XPOLE RATE
C                       Y POLE RATE UT1-TAI AND UT1-TAI RATE to rlabel. Write
c                       ibin, prevals and rlabels into correct positions in
c                       islot1, aprval, alabel arrays. McCLUSKY 940429
c               FILLS2: Add earth orientation parameters to islot2 array. McCLUSKY 940429
c               QUERYP: Check to see if eop partials available. (only available
c                       if orbit partials available). McCLUSKY 940429
c               QUERYP: Fix bug in check for presence of orbit/eop partials.  King 940516
c version 9.16  CFMRG, APRIOR, FILLS1: Add a priori values for zenith delays.  King 940609
c               Makefile: Move ranlib to be executed only once.  Herring/King 940622
c version 9.17  CFMRG: Initialize common arrays with assignments, not data statement,
c                       and fix formats, to avoid RISC-machine compiler errors.  King 94120
c               GETCFS:  Fix formats for XL compiler on RISC machine.  King 950103
c               REORDR:  Minor compile warning.  King 950403
c               Declare variables explicitly:  CFMRG, CVERSN, FILLS1, FILLS2, QUERYP. King
c version 9.18  APRIOR, GETCFS, FILLS1:  Add 6 additional orbit parameters.    King 950515/
c               CFMRG, GETCFS, FILLS1: Pass C-file labels and slots to FILLS1.  King 950517
c               APRIOR: Temporary code get shift zenith and clock a priori values.  King 95

c version 9.20  New C- and M-file formats:  CFMRG, GETCFS, FILLS1, FILLS2, QUERYP.  King 95
c               Make doubly dimensioned variables except nfile and fname MAXSIT x MAXTFL ra
c                   than MAXNET X MAXTFL since only these two include all stations in the n
c                   Others (e.g., cfname, lstawg,t) correspond to C-file session records an
c                   only the stations in the session.  Caution: The latter may have been a
c                   move to overcome a coding problem  in the past.  Changes to CFMRG, REOR
c                   King 950526.
c               Makefile:  Name all .a files the same.   King 950608
c version 9.21  CFMRG: Add decimation factor to record 2 of M-file.  King 950616
c version 9.22  GETCFS:  Call to READC2 for additional variable (offarp).  King 950717
c               FILLS2:  Remove multiple declaration.  Sanli/King 950728
c version 9.23  FILLS1:  Add trap for zero orbital elements.  King 951005
c               FILLS1:  Fix bug in above trap.  King 951026
c version 9.24  CFMRG, APRIOR, GETCFS:  Fix calculation of number of orbital parameters.  K
c version 9.25  QUERYB:  Fixed minor format problem for DEC.  Sanli/King 951111
c version 9.26  CFMRG: Fix bug in writing station elevation cutoffs to M-file header.  King
c               GETCFS: Fix bug causing stop if no orbit partials.  Calais/King 951129
c version 9.27  Makefile:  Variable ranlib to make compatible with Solaris 2.  Fang/King 95
c version 9.28  Added report_stat calls to all routines giving better error handling: APRIOR
c                  CFMRG, CVERSN, FILLS1, FIllS2, GETCFS, QUERYB, QUERYP, and REORDR. McClusky
c version 9.29  Changes for new kf/gamit structure and Makefiles from unimake.  King 960718
c                   All routines:  removed trailing blanks and replaced lib/*.fti 
c                   includes by ../includes/*.h. 
c version 9.30  CFMRG: Skip executation if GAMIT.fatal exists.  King 960828 
c version 9.31  Changes to APRIOR, CFMRG, FILLS1, FILLS2, for new atmospheric
c               gradient parameters. Modified the common block MERCOM, adding
c               the real*8 variable aprgrd(maxtfl,maxnet,2). McClusky 960926 
c version 9.32  FILLS1:  Trap and fix mixed single/dual-freq C-files; fix
c                  error in report_stat call for mtpart check.   King 961123
c version 9.33  G77 compiler changes:
c                 QUERYB,FILLS1,APRIOR: Initialised previously uninitialised variables,
c                 removed unused variables, and explicitly typed all unexplicitly typed 
c                 variables. McClusky 970110   
c               Fix report_stat call for duplicate c-file names.  King 970306
c version 9.34  QUERYP: Fix report_stat call for 
c version 9.35  GETCFS: Add norb=12 case (BERN2 model).   King 980708
c version 9.36  APRIOR, CFMRG, FILLS1, FILLS2, QUERYP:  Remove satellite clocks, use labels
c                 from C-files when available, add SV antenna offsets. King 980905
c               GETCFS: Add 'norb' to calling argument for readc2; comment out
c                 local calculation of 'norb'.  King 980911    
c version 9.40  CFMRG, FILLS1, FILLS2, QUERYP:  Allow multiple gradient parameters.  King 980929
c               QUERYP:  Fix bug in dimensioning of ngrad.  Herring/King 981021
c version 9.41  CHECK_LIMITS:  Fix bug in report_stat call.  King 990602
c version 9.42  FILLS1: Fix bug to make L2 biases always estimated when mixed dual/single
c                   frequency observations.  King 990705
c version 9.43  CFMRG, FILLS1, GETCFS: Define L2 biases even when all L1-obs in order to 
c                   make SOLVE controls more flexible.  King 990907   
c version 9.44  Remove superfluous writing of sites to status file.  King 010620  
c version 9.45  GETCFS: Change format for SV list from 18 to 30.  King 021204
c version 9.46  QUERYB: Replace function upperc to subroutine uppers.  King 031027
c               CFMRG: Fix length of 'iscr' in assinging file names.  King 031027       
c version 9.47  FILLS1, FILLS2: Include average zenith delay in slot 4 along with 
c                 multiple zenith delays in slot 28.  King 040302
c version 9.48  GETCFS: New C-file format.  King 050201
c               CFMRG, FILLS1: Make c-files lowercase, site names uppercase. King 050214
c version 9.50  FILLS1: Fix bug in counting zenith delay parameters when numzen=1.  King 060206
c               CVERSN: No longer print to screen (status msg enough). King 061112
c version 9.51  FILLS1: Change slots, removing implicit (no longer supported) and
c                appropriating these to expand L1 and L2 explicit to 45 SVs.  King 061128
c               CVERSN: Removed extra blank lines written to screen. King 070124
c version 9.52  CVERSN: Remove call to clrscr; remove 'iun' argument for lversn..  King 070416
c               APRIOR, CFMRG, CVERSN, FILLS1, FILLS2, QUERYB, QUERYP, REORDR: Remove unused calling arguments.  King 070906   
c version 9.53  QUERYP: Add trap for #sites times #zenith-delays < 2500 (limit of 'solve' pointers). King 100129
c version 9.54  GETCFS: Arguments added to readc2.  King 100827      
c version 9.55  GETCFS: Arguments added to readc2.  King 130115
c version 9.56  GETCFS: Arguments added to readc2.  King 140327
c               GETCFS: Arguments added to readc2.  King 140401
c version 9.57  GETCFS: Arguments added to readc2.  King 141004
c version 9.60  GETCFS: Modify c-file format for GNSS. King 141206 
c version 9.61  GETCFS: Fix cmfrg.out format to allow 32, not 30 SVs. King 170328
c               FILLS1: Allow more than 99 zenith delays per site. King 17030
c version 9.62  FILLS1: Remove logic for implicit biases and multi-sessions (no longer supported
c                 in other modules); add slots for ECOM2 radiation-pressure parameters.  King 190429 
c               CFMRG, APRIOR. FILLS1, FILLS2, QUERYB, QUERYP,remove REORDR; Makefile: Remove multi-session 
c                 logic, change all'maxnet' to 'maxsit', dimpar.h.  King 190502 
c               FILLS1 Added site numbers back into L1 and L2 islot numbers (all sites had the same numbers
c                 for the L1 and L2-L1 islot values.  Herring 190615.
c version 9.63  GETCFS, WRITC. Updated for L1/L2 Satellite PCO values in c-file
c                 svantdx(2,3,nsat). Herring 200126. Added antdaz. Herring 200205. 
c               GETCFS: Updated 32I to 50I to allow for 35 Beidou satellites TAH 200618.
c version 9.64  FILLS1: Reassign islot1 values for biases, atm, and grads to allow for 45 SVs.  King 210619 
c               

      RETURN
      END
