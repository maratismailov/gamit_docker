      subroutine CVERSN(VERS)

      implicit none

      character*10 getmac,mach
      CHARACTER*80 VERS
      integer nblen

      mach = getmac(1)
      write (vers,5) mach(1:nblen(mach))
    5 format ('CLEAN v. 9.89 2021/04/09 21:35 UTC (',a,')')

      WRITE(6,10) vers
   10 FORMAT(a80,//)

      WRITE(6,20)
   20 FORMAT (
     .   1X,'This program is for viewing and editing GAMIT C-files.',//,
     .   1X,'New features include:                                 ',/)
c    .   5X,'PART replaces JUMP.  Without brackets, same as old PART.',/,
c    .   5X,'SAVE and ABORT are now next to each other.',/,
c    .   5X,'Plot marginal points as open circles.',/,
c    .   5X,'REWEIGHT is now available to undo an UNWEIGHT.',/,
c    .   5X,'NOTE will record a suspect epoch in the V-file.',/,
c    .   5X,'Display standard deviation s = RMS /(n-1).',/,
c    .   5X,'Display elevation angles (EL).',/,
c    .   5X,'Bias parameters needed at points with vertical bars.',/,
c    .   5X,'BIAS will flag a point as needing a bias parameter.',/,
c    .   5X,'VISIBLE/INVISIBLE toggles plotting of marginal points.',/,
c    .   5X,'Read new format C-files.',/,
c    .   5X,'Read old format C-files.',/,
c    .   5X,'Read new format INT4 C-files.',/,
c    .   5X,'Plot post-fit residuals via M-file.',/,
c    .   5X,'ALL and PART (without brackets) are now toggles.',/,
c    .   5X,'STOP is out of the way.',/,
c    .   5X,'RMS is still wrong.',/,
c    .   5X,'I1 and I2 ionospheric path delays.',/,
c    .   5X,'Testing multipath.',/,
c    .   5X,'Conjuction of error flags correct.',/,
c    .   5X,'Eliminate 7s for unweighted points.',/,
c    .   5X,'Make sky plots.',/,
c    .   5X,'Title bug is fixed.',/,
c    .   5X,'Half-cycle slips allowed for iabs(lambda) .eq. 2.',/,
c    .   5X,'Sign of residuals to polynomial is fixed.',/,
c    .   5X,'Time label is fixed.',/,
c    .   5X,'Read X-file correctly.',/,
c    .   5X,'Post-fit mode no longer garbles data.',/,
c    .   5X,'Pre-fit mode no longer reads an M-file.',/,
c    .   5X,'Choice of pre- or post- fit residuals from M-file.',/,
c    .   5X,'Time labels are now from start of scenario.',/,
c    .   5X,'Quarter cycle slips now allowed in manual MOVE.',/,
c    .   5X,'FILE is now supported through XYPLOTC and SKYMAP',/,
c    .   5X,'STOP now allows a temporary suspension of editing.',/,
c    .   5X,'Quarter cycle slips now actually saved.',/,
c    .   5X,'Optimized for 560 CPU',/,
c  version 1.15  Mods: cview,editor,polyft -- MHM for KF
c    .   5X,'Warn if polynomial fit fails.',/,
c    .   5X,'MOVE 10 is fixed.',/,
c    .   5X,'Selection of subset of C-files on M-file should work.',/,
c    .   5X,'PATCH now works on LC and LG.',/,
c  version 1.17  New dimpar: MAXSAT=12
c    .   5X,'MAXSAT=12',/,
c  version 1.18  New dimpar: MAXSAT=12
c        changes to EDITOR,GETSER,COMBO,CVIEW,READC
c        FEXIST replaced by library routine FCHECK
c    .   5X,'MAXSAT=12',/,
c    .   5X,'Case insensititve.',/,
c    .   5X,'Plot station clock corrections with CL.',/,
c    .   5X,'Units and zero point labelled.',/,
c    .   5X,'Polynomial coefficients are displayed.',/,
c  version 1.19
c    .   5X,'Post-fit pseudoranges are now obs. - calc.',/,
c    .   5X,'Subset of C-files on M-file now works.',/,
c    .   5X,'Increased resolution of clock (CL) observable .',/,
c  version 1.20
c    .   5X,'Case insensitive: should work under SR10 .',/,
c  version 1.21
c  version 1.22
c  version 1.23
c  version 1.24
c    .   5X,'Handles RINEX.',/,
c    .   5X,'Testing time tags.',/,
c    .   5X,'MOVIE and STACK feature.',/,
c    .   5X,'PART replaced by SPAN.',/,
c    .   5X,'Time tags correct for marginal points.',/,
c    .   5X,'Label output C-files with your name.',/,
c  version 1.25
c    .   5X,'Minimum vertical range now +/- 0.25 cylces.',/,
c    .   5X,'Spectral estimation (SPECTRUM).',/,
c    .   5X,'BIAS flags now larger.',/,
c    .   5X,'X-files now read without crashing.',/,
c    .   5X,'SMOOTH now works most of the time.',/,
c    .   5X,'SEEK button now under ABORT button.',/,
c    .   5X,'ABORT no longer requires double click.',/,
c    .   5X,'SAVE requires double click.',/,
c    .   5X,'SEEK should be faster.',/,
c    .   5X,'Subsets of M-files now work.',/,
c version 1.26
c Changed CVIEW to allow periodic update of C-files without
c   having to exit program - Shimon Wdowinski 10/11/90
c     .   5x,'New C-file format',/,
c     .   5x,'Can update C-files periodically',/,
c     .   5x,'Correct differencing of CL',/,
c     .   5x,'BIAS button changed: Left to remove, Right to add.',/,
c version 1.27
c     .   5x,'BIAS button handles spans of data.',/,
c     .   5x,'P+B = PATCH + BIAS.',/,
c     .   5x,'FIND to seek out and go to chan/site pair with gap.',/,
c     .   5x,'PLOT after stopping MOVIE will SEEK backwards.',/,
c *** Last changes from MIT up to here
c version 1.28
c  Fixed bug in POLYFT for case of 1 observation - YB 11/20/90
c  Added kk0 and kk1 arrays to cview.fti - YB 11/26/90
c version 1.29
c  Cominclude '../includes/errflg.h'bd = 8) - YB 12/11/90
c  Modify CVIEW to check if C-files exist - YB 12/11/90
c version 1.30
c  Modify FIND in EDITOR to search for a bias, and if there
c              is no bias look for a gap - SW 12/19/90
c version 1.31
c  Broke EDITOR into 2 files (editor.ftn, editor_subs.ftn)
c  Fixed bug in editor_plotts to show a bias in last point
c version 1.32 MAXNET=30, MAXSAT=18 in dimpar.fti
c              Yehuda Bock 12/23/90
c version 1.33 MAXEPC=1500 in makex.fti
c              Yehuda Bock 12/25/90
c version 1.34 Change dimpar.fti, DATTYP=5
c              Yehuda Bock 2/3/91
c version 1.35 In editor.ftn, fix bracket on epoch 1 bug
c              Yehuda Bock 2/14/91
c version 1.36 broke EDITOR to 7 files and many subroutines
c              replaced NOTE with ELIM that eliminates all bias and save
c                                                          at the second clic.
c              after FIND use PLOT (instead of SEEK in the past)
c              when the cursor in the data-area, a right bottom clic preform
c                                the following: FIND, PATCH, SAVE, and PLOT
c              modified FIND to search bias or gap or bias+gap at any combination.
c              Shimon Wdowinski 2/18/90
c
c      WRITE(6,21)
c   21 FORMAT (
c     .   5x,'************************************************',/,
c     .   5x,' magic right click in data area           ',/,
c     .   5x,' ELIM eliminates all bias in span         ',/,
c     .   5x,' improved FIND key                        ',/,
c     .   5x,' can delete old cfiles                    ',/,
c     .   5x,'                                          ',/,
c     .   5x,' !!!!!! Important !!!!!!!!!!!!!!!!        ',/,
c     .   5x,'                                          ',/,
c     .   5x,' Each of the MOUSE keys has a different functions ',/,
c     .   5x,' Left   - puts brackets - as in the old version   ',/,
c     .   5x,' Middle - changes menus SEEK or ABORT             ',/,
c     .   5x,' Right  - starts the FIND, PATCH sequence and     ',/,
c     .   5x,'          waits to your response:                 ',/,
c     .   5x,'          Right  - saves the patch and replot     ',/,
c     .   5x,'          Middle - aborts and replot              ',/,
c     .   5x,'          Left   - undo                           ',/,
c     .   /)
c      WRITE(6,22)
c   22 FORMAT (
c     .   5x,'************************************************  ',/,
c     .   5x,'  SLIP - search,find,patch, and save a slip      ',/,
c     .   5x,'         Left   key - only search the next slip  ',/,
c     .   5x,'         Middle key - fix slips in this plot     ',/,
c     .   5x,'         Right  key - fix slips in all plots     ',/,
c     .   5x,'  ELIM - can define a span to eliminate bias     ',/,
c     .   5x,'  FIND - improved to use the RIGHT click         ',/,
c     .   /)
      WRITE(6,23)
   23 FORMAT (
     .   5x,'  SPAN[*] includes extra 5 epochs on each side   ',/,
     .   5x,'  PUSHED BIAS is represented by a triangular     ',/,
     .   5x,'                                                 ',/,
     .   5x,'  HIDE - shows hidden biases                     ',/,
     .   5x,'  DT   - a new data type, like LC, LG, ...       ',/,
     .   5x,'         shows the 4 1-WAY series that           ',/,
     .   5x,'                contribute to the displayed data.',/,
     .   5x,'  Allows 0.5 LC PATCH in double-difference       ',/,
     .   5x,'         if total slip <= 1.0 !!                 ',/,
     .   5x,'  Allows 0.25 MOVE of L1 and L2 for codeless data',/,
     .   /)



c version 1.37 Patch in NICER to avoid CCVIEW crash
c         YB 3/10/91
c version 1.38 debug in CVIEW,EDITOR, and EDITOR_SUBS
c version 1.39 merge SIO and MIT changes
c     ALLAN:          calculate Allan Standard Deviation.
c     ED_BIAS:        allow removal of more than one bias in span
c     EDITOR_PLOTAL:  new routine to plot allan
c     NICER:          even more robust crash prevention.
c     CVIEW:          allow correct editing of a subset of C-files
c                     from the M-file list.  This has been tested
c                     and works as it should.
c                     Ask if old C-files should be deleted.
c                     NOTE: we have a potentially dangerous
c                     situation here, in that multiple "updates"
c                     WILL NOT work if when idelc = 0.
c     WRITEC:         delete old C-file if idelc = 0, per SIO
c                     Now work in place, calling new library routine irename.
c     EDITOR_PLOTTS:  Show vertical bar on last point if it is a BF.
c     ED_DATA:        Get difference calc. correct.
c     ED_MNU0:        Allow ABORT after automatic patch.
c                     Only 1 click to SAVE.
c                     Exchange positions of SAVE and ABORT.

c version 1.40 final clean up
c     GETSER          change calling sequence to return yjump and ystep,
c                     and the epochs at which they occur.
c version 2.1         minor bug on putting L1 in last slot: EDITOR and ED_MNU1.F
c                     Kurt July 8

c version 8.2 debug in GETFIL, READM, FINISH
c                     removed the 'call readm1( )' from GETFIL. The value of
c                     'nsessn' is evaluated in READM( ).
c                     Shimon Wdowinski (incorporated by YB 7/9/91)
c             SCAN:   Shimon: 7/9/91 - version 8.11
c                     add 2 types of LC-RMS, 1)sum of RMS between flags
c                     2)sum of RMS excluding big gaps
c                     SCANRMS.FTN include a message about unix sorting
c                     SCANNER.FTN include a header to the output file 'scan.rms'
c                     add a call to a new subroutine GET_RMS (in DR_PTCH.FTN)
c                     DR_PTCH.FTN add a new subroutine GET_RMS
c              DBLCLN: Initialize "      more = .true."  before calling GETFIL
c                     Shimon 7/7/91
c version 8.3         Remove extra 'endif' in SCAN.   rwk 911102
c                     WRITEC: Remove 'status' in V-file open; will not work
c                       on HP.  Donnellan/King  920226.
c                     Scripps changes (added at MIT 920327):
c version 8.12 CVIEW: Add new features: LIST and 1-WAY, modified CVIEW.fti !!!!!!
c                     add new file to Makefile, ED_LIST.FTN that contains the subroutines: READ_LIST and ED_LIST.
c                     modified CVIEW.FTN by adding a call to READ_LIST
c                     modified EDITOR_SUBS, ED_MENU0, and ED_MENUE1 to allow to new features
c                     Shimon Wdowinski  7/23/91
c             SINCLN and DBLCLN: Improved the slip-searching and cleaning algorithms
c                     changes in the files: DR_PTCH, DR_SRCH, DR_UNWT, DR_SIN1, DR_SIN2
c                     Changed output format of SCAN.DD, and the message
c                       at the end of the program.
c                     Shimon Wdowinski  7/23/91
c version 8.13 changed the list reading format in ED_LIST.FTN.
c                     Shimon Wdowinski  7/25/91
c version 8.14 changed the list reading format in ED_LIST.FTN.
c                     Shimon Wdowinski  7/25/91
c version 8.14 changed the list reading format in ED_LIST.FTN.
c                     Shimon Wdowinski  7/29/91
c              SCAN:  Shimon: 7/31/91 - version 8.13
c                     Replaced the line "write (iv,*) vers" after opening file (iv).
c version 8.15 Patch for option 3 problem, introduce jstat variable in
c                     WRITEC, FINISH, CVIEW
c                     Yehuda Bock 7/31/91
c version 8.16 Add W* - a flatening version of WL (Wide Lane)
c                     COMBO, ED_MNU1, EDITOR
c                     Shimon Wdowinski 9/30/91
c             DBLCLN: Added a slip search in WL (DR_SRCH); improved the speed of
c                       double-diff. calculation (GETSER1)  Shimon Wdowinski  9/30/91
c                       Add 3 checks to confirm a slip (DRV_DBL1, DR_SRCH)
c                       Shimon Wdowinski  10/04/91
c             SCAN  : Shimon: 9/30/91 - version 8.14- Add a slip search in WL (DR_SRCH).
c             SINCLN: Add 3 checks to confirm a slip (DRV_SIN0,DRV_SIN1, DR_SRCH)
c                     Shimon Wdowinski  10/04/91
c                     Fixed a bug in the new test of confirming a slip (DRV_SIN1, DRV_SIN2)
c                     Shimon Wdowinski  10/08/91
c                     Modified the slip detection algorithm, in particular the WL
c                     (DR_SRCH, DR_PTCH, DRV_SIN1, DRV_SIN2)
c                     Fixed a bug in the WL test. (DR_SRCH)  Shimon Wdowinski  10/14/91
c             SCAN:   Shimon: 10/17/91 - version 8.15
c                       Modified the silp confirmation routins (DRV_SCN, DR_SRCH)
c                       Changed the 3 types of LC-RMS to quick, full, and total
c                       (SCANNER, DR_PTCH)
c version 8.17 Replaced button CANCEL, which is not very useful with SAVE
c                     rearranged the order of commands in the upper panel.
c                     ED_MNU0, EDITOR_SUBS
c                     Shimon Wdowinski 10/30/91
c version 8.18 Fixed multi-session read from M-file problem in GETFIL
c                     Yehuda Bock 12/19/91
c version 8.19 Added FILE and SEEK-MODE buttons to lower menu.
c                     The SEEK-MODE determines how to search the next observation
c                          SEEK-FORW - forward seek [default]
c                          SEEK-BACK - backward seek.
c                          SEEK-LIST - forward search along a list.
c                     The SEEK-MODE affects SEEK, MOVIE, and middle button.
c                     changed ED_MNU0, ED_MENU!, EDITOR_SUBS, NEXTOB,
c                     Shimon Wdowinski 12/19/91
c version 8.20 Changed the names of the default lists, compatble with SORTER. (ED_LIST)
c                     Shimon Wdowinski 12/22/91
c version 8.21 Added wide-lane patch in PATCH and ED_PTCH, improved list massages (ED_LIST)
c                     Shimon Wdowinski 2/4/92
c version 8.22 WRITEC: do not stop if problem with v-file output
c                     Yehuda Bock 2/23/92

c version 8.4  Incorporated SIO changes 8.11 - 8.22 for EDITOR, ED_MNU1, GETFIL
c              FINISH, CIVEW, EDITOR_SUBS, COMBO, ED_LIST (new), NEXTOB, ED_PTCH,
c              CVIEW.FTI, Makefile.
c              Merged SIO/MIT changes in WRITEC
c              August changes to GETSET and READM already in.
c              Makefile:  Add ED_PTCH for dblcln.
c              Additional MIT changes:
c              READX: Add xtype,skd to call of XHDRED.
c              ED_MBLS, ED_SLIP, DRV_SIN1, DRV_DBL1 : remove redundant declarations.
c              GJ, POLYFT, READX, PATCH, SCAN : Explicitly declare variables.
c                 King 3/30/92
c              GETFIL, CVIEW, DBLCLN, SCAN, SCANDD, SCANFLG, SCANRMS - remove idelc from
c                 calling sequence (not used)    King 4/22/92
c              Incorporated SIO changes for SINCLN, DBLCLN, and SCAN  (7/23/91-2/23/92):
c                    WRITEC, GETFIL, DR_SRCH, DR_PTCH, SCANNER, DRV_SCN
c                    DRV_SIN0, DRV_SIN1, DRV_SIN2, DRV_DBL0, DRV_DBL1 - King 4/22/92
c              Merged all version files into this one.  Changes to Makefile,
c                    CVERSN, DBLCLN, SINCLN, SCANDD, SCANRMS, RMBIAS.  King 4/22/92
c              New C-file format for skd and ircint:  READC, WRITEC.   King 4/23/92
c              Add session number to C/X-file format: READC, WRITEC, READX.  King 4/29/92
c              Removed unused variables identified by Sun compiler:  GETFIL,
c                  ED_MLBS (MULTIBIAS), ED_MNU0 (ED_MENUS), EDITOR_SUBS (EDITOR_PLOTAL)
c                  ED_SLIP (AUTOSLIP, FINDSLIP, CHECK_SLIP, NOISE_LEV)
c                  NEXTOB (NEXTOBSERV), ED_MNU3, FINISH, READX, RINEX
c                  GSUBS_X (GDUMP), ED_LIST (READ_LIST, ED_LIST), DBLCLN
c                  DRV_DBL1 (IDEN_SLIPS, SORT_SLIPS, PATCH_SLIPS, DO_PATCH0, DO_PATCH1, DO_PATCH2)
c                  DRV_DBLL0 (DRIVER),
c                  DR_SRCH (SEARCH_SLIP, FIND5, GET_EPS, CHECK_SLIP, CONFIRM_SLIP, CHECK_GLOBAL, CHECK_WL)
c                  DR_PTCH (DR_PATCH, UNWT_LEFT, UNWT_RIGHT, CAL_RMS, GET_RMS)
c                  DRV_SIN0 ( DRIVER1, DRIVER2 ), DRV_SIN1 (FINDSLIP1)
c                  DRV_SIN2 ( IDEN_SLIPS, PATCH_SLIPS, DO_PATCH1, DO_PATCH2 )
c                  GETSER1 (GETSER1, GERSER2, GERSER4, MINMAX ), GMSG_DUMB
c                  DR_UNWT ( UNWT_SNGL, UNWT_TRBL1, UNWT_ALL )
c                  DRV_SCN0 ( DRIVER0, DRIVER1, DRIVER2 ), DRV_SCN (FINDSLIP)
c                  SCANNER (SCANNER, SCANNERDD ), SCANDD

c version 9.11  Modified CVIEW and SINCLN to treat cfiles with
c               data collection intervals > interval between epoches.
c               CVIEW.FTI, READC, EDITOR_SUBS, DR_UNWT, DR_PTCH
c                 Shimon Wdowinski 5/31/92
c version 9.12  Fix typos and omissions in above:  EDITOR, EDITOR_SUBS, READX, WRITEC
c                 King  6/5/92
c               Add SHOWRMS to library and Makefile   Fang/King 6/5/92
c version 9.13  READX, RINEX :  Initialize epoch-skip variables.  King 7/6/92
c version 9.14  Makefile :  Remove obsolete GETSERS DRV_SCN0 and DRV_SCN1, and add
c                           DRV_SCN (FINDSLIP) to /stdrel.
c               Remove obsolete Fortran files from /stdrel:
c                  DRV_SCN0, DRV_SCN1, DS2HMS, DVERSN, FEXIST, GETSERD, GETSERS,
c                  GSUBS, HOLD, ICLARG, IRENAME, MINMAX, NBLEN, NEWCHR,
c                  ONE_PTCH, RVERSN, SCAN, SCANFLG, SEEKOB, SIN_PTCH, SVERSN,
c                  TNEXTOB, TSMOOTH.
c               DRV_SCN (FINDSLIP), SCANRMS:  Remove unused variables.
c                   King 92/7/16
c               EDITOR_SUBS, GETSER : Add diagonal 'flags' to the bias flag to
c                   show which satellites/stations have biases.  Oral/King 92/7/16.
c               Makefile : new routine GET_COLOR.RO.C to use color display.
c                   Oral/King  92/7/16.
c               CVIEW.FTI : Change order of variables in common  /STIME/ to avoid
c                   performance degradation. King  92/7/17.
c               EDITOR_SUBS: Remove unused variables and correct non-executed DO loop.
c                   King  92/7/17

c version 9.15  (MIT)
c               SINCLN : Optional 3rd command-line argument to delete input file.
c               DBLCLN : Change comment to clarify input file delete.  King 92/9/30
c
c version 9.15  (SIO)  : Debuged the problem with the diagonal bias flags on SUN
c                        (subroutine editor_plotts in EDITOR_SUBS.F).
c                        Shimon Wdowinski  8/19/92

c version 9.16  Added 2 new features to view hidden data:
c                HIDE - shows hidden biases
c                DT   - a new data type, like LC, LG, ...
c                       shows the 4 1-WAY series that
c                              contribute to the displayed data
c
c                Split the file editor_subs.f into 2 files:  ed_subs1.f, ed_subs2.f
c               (files Makefile, CVIEW.FTI !!!!
c                     EDITOR, ED_SUBS1, ED_SUBS2,
c                     ED_MENU0, ED_MENU1, ED_MENU2).
c                   Shimon Wdowinski  9/24/92

c version 9.17  Modified the code to allow:
c                 0.5 LC PATCH in double-difference
c                      if total slip <= 2.0
c                      (applies only to CVIEW, DBLCLN)
c                 0.25 MOVE of L1 and L2 for codeless data.
c                 (files ED_MNU0, ED_MOVE, ED_PTCH, DR_PTCH, DR_PTCH1, PATCH).
c                   Shimon Wdowinski  10/14/92

c version 9.18  RINEX: set epochs correctly, allowing holes.
c               GETFIL : allow RINEX file name starting with m or c.
c                   Feigl/King  10/16-28/92
c               DR_PTCH : Fix missing comma in calls to PATCH.  King 10/16/92
c               ED_SUBS2: Declare and dimension kpl,plt in editor_plotdt.  King 10/18/92
c               RINEX, ED_SUBS2: Remove unused variables.   King 10/20/92
c               Fixed a bug in ED_SUBS2 that appeared only in the SUN CVIEW's
c                   data presentation.   Shimon Wdowinski  10/20/92

c version 9.19  added 5 data on each side of a SPAN[*] data representation.
c               reduced the L1 half cycle patching from 1.5 to 0.5.
C               eliminated the 0.5 L1 patch from SCANDD
c               added pushed bias display
c                   EDITOR, PATCH, DR_PTCH, ED_SUBS2
c                   Shimon Wdowinski  11/20/92
c               Restore MIT 9.18 changes to CVERSN, ED_SUBS2.   King 12/3/92
c               Implement color display and slight narrowing of window in CSUBS.C
c                     Oral 2/19/92; King 12/3/92

c               ED_SUBS2:  Fix bug causing no SPAN display with mixed samples.
c                     Oral/King 12/16/92.
c               GET_COLOR_RO.C:  Remove screen listing of loaded colors.
c                     Oral/King 12/16/92
c               SHOWRMS:  Enhanced features.  Oral/King 12/16/92
c               SORTV (script):  Enhancements.  Feigl/Oral/King 12/16/92

c version 9.20  CVIEW
c                   fixed bug in SPAN[] mixed sampling rate presentation
c                   expanded the margin of SPAN from 5 points to 5%
C                   fixed bug in DT data presentation
c               SCAN
c                   modified the GET_RMS algorithm to include:
c                       mixed sampling rate
c                       pushed bias
c                   modified the output of scandd to include bias information
c               EDITOR, ED_SUBS2, DR_PTCH, SCANNER, DRV_SCN
c                   Shimon Wdowinski  12/17/92
c               Fix typo in SCANNER and declare variables explicitly in
c                   EDITOR, ED_SUBS2.
c               Push biases properly for vscan.out=> vexpty.ddd.
c                   SCANNER.    Shimon Wdowinski/King 12/22/92
c               SHOWRMS: Enhancement, standardization.  Fang/Oral/King 12/22/92
c               READC: Output up to 24 satellites.  Bock/Fang 01/12/93
c         9.21  READC,GETFIL: Small modifications.  Bock 01/26/93
c               GETFIL: More small bugs.   Bock  01/27/93
c         9.22  Modify Apollo Makefile to include -save compiler option
c               (fixes the menu problem associated with find) Bock 02/06/93
c               CVIEW (MIT v. 9.21)
c                   ED_SUBS2: Temporary changes at MIT in parallel to SIO development:
c                       Fixed bug in display of unmatched sampling data when first
c                       (real) data points don't coincide. Additional changes needed
c                       to calls.
c                   PATCH:  Changes to make LGWL patch to work
c                   ED_PTCH: Changed to use more numbers of epochs so that for sparsely
c                       sampled data there are enough epochs to get the number of
c                       data points.  Needs some more mods to change calling arguments
c                       to get more information.  Current work-around could be improved.
c                       (Corresponding changes needed to scan routines.)
c                   Makefile:  Changed 03 (zero3) to O3 (ohh3).
c                       Herring 1/28/93
c                   ED_SUBS2: Fixed display bug when first satellite has no data.
c                       Herring/King 2/4/93
c                   READX: Add antcod to XHDRED call.  King 930212
c               GET_COLOR_RO.C (Sun only): Rewrite to avoid system crashes.  Oral/King 930308
c               GETFIL:  Allow second argument for day-of-year in multisession
c                       M-files.   Murray/King 3/25/93
c               PATCH : Comment out half-cycle patches (too dangerous).  King 930409
c               Makefile: Remove test program tgsubs from Apollo Makefile, and
c                      Apollo and Sun stdrel directories.   King 930413
c         9.23  GETFIL : Reinstate some old SIO bug fixes. Bock 930429
c         9.24  ED_SUBS2: Incorporated Burch change, 3000 => maxepc Bock 930514
c                         Copied to SUN version, Bock 930706
c         9.24  PATCH: Improve handling of large gaps; don't add half-cycles
c                         to L1.  Herring (930519)/King 930910
c               SINCLN: Allow 1 command-line argument (write in place) without
c                       killing input and output.  Herring/King 930910
c               POINTR, ADDADJ:  Debug and stop. King 930913.
c               ADDADJ, READC, WRITEC : Correct residuals for multiple-zenith-delay
c                       parameter adjustments.   King 930917
c               POINTR: Remove debug and stop.   King 930925
c               READC, SHOWRMS:  Increase format to 32 satellites.   King 931018
c               ADDADJ : Fix code for multiple-zenith delays.   King 931019
c               READC: Set default zenith model to piecewise linear.  King 931020
c          9.25 READC, WRITEC, ADDADJ, CVIEW.FTI: Fix bugs in passing zenith-
c                     delay quantities (put in new common).   King 931213
c               READC: Remove debug.   King 931229
c               Makefile: Add explicit dependencies.  Fang 940215
c               ADDADJ, READC: Define zenmod for numzen=1.   King 940325
c               ADDADJ: Fix bug for nsite*numzen >99.  King 940328
c               Makefile: Remove forced object making for C routines.  Fang/King 940412
c               ED_LIST: Set # in list for cview.list.  King 940413
c          9.26 CSUBS:  Duplicate code to allow HP/Sun compatibilty.  Fang/King 940428
c               RINEX:  Change call to lib/rrinex.f.   King 940505
c               Makefile: Change HP compiler flags from +e +E1 to +U77.  Fang/King 940505
c               CSUBS:  Use 'sleep' vice 'usleep' on both Sun and HP.  Feng/King 940606
c               Makefile:  Shift ranlib to be executed only once per module.  Herring/King 940622
c          9.27 ADDADJ, POINTR: Add earth-rotation parameters to residual correction.  McClusky/King 940715
c          9.28 Change vfile opens from 'append' to 'unknown' to avoid non-ANSI standard:
c               ED_MNU3, ED_SUBS1, WRITEC.  King  941205
c               SHOWRMS: Fix format in reading values.  Vigny/King 941214
c               ED_SUBS1 , SCANNER :  More changes for XL/RISC compiler.  King 950103
c          9.29 Change vfile opens to include access = 'append' this is ANSI standard:
c               ED_MNU3, ED_SUBS1, WRITEC.  McClusky  950405
c               ED_SUBS2: Modified to plot low elevation data as rectangles (squares) when
c               marginals are requested. McClusky 950412
c               EDITOR: Modified to correctly scale and plot satellite elevation
c               angle data. McClusky 950412
c               CSUBS: Added new subroutine called GRECT and GRECT_, which
c               plots open rectangles. McClusky 950412
c               WRITEV: New routine (similar to writec) to write out vcview.out file, without
c               writing any cfiles. McClusky 950413
c               FINISH: Modified routine to allow 5th finishing option. Just write
c               out vcview.out log file. McClusky 950413
c               FINISH, CVIEW: Modified to correctly write out vcview file when option 5
c               is selected. McClusky 950421
c
c          9.30 New C- and M-file formats:  READC, WRITEC, WRITEV, READM, FINISH.  King 950524/950607
c               New logic to handle pointers for postfit residuals:  cview.fti (common /postft)
c                  Makefile, ADDADJ, POINTERS (new), READC, WRITEC, WRITEV.  Remove old
c                  function routine POINTR.   King 960606
c          9.31 READM:  Add decimation factor to M-file record 2.   King 950616
c               READC, WRITEC, WRITEV:  Fix definition of maxpar.  King 950622
c               WRITEV:  Remove extra calling argument for READC1.  King 950711
c               Makefile, FINISH:  Change WRITEV to WRITE_VCVIEW to avoid conflict with HP and
c                   Sun system routines.   Herring/King 950712
c               ADDADJ: Remove duplication declaration and unused variable.  King 950712
c          9.32 READC, WRITEC, WRITE_VCVIEW:  Add antenna reference point offset to C-file
c                    record 2.  King 950717
c          9.33 Makefile: Add lines for DEC.  Sanli/King 950719
c               ED_BIAS, ED_MNU0: Add start and stop epochs to calling arguments (previously
c                    undefined in ED_BIAS).   Sanli/King 950719
c               ED_MNU1 (ED_PRINT): Comment out lines to avoid compiler warning about
c                    inaccessible statements (ED_PRINT for dump doesn't work).  Sanli/King 950719
c               DR_UNWT(UNWT_TRBL): Bug in igap.   Sanli/King 950719
c               ED_MNU0:  Fix overrun comment line.  Sanli/King 950719
c               ED_SUBS1(EDITOR_MENUS):  Add comment about compiler warning.  Sanli/King 950719
c               SCANNER(GET_RMS): Fix misnamed next=>inext.  Sanli/King 950720
c               DR_UNWT(UNWT_TRBL): Bug: undefined igap--not yet fixed.   Sanli/King 950719
c               READC: Stop if number of zenith delays exceeds dimensions.  King 950721
c          9.34 DR_UNWT(UNWT_TRBL): Fix use of undefined igap and igood.  Shimon/King 950725
c               SCANNER(GET_RMS):  Correct fix to 'inext' bug .  Shimon/King 950725/950728
c          9.35 GETFIL: Fix too-short argument list for READM3 (segmantation fault when
c                  selecting only some C-files). McClusky/King 950809
c               GETFIL:  Remove redundant declarations.  King 950818
c          9.36 GETFIL, READC: Change site index for calls to READX and RINEX from 'i' (m-file index)
c                  to 'isite' (index of sites actually used; and add 'i' (m-file index) to call to
c                  READC so that postfit o-c's are correct.  King 950825
c          9.37 ADDADJ, POINTERS, READC:  Rewrite logic to avoid looping over all parameters
c                  at each epoch--significant time-savings in reading C-files in postfit
c                  mode.  (Still need to remove 'parptr' from cview.fti.) Herring/King 950830
c               SHOWRMS:  Fix computation bug for 'total' or 'full' values.  Fang/King 950831
c               ADDADJ:  Fix bug in postfits from atmospheric parameters.  King 950902
c     include '../includes/cview.h'parptr from common /postfit/ (no longer used. King 950904
c               SCANNER:  Remove printing of explcit plus signs in vscan.out to avoid
c                   problems with sortv on HP.  Herring/King 9508904
c          9.38 READM: Fix bug in call of readm3. McClusky 950908
c               GETFIL: Now read dummy argument elvcut_dum not elvcut. McClusky 950908
c               READC: Now elvcut_cfiles and elvcut read into an array. McClusky 950908
c               READM: Now read dummy argument elvcut_dum not elvcut. McClusky 950908
c               FINISH: Now read dummy argument elvcut_dum read not elvcut. McClusky 950908
c               WRITE_VCVIEW: Now elvcut read into an array. McClusky 950908
c               WRITEC: Now elvcut read into an array. McClusky 950908
c               CVIEW.FTI: Added /cutoff/elvcut as a common. McClusky 950908
c               READC: Added code to check solve elevation cutoff against cfile low elevation flags,
c               and change the loel flag if solve cutoff higher. Bias pushing included. McClusky 950908
c               WRITE_VCVIEW: Modified to work with new version of addadj McClusky 950908
c               WRITEC: Modified to work with new version of addadj. McClusky 950908
c               FINISH: Modified to pass mfile index into WRITEC and WRITE_VCVIEW. McClusky 950908
c               WRITEC: Added code to check solve elevation cutoff against cfile low elevation flags,
c               and change the loel flag if solve cutoff higher. McClusky 950908
c               READC: Put in a kludge fix to convert cfile elevation cutoff written
c               on the cfile header from radians to degrees. McClusky 950912
c               WRITEC: Same as for readc above. McClusky 950912
c               WRITE_VCVIEW: Same as for readc above + some cosmetic changes. McClusky 950912
c               FINISH: Some cosmetic changes to screen output for vcview.out. McClusky 950912
c               READM:  Remove duplicate declarations.  King 950923
c         9.39  SINCLN, RMBIAS: Add extra argument to READC, WRITEC calls.  King 950926
c               WRITE_VCVIEW: Correctly identify subroutine on stop.  King 950926
c               Correct spelling of 'flagged', 'weight', and 'channel' in 9 routines:
c                 DR_CLNR, DR_PTCH, DR_UNWT, DRV_DBL1, DRV_SIN1, ONE_PATCH, SIN_PTCH,
c                 DRV_SIN2, DR_CLNR1.  King 950926
c               READC: Removed a kludge to convert cfile elevation cutoff written
c                 on the cfile header from radians to degrees. McClusky 951002
c          9.40 READX:  Add receiver and antenna ids to lib/xhdred call.  King 951019
c          9.41 Makefile: Variable ranlib to make compatible with Solaris 2.  Fang/King 951208
c               SCANRMS, SCANDD: Update message at end on sorting.  King 951229
c          9.42 report_stat changes to routines; addadj, combo, cversn, cview, dr_srch,
c               drv_scn, ed_list, editor, getfil, nextob, polyft, readc, readm,
c               readx, rinex, scandd, scanner, vcopy, write_vcview, writec. McClusky 960229
c               READM2 (lib):  Fix bug in defining 'nsat'.  McClusky/King 960301
c          9.43 Makefile, SCANNER: Remove obsolete SCANRMS, fix print of full rms in
c                  vscan.out.   King 960522
c          9.44 RINEX: Change observable dimension from 6 to maxdat.  King 960704
c               CVERSN: Remove long comment line (DEC complains). Bock 960609/960710   
c          9.45 Changes for new kf/gamit structure and Makefiles from unimake.  King 960718
c                   All routines:  removed trailing blanks and replaced lib/*.fti 
c                   includes by ../includes/*.h.    
c          9.46 SCANDD: Stop if GAMIT.fatal exists.  King 960912 
c          9.47 CSUBS_HP:  Remove underscore for gclip.  King 960930 
c               Remove GSUBS, GSUBS_CGI, GSUBS_GPR.  King 960930
c          9.48 READC: Modified call to POINTERS and ADDADJ to pass grad_index
c               Also added call to ATM_GRAD_PART to compute grad partials. McClusky 961002
c               WRITEC: Same mods as to READC. McClusky 961002
c               WRITE_VCVIEW: Same modes as to READC. McClusky 961002
c               PONITERS: Modified logic to compute the grad parameter indices. McClusky 961002
c               ADDADJ: Modified to add atm grad adjustmets to postfit residuals. McClusky 961002 
c          9.49 DR_PTCH, EDITOR, ED_MOVE, ED_SAVE, READC, RINEX< WRITE_VCVIEW:  Call new /comlib 
c                routines for integer*2 functions (necessary for gcc compiler).  Trengoning/King 961017
c               ATM_GRAD_PART: removed duplicate implicit none statement  Tregoning 961018
c          9.50 READC, RINEX New routines min02 and max02 require that both input arguments are
c               integer*2 this has been fixed throughout clean now. In the future when we do
c               kinematic surveys that have over 32,768 epochs many i*2 variables refering to
c               epoch number will have to be changed to I*4. McClusky 961203
c          9.51 G77 compiler changes:
c               GSUBS_X: Changed call to the library function fseek to the new wrapper 
c                 subroutine call to fseekg. McClusky 970113
c               PATCH:  Changed iabs call to generic abs call....
c               DR_UNWT, DR_SRCH, DRV_SIN2, DRV_DBL1, POINTERS, ED_SUBS2, RINEX, GJ, GETFIL, 
c                   ED_MLBS, EDITOR:  Initialised previously uninitialised variables, removed unused 
c                   variables, and explicitly typed all unexplicitly typed variables. McClusky 970110
c          9.52 GETFIL:  Allow exit after failing to find a file.   King 970129
c               DR_PTCH, ONE_PATCH, PATCH, SIN_PTCH: Change name of sb 'scan' to 'scanlgood' to 
c                   avoid conflict with Fortran 90 intrinsic.  McClusky/King 970321
c          9.53 ED_FIND, ED_LIST, ED_MLBS, ED_PTCH, ED_SLIP, NEXTOB:  Avoid splitting Holleriths
c                   over several lines to satisfy the picky DEC compiler.  King 970606
c               SHOWRMS:  More Hollerith splitting.   King 970807
c          9.54 CVIEW, DBLCLN, FINISH, GETFIL, READM, SCAN, SCANDD: Increase m-file name to 16 
c                 characters to be consistent with CFMRG and SOLVE and to allow post-fit editing.   
c                 King 970809
c          9.55 Changes to remove duplicate subroutine names and thus allow saving of object 
c                 modules for DEC OSF4 Make.  King 970918
c                 Makefile.generic: add CHECK_WL
c                 Make CHECK_WL a separate file; remove from DR_SRCH and ED_PTCH.  Use
c                     the latter version, modified by Herring for decimated data.
c                 Rename check_slip within DR_SRCH to dr_check_slip, and within 
c                     ED_SLIP to ed_check_slip (routines are very different).
c                 Rename gmsg within GMSG_DUMB to gmsgd, called by non-CVIEW routines DRV_DBL0, 
c                     DRV_SIN0, DRV_SCN0 to avoid conflict with gmsg within gsubs_x, which is called
c                     by the CVIEW routines to write to an X-window.  In the shared routines
c                     PATCH (patch, erscan, patch_lgwl) and NEXTOB (seeklist), decide which to 
c                     call by checking the name of the executing program.  
c                 Rename iden_slip within DRV_SIN2 and DRV_DBL1 to iden_slip_sin and iden_slip_dbl,
c                     respectively, and mod calls DRV_SIN0 and DRV_DBL0 accordingly.   Ditto with
c                     sort_slips, patch_slips, do_patch1, and do_patch2.
c                 Rename findslip within ED_SLIP to findslip_ed and within DRV_SCN and SCANNER to 
c                     findslip_scn, and mod calls in ED_SLIP and DRV_SCN0.  
c                 Remove unused routines DRIVER, DRIVERA, DRIVERS, DR_CLNR, DR_CLNR1.
c                 Makefile: add (new) CHECK_WL to cview_lib, also gmsg_dumb; add X11LIB
c                     to dblcln, sincln, and scandd since gsubs_x (gmsg) is present in (but not)
c                     used by) PATCH and/or NEXTOB called by these programs.  
c          9.56   RINEX: Changed incorrect integer calling argument types ircvr and iant
c                     to character*20 rcvnum and antnum. McClusky Herring 980115.  
c                 POLYFT:  Replace i(i-1)/2 calculations by function call to avoid HP compiler bug. 
c                     Herring/King 980121        
c                 RINEX:  Write times in GPST, not UTC, using modified lib/wtime.f.  McClusky/King 980123  
c                 AUTOSLIP: Wrong number of arguments in call to patch_save.  King 980227
c                 AUTOSLIP, ED_CHECK_SLIP, FINDSLIP_ED, NOISE_LEV:   Three variables declared single 
c                    precision in called routine but double precision in caller.  Change variables to
c                    double precision and intrinsics using them to generics.  King 980227 
c                 GETFIL, READC, RINEX, WRITEC, WRITE_VCVIEW: Fix calls to report_stat.  King 980228
c                 GET_EPS:  Remove unused last argument.  King 980228
c          9.57   GETFIL:  Allow reading of a 3rd command-line argument for c-file version
c                    (useful for cview and scandd in postfit editing mode).  King 980304  
c          9.58   READC, WRITE_VCVIEW, WRITEC:  Add 'norb' to calling argument for lib/readc2,writc2.  King 980911
c                 ADDADJ, POINTERS, READC, WRITEC, WRITE_VCVIEW: Change 'svclk_index' to 
c                     'svant_index' and add logic for adjustment.  King 980911/980914 
c                 READC, WRITEC, WRITE_VCVIEW: Add svclock terms to calling arguments for lib/readc5,writc5.  King 981002
c                 cview.h, ADDADJ, GETFIL, FINISH, READC, READM, WRITEC, WRITE_VCVIEW:  
c                    Add adjustments for multiple gradient parameters.  King 981002
c          9.59   CSUBS_SUN CSUBS_SUN, CSUBS_HP, CSUBS_DEC: Modified so that a smaller text font
c                     is used on resolution screens (like laptops). McClusky 981026 
c                 WRITEC:  Fix calling argument for writc2.   King 981118 
c                 Remove implicit statement to statisfy SGI compiler, and add implicit none;
c                    all routines compiled ok previously with explicit compiler switch:
c                    FINISH, GETSER.  Morgan/King 981231
c                 Remove unused DS2HMS and DAYJUL.   King 990225
c                 WRITE_VCVIEW: Fix dimensions of svant_index.   King 990225
c                 SHOWRMS: Fix undefined loop variable.  Fang 990308   
c                 WRITE_VCVIEW: Avoid keyword 'access' in 'open' for compatibility with IBM 
c                     and f90).  Fang/King 990324   
c          9.60   SHOWRMS:  Fix loop bug.  Fang 990429
c                 GETFIL: Initialize lclcfv.  Espinoza/King 990528
c          9.61   RINEX: Fix timcon call for y2k.  King 991012
c          9.62   CVIEW.H:  Change name of common STIME to TIMES to avoid conflict with
c                    IRIX system routine.  King 000816
c                 DBLCLN, DR_PTCH, GET_SER1, RMBIAS:  Change to implicit none.  King 000816   
c                 CSUBS_SOL: Declare void to avoid (name argument not assigned) to avoid 
c                     IRIX warning.  King 000816
c                 CVIEW_ICON: Add 'const' to cview_icon_bits declaration to satisfy IRIX.  King 000816
c          9.63   RINEX:  Changes for RINEX version 2.10:  version and interval now real,
c                      extra decimal place for time of first observation.  King 010301
c                 RINEX: Change calling argument to RRINEX to return type-of-SV flag; call
c                      GPS_ONLY to remove non-GPS SVs.  King 010313     
c          9.64   BOUND, ED_SUBS2, READC,  READX, RINEX: Add 'int' to intrinsic calls to avoid 
c                     type mismatch (f90).  Herring 010610   
c                 ED_MNU3: Add trap for undefined variable (relic code?).  King 010731
c          9.65   ED_SUBS1:  Fix statement extended beyond col 72 (harmless). King 031027
c                 SHOWRMS: Fix line-length typo in comment.  King 031027   
c          9.66   CSUBS: Reduce window size.  Herring 031218 
c                 EDITOR, ED_SUBS1: Add y-axis average to display.  Herring 03122
c                 ADDADJ, POINTERS, READC, WRITEC: Combine average and multiple zenith delays.  King 040628 
c          9.67   READC, READX,WRITEC, WRITE_VCVIEW: Change format of C-file record 2. King 050201
c          9.68   GSUBS_X: Made thsi routine compiler dependant to avoid missing gfortran fseek. McClusky  060823
c                 READC, READM: Fix bugs in calling arguments for readc2 and readm3.  McClusky/King 060824
c          9.69   MULTIP ED_PTCH: Reordered subroutines within these routines to work around a bug
c                 in the gcc3.4.0 compiler under intel OSX 10.4.8. McClusky 070202 
c          9.70   RINEX: Remove unneeded variables from call to /lib/rrinex; removed 
c                  unneeded variable from call to /lib/gps_only. King 070416
c          9.71   ADDADJ, COMBO   : Add dummy statement for unused calling arguments. King 070906 
c                 ED_SUBS1: Fix obsolete do-loop construct.  King 070910  
c                 Makefile: Override -Wunused option for gfortran (too many problems).  King 070910
c                 SHOWRMS: Fix potential infinite loop on read.  King 070910. 
c          9.72   CVIEW: Remove unused format statement. King 080306   
c                 READX: Fix calling arguments for lib/xhdred.  King 080509
c                 SHOWRMS: Add iostat to file read.  King 080804/081006
c          9.73   CSUBS: Fix to work with high-aspect-ratio screens. Herring 080925
c                 ED_SUBS1:  Remove non-interger do loops to conform with f90.  Morgan 081203   
c                 CSUBS: Increase the height tolerance a little more to avoid hiding buttoms. Herring/King 090508
c                 ED_SUBS1: Fix undeclared 'icount'.  Herring 091001
c          9.74   READC, WRITEC, WRITE_VCVIEW: Add variable to C-file calling routines.  King 100827
c          9.75   READC, WRITEC, WRITE_VCVIEW: Add variables to C-file calling routines. King 130115
c          9.76   ED_MNU3, WRITEC, WRITE_VCVIEW : Change 'access' to 'postion' in open statement (required for AIX compiler). King 131121
c          9.77   READC, WRITC,  WRITE_VCVIEW: Add variables to C-file calling routines. King 140327
c                 READC, WRITC,  WRITE_VCVIEW: Add variables to C-file calling routines. King 140401
c          9.78   READC, WRITC, WRITE_VCVIEW, READX: Variables added for GNSS. King 141004
c          9.80   READC, READX, WRITC, WRITE_VCVIEW: Change c-file format to support GNSS. King 150109 
c                 RINEX: Hard-wire GPS in call to lib/get_gnss (need to fix later). King 150327
c          9.81   GETFIL, RINEX: Add terminal input for GNSS requested from RINEX.  King 151006  
c                 RINEX: Increase dimension of iobtypx to allow backup L1, L2 pseudoranges.  King 151028
c          9.82   READC: Fix argument mismatch for readc5.  King 151116
c          9.83   Makefile: Remove 42 routines used only by sincln, dblcln, and rmbias, no longer supported. King 151229
c                 CVIEW, CHECK_WL, COMBO, EDITOR, ED_SLIP, ED_PTCH, GETFIL, GETSER, PATCH,  READC, READX, 
c                   RINEX, SCANNER, SETL_FREQS(new), WRITEC, WRITE_VCVIEW, Makefile, includes/cview.h:  Make 
c                   the frequencies SV-dependent arrays, set from the observation files, not hard-wired 
c                   for GPS.  King 151230
c                 READX: Set the gnss code from the x-file satnam.  King 160513
c                 RINEX: Fix bug in reading RINEX observable types.  King 160513
c          9.84   READX, RINEX: Add frequencies for IRNSS.  King 160823 
c          9.85   READX, RINEX, SET_FREQS: Fix naming of GNSS frequencies.  King 170605      
c                 RINEX: Change thge frequency assignments to match the libraries/include/freq_def.h
c                   variable names; fix dimension bug.  Herring 170726  
c                 RINEX: Fix bug in setting frequencies.  King/Herring 170727
c                 READX: Fix missing calling argument for lib/xhdred. Herring/King 180309
c                 READX: Change obsolete 'skd' to 'gnss' in calls to lib/xhdred. King 180322
c          9.86   FINISH, GETFIL, READM: Remove multi-session variables. King 190503
c                 READC, READM, WRITEC, WRITE_VCVIEW, includes/cview.h: Change obsolete 'maxcfl' to 'maxsit'.  King 190524 
c          9.87   READC, WRITE_VCVIEW, WRITEC: Updated for L1/L2 Satellite PCO values in c-file
c                   svantdx(2,3,nsat). Herring 20016.  Added antdaz. Herring 200205.
c                 SHOWRMS: Updated 32I to 50I to allow for 35 Beidou satellites TAH 200618
c          9.88   CSUBS: Added stdlib.h to list of includes to accommodate calls of exit and malloc, and
c                        explicitly added return type "void" to all function definitions except get_colors ("int");
c                 CVIEW_ICON: Changed variable type to "const char" from "unsigned char" in declaration of cview_icon_bits;
c                 GET_COLOR_RO: Added stdlib.h to list of includes to accommodate call of exit, and
c                               explicitly added return type "int" to function definition. MAF 20201201
c          9.89   CSUBS: Changed declaration of get_colors() to preamble from call within gsetup() function
c                        and changed output type of gsetup() function to "int" from "void". MAF 20210409
c  
      return
      end





