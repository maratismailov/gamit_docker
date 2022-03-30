Copyright (c) Massachusetts Institute of Technology,1992. All rights reserved.
      SUBROUTINE CVERSN(VERS)

      implicit none

      CHARACTER*10 GETMAC,MACHIN
      CHARACTER*40 VERS
      character*45 libver
      character*256 message

c     return non-blank length of string
      integer nblen

      MACHIN = GETMAC(1)
      WRITE (VERS,5) MACHIN(1:nblen(machin))
    5 format ('9.67 of 2020/06/18 14:42 UTC (',a,')')
c   5 format ('%I% of %E% %U% (',a,')')
c     The above line contains magic SCCS keywords which will be
c     automatically updated by delta and get.
c     Please do not change them.

      call lversn(LIBVER)

      write(*,'(a)')' '
      WRITE(message,10) vers,libver
   10 FORMAT('Started CTOX Version ',a34,' Library version ',a45)
      call report_stat('STATUS','CTOX','cversn',' ',message,0)

c   pre-sccs:     Written 87/4/13 by R. King
c   Version numbers prior to 92/2/11 refer to CTOX; afterward to CVERSN
C   MOinclude '../includes/dimpar.h'
C        Add pseudo-range and provide for reading old X-files to get it-rwk 880607
C     include '../includes/dimpar.h'hese changes in 6.4  -- rwk 880726
c   version 6.6  Output non-zero flags in XDTRIT -- 880726 mhm
c   version 6.7  Ignore extraneous flag=2's due to problems in
c                cfedt -- 880726 mhm
c   version 7.1  Read a D-file to get the X-file names.  Decimate the output if
c                  requested.  --900105 rwk
c   version 7.2  Modify READD to handle the J-file name. -- 900105 rwk
c   version 8.1  Write X-file header with latflag and lonflag
C   version 8.2  Clean up print outs.
c   version 8.3  normal point
c   version 8.4  Modify CDTRED for correct ASCII output of o's and omc's
c                Lowercase dimpar.fti in CDTRED, CHDRED, XDTRIT, XHDRIT
c   version 8.41 New Makefile, cpu mathlib_sr10, increase
c                maximum stacksize for NMPNT.FTN -- YB 5/28/91
c   version 8.5  CHDRED & XHDRIT - increase antenna heights by 1 decimal place
c                YB/KS 7/10/91
c   version 8.51 XHDRIT - add significant place for radius (0.1 mm) -- YB 8/25/91
c   version 8.52 CHDRED - comment out much of screen output -- YB 10/6/91
c   version 8.53 (9.1) Expanded X- and C-file format (for kinematic). Changes to CDTRED
c                  and new routine ROLDC4 -- JFG & YB 11/11/91
c                CHDRED - restore print output for dump mode, but write
c                  only the title and site lines in copy mode.  Change
c                  'phi' to 'psi' in print - king  911221

C   Version 8.1   Create CVERSN and install 91/12/21 version as sccs.
c
c           8.2   Make CTOX call CVERSN to check library version. King
c                 CDTRED, XDTRIT, XHDRIT :  Correct call to DS2HMS.   Bock
c                 CHDRED:  Add logical variable to specify dump mode for print.  King
c                 Write original data interval, if different, on the X-file.
c                    Pass EXTRA from CTOX to XHDRIT.  King 920211
c                 Re-started (with loss of history) sccs files.  King 920409.
c                 Makefile:  removed TIMINC (in library).  King 920409.
c                 CTOX:  Fix order of calling arguments for CHDRED.  King 920409.
c                 New C-file format for skd, ircint:  CTOX, CHDRED, XHDRIT.
c                    Change static kinflg from 'X' to 'S'.
c                    Remove ROLDC4,RNEWC4 from Makefile.  King 920423
c                 Eliminate unused variables:  BEAT, CDTRED, CHDRED, CTOX
c                    GROUPN, NMPNT, READD, XDTRIT, XHDRIT.   King 920424
c                 Add session number to X-, C-file formats.  Changes to
c                    CTOX, CHDRED, XHDRIT.    King 920429

c   Version 9.1   Changes copied to SIO by YB.  End of sccs.   King 920506
c           9.11  XHDRIT:  Fix format for 'SATELLITES' line.  King 920507.
c           9.12  CTOX  :  Fix print for option 1.  Bock 920513
c                 CHDRED:  Remove warning about old format.  King 920513
c                 XHDRIT:  write X-file header output to correct unit
c                            correct output for static expanded format case
c                            Bock 920513
c                 CDTRED:  Add delay rate to printout.  King 920514
c           9.13  READD :  Read L-file from D-file.  King 920604
c                 CTOX, XDTRIT: Output epoch by epoch coordinates properly
c                          Bock 920606
c           9.14  CTOX :   Fix format for dump file.   Feigl/King 920716
c                 Makefile : Remove extranesous routines from /stdrel to match
c                          Apollo and Sun /sccs.    King 920717
c           9.15  READD:   Change integer reads to free format.   King 920915
c           9.16  CTOX, CDTRED: Add capability to write a file of azimuth, elevation,
c                          IER, and ISNR.    King 930226
c           9.17  XHDRIT: Initialize header to avoid nulls, and call new library
c                         routine sub_char to remove nulls from C-file header.
c                         Herring/King  930311
c                 CDTRED: Initialize entire data array.  King 930311
c                 CVERSN: Correct comment on sub_char.  King 930316
c           9.18  READD : Change call to ORDSIT to match FIXDRV version moved
c                         to /lib.    King 930420
c                 GROUPN: Remove unused variable.
c                 Makefile and directory: Remove unused XDTRED.  King 930420
c           9.19  XHDRIT: Add year to NSNPRN call.   King 931019
c                 CTOX, CHDRED: Increase formats to 32 satellites.  King 931018
c                 XDTRIT, GROUPN: Remove unused variables.  King 931018
c           9.20  CTOX, CDTRED: Allow selection of epochs, channels,
c                         and partials.  King(Canberra) 931208
c           9.21  XHDRIT: Change NSNPRN call to subroutine.  King 940127
c                 Makefile: Change HP compiler switches from +e +E1 to +U77.  Fang/King 940506
c                 Makefile: Shift ranlib to be executed only once.  Herring/King 940623
c           9.22  CTOX, CHDRED, XHDRIT: Write output X-files in GPST.   King 940726
c                 CDTRED:  Initialize partials to zero before reading to avoid carry-over
c                      values when a satellite sets.   McClusky/King  940726
c           9.23  CTOX: Batch mode implemented.   King 941006
c                 XHDRIT: Add space and 'CTOX' to text written in header.  King 941006
c           9.24  XHDRIT: Declare undeclared literal.  King 950103
c                 CHDRED: Fix format to satisfy XL compiler on RISC machine.  King 950103
c                 CDTRED: Change subroutine name in error message.  King 950105
c           9.25  XHDRIT: Replace call to NSNPRN by call to SVNAV_READ.  King 950403
c                 XHDRIT: Add hrs,min to SVNAV_READ call.  King 950405
c                 XHDRIT: fix mismatch in ihour, and min in SVNAV_READ call.  McClusky 950510
c                 Declare variables explicitly:  CTOX, BEAT, CDTRED, CHDRED, DIAGNO, GJ
c                     GROUPN, NMPNT, OPNCHK, POLYFT, READD, XHDRIT. King 950511

c           9.30  New C-file format: Changes to CTOX, CHDRED, and CDTRED  for new arguments
c                     to/from library read routines (all except readc3).  Add new routine
c                     CTOX92 to translate old-style C-files, containing explicit copies of
c                     the old chdred, cdtred, and library read routines.  Also add new
c                     routine CTOC92 to translate old-style to new-style C-files.  Add
c                     to Makefile.
c                 Makefile:  Make all .a file names the same.  King 950608
c                 CHDRED:  Fix spelling of 'srpmod'.  King 950613
c                 CHDRED: Fixed bug in format descriptor. McClusky 950701
c                 CTOC92: Fixed bug in format descriptor. McClusky 950701
c           9.31  CHDRED, CTOX, GROUPN, PLOTC:  Add antenna reference point offsets to new C-file
c                     format and calls to CHDRED (but not yet to XHDRIT).  King 950717
c           9.32  CDTRED: Correct bug in condition of bad jdoy from READC4.  Sanli/King 950719
c                 CTOX, GROUPN: Correct spelling of 'iscrn'.  Sanli/King 950719
c                 READD: Remove unreachable statements and fix error condition branch for
c                      reading D-file.  Sanli/King 950719
c                 CHDRED, CTOC92:  Change variable names itide/iut1pol=>ietide/isptide
c                      for conformity with doc, lib, model, and solve.  King 950720
c           9.33  CDTRED: Fix mismatch in (isnr/amp) calling arguments for READC5.  King 950725
c                 NMPNT: Remove redundant declaration.  Sanli/King 950728
c                 CTOX92(CDTRED92): Correct bug condition in bad jdoy from READC4.
c                     Sanli/King 950728
c                 CTOX92: Fix spelling 'iscren'=> 'iscrn'.   Sanli/King 950728
c                 BEAT: Fix compiler complaint on declaration order.  Sanli/King 950802
c           9.34  CHDRED: Fix printout of (unused) slip quantities.  Sanli/King 950807
c                 CTOX92: Comment unused variable.  King 950807
c           9.35  CTOX: Variable ranlib to make compatible with Solaris 2.  King 951208
c                 Makefile: Remove extra spaces after DORANLIB at bottom.  Fang/King 951221
c           9.36  Added report_status error handling routines too:  CDTRED, CHDRED, CTOX
c                    CVERSN, GJ, GROUPN, JULDAY, NMPNT, POLYFT, READD, XDTRIT. McClusky 960223
c                 CTOX, GROUPN: Disallow normal pointing in CTOX (can use AUTCLN to do
c                    it more rigorously; also, a code bug in GROUPN call to XHDRIT. King 960226
c                 GJ: Change order of declarations (DEC compile bug) Tregoning 960527
c           9.37  Changes for new kf/gamit structure and Makefiles from unimake.  King 960718
c                      All routines:  removed trailing blanks and replaced lib/*.fti
c                      includes by ../includes/*.h.
c                 Removed unused OPNCHK.  King 960719
c                 CTOX92: Declare calling argument explicitly.  King 960815
c           9.38  G77 compiler changes:  CDTRED, CTOX, CTOX92, DIAGNO, RADDEG
c                    XHDRIT,RADDEG,DIAGNO: Initialised previously uninitialised variables,
c                    removed unused variables, and explicitly typed all unexplicitly 
c                    typed variables. McClusky 970113
c                    READD, NMPNT, GROUPN, CDTRED.  McClusky  970121  
c           9.39  READD: Replace call to lib/ordsit with call to lib/sort_string.   King 970307
c           9.40  Makefile:  Move READD to /lib.   King 970321
c           9.41  CDTRED:  Fix passing of data_flag from readc5, and print in octal.   Herring/King 970831
c           9.42  CDTRED: Fix call to report_stat.   King 980227
c           9.43  XHDRIT: Add SV antenna offsets to svnav_read call.  King 980615
c           9.44  CTOX, CTOC92, CTOX92: Replace sb gamit/lib/pickfn.f with function libraries/comlib/pickfn.c. King 980720
c                 CTOX, CTOC92, CTOX92: Add length restrictor after pickfn call.  Herring/King 980905
c           9.45  CHDRED: Change print format for more parameters.  King 980907  
c                 CTOX: Add trap for input/output names the same.  King 980907  
c                 CHDRED: Add 'norb' to lib/readc2 call. King 980911
c                 CDTRED: Add 'svclock' to lib/readc5 call.  King 981002
c                 CTOC95: Add 'norb' to writc2 and 'svclock' to writc5 calls.  King 980911/981002
c                 CTOC95: Replace imbedded readc2.f with 1995 version and add 1995 readc5.  King 981007
c                 CDTRED, CTOX: Add ampl1,ampl2 to dump  King 981208
c                 CTOX: Call pickfn with only lowercase 'c' and 'd' to avoid DPH files.  King 990225 
c           9.46  CTOX:  Fix length of wildcard for c-, d-file search; fix read of iterm input for 
c                      C-file character.  Fang/King 990414
c           9.47  XHDRIT:  Changes for 4-digit years.    King 990805
c           9.48  CDTRED, CTOX92: Correct undefined ihr,min,sec.  King 010103
c                 CDTRED: Initialize a more arrays to avoid ugly prinout.  King 010112
c                 CTOX, CDTRED: Fix to handle decimation in dumping (worked for x-files only).  King 010216
c           9.49  CTOX: Fix bug in calculating stop when dumping all epochs.  King 030507
c                 CDTRED: Change name of modkin.h variable from kradk to krad.  King 030509
c                 CTOX, CTOC92, CTOC95, CDTRED, XDTRIT: Make latr,lonr,radius local (no longer in modkin.h).  King 030521
c                 Makefile:  Remove ctox92 to avoid making updates.  King 030521    
c                 CTOX: Fix extra comma in cdtred call. Herring 030815
c                 XHDRIT: Fix number of blanks assigned to format.  King031027   
c           9.50  CTOX, NMPNT: Fix missing calling arguments for XDTRIT (longstanding bug).  Herring/King 040416 
c           9.51  CHDRED: New c-file format (record 2).  King 050201
c                 CTOX, XHDRIT: Remove 'extra' variable from calling sequence. King 050201
c                 Makefile remove obsolute ctox92, ctoc92, ctox95x.  King 050201 
c                 CHDRED: Rename 'ichan', 'ischan' to 'nsat','isprn' to match rest of GAMIT. King 050202     
c                 CHDRED: Fix format for SV ant print.  King 050225
c           9.52  CTOX, CDTRED, XDRRIT: Change modkin.h to model.h.  King 060823
c                 CDTRED, CTOX:  Remove declation of model.h variables 'pres' and 'temp'.  King 060829
c           9.53  CVERSN: Remove 'iun' argument for lversn.  King 070416
c                 CDTRED, NMPN, XDTRIT: Remove unused labels.  King 070907
c                 CTOX, CVERSN: Remove unused 'iun' argument.  King 070907  
c                 CTOX, NMPNT, XDTRIT: Removed unused calling argument.  King 070910 
c                 CHDRED: Remove maxpar dimension of 'preval' in favor of maxprm.  King 080117 
c                 CDTRED: Increase format width for printing receiver clock value.  King 080123
c           9.54  CTOX: Change 'readd' to readdf' for rename in /lib.  King 090415 
c           9.55  CTOX, CDTRED, NMPNT:  Changes to accomodate new model.h.  King 100212 
c           9.56  CHDRED, CDTRED: Add new variables to C-file records 2 and 4. King 100827
c                 CTOX: Remove declaration of 'swver', not in model.h. King 100908
c                 CTOX, CDTRED, CHDRED: Fix bug and confusion in assigning logical unit numbers;
c                   remove 'iup' debug feature.  King 1009014
c                 CHDRED: Add back SV antenna offsets.  King 101111
c           9.57  CHDRED: Add new variables to C-file record 2. King 130115
c           9.58  CDTRED: Fix bug in testing for no data. Herring 140325
c           9.59  CHDRED: Add new variables to C-file record 2.  King 140327
c                 CHDRED: Add new variables to C-file record 2.  King 140327
c           9.60  CTOX, CDTRED, CHDRED, XDTRIT, XHDRIT:  Modifiy C-file format and put 
c                   variables in common in model.h.    King 141206 
c           9.61  XHDRIT: Add start/stop times to svnav_read call. King 150520
c           9.62  XDTRIT: Fix epoch number in x-file when decimating.  King 160420
c                 CTOX, CDTRED,CHDRED, XDTRIT, XHDRIT, Makex:  Fix variables for writing x-file header;
c                   no longer write the kinematic information; remove routines for normal points.  King 160420/160427
c                 CDTRED: Fix formats for printing jdobs and partials in dump.  King 160427 
c                 CDTRED: Add two decimal places to the 'Delay' printout. King 160822
c           9.63: Fix incorrect declaration of 'l' coordinate variables.  Herring 180309      
c           9.64  CTOX, CDTRED, CHRED: Add new common includes/units.h to provide unit numbers
c                   removed from includes/model.h. King 180320
c                 CHDRED: Temporarily declare obsolete 'skd' from C-file locally since removed
c                   from includes/model.h.  King 180322
c           9.65  CTOX: Change 'maxnet' to 'maxsit' since multisesson no longer supported. King 190503
c           9.66  XHDRIT: Added dummy for antenna power in svnav_read call.  Herring 190702
c           9.67  CHDRED (different from other CHDRED routines):  Updated for L1/L2 Satellite PCO values in c-file
c                   svantdx(2,3,nsat) 200126, Added antdaz. Herring 2020/02/05.
c                 CTOX: Updated 32I to 50I to allow for 35 Beidou satellites TAH 200618
c
      RETURN
      END

