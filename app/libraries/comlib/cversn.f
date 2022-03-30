      subroutine cversn
                                                         

      implicit none
c                Documentation but not yet called by anything.

c comment this out until called to avoid compiler warning -rwk 070405
c    5 format ('2.32 of 2021/01/13 20:27 UTC (',a,')')

c version 0.0    First version with cversn documentation 
c         0.1    Changes to allow compilation under G77
c                CVERSN: New version routine
c                FSEEKG_G77 FSEEKG_SUN: Wrapper subroutines for the function and intrinsic
c                  subroutine versions of the system fseek routine.  McClusky 960109 
c                FMPRUNPROG: Changed fseek function call to fseekg subroutine call.  McClusky 960109
c                MULTIREAD, FMPNEXTMASK: Initialised uninitialised variables found
c                   by G77.  McClusky 960109
c         0.2    READ_LINE: Change dimension of 'value' from 1 to 8 to allow bounds checking
c                   King  970111  
c         0.3    Makefile; FSEEKG_SUN, FSEEKG_G77: wrappers for function (Sun,HP,etc) or subroutine (G77);
c                   FREEG_G77, FREEG_SUN.C, MALLOCG_G77, MALLOCG_SUN:  Call as Fortran function (Sun etc),
c                   or C subroutine (G77).  McClusky 970325 
c                FMPRENAME_SUN, FMPRENAME_G77: wrappers for RENAME Sun/HP/etc F77 library function or G77 
c                   intrinsic subroutine.  McClusky  970325   
c                FMPSHORTNAME_SUN, FMPSHORTNAME_G77: wrapper for GETCWD Sun/HP subroutine or G77 intrinsic
c                   function.  McClusky 970325   (G77 will soon change to a subroutine.)
c         0.4    SYSTIME:  Declare IDATE external to satisfy DEC.  King 970805  
c                OPEN_LU:  Add check for available unit number.  King 970828 
c         0.5    REPORT_STAT: Add ORBITS/yawtab and UTILS/autecl to list for GAMIT.ext.  King 980116
c         0.6    Makefile, FCEM: Remove fcem.f (not used in gamit or kf).  Herring/King 980302
c                REPORT_STAT: Change id for yawtab and autecl to just program name.  King 980313
c         0.7    REPORT_STAT:  Fix bug in checking AUTECL for assigning output file to GAMIT.ext.  King 980501 
c         0.8    Makefile, GETDIR, PICKFN, IRENAME: C routines to replace Fortran routines that use
c                  Fortran system calls that are non-standard.  Fang/McClusky/King 980720
c                  Move getdir.f and pickfn.f from gamit/lib for use with Sun OS/4 (missing 
c                    fnmatch used by the C routines.  King 980721
c                Makefile, LJUST, LUNIT: Move from gamit/lib to be accessible by getdir,pickfn. King 980721
c                GETDIR: Fix bug in C routine when used with pathnames.  Fang 980805
c                GETDIR: Cast NULL before assignment at line 65 to satisfy gcc compiler.  Tari/King 980807
c                PROPER_RUNST: If help file not found, tell the user what to do.  Herring 980813
c                GETDIR, PICKFN: Fix trap for exceeding # files.  Fang 980820
c                FERROR, Makefile: Move from gamit/lib.  Morgan/King 980901
c                GETMAC, IRENAME:  Add 'void' statement (some arguments dumy) to avoid compiler
c                   warning (IRIX on SGI).  Morgan/King 980901
c                PICKFN_G77: Removed calls to ftell and fseek. Were screwing up LINUX version. McClusky 980914
c                GETDIR_G77: Removed calls to ftell and fseek. Were screwing up LINUX version. McClusky 980914
c                PICKFN:  Replaced C versions with new code to handle piped input problems.  Fang/King 990122
c                FMPSHORTNAME: Declare type of function getcwd.  King 990128
c         0.9    RAY, Makefile:  New routine for dirunal pole, UT1.  Herring/King 990319
c                PICKFN_xx.C: Add statement to remove trailing blanks from filename.  Fang/King 990324
c         1.0    GETDIR_SUN.F, PICKFN_SUN.F: Replace fcheck (in gamit/lib, not accessible) by inquire.  King 990709
c         1.1    REPORT_STAT: Add OCTTAB to list of 'GAMIT' programs.  King 000225 
c         1.2    New routines for SGI from Lada Dimitrova via Steve Shapiro:  BIT_UTIL, FREEG, FMPOPEN,  FMPSHORTNAME,
c                    INKEY. Link to sol/hp:  CAND, COR, FMPRENAME, FSEEKG, GETDIR, GETMAC, IRENAME, MALLOCG, PICKFN  King 000814  
c                Makefile:  Make FMPOPEN system-dependent.   King 000814
c                GETDIR, PICKFN:  Declare void (name argument not assigned) to avoid IRIX warning.  King 000816
c         1.3    GETDIR: Increase number of files to 500.  McClusky/King 010604  
c         1.4    RCPAR: Detect if getarg starts at 0 or 1 for first argument, to be compatible 
c                   with F90.  Herring 010610    
c                FSEEKG, BIT_UTIL:  Link _ibm versions.  Matheussen/King 010731 
c                MALLOCG,
c         1.5    SYSTIME: Added a new routine systime_ifc.f since ifc has different implementation / calling sequence 
c                for function idate than compilers like g77 etc. McClusky 031022
c                MULTIREAD: changed declaration of ivalue from i*2(1) to i*2(*) for campatability with ifc. McClusky 031022
c                REPORT_STAT: Add trap and message for iostat 101 with g77.  King 031111
c         1.6    SYSTIME: Removed architecture dependent systime routines since ifort has now fixed calling sequence McClusky 040218
c                FERROR: Added a new routine ferror_intel.f since ifort now has different name for function perror
c                        called for_perror. all other compilers like g77 etc use perror (go figure!). McClusky 040218
c                 ALL: _ifc routines moved to _intel naming convention. McClusky/King 040218
c         1.7    GET_COMMAND: Changed return for 'ALL' selection from 999 to 999999 Herring 040309 
c         1.8    Makefile, MEMASSIGN (new): Replace malloc and free routines by new routine to handle 64-bit addressing. Herring 040801
c                GETDIR_G77: Add include sys/types.h for Mas OSX.  Celli/King 050309
c                INKEY_G77:  Substitute POSIX-standard 'termios' for 'termio' (no longer link to inkey_sol.c). Celli/King 050309
c         1.9    INKEY:  Use SGI version for g77 and intel (ok for MacOSX). Cp to g77, link sgi and intel. Herring/King 050314
c         2.0    REPORT_STAT: Change to write status/warning/fatal messages for MAKEXP, MAKEJ, MAKEX, and FIXDRV 
c                   into the GAMIT files, not ones for each module.  King 060628
c         2.1    ADDED _gftn system dependent source identifier to allow gfortran compiler to be selected. McClusky 060818
c                FMPRUNPROG: system routine external statment commented out since this is an internal intrinsic. McClusky 060818
c                SYSTIME: system routine external statment commented out since this is an internal intrinsic. McClusky 060818
c                GET_COMMAND: renamed to GET_CMD to avoid conflict with F95 intrinsic finction get_command. McClusky 060818
c                FMPRUNPROG:  made system dependant routine to allow fseek to be commented in gftn version. McClusky 060823
c                REPORT_STAT: Add 'GRDTAB' to list of file names changed to 'GAMIT'.  King 060909
c         2.2    Makefile, SYSTIME: Make again OS-dependent since HP uses an external for 'idate'.  Herring/King 061228
c                FMPRUNPROG: Make HP version different since it uses an external for 'system'.  Herring/King 061228
c                FMPRUNPROG: Make Intel ver linked to HP since it needs external for 'system'.  McClusky 070118
c         2.21   REPORT_STAT: Add MSTINF/MSTINF2 to programs writing to 'GAMIT' status, warning, fatal. King 070130
c                Makefile, EXECUTE:  Moved 'execute' from kf/gen_util since it is OS-dependent.  Herring/King 070201
c         2.22   FMPINITMASK, IFBRK, LOGLU, FMPOPEN:  Add dummy statements of avoid a gfortran warning. King 070416
c                FMPOPEN: Removed dummy statement to avoid segmentation violation. King 070420 
c         2.23   REPORT_STAT: Expand search of available units (needed for 99 stations).  McClusky/Herring/King 080124  
c                FMPPURGE: Fix undeclared variable. Herring 080207
c         2.24   Makefile, GEOD_to_XYZ, XYZ_to_GEOD, XYZ_to_NEU: move routines from kf/gen_util so that
c                  they are accessible by GAMIT.  Also add /includes directory to /libraries, and move
c                  const_param.h from kf/includes (put link into kf/includes).  King 080307  
c                PICKFN_SOL, PICKFN_HP: Modified to allow more files. Herring 080512 
c                PICKFN_SOL: Change // comment syntax not supported on old Solaris compiler.  Herring 080717
c                Makefile, DECYRS_TO_YDHMS, DS2HMS, LEAPYR : Move decyrs_to_ydhms, ds2hms, leapyr from gamit/lib.  King 081118
c         2.25   REPORT_STAT: Remove restriction on only 6-character program names.  Herring 091228
c         2.26   Makefile, new PMUT1_OCEANS: IERS2010 short-period EOP model. Herring 110210
c                REPORT_STAT: Add NGSTOT to modules for message written to GAMIT [status/warning/fatal]. King 110416
c         2.27   Makefile, SD_COMP SD_COMBBD (moved from gamit/lib), TIDE_ANGLES (separated from SD_COMP), GIPSON (new):
c                   New code for Gipson short-period EOP model. Wang/King  151022
c         2.28   REPORT_STAT: Add ORBFIT to modules for message written to GAMIT [status/warning/fatal]. King 160812
c         2.29   Added MEAN_POLE to compute mean or secular pole position for pole-tide model Herring 200219.
c         2.30   NAME_TO_BLK: Added GLONASS M+, K1A and K1B body types. TAH 200603.
c         2.31   GETDIR_{G77,HP,SOL}: Added fnmatch.h to list of includes to accommodate call of fnmatch;
c                INKEY_{G77,HP,SOL,SUN}: Added unistd.h to list of includes to accommodate calls of isatty and read;
c                IRENAME_{HP,SOL}: Added string.h to list of includes to accommodate call of strtok;
c                PICKFN_{HP,SOL}: Added fnmatch.h, ctype.h and stdlib.h to list of includes to accommodate
c                                 calls of fnmatch, isdigit and atoi, respectively. MAF 20201130
c         2.32   NAME_TO_BLK: BEIDOU-3I body type added (not coded in earthradTUM.f yet),
c                replacing dummy placeholder "BEI NEXT 07" (ID 37). MAF 20210113
c
       return
       end
