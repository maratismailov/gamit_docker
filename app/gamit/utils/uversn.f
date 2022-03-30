Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego, 1994. All rights reserved.

      Subroutine UVERSN(iun,vers)

      character*40 vers
      character*10 getmac,machin
      character*45 libver
      integer iun
c     function to return non-blank length of string
      integer nblen
             
c     get the machine name
      machin = getmac(1) 

c     write the machine name and utilities version into 'vers'
      write (vers,5) machin(1:nblen(machin))    
c     ** update the following each time you change the program
    5 format ('10.06 of 2021/04/15 14:36 UTC (',a,')')
          
c     get the library version
      call lversn(libver)
           
c     wrrite the versions to the screen or log file
      write(iun,*) ' '
      write(iun,*) 'UTILITIES v. '//vers
      write(iun,*) 'Library v. '//libver


c   Old History to be filled in gradually.

c   9.11  XTORX:  Merge Oral/Herring changes of Jul 92 with Feigl changes
c                          of Jan 93.   King 930126.
c   9.12  XTORX:  Comment out code for bit-2 of loss-of-lock indicator.
c                 Write out only weighted data.   Fix bugs.  Herring 930128
c   9.13  XTORX:  Add code to calculate RINEX antenna heights from either
c                 the X-file values or station.info.   King 930212/930215
c         COUNTX, PLOTK : Add antcod to XHDRED call.   King 930212
c         XTORX:  Remove unused variables.   King 930219
c   9.14  XTORX:  Modified two comment statements in RINEX header Bock 930308
c   9.15  CCX, Makefile:  Fang/Oral program added from gs/contributions/oral.
c                       King 930407
c         XTORX: Add Turbo-Rogue and fix call to ORDSIT.  Fang  930502
c   9.15  Remove HPRN from MIT active and stdrel.
c         CCX: Remove comment on # sats restriction.   King 931018
c   9.16  ARGO2FIC: Allow more read errors and fix format.  Oral/King 931220
c   9.17  XTORX: Fix bad JULDAY bug with empty file.  King 940124
c         EMNTAB: New program to write GAMIT Sun, Moon, nutation tables from
c             a PEP N-body file.  King 940125/940204
c         Makefile: Change HP compiler switches from +e +E1 to +U77.  Fang/King 940506
c   9.18  XTORX: Add calling argument to library RSTNFO.   King 940511
c         COUNTX: Increase length of filename to allow full path.  King 940526
c   9.19  Add DECYR_ATM and POLY01 for Oral scripts; change Makefile.  King 940623
c   9.20  PLOTK:  Change to use GPST vice UTC.   King 940728
c         XTORX:  Fix iyr=0 bug when last record has no obs.  Fang/King 940818
c   9.21  RX2APR: New program to convert RINEX header positions into globk
c           .apr files.   King  941006
c         XTORX: Add LEICA information.  Bock/Genrich 941108
c   9.22  XTORX: Fix bug with UTC=>GPST conversions; replace week,sow with jd,jdoy,sod
c            throughout.   King 941114
c         DECYR_ATM:  Fix minor compile warning in write statement.  King 941201
c         XTORX: Fix minor logic problem caught by RISC compiler.   King 941201
c         ARGO2FIC, CCX, WBSLFILT, XTORX :  Fix format problems caught by RISC compiler.  King 941205
c         MERGE_RINEX: New program to merge multiple rinex files into single rinex file
c                      spanning one day 00:00:00.000 23:59:59.999. McClusky 950201
c         XTORX: Add Leica receiver; reimplement changes of 941114 (lost?).  King 950210
c         CCX:   Remove illogical return at line 465.  King 950331
c         MREGE_RINEX: Fix minor HP compiler complaints.  King 950331
c   9.23  AUTECL: New utility to create autcln.cmd edit lines to remove eclipse and
c                 post eclipse attitude recovery data.  McClusky 950404
c         AUTEDT: New utility to create autcln.cmd edit lines to remove data
c                 flagged during a cview session.  McClusky 950407
c         CCX: Increase dimensions of SV arrays to 99 to avoid dependency on svnav.dat. Oral 950212.
c                 Merge with MIT changes of 941205 and 950331.   King 950410.
c   9.24  CCX: Allow for GPST X-files.  King 950418
c         CCX: Don't let data interval be less than interval of input x-file.  King 950420
c         AUTEDT: fixed a bug that caused grief when site nmaes were numbers McClusky 950517
c         MAKEFILE: Added erptab to makefile.  Tregoning 950525
c         ERPTAB: Add 1 significant fig. to ut1. pole. tables and
c                 convert either UT1-AT or TAI-UTC into ut1. table  Tregoning 950525
c         AUTECL: Fix bug that missed reading hours and mins from session.info. McClusky 950525
c    9.25 CCX: Fix problem reading session.info.   King 950526
c    9.26 EMNTAB: Allow rotation to J2000.   King 950616
c         ERPTAB: Increase parameter dimoension from 1000 lines to 2000  Tregoning 950622
c    9.27 ARGO2RX: Fix format statements (DEC compiler complained).  Sanli/King 950719
c         CCX: Cannot have null strings on DEC and IBM.  Sanli/King 950719
c         EMNTAB: Remove blanks at end of too-long record. Sanli/King 950719
c         PULL_BIAS: Open statements must have 'status' for DEC.  Sanli/King 950719
c         XTORX:  Remove unused statements.  Sanli/King 950719
c ------- (SIO 9.22)
c         XTORX: Fix time conversion for GPST X-files & update version number.
c                Hard wire day number from X-file for RINEX file.
c                (year boundary is still a potential problem).
c          Bock 950211
c         XTORX: Add LEICA information, again.
c         Fix time conversion for UTC X-files near day boundary.
c --------  Bock 950213
c         XTORX:  Replace MIT version with SIO version, after merging 950719 changes.
c           (Save MIT version with doy,sow in ../obsolete/xtorx.f.mit950719.)
c         AUTECL:  Add status=unknown for autcln.cmd open.  Sanli/King 950728
c         AUTEDT:  Add blank to zero-length string.  Sanli/King 950728
c         MERGE_RINEX: Remove redundant declarations.  Sanli/King 950728
c         RX2APR: Add missing unit number in screen write.  Sanli/King 950728
c         ARGO2RX: Change double quotes to single for DEC.  Sanli/King 950728
c         XTORX: Modified to write correct rinex name for doy numbers < 100. McClusky 950728
c         RX2APR: Missing declaration of line_count added to editxyz McClusky 950728
c         ARGO2FIC:  Rearrange declarations to satisfy DEC compiler.  Sanli/King 950802
c   9.28  PLOTK: Add words to stops.  King 950812
c         XTORX: Don't allow overwriting of RINEX files--too dangerous.  King 950823
c         SOLN_INFO: fixed bug in eop a priori constraints  Tregoning 950901
c         WRTSNXHD: Cleaned up FILE/COMMENT block  Tregoning 950901
c         RDHHED: Fixed declaration of date variable to be integer  Tregoning 950901
c         POLY01: Fixed the way the unc_mean was calculated. King McClusky 950904
c         Makefile, SVDLIN SVDSIG: Added two new routines from /sr/contributions/oral/progs. McClusky 950904
c         Makefile: Merge SIO and MIT changes.  King 950907
c         WRTANT_INFO: Output only one default SITE/GPS_PHASE_CENTER record per antenna
c                      type (sinex convention)   Tregoning 950909
c         WRTCOV_APR: Change a priori covariance block name from CORR to COVA  Tregoning 950909
c         RDHHED: Change read of hfile header to search for "session:" line  Tregoning 950909
c         COMP_TIME: fix bug in computation for cfiles > 24hrs  Tregoning 950909
c         HEADERS: Fix dimension bug in declaration of est_type   Tregoning 950911
c         ERPTAB: Increase maximum file name length  Bock/Fang 950912
c         WRTCOV_APR, WRTSOL_APR: fixed bug in converting llr to xyz covariances  Tregoning 950914
c         SOLN_INFO, WRTSNXHD, WRTSOL_APR : Remove extraneous commas in format descriptors.  King 950922
c         WRTCOV:  Remove unused variable.  King 950922
c         OPEN_SNX: Minor changes to FILE/REF help info  Tregoning 950922
c         Makefile: Change htosnx.d to htosnx.a for /test version.  King 950925
c         ERPTAB:  Increase dimensions to handle 15-yr, 1-day interval USNO tables.  King 950925
c         WRTSNXHD:  Change format statement in FILE/REF block to match new hfile format  Tregoning 950925
c         WRTSNXHD:  Changed it again - for the last time hopefully!  Tregoning 951002
c         RD_COV:  Fixed bug in handling of velocity parameters (not available in hfiles)  Tregoning 951002
c         SOLN_INFO: Fix bug in eop indices.  Tregoning 950925 (resent 951003)
c         MERGE_RINEX: fix bug to read file when there are 7 observables, allow decimation
c                      of input rinex files    Tregoning 951003
c         MERGE_RINEX: change 'err=' to 'end=' in read statements (for HP)  Tregoning 951003
c         MERGE_RINEX: add 5th argument - to limit # observables (and hence rinex file
c                      size). 5th argument is optional  Tregoning 951005
c         MERGE_RINEX: fix bug when duplicate time lines in rinex file  Tregoning 951012
c    9.29 MERGE_RINEX, ARGO2RX, XTORX:  Pass receiver and antenna serial numbers as
c            20-character strings, not integers.  King 951019
c         COUNTX, PLOTK, XTORX: Add receiver and antenna ids to lib/xhdred call.  King 951019
c    9.30 CONVERT_ANTMOD: New routine to manipulate antenna VPC tables. McClusky 951023
c         MERGE_RINEX: fix bug in handling incomplete rinex files   Tregoning 951023
c    9.31 XTORX:  Remove duplicate declarations.   King  9501103
c    9.32 COMP_TIME, RADDEG:  Make argument types consistent for intrinsic calls
c            (caught by DEC compiler)
c         SOLN_INFO:  Declare array dimensions before used.  Sanli/King 951111
c    9.33 Makefile: Variable ranlib for compatibility with Solaris 2.  Fang/King 951208
c         MERGE_RINEX: Avoid copying nulls; clarify argument in documentation.  King 951228
c    9.34 MERGE_RINEX: allow reading of rinex files greater than 2880 epochs long. McClusky 960121
c         MERGE_RINEX: Fix DEC compile bug in mod statement.  Tregoning 960528
c         AUTECL: Rewrite free format char read ( for DEC).  Tregoning 960528
c         AUTECL: Remove unused variables.  King 960529
c         CONVERT_ANTMOD: Define IPRNT variable.  Bock 960618   
c    9.35 Changes for new kf/gamit structure and Makefiles from unimake.  King 960718
c             All routines:  removed trailing blanks and replaced lib/*.fti 
c             includes by ../includes/*.h. 
c         RX2APR: Fix truncation of variable name.  King 960816 
c         EMNTAB:  Use variable rather than constant in reading kbd to circumvent
c             Solaris 2.5 compiler bug.  Murray/King 960827   
c    9.36 MERGE_RINEX:  Correct bug in setting start time for RINEX files.  Tregoning/King 961015  
c         AUTEDT, WBSLFILT:  Fix format statements for gcc compiler.  Tregoning/King 961017
c         RX2APR: Use 'eqv' not 'eq' in logical comparison, for gcc compiler.  Tregoning/King 961017
c         Makefile:  Add in svdsig (omitted by mistake in Makefile.generic).  King 961108
c    9.37 XTORX:  Add calls to read_rcvant to get full names; add 'xtorx.out' file, 
c               restructure some.  King 970113  
c    9.38 G77 compiler changes: ADDVEC, AUTEDT, AUTECL, ARGO2FIC, CCX, ERPTAB, EMNTAB
c             RD_COV  WBSLFILT, XTORX, RDHHED, RDSIT_INFO, READ_APR, SOLN_INFO: Initialised 
c             previously uninitialised variables, removed unused variables, and explicitly 
c             typed all unexplicitly typed variables. McClusky 970121
c         XTORX:  Fix bug in checking dattyp.  McClusky/King 970121
c    9.39 XTORX: Set RINEX version number to '2' since END OF HEADER written.  King 970303
c         XTORX: Replace call to lib/ordsit with call to lib/sort_string.   King 970307  
c         Makefile: Move XTORX to /makex directory.  King 970321
c         CCX: Rename function itime=>ktime to avoid conflict with Fortran system call. King 970321
c         POLY01: Fix calculation of wrms.  Demir/King 970331
c    9.40 PLOTK: Change 'nbatch'=> 'clkft' in call to lib/clkera for clarity.  King 970501  
c         PLOTK: Add interactive query for order of polynomials for jump and no-jump fits,
c              and add calling arguments to clkera.  King 970522
c    9.41 AUTEDT: Routine moved to CVEDT
c         CVEDT: Routine autedt moved to cvedt, and output now written to std output not
c                the autcln.cmd file. McClusky 970618 
c         AUTEDT: Routine moved to cvedt
c         AUTECL: Output now written to std output not the autcln.cmd file. McClusky 970618
c         RX2APR: Initialised variable "rs". Necessary for the G77 compiler. McClusky 970722  
c    9.42 Makefile, VEXCLUDE,CHANGE_EXCL, READ_Q: New utility (vexclude) to remove excluded 
c             sites/satellites in solve run from vscan.out residuals for cview  Tregoning/King 970831
c         CONVERT_ANTMOD, RD_COV: Break Hollerith string to avoid DEC warning.  King 970920
c         PROFILE: New routine to compute velocity profiles from getrel output. McClusky 970918
c         Makefile: Added routine profile to Makefile
c    9.43 POLY01R: Similar to routine poly01, now returns residuals to fit when requested. McClusky 970926 
c         Makefile: Added routine poly01r to Makefile 
c    9.44 ARGO2FIC:  Fix output FICA file name (ssssyddd.fic vice sssswwwwy.ddd).  King 971024
c         PROFILE: Fixed bug in output of map projections of profile and extents.  McClusky 971205
c         ARGO2RX: Rework to run with (new) sh_argo2rx using an input stream, and to handle multiple
c                  RINEX output files.   King 971217
c         MERGE_RINEX: Remove debug and reduce number of lines to screen.  King 971217
c         ARGO2RX: Fix bugs in opening antmod.dat and in calls to report_stat (segmentation
c                  fault).   King 971226   
c         EMERGE, Makefile: Add SIO program to merge RINEX nav files.  Fang/King 971231
c         AUTECL: Set module for report_stat to 'GAMIT' since called in batch run; skip execution
c                 if a previous step has failed.  King 980106/980116 
c         ARGO2RX: Remove sb wtime in favor of library version; call with GPST rather than UTC.  King 980123
c    9.45 AUTECL:  Don't delete Block IIRs during eclipse.  Tregoning/King 980204
c    9.46 CONVERTC: New routine that does lots of different coordinate conversions. McClusky 980204
c         MERGE_RINEX:  Fix decimation so that input interval not required and also so that
c           times up to 5s from an even 10s will get picked up in decimation.   King 980211
c         ARGO2RX:  Fix bug in passing station name line (failure to initialize string
c           results in header line missing the line id.  King 980211  
c         EMERGE: Fix spelling of 'unknown' in open.  Morgan/King 980211  
c         MERGE_RINEX: Put stop for 0, warning for 1 input file.  King 980216
c         ARGO2RX: Add code to read a globk apr file for approx coordinates.  King 980223    
c         ARGO2RX: Fix nwave2 for MiniMac.  Morgan/King 980225  
c         ARGO2RX: Initialize error counter.  King 980227  
c         AUTECL: Fix bug in call to report_stat.  King 980227 
c         CONVERTC: Initialize error flag.  King 980227 
c         EMERGE: Initialize header counter redundantly (G77 pickiness).  King 980227
c         PROFILE: Replace dsind/dcosd not supported by G77.  McClusky/King 980227
c         AUTECL:  Initialize variables to prevent G77 warnings.  King 980303   
c         MERGE_RINEX:  Clarify comments regarding the max_obs entry.   King 980310
c    9.47 AUTECL: Change report_stat id to 'AUTECL','utils/autecl'.  King 980313
c    9.48 CONVERTC: Fixed convertc to read an globk .apr file. McClusky 980507 
c    9.49 COUNTX:  Replace calls to gamit/lib sb pickfn.f with libraries/comlib function pickfn.c; fix
c             error checking on reads.  King 980720 
c         COUNTX: Correct pickfn call.  Herring/King 980905
c         VECTOR_STAT, Makefile : New routine added need by sh_histogram. McClusky/Herring 990216
c    9.50 POLY01R: Allow sigma editing.  Herring 990304
c         VECTOR_STAT: Fix minor format error.  Fang/King 990305  
c         WRTANT_INFO: Do not use comma as continuation character 
c              (IBM and f90 restriction).  Fang/King 990324
c    9.51 UPD_STNFO: New routine to add an entry to station.info from the RINEX header.  King 990421
c         UPD_STNFO: Allow full pathname for RINEX file.  King 990426
c    9.52 CONVERTC: Fixed bug in conversion of LLH to XYX. Updated help. McClusky 990511
c         MERGE_RINEX: change time rounding from 10 sec to 1 sec when thinning files  Tregoning 990608
c         UPD_STNFO: Add error checking, consistent casefolding, other checks.  Herring 99062
c         UPD_STNFO: Trap too-large length of input format and non-zero ioerr on line decode. King 990805
c         Changes for 4-digit years: ARGO2FIC, ARGO2RX, AUTCLN, CCS, CONV_DATE, J2SAT, 
c           MERGE_RINEX, RDHHED, UPD_STNFO.  King 990805
c         UPD_STNFO: Fix bug in writing out header line.  King 990816
c         Remove: CCX, HTOSNX, J2SAT, TRMM (change Makefile).  King 990817   
c         AUTECL: Read session.info using lib/rsesfo; allow search only on day number.   King 990818
c         ERPTAB:  Y2K bug fixed.  King 990826
c         UPD_STNFO: Fix format for warning of non-unique station.  King 990830
c         EMNTAB: Increase size of arrays.  King 990907
c         UPD_STNFO:  Fix check of whether rcvr firmware is a number.  Herring 990912
c    9.53 CONVERT_ANTMOD:  Old routine for GAMIT and AIUB conversion replaced by new one for
c            NGS tables conversions, to be generalized later.  King 991129 
c    9.54 EMERGE: Fixed to generate merged rinex navigation containing unique entries files. McClusky 000104
c         OCTTAB, Makefile: Create a GAMIT ocean tide (u-) file from Scherneck station table.  King 000217/000224 
c    9.55 ERPTAB: Remove call to check_y2k. Added code to erptab. Report_stat call too inefficient. McClusky 000218
c         ARGO2RX: Fix y2k problem.  King 000321
c         ARGO2NAV: New program to translate ARGO orb file to RINEX n file.  King 000324
c         ARGO2RX: Fix screen display of first epoch (no effect on files).  King 000327
c         OCTTAB: Don't use 'cosd' as intrinsic (not available with g77).  King 000403
c    9.56 DAF2DA, Makefile: Create a GAMIT ocean-tide grid table from Scherneck grid ('daf') files.  King 000427
c         DAFEXPAND, Makefile: Expand a Scherneck 1x1 deg grid to 0.5x0.5 deg.  King 000718
c         OCTTAB: Add byte-swapping of binary grid file for Linux/DEC.  Herring 00809
c         OCTTAB: Change site dimensions from maxsit to maxnet to match lib/readdd.  King 001102
c    9.57 Makefile, CRX2RNX RNX2CRX:  Add Hatanaka compression routines.  King 001206
c    9.58 ARGO2FIC: Fix bug from undefined read-error variable.  King 010103 
c         POLY01R, READQ: Change logic to avoid technically undefined variables.  King 010103
c         ADDVEC:  Initialize variables.  King 010116
c         ARGO2FIC: Remove unused variable.  King 010116 
c    9.59 ARGO2RX, ARGO2NAV, MERGE_RINEX, UPD_STNFO: Changes for RINEX version 2.10:  version and 
c            interval now real, extra decimal place for time of first observation.  King 010301
c         MERGE_RINEX:  Changes for RINEX version 2.10: decimal version number and sampling interval. King 010515
c    9.60 UPD_STNFO: Fix bug in reading pure-number swver.  King 010515  
c         FIXSST: New program to change header info for SST translated incorrectly as SSE.  King 010515   
c         OCTTAB:  Avoid null string.  Tregoning/King 010518   
c         FIC2NAV, Makefile: New program to translated FICA Block 9 to RINEX navigation file.  King 010615
c         POLY01, VECTOR_STAT: Call rcpar rather than getarg to avoid getarg conflict with HP UX-11.  Herring 010706 
c         FIXX:   New program to add the missing LAMBDA values to old x-files.  King 010710
c         FIXX:   Add trap for binary characters in first line of files written from old CTOX.  King 010712 
c         AUTECL, CONVERTC, CVEDT, DECYR_ATM, FIC2NAV, HISTGM, POLY01, PROFILE:   Use rcpar rather than 
c            getarg to avoid problem with HP UX-11.  Herring/King 010718
c         RX2ARR: Remove local rcpar routine.  King 010718   
c         FIC2NAV: Change calling argument to include yr,doy, not output RINEX name, to make sh_fic2nav easier.   King 010723 
c         OCTTAB:  Change order of variables in 'int' expressions to satisfy Leahy compiler.  Morgan/King 010731 
c         FIXX:  Modify to do only the binary-character removal, not LAMBDA fix unless set internal control.  King 010828
c         FIXX:  Modify to remove additional binary characters in the 'Receiver' line of TI x-files. King 011008
c    9.61 ARGO2RX: Allow new-style station.info; change 'swver' from R*8 to R*4 to match rest of GAMIT.  King 020312
c    9.62 OCTTAB: Form ufile name from ifile rather than dfile name (consistent with model/open.f, but neither
c            works when the i-file has 6th character 'a' or is missing). Tregoning/King 020220/020327
c    9.63 OCTTAB: LF95 compiler need extea brackets in expression (2.*-alat1) Ie. (2.*(-alat1)). McClusky 020417
c         POLY01R: Fix problem with zero line passed to sh_baseline.  Herring/King 020418
c         EMERGE:  Change format for reading file names from * to a80 to satisfy DEC.  Wang/King 020515  
c         FIXSST: Add blanking out of firmware version to header changes for mistranslated Trimble SST.  Feigl/King 020617 
c         OCTTAB: Add NAO ocean tide model.  King 020618
c         FIXX: Add some more binary-character recognition code for December 1986 x-files.  King 020628
c    9.64 OCTTAB: Allow reading of an apr, as well as an l-file.  King 020807  
c    9.65 OCTTAB: Fix bug in 020618 changes for NAO model (missing long-period terms in grid for OSO model).  King 020923
c    9.66 Makefile, SITE_ID: Move GEOD_XYZ to /lib and rename GEOXYZ to be called by SOLVE and TFORM. King 021002 
c    9.67 POLY01R:  Fix integer/real problem in calculating wrms.  King 021126
c    9.68 ARGO2RX: Add 'span' to arguments for rstnfo2.  King 030107 
c    9.69 CHECK_LFILE, Makefile: New routine to be called from sh_gamit to get L-file type (old or apr).  King 030215 
c    9.70 NAO2DA, OCTTAB, Makefile: New routine and modifications to create an ocean tide grid file from NAO files.  
c            Matsumoto/King 030305    
c    9.71 CONVERT_ANTPCV, Makefile: New routine to convert between NGS, GAMIT, and ANTEX format antenna phase center
c            models; calls new library routines read_antmod and read_antex.  King 030417
c         ARGO2RX: Remove obsolete 'icall' from call to hisub.  King 030418
c    9.72 OCTTAB: Change calling arguments for lread  King 030519
c         OCTTAB: Set l-file type.  King 030605
c    9.73 CONVERT_ANTPCV: Add conversion PCV model designators 'IGS01' between GAMIT and ANTEX formats.  King 031015
c         ARGO2RX: Fix bug in reading date from ARGO header (line too long).  King 031027
c         AUTECL: Fix string-length in 'eclipse_out' assignment (ifc warning). King 031027 
c         FIC2NAV: Fix string-length in assigning 2-digit yr to RINEX file name (ifc warning). King 031027 
c         PLOTK: Fix string length in assigning x-file name (ifc warning). King 031027
c         PULL_BIAS: Fix string length in assigning input file name.  King 031027
c         SVDLIN: Fix line truncation resulting in single precision. King 031027
c         UPD_STNFO: Fix string length in assigning firmware codes (ifc warning). King 031027
c         WBSLFILT: Remove superfluous characters beyond line 72.  King 031027
c    9.74 CONVERT_ANTPCV: Add 'sinex_code' to call for read_antex_head.  King 040119
c         VEL2HIST, CHICURVE: New programs for (new) sh_velhist to plot a histogram of
c           velocity residuals with the theoretical chi distribution superimposed. McClusky/King 040217
c            freedom and plot this with a histor
c         VEL2HIST: Add an argument allowing scaling of the sigmas.  040317  
c    9.75 ATMTOASC, Makefile: new utility to extract a time series of atmospheric pressure loading
c            displacements from a binary 2.5 x 2.5 grid file, one site at a time. Tregoning 040107
c         GRDTAB, BILIN, INTERP_ATM: new program (+ subroutines) to interpolate 
c            a global gridded binary file to extract atmospheric loading displacements
c            for use in MODEL. Tregoning 040107
c         GEOC_TO_GEOD: moved out of OCTTAB and made into an independent subroutine. Tregoning/King 040107/040331
c         METUTIL: New program to read a z- (met print) file from MODEL and an o-file from SOLVE and
c            write out estimated values of dry and wet zenith delays and precipitable water.  King 040625 
c         POLY01: Increase significant digits in header.  Herring 040714
c    9.76 Makefile, BILIN: Attach bilin.f to interp_atm.f to avoid single/double precision confusion
c           in other uses.  King 040903
c    9.77 ARGO2RX, CONVERT_ANTPCV: Add 'pcncod' argument to read_rcvant call.  King 040930
c         MERGE_RINEX; Fix to allow #SVs > 12 (extra time line); change check_2yk to fix_y2k.  King 041011
c    9.78 CONVERT_ANTPCV: Add iostat check reading NGS header file.  King 041222
c            Makefile: remove obsolete CONVERT_ANTMOD.  
c         METUTIL: Complete rewrite to handle both o- and SINEX file estimates, and both z- and
c            RINEX met file values for pressure and temperature; change name and format of output file. Tregoning/King 041229
c    9.79 ATM_COMP, GRDTAP, INTERP_ATM, Makefile: Changes to use a header on the ascii and binary grid files. Tregoning 050104  
c         DCBTAB: Modify for AIUB format change 4 January 2005.  King 050104
c         INTERP_ATM: Fix problem in computing day of year from grid file.  Tregoning  050110 
c         INTERP_ATM, INTERP_IMF, BILIN, Makefile:  Change name of BILIN to BILIN4 to make clear distinction 
c              from real*8 version in /lib/get_antpcv.  King 050112
c         ATM_COMP, GRDTAB: Change version in binary file header to integer*2 to allow for more information.  Tregoning 050114
c         ATMTOASC, ATM_COMP, GRDTAB: Fix formats for header reads and writes.  King 050114
c    9.80 COUNTX, PLOTK: Remove 'extra' variable from call to lib/xhdred.  King 050129
c         CONVERT_ANTPCV: Add antenna serial number (or SV PRN) to calls to read_antex.  King 050208
c         GRDTAB: Move grid file indices print to after byte-swapping.  Tregoning 050209
c         NUTTAB, Makefile: New routine to write the nutation table from the MHB_2000 series. King 050307.
c         DCBTAB: Fix logic to read the header more robustly.  King 050307
c         GRDTAB: Fix byte-swapping for new header format.  King 050309
c         ATMCOMP:  Fix byte-swapping for new header format; rearrange logic some.  King 050311
c    9.81 ATMCOMP, ATMTOASC: Fix syntax error in write statements.  King 050314
c         AUTECL: Fix module name in report_stat call.  Saba/King 050404
c         INTERP_VMF: Fix bug in finding site, and trap uppercase entry in file.  Tregoning 050406/050407
c         DAF2DA: Allow a different header format; more lat/lon values; uppercase keys for
c            long-period tides in Scherneck daf files; 8-character model ID on first record. King 050411 
c         OCTTAB: Read 8-character tide model from grid or stations.oct header and write into
c            the u-file; allow stations.oct to be missing or blank.  King 050413/050419  
c         GRDTAB, INTERP_ATM: fix bug in reading beyond end of atm grid file.  Tregoning 050418
c         VEL2HIST: Add vel2hist.out print file to better examine outliners; fix some 
c            inconsistencies in applying sigma scalilng; remove _GPS-only criterion.  King 050502
c    9.82 ARGO2RX, CONVERT_ANTPCV:  Add radome to read_rcvant and hisubcalls.  King 050618
c    9.83 ARGO2RX:  Add warnings variable to hisub call.  King 050719
c    9.84 CRX2RNX, RNX2CRX:  Replaced Dec 2000 version of Hatanaka compression routines
c           with new versions 2.4.1.  King 050816 
c         PROFILE: Correct bug in using correlation coefficient.  Vernant 050831
c    9.85 OCTTAB: Fix incorrect declaration of 'tol', passed to get_station.  Shimada 050906  
c         GRDTAB, INTERP_ATMTIDE, Makefile: New code to interpolate grid for atm tides. Tregoning 050718/050906     
c         INTERP_VMF: Trap to catch case when time requested runs beyond the end  of the vmf input file  Tregoning 050721/050906
c         ATMTOASC: Fix bug due to addition of second header line.  Tregoning 050908
c         AUTECL: Change to formatted read of arcout file to allow asterisk.  King 050929
c         IMFTOASC, Makefile: New program to read binary IMF grid file and print values.  Tregoning 050921 
c         INTERP_VMF: Allow upper or lower case, and 4-character or 8-character site names. Tregoning 050929
c         INTERP_ATMTIDES: Remove unused variables.  King 051029
c         IMFTOASC: Remove extranious comma.  Saba/King 051104
c         ATMTOASC: Increase length of atm file name to 15 characters.
c         ATMTOASC: Fix bug in swap_byte call for version  Tregoning 051110
c         ATMCOMP: Allow conversion of a file extending into the next year (369 entries). Tregoning/King 051110
c         ATMCOMP: fixed bug in writing grid size into binary grid file  Tregoning 051114 
c    9.86 GRDTAB, INTERP_VMF: added VMF1 mapping function Tregoning 060702
c         METUTIL: use a priori ZTD from zfile. Output press and temp. Tregoning 060702 
c         GRDTAB, INTERP_VMF, INTERP_VMF1GRD (new):  add reading of VMF1 coefficients from a grid  Tregoning 060702
c    9.87 AUTECL: Change arguements in call to rsesfo (see lib/uversn.f).  King 060815/060829
c         RX2APR: tabs is format line were a problem for gfortran (removed). McClusky 060818
c         CONVERT_ANTPCV: Problem with radome when reading NGS.  King 061011
c         Makefile; remove: ATMTIDE, ATMTOASC, ATMCOMP, DAF2DA, DAFEXPAND,
c             GRDTAB, DAF2DA, INTERP_ATM, INTERP_ATMTIDE, INTERP_IMF, INTERP_VMF,
c             OCTTAB VMFTOASC  (in grdtab or obsolete).  King 061103
c    9.88 AUTECL: Fix calling argument for rsesfo so that only day-of-year is checked (see /lib/lversn.f);
c             fix other logic.  King 061208
c         CONVERT_ANTPCV: Fix bugs in converting a specified antenna or all antennas from an input
c             NGS file containing many antennas. Herring/King 070115   
c         CRX2RNX, RNX2CRX: Update with latest version from Hatanaka.  King 070122
c         METUTIL: Use ZHD rather than pressure if a priori values from z-file rather than 
c             a RINEX met file.  King 070124
c         METUTIL: Use linear interpolation rather than closest-value for pressure/ZHD from
c             RINEX/z-file.   King 070129
c         VEL2HIST: Allow comments in header lines of vel file.  King 070215
c         METUTIL:  Fix bug in adding a value at end for interpolation; fix bug in returning
c            temperature values from RINEX files in Celsius rather than Kelvin. King 070216
c         UVERSN: Remove print unit number from 'lversn' call.  King 070416
c         METUTIL: Fix header units for output met file.  King 070611
c         ARGO2FIC, ARGO2NAV, CML, MERGE_RINEX, PLOTK :  Remove unused statement label.  King 070910  
c         ARGO2FIC:  Add dummy statement for unused calling argument. King 070910    
c         UPD_STNFO: Remove unused calling argument.  King 070910
c         RX2APR:  Fix format statement. King 070910
c         VEL2HIST: Fix syntax errors. King 070911
c         ERPTAB: Minor change in report_stat report.  King 070912  
c         VEL2HIST: Allow setting a minimum as well as maximum sigma, and a geographic box. King 071019
c         VEL2HIST: Fix extraneous declaration of 'inbox'. Casula/King 071106  
c         EMNTAB:  Add JPL DE capability; clean up interfaces and fix for Linux; Makefile to add temporarily TESTEPH. King 071231
c         VEL2HIST: Change variable names and input so that horizontal magnitude is 'h'/'H' (not 'v'/V')
c            and vertical is 'u'/'U' (not 'h'/'H').   King 080109 
c         CONVERTC:  Add Google Earth KML format output.   McClusky 060117/080204    
c         INTERP_LOG_TAB, Makefile: Interpolate log functions and velocities from tsview.  McClusky 080204
c         EMNTAB, TESTEPH: Fix missing declarations in JPL routines. King 080207
c         INTERP_LOG_TAB:  Change 'return' to 'stop' to satisfy HP compiler.   King 080207
c         CONVERTC: Convert longitudes to +-180 for KML output.  King 080227
c         VEL2HIST: Fix undefined horizontal magnitude for print when E, N,or U selected for histogram. King 080310
c   9.89  COUNTX, PLOTK: Change antenna offset arguments for lib/xhdred.f.  King 080509
c         ARGO2RX: Remove calls to old-style hisub and rstnfo; change calls of hisub2 and 
c             rstnfo2 to hisub/rsntfo. King 080513  
c         METUTIL: Remove obsolete comments re model.h common not used here.  King 080513
c         Makefile: Remove CONVERTC (new version in /tform).  King/McClusky 080522
c         CONVERT_ANTPCV:  Add antenna name to 'not found' fatal message.  King 080929
c         EMNTAB, TESTEPH, Makefile:  Remove subroutine const from emntab (not used); remove
c            testeph from Makefile (not yet correct).  King 081002
c         CRX2RNX, RNX2CRX: Updated to version 4.0.3 from Hatanaka ftp site.  King 081008
c         METUTIL: Fix units for Temp in output met file.  King 081122
c         Makefile, CONFPOL CONVEULER EULERVEL NE2AZEL: Added these programs from the euler_utiles directory. McClusky 081122
c         AUTECL: Fix logic in do-loop.  Morgan 081205 
c         Add GLIST2CMD to generate a use_site list from a glist output.  King 090102
c         EMNTAB: Fix bugs in writing PEP tables.   King 090129
c         GRW, Makefile:  Add utility for easy creating of sig_neu commands for editing.  King 090424
c         Makefile: remove grw, glist2cmd, org2vel (moved to kf/utils).  King 090603
c   9.90  CRX2RNX, RNX2CRX: Updated to version 4.0.4 from Hatanaka ftp site.  New version avoids conflict
c            of getline with gcc 4.4, traps LF in last line, and increases MAXTYPE from 30 to 50. Hatanaka/King 090707 
c         METUTIL: Trap empty z-file from no-obs case.  Murray/King 090816
c         MERGE_RINEX: Increase number of files allowed to 150; remove check on end time to avoid losing
c           the last epoch (logic may still need fixing).  Lau/King 090820
c         CHECK_SITEID, Makefile:  Add utility to compare 4-character IDs.  King 090821
c         METUTIL: Add gradients to output file.  Murray/King 090827.
c         NUTTAB, Makefile: Remove this program since it's duplicated by kf/utils/nutIAU2000.  King 091221
c         EMNTAB: Change to read the ascii export version of the PEP nbody file, rather than binary. King 091218
c         MAKE_STNFOLIST, Makefile: Add utility to create a site list for sh_upd_stnfo/mstinf from sites.defaults
c            and/or files in the [experiment]/rinex directory.  King 091230/100114  
c  9.91   ARGO2RX: Remove 'trkcod' from call to rstnfo; remove call to obsolete check_oldstnfo  King 100209/100216
c         CONVEULER: Added output lines in GLOBK PLATE - XYZ and LLM formats.... McClusky 100319    
c         POLY01R: Increase format size for printing the mean.   McClusky 100129/100324
c         CONFPOL:  Modify formats.  McClusky 090812/100324 
c         CONVERT_ANTPCV: Add command-line argument for absolute or relative and make absolute default.  King 100607
c         CHECK_SITEID: Convert to command-line input, increase precision, make source ID more flexible, add
c           tolerance for mismatch, fix bugs.  King 100819.  
c         ARGO2RX: Add rcvers, rcvrsn, antsn to rstnfo calls. King 100908
c         Makefile: Remove obsolete upd_stnfo. King 100908     
c         NETSEL, Makefile: Move netsel.f and netsel.h from kf/utils to gamit/utils. King 101027
c         INTERP_LOG_TAB: Use const_param.h from kf/includes rather than gamit/includes.  King 101027
c         Makefile: Remove break_to_rename (moved to kf/utils).  King 101112
c         DCBTAB: Fix test on reasonable year to extend beyond 2010 (to 2050). King 110304
c         PROFILE: Change the reads to free-format. Simon 110331
c         DCBTAB: Add new header format.  King 120824
c 9.92    DCBTAB: Add yet another new header format.  King  121011
c 9.93    METUTIL: Increase maxmet (times on z-file) from 3000 to 10000 and add check. King 130715
c 9.94    CHECK_SITEID: Fix typos in call to read_cart and write statements.  King 131011
c         CHANGE_EXCL, VEXCLUDE:  Fix missing field separate for $ in format statements to satisfy the IBM AIX compiler. King 131011
c         Makefile.generic: Add flags for IBM. King 131011
c 9.95    AUTECL: Added Block IIA to test for post-ellipse data edit. Herring 131222
c         DCBTAB: Fix date problem at end-of-year. Shimada/King 140108 
c 9.96    AUTECL: Change to read the binary y-file rather than arcout.DDD.  King 140123
c         CONVERT_JPLYBIAS (new), Makefile: Convert the times in a JPL yaw_bias_table.  King 140217
c         Makefile:  Remove netsel (duplicate of /kf/utils version). King 140418
c         GEOC_TO_GEOD: Change 'octtab' to 'grdtab' in report_stat message. King 140630 
c 9.97    COUNTX, PLOTK: Change variable names to match those used in other modules. King 141129 
c         MERGE_RINEX: Allow for more than 10 observables.  King 150114
c         VEL2HIST: Trap zero residual vector.  King 150226
c 9.98    ERPTAB: Add a decimal place to the output for ut1 and pole; change position of partial-row value. King 150422/150429
c 9.99    DCBTAB2 (new), Makefile:  Create Version 2 dcb.dat files.  King 150529
c         DCBTAB, DCBTAB2: Change the names of the output files to dcb.dat.newgps and dcb.dat.newgnss. King 150803
c         GEOC_TO_GEOD: Put the right information into the report_stat call. King 150910
c         VEL2HIST: Allow old or new vel-file line length.   King 151008
c         METUTIL: Trap bad radius value.  King 151024
c         FIC2NAV: Match 'reade' calling arguments for new code. King 151119
c         CONVERT_ANTPCV: Add 'gnss' to calling arguments for lib/read_antex. King 151229
c         BSX2DCBTAB, Makefile: New routine to translate DCBs from SINEX to dcb.dat.  King 151231
c         DCBTAB2: Change so that no-long-active SVs do not get their stop years changed to 2100. King 160104
c         DCBTAB2: Add ability to read new CODE header and fix bug in setting stop/start dates for incremented files. King 160121
c         DCBTAB, DCBTAB2: Set all of the file names from the command line. King 160908
c         CHECK_SITEID: Fix glist read for new format.  King 161027
c         COUNTX, PLOTK: Add observable types for call to lib/xhdred.  King 170525
c         METUTIL: Corrected 5th decimal place in value of pi. King 170603  
c 10.00   MERGE_RINEX: Fix declaration of rxobtyp (c*3 now not c*2); add trap for RINEX 3. King 170720/170802
c         METUTIL: Add 'include ../includes/dimpar.h' for 'maxatm' instead of declaring is locally. King 170901 
c         STNFO_CONT, Makefile: New program to enforce continuity in station.info entries when made from
c            sparse individual RINEX headers.  King 170921
c 10.01   COUNTX, PLOTK: Change obsolete 'skd' to 'gnss' in calls to lib/xhdred.  King 180322 
c         ARGO2RX: Remove obsolete 'skd' from call to lib/rstnfo. King 180322 
c         AUTECL: Fix check on current version (1061, not 1051). King 180518
c         AUTECL: Fix to read header correctly for 1061 format. King 180521
c 10.02   ATX2SVNAV: Fix minor format errors.  S.H. Park/King 180827 
c         METUTIL: Modified sb lininterp to handle a vector rather than a matrix (see model/lininterp.) King 180827 
c         DCBTAB2: Fixed bug in setting year with one form of CODE header.  King 180904
c 10.03   SCAN_RINEX: Now utilty to report observables in RINEX 2/3 files. King/Herring 200512
c         CONVERT_ANTPCV: Updated atxfrq declaration to be consistent with read_antex call. Herring 200512.
c         SCAN_RINEX: Fixed error in Rinex3 handling of event flags.  Removed goto logic (as always 
c           was not needed).  TAH 201113.
c 10.04   SCAN_RINEX: Remove comma before output string on line 1361 and output variable(s) on line 1479. MAF 20201201
c 10.05   DCBTAB2: Added functionality for GLONASS DCBs. MAF 20201201
c 10.06   SCAN_RINEX: Updated warning for RINEX version (allow up to 3.05). Floyd 20210415

      return
      end


