      SUBROUTINE TVERSN

      character*10 getmac,machin
      CHARACTER*40 VERS
      integer*4 nblen
C.....
      MACHIN = GETMAC(1)
      WRITE (VERS,5) MACHIN(1:nblen(machin))
    5 format ('3.14 of 2018/8/27 09:05:00 (',a,')')

      WRITE(6,10) vers
   10 FORMAT(10x,'Program TFORM ',a40,//)
C
      RETURN
C
C version 1.2     First release -- MHM 880504
C version 1.3     Correct dimension bug in ESTREF, Semi-implement
C                 Cylindrical coordinates in GETCYL (no PUTCYL yet),
C                 Implement batch mode for local coordinates in
C                 GETLOC,PUTLOC,GETCOR,FORSYS,TFORM -- MHM 880911
c ???  updates between 9/88 and 10/92 ?

c  version 1.4
c                 PUTLOC : Correct batch-mode local to cartesian.
c                    [Whose code? Date?  Implemented on MIT Apollo by king 921005]
c                 PUTLOC, GETLOC : Move radcon from parameter statement to avoid
c                     intrinsic-in-parameter-statement problem with some compilers.
c                    [Whose code? Date?  Implemented on MIT Apollo by king 921005]
c                 GETTFR: Fix pathname for tform.dat to use gu/tables/tform.data
c                 correctly for any home path.  Herring 921002
c                 GETTFR: For now, point to just tform.dat (use links). King 921005
c version 1.5
c                 Change to use gdetic.dat rather than hard-wired code for geodetic
c                     coordinates.  Changes to TFORM, FORSYS, GETCOR, GETGEO, GETLOC,
c                     PUTGEO, PUTLOC, Makefile; new routine GDATUM (from /model).
c                     (Need eventually to put rationalize /model, /lib, /solve, and /tform
c                     versions with two library routines.)
C                 GEOXYZ : Skip DATAN2 if x and y are zero.  King 921006
c                 TFORM : Clarify options message.  King 921006
c                 Remove unused variables:  GETGEO, GETTRF, PUTCYL, PUTGEO
c                     PUTSPH, TFORM.
c version 1.6     Allow use without gdetic.dat linked:  GETGEO, GETLOC, PUTGEO,
c                     PUTLOC, GDATUM, FORSYS.    King 921127
c                 Remove debug from PUTGEO.  King 921207
c version 1.7     Correct DMSDEG/DEGDMS for negative angles
c                 (Note: should be changed to read N/S and E/W)
c                 Change SITNAM to A16: GETCAR, GETCOR, GETCYL, GETGEO, GETSPH,
c                                       PUTCAR, PUTCYL, PUTGEO, PUTLOC, PUTSPH,
c                                       TFORM, ESTREF
c                 Change include statement: ESTIMT, ESTREF, FORSYS, GETCAR,
c                                           GETCOR, GETCYL, GETSPH, PUTCAR
c                                           PUTCYL, PUTSPH, GETTFR
c                 PUTSPH: Output file in l-file format
c                 Bock 921228
c version 1.71    Consolidate Bock/Feigl changes for L-file read/write:
c                 FORSYS, GETSPH, PUTSPH
c                 Bock/Feigl 930122
c                 FORSYS: Updated in MIT version from SIO 920122 (missed by
c                     King)  930210
c version 1.72    PUTSPH: Fix site-name length.  Ferhat/King 940412
c                 Makefile: Change HP compiler switches from +e +E1 to +U77.  Fang/King 940506
c                 Makefile: Shift ranlib to execute only once.   Herring/King 940624
c version 1.73    GETCAR, GETCYL, GETGEO, GETLOC, GETSPH:  Avoid branching into
c                     if/then (RISC compiler complains).  King 941201
c                 GETLOC:  Finish branch correction (HP compiler).  King 950331
c                 Declare integer variables explicitly:  DEGDMS,DMSDEG, ESTIMT,  ESTREF, GEOXYZ,
c                    GETCAR, GETCOR, GETCYL, GETGEO, GETLOC, GETSPH, PUTCAR, PUTCYL, PUTGEO
c                    PUTLOC, PUTSPH, SYMINV, TF2REF
c version 1.74    ESTIMT, ESTREF: Define output print unit.  Sanli/King  950719
c                 FORSYS, GETTRF: Status required for open with DEC.  Sanli/King 950719
c version 1.75    Makefile: Variable ranlib for compatibility with Solaris 2.  Fang/King 951208
c                 GETTFR:  Enlarge format field for translations.  Ferhat/King 960208
c                 GETGEO:  Allow L-file format for geodetic coordinates.  King 960619 
c version 1.76    Changes for new kf/gamit structure and Makefiles from unimake.  King 960718
c                      All routines:  removed trailing blanks and replaced lib/*.fti 
c                      includes by ../includes/*.h.  
c                 GETGEO: Declare 'latflag', 'lonflag' character*1.  King 960815
c                 GETCAR, GETCOR, GETCYL, GETGEO, GETTFR, GETLOC, GETSPH, ESTIMT, PUTCAR, PUTCYL
c                      PUTLOC, PUTSPH : Declare calling arguments(and others) explicitly.  KIng 960816
c                 GETGEO, GETSPH:  Add extra instructions on l-file input.  King 960904
c version 1.77    GETSPH: Trap and bypass illegal l-file lines.  King 960912
c                 FORSYS: Allow only spherical coordinates into L-file.  King 960912 
c                 ../includes/tform.h:  Increase dimensions to 2000.  King 961008
c version 1.78    G77 compiler changes:
c                 SYMINV, PUTSPH, GETLOC, ESTREF, TFORM, TF2REF, GETSPH, GEOXYZ, DEGDMS: 
c                    Initialised previously uninitialised variables, removed unused variables, 
c                    and explicitly typed all unexplicitly typed variables. McClusky 970113
c version 1.79    GETCAR, GETCYL, GETGEO, GETLOC, GETSPH: Remove explicit limits of do loops 
c                    (1000 stations) in favor of nsdim in tform.h
c                 DEGDMS, PUTSPH, PUTGEO:  Fix problem with latitudes within one south degree 
c                    and longitudes within one degree west (i.e. negative values).  Morgan 980512 
c                 PUTGEO: FIx bug in defining sign of lat and long   Tregoning 980604  
c                 GEOXYZ: Fix bug for southern hemisphere.  Fang 990102      
c                 GETLOC: Fixed bug when dealing with local coordinates.  Bock 990108
c                 PUTSPH, Makefile: Removed CARSPH as individual subroutine.  Bock 990108
c                 GETGEO: Fix output format.  MKing, RKing 990819      
c version 1.80    DMSDEG, FORSYS, GDATUM:  Make implicit none.  King 000816
c                 GETSPH: Fix sign change so works for lat,long between 0 and 1.  GeMaorong/King 001004
c                 PUTGEO: Fix extra space in lat print for case of positive long and negative lat.  King 001102
c version 1.81    GETSPH: Add 001004 fix for file-read (done origionally only for terminal input).  Bouin/King 011105  
c version 1.82    GETSPH, GETGEO:  Add warning for -0 lat or long entries.  King 020404  
c                 GETSPH: Correct rwk blunder in leaving duplicate code.  Herring/King 020503
c version 1.83    Makefile, GETDMS:  Move getdms to /lib for use in /makex.  King 020923   
c version 2.00    Makefile, GETGEO, PUTGEO, GETLOC, PUTLOC:  Replace GDATUM and local GEOXYZ by library 
c                    routines read_gdatum and geoxyz.  King 021002  
c                 Makefile:  Move dmsdeg.f to /lib.  King 021003   
c                 TFORM, CONV_XYZ2GEO, FORSYS,GETGEO, GETCAR, GETSNAME, GETSPH, GETLOC, PUTGEO, PUTCAR, PUTSPH,
c                    Makefile:  Major changes to read new-format gdetic.dat, compute NGS ITRF-to-NAD83
c                    transformations, new routine for command-line conversions, new restrictions on
c                    file formats.  King 021022    
c version 2.01    GETSPH:  Allow reading of non-standard L-file lines that may be comments.  King 021218 
c version 2.02    PUTGEO:  Fix bug in outputing site names to file when deg/min/sec.  King 030328
c version 2.03    GETGEO, GETSPH, PUTGEO:  Fix calliing arguments for new verisons of lib/dmsdeg.f
c                    and lib/degdms.f   King 030521 
c version 2.04    GETGEO, GETSPH: Fix bug in converting S lat and W lon.  King 030605   
c                 PUTGEO, PUTSPH: Fix bug in converting deg to dms.  King 030623
c version 2.06    PUTGEO, PUTSPH: Fix 1-column format error in radius.  King 030717   
c                 GETSPH: Fix bug in converting decimal input via the terminal.  King 030821
c version 2.07    GETSPH: Add numsit to warning about comment lines.  King 031209
c version 2.08    GETCOR, GETGEO, GETSPH, TFORM : Fix l-file format, fix problem reading geodetic 
c                      coordinates from a file.  King 050825  
c                 GETGEO, PUTGEO: Allow NOD83 defaults if gdetic.dat missing. Agnew/King 050826
c                 GETCOR, GETCYL: Remove unused variables.  King 051029
c version 2.09    PUTSPH: Fix bug in output format for decimal latitude.  King 060703
c version 2.10    ESTIMT, ESTREF, GETLOC, GETTRF : Remove unused statement label.  King 070910    
c                 FORSYS, GETCOR, GETCYL, GETLOC, PUTCYL, PUTLOC, TFORM: Remove unused calling argument.  King 070910 
c version 3.0     CONVERTC: New routine to do generalized conversions for GAMIT, GLOBK, and Google Earth
c                 files (lfile, apr, vel, glist, kml) from command-line input. Changes to Makefile. New
c                 routines: CONVERTC_HELP, GET_FILETYPE, READ_GLISTF, READ_LFILE, READ_VELF, WRITE_KMLF, 
c                  WRITE_LFILE, WRITE_VELF.   King 080304
c                 CONVERTC: WRITE_KML: Write correct kml header, footer and lon/lat coord order for KML output. McClusky 080521
c                 CONVERTC: Add site name to full-precision geodetic output.  King 080612
c version 3.10    CONVERTC: Add 'GEP' as 'GEO' output option with longitudes 0-360.  King 080903
c                 CONVERTC: Fix GEP option for screen output.  King 080910/080911
c                 CONVERTC: Fix bug in outputting spherical L-file.  King 080913
c                 CONVERTC_HELP: Clarify that only GEO, GEP, and XYZ allowed as screen output. King 090421
c                 CONVERTC_HELP: Fix typo GEOP --> GEP.  King 100514
c version 3.11    CONVERTC, READ_GLISTF, READ_LFILE, READ_VELF, WRITE_KMLF, WRITE_VELF: 
c                   Use const_param.h from kf/inclues instead of gamit/includes (removed) King 101027
c version 3.12    CONVERTC, GET_FILETYPE, READ_LFILE: Get the epoch for an old-style L-file from the first line
c                   and check each site line for reasonableness; allow input GEO or GEP files; reverse
c                   order (lon/lat rather than lat/lon) for GEO/GEP files.  King 110701
c                 CONVERTC: Add error messages for apr-file reads.  King 111102
c                 GET_FILETYPE, READ_VELF: Fix problems reading vel files from velrot or with longer lines.  King 130531
c version 3.13    CONVERTC, READ_GLIST, READ_LFILE, READ_VELF, WRITE_VELF, WRITE_KMLF: Change const_param.h reference
c                   from /kf to /libraries. King 151008
c                 PUT_GEO: Fix extraneous character warning. King 151008
c version 3.13    CONVERTC: Update output format for VEL files.  King 180111
c version 3.14    GETGEO: Fix format syntax. S.H. Park/King 180827

      END


