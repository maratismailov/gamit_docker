      subroutine museum(version)
c
c     Welcome to visit the museum of SOLVEM.
c     Museum stores all history of the evolution of SOLVEM.

      character*16 version

      version = '1.35 beta May 30 2001'

c
c------------------------------------------------------------------------------------
c   version    |  time   | designer   |       description 
c------------------------------------------------------------------------------------
c   0.00       1992.5     D.N.Dong     finish the prototype of SOLVEM 
c   0.10 K     1992.8     K.Fiegl    1. repair english spelling and expression
c                                    2. modify i/o format
c                                    3. modify hwrite to create digestable h-file
c   0.11 TEST  09/12/92   Dong       1. restore the auxiliary parameter approach in
c                                       hirizontal angle observable (FILLN,RESIDL,READRV,LSWGHT)
c                                    2. fix priori coordinate file format and frame:
c                                       always 8-c site name + 12-c full name
c                                       always geocentric frame (GAMIT L-file format)
c                                    3. accept 3-D baseline vector variable from l-file 
c   0.12 TEST  09/30/92   Dong       1. accept full covariance matrix input from h-file
c                                    2. accept sequential input file list(FRMNOR_SEQE)
c                                    3. use FONDA-style apriori file only
c   0.13 TEST  10/08/92   Dong       1. accept fonda-output h-file data 
c                         Feigl      2. format modification
c   0.14 TEST  10/20/92   Dong       1. rewrite earthquake correction subroutines
c                                       (CALC_RES,ADD_QUAKE_C,CHK_QUAKE_LIST,GET_QUAKE_LIST)
c   0.15 TEST  11/07/92   Dong       1. chnadge net-file format with globk style
c                                    2. ouput modified network file with uncorrelated coordinates.
c                                    3. set up by-pass (if not output h-file, the running
c                                       time will be at least 3 times fast)
c                                    4. fix a bug for the sites undergoing multiple earthquakes
c   0.16 TEST  11/21/92   Dong       1. add network assign constraints
c                                       (EXERT_LINK,PUT_ASSIGN)
c   0.17 TEST  01/17/93   Dong       1. rewrite inner coordinate solution package to save space
c                                       (add INNER_SLN, remove INERCO, modify PUT_ASSIGN,SOLVEM)
c                                    2. rewrite outer coordinate solution package to save space
c                                       (add OUTER_SLN, remove OUTECO, modify PUT_ASSIGN,SOLVEM)
c									 3. allow different data types in sequential input mode
c                                       (add READ_OBS_BLOCK, modify FRMNOR_SEQE)
c                                    4. replace huge array aprm by a small aprm and two index arrays
c                                       (add REORDR_INDX,LAXB_INDX,LATWA_INDX,FILL_APR_MTX,PUT_APR_WGHT, 
c                                        remove PUT_APR_COV, rewrite FILLN_FULL, modify INITN,SOLVEM).
c                                       Now for 400 sites, the size of SOLVEM is only 30 Mb.
c                                    5. remove MODLCO temporary (need to be rewriten).
c                                    6. calculate cumulative chi squares
c                                       (FILLN,FILLN_BLOCK,FILLN_FULL,FRMNOR,FRMNOR_SEQE)
c                                    7. add Cartesian velocity observable
c                                       (FILLN_FULL)
c                                    8. modify output format (RESIDL,SOLVEM)
c 1.0 beta test  03/01/93  Dong      1. copy to stdrel.
c 1.01 beta test 03/07/93  Dong      1. add driving file name output(SOLVEM,WTHEAD,HWRITE)
c                03/11/93  Dong      1. correct fact in getpri.f (v unit in netfile is m/y)
c 1.02 beta test 03/12/93  Mark,Dong 1. remove unit bug in filln_full.f
c                03/20/93  Dong      1. correct episodic event index error (FILLN_FULL)
c                                    2. output chi2 and dchi2 to output file (READRV,PUT_LINK,PUT_ASSIGN)
c                                    3. modify event file format (FRMNOR,FRMNOR_SEQE,READ_OBS_BLOCK)
c                                    4. modify counting (REMEDY)
c 1.03 beta test 03/27/93  Dong      1. restore outer coordinate solution for subset of sites
c                                       (READRV,PUT_ASSIGN,OUTER_SLN,GET_IOM_LIST)
c                                    2. restore inner coordinate solution for subset of sites
c                                       (SOLVEM,PUT_ASSIGN,INNER_SLN,GET_IOM_LIST,READRV)
c                                    3. change format for char_4_8 (GETSIT)
c 1.04 beta test 03/30/93  Dong      1. calculate residuals for data with full covariance matrix
c                                       (SOLVEM,READRV,RED_LIST,CALC_OMC_FULL,GET_IOM_LIST,GET_RES_SGL,
c                                        GET_RES_FULL,PUT_ASSIGN)
c 1.05 beta test 04/07/93  Dong      1. calculate strain rate by strain_list file
c                                       (READRV,CALSTR,DELAUNAY,STRAIN)
c                                    2. allow to assign apriori uncertainties to episodic parameters
c                                       (GET_QUAKE_LIST,READRV,SOLVEM)
c                                    3. add: link episodic parameter (CHK_EPI_LIST,READRV,PUT_LINK)
c                                    4. full site name output to mapping file (OUTPUT)
c 1.06 beta test 04/28/93  Dong      1. restore the it=24 unit=mm/year (FILLN_BLOCK)
c                                    2. modify output format (PUT_LINK)
c                                    3. allow soft link constraint (PUT_LINK,PUT_ASSIGN,EXERT_LINK)
c                                    4. remove unit comfusion (FILLN_BLOCK,CALC_OMC_BLOCK,GETOMC,RESIDL,
c                                         GET_RES_SGL)
c                                    5. restore model coordinate solution (MODEL_SLN,READRV,PUT_ASSIGN,
c                                       GET_IOM_LIST)
c                05/01/93  kurt      6. modify output format(CALSTR)
c 1.07 beta test 05/03/93  Dong      1. remove singularity of ib = 0 (RISIDL)
c                                    2. fix cria unit comfusion (FILLN,GETOMC)
c                                    3. add it=28 (FILLN,FILLN_BLOCK,CALC_OMC_BLOCK,GETOMC)
c                                    4. remove unused array (CALSTR,FRMNOR,FRMNOR_SEQE,CPSOLN,OUTPUT,GEOVCV,
c                                       RESIDL,INNER_SLN,MODEL_SLN,INITL,GEOCNT,READ_OBS_BLOCK)
c                                    5. recompile whole package.
c 1.08 beta test 05/13/93  Dong      1. add episodic parameter calculation to block observable
c                                       (CALC_OMC_BLOCK,FILLN_BLOCK,GET_GK_ADJ)
c                                    2. link c,v in geodetic frame (PUT_LINK)
c                                    3. modify array dimension (CALSTR,DELAUNAY)
c                06/27/93  Dong      4. remove index bug in filln_block for episodic parameters
c 1.09 beta test 07/05/93  Dong      1. modify mapping file format (OUTPUT)
c 1.10 beta test 07/11/93  Dong      1. put refraction correction file explicitly (GET_REFRACTN,SOLVEM,READRV)
c                                    2. remove original refraction quasi-observation processing (CALC_OMC_BLOCK,
c                                       CALC_OMC_FULL,CALC_RES,FILLN,FILLN_BLOCK,FILLN_FULL,FRMNOR,FRMNOR_SEQE,
c                                       GET_RES_FULL,GET_RES_SGL,GETOMC,OBSRED,PRECHK,PRECHK_FULL,READ_OBS_BLOCK,
c                                       READ_OBS_FULL,RES_LIST,RESIDL,SOLVEM)
c 1.11 beta test 07/13/93  Dong      1. put frame alignment parameters in it = 31 data (FILLN_FULL,CALC_OMC_FULL,
c                                       LSWGHT,RESIDL,GET_RES_FULL,GET_RESI_SGL,SOLVEM,FRMNOR_SEQE,CPSOLN,PRECHK_FULL)
c 1.12 beta test 11/08/93  Kurt      1. add gamma rate file (READRV,CALSTR)
c 1.13 beta test 04/04/94  Dong      1. fix index bug for o-file type 25 data (FILLN_BLOCK)
c                                    2. calculate residulas handling episodic parameters for type 25 observable
c                                       but the others are still need to be modified (CALC_OMC_BLOCK).
c 1.14 beta test 06/14/94  Dong      1. fix index bug for strain rate calculation (CALSTR) 
c                06/22/94  King      1. Change Makefile to execute ranlib only once.  
c 1.15 beta test 12/15/94  Dong      1. fix bug for in sigma of gamma (STRAIN)
c 1.16 beta test 06/05/95  Dong      1. remove pmode, fcode option (READRV)  (temporary version)
c 1.17 beta test 07/03/95  Ferhat    1. add a sketch file and gmt script in driving file (READRV, GMT_FILE, GMT_SKETCH)
c                                    2. Modify CALSTR to have correct triangle name with 8-char 
c                                    3. Modify RESIDL to check problem with horizontal angle equals to zero
c 1.18 beta test 02/11/97  Ferhat    1. change calcstr.f pb woth not >o definite
c                           
c 1.19 beta test 11/29/99  MKing     1. Modify output.f to print southern hemisphere output correctly
c                                    2. Remove format mod to display std dev correctly in output file
c 1.20 beta test 08/12/99  MKing     1. Add '<horizon. dir.>' as option to wthead.f
c                                    2. Introduce maxfix to define number of fixed sites allowed (solvem.fti)
c                                    3. Correctly write entire fixed site list to output file
c                                    4. Move check for from-site=to-site from RECOD to PRECHK
c                                    5. Output site name rather than number in check
c                                    6. Add headings for solution screen dump 
c                                    7. Output site name for chi^2 update
c                                    8. Clean up GMT_FILE & GMT_SKETCH screen output
c                                        Current status of the GMT subroutines is that they do not produce legal GMT commands
c                                    9. Sequential obs files not read yet - clean up idatf assign in SOLVEM
c                                   10. Change reading/writing of the coordinate file to accept velocity in m/yr
c                                        READRV, OUTPUT, GETPRI are changed. Internal workings remain the same!
c                                   11. Change the error limit to output velocity in imnetf (OUTPUT) to be compatible
c                                        with the units (m/yr or mm/yr)
c                                   12. Print site name when identifying an ill site 
c                                   13. More output clean-up (RESIDL)
c                                   14. Rename OMC to adjust and prefit to O-C in residual output file
c                                   15. Print m/yr in mapping file when flag is set
c                                   16. Change format for the velocity components of the coordinate file OUTPUT, GETPRI
c                                   17. Read in and write out the full names of the sites
c                                   18. O-C test against user-specified criteria incorrectly implemented when O-C is in the
c                                        LH quadrant. Ie, especially when it is close the 360 degrees
c                                   19. Correct apriori-coordinate computation at standard epoch, due to incorrect units
c                                        Used in the multiplication for the velocity
c                                   20. Fix output in events file to show site name and more characters
c                                        Also unfix some of the changes to the residual file
c                                   21. Mapping file now outputs adjusted coordinates, rather than input 
c 
c 1.20 beta test 08/12/99  MKing     1. changes to output.f in solvem to correct error limits
c                                    2. Correct computation of residuals for angles and azimuths
c                                        Previously they were incorrectly computed at the reference epoch
c                                    3. Add uncorrelated version of output coordinate (FONDA format) file
c                                    4. Add residual/error term in the residual file
c                          VKotzev   5. Numerous small fixes to allow compile on Linux
c                          MKing     6. Reinstate previous (incorrect) removal of in_list in solvem.f - this is used in the hfile import
c                                    7. Fix output when coords are in Western hemisphere
c                                    8. Sequential obs files read when reading hfiles - reinstate one line in solvem.f 
c 1.22 beta test                     1. add conversion to uppercase for apr constrains in readrv.f
c                                    2. add v/s term to get_res_full.f to match residl.f
c                                    3. Fix scaling of results in get_res_full.f
c 1.33           17/05/01  Mking     1. change read_obs_block_full  also get_res_full to accept the right itp=34 vcv 
c 1.35           30/05/01  MKing     1. Tidy up the output when adding constraints. Fix up manual
      return
      end
