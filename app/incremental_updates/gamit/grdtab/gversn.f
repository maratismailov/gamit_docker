Copyright (c) Massachusetts Institute of Technology and the University of
California, San Diego, 2006. All rights reserved.

      Subroutine GVERSN( vers )

      character*40 vers
      character*10 getmac,machin
      character*45 libver
c     function to return non-blank length of string
      integer nblen
             
c     get the machine name
      machin = getmac(1) 

c     write the machine name and utilities version into 'vers'
      write (vers,5) machin(1:nblen(machin))    
c     ** update the following each time you change the program
    5 format ('1.24 of 2020/11/11 23:45 UTC (',a,')')
          
c     get the library version
      call lversn(libver)
           
c     write the versions to the screen or log file
c      write(iun,*) ' '
c      write(iun,*) 'GRIDS v. '//vers
c      write(iun,*) 'Library v. '//libver


c   Old History to be filled in gradually.

c < 1.0   This directory created by removing from /utils all of the programs and
c         subroutines pertaining to station-list and grid files used for models
c         of ocean and atmospheric loading and met data. 
c   1.0   GRDTAB; New version merging OCTTAB for ocean loading (R. King) 
c         with GRDTAB (P. Tregoning) for atmopheric loading and mapping functions.  
c            King/Tregoning 060812-060913
c         GET_OTL_LIST: Fix mis-assignment of read lon/lat.  Shimada 060914
c         RD_ATML_GRID, GET_ATML_GRID : Fix bugs in reading/converting headers. Shimada 060915
c         GET_OTL_GRID: Fix bug in indexing for S-to-N grids.  King 060920
c         VMFTOASC: byte-swap orography for ZHD values   Tregoning 060926
c         GRDTAB; Add trap for case of list-file mismatch and no grid file.  Tregoning/King 060928
c         FIXATMLG: Add some more printout and end more gracefully on EOF.  King 061003
c   1.01  UPDATE_ATMLG, Makefile: Replaces functions of FIXATMLG and EXTEND_ATM.  King 061013
c         ATMTOASC: Generalize for all versions of grid files.  King 061013
c   1.02  RD_OTL_GRID: Change default modules from 'CSRx  CE' and 'NAO99bCE' to
c            'CSR4   E' and 'NAO99b E'.  King 061107
c         RD_OTL_GRID: Fix logic for assigning name for FES2004.  King 061108  
c         GET_MAP_GRID: Fix dimensioning of 'map'.  Shimada/King 061108 
c         GVERSN: Remove print to screen (status msg enough). King 061112
c   1.03  Makefile, REDUCE_ZHD: new routine to replace correct_zhd  Tregoning 061213
c         GET_MAP_GRID: fixed interpolation of ray-traced ZHD
c           values. Also now calls saaszd, wpress, ffun.   Tregoning 061213
c         ATMTOASC, GET_MAP_GRID, RD_MAP_LIST, VMFTOASC:  Remove extra commas in write statements.   Saba/Herring/King 061229
c         GET_ATML_GRID: Fix variable type in test of 'yr_start'.  Saba/King 061229      
c         REDUCE_ZHD: Remove unused variables. King 061229 
c         GRDTAB: Fix name of l-file unit to match new model.h.  King 070124 
c         REDUCE_ZHD: Substitute explicit code from model/atmdel instead of old routines
c            saaszd, wpress, and ffun.  King 070124
c         GET_MAP_GRID, REDUCE_ZHD: Remove unused variables.  King 070130
c   1.04  GET_OTL_GRID: Trap erroneous values at South Pole for CSR OTL grids.  King 070315
c         GRDTAB, GET_MAP_LIST, RD_MAP_LIST: Allow three versions of VMF1 file, and add 
c            a 9th variable (orthometric height) for version 1.1; fix bugs.  King 070404  
c         GRDTAB: Fix bug in checking existence of l-file. King 070404
c   1.05  RD_ATML_GRID: Add trap for old versions.  King 070411
c         UPDATE_ATMLG: Change command-line help info to reflect new program name. King 070413
c         GVERSN: Remove print unit number from 'lversn' call. King 070416    
c   1.06  GET_ATML_GRID, GET_MAP_GRID, GET_MAP_LIST: Write the results (list or grid)
c           and model into grdtab.out.  King 070420
c         GRDTAB, GVERSN: Write the program version to grdtab.out rather than the screen. King 070420
c         GET_OTL_LIST: Fix bug in printing model and coordinates to grdtab.out. King 070421
c         ATMTOASC: Enlarge string size for input atml file.  King 070510
c   1.07  GET_ATML_GRID: Change misleading comment about record counts.  King 070718
c         ATMTOASC: Fix problems with calculating in time.  King 070718
c         UPDATE_ATMLG:  Allow adding a short incremental file to an existing file to extend
c           the grid, used for regular updates from U Luxembourg.  King 070718
c         UPDATE_VMFG, Makefile: New program to update a VMF1 grid from either an incremental
c           grid or ascii files from Vienna.  King 070803
c         Makefile: Add VMFASC2BIN.  King 070806  
c         GET_ATL_LIST, GET_MET_LIST:  Add dummy statement to avoid compiler warning.  King 070907   
c         GET_OTL_GRID: Add dummy statement to avoid compiler warning of unused 'n'.  King 070907
c         GVERSN:  Remove unused 'iun' from calling argument.  King 070911
c         UPDATE_VMFG, VMFTOASC:  Fix input file size and format for open error.  King 070920
c         UPDATE_VMFG: Fix dimensioning problem for # files.  King 070920
c         UPDATE_VMFG: Change sb 'diffdays' to use comlib/yds_to_jd instead of gamit/lib/julday
c           so that days 366 and 367 work.  King 071017
c         UPDATE_ATMLG, UPDATE_VMFG, VMFASC2BIN, VMFTOASC: Remove function nydays to /lib. King 071023
c         RD_ATML_GRID:  Convert start/stop year to 4-digits.   King 071024
c    1.08 GET_OTL_GRID:  Fix bug when site longitude between 0 and the next grid point. King 080111
c         GRDTAB, Makefile, grdtab.h, new routines GEOD_TO_XYZ, XYZ_TO_GEOD, HARPOS_HEAD: Add a second
c           mode to compute OTL from command-line arguments and write an output file in 
c           HARPOS format. King 080116   
c         UPDATE_VMFG: Allow omitting an input file (start with add-file).  King 080201     
c         VMFTOASC: Change comments and a few print statements.  King 080201  
c         ATMTIDE: Fix length of character string.  King 080207
c         GRDTAB; Open eq_rename file (mainly to avoid a warning in /lib/lread.f). King 080522
c    1.09 GRDTAB; Change 'readd' to 'readdf' to accommdate name change in gamit/lib.  King 090415    
c         VMFASC2BIN: Remove extraneous common in read statement.  King 100129 
c         GET_ATL_LIST: GET_ATL_GRID, GET_ATML_LIST, GET_ATM_GRID, GET_MAP_LIST, GET_MAP_GRID,
c            GET_OTL_LIST, GET_OTL_GRID, GET_MET_LIST, GET_MET_GRID: Add ../includes/dimpar.h 
c            to get 'maxsat' for model.h.  King 100209
c    1.10 ATMFILT, ATM_BUTTER, FILTER_ATML, WRITEFILTNODE, Makefile: Add program to 
c            remove the S1, S2 tides from an ATML file using a Butterworth filter. 
c             Tregoning/King 100828
c    1.11 ATMFILT, READ_ATML_HEADERS, remove WRITEFILTNODE, Makefile: Modify to handle part-year 
c           files; reduce previous-year values used for filtering from 50 days to 13 days (as 
c           as already the case for the end of the span).  King 100831
c         AMTFILT: Fix bug in stop time for filtered file. King 101014
c         ATMFILT: Enhance the documentation message for users.  King 110223
c         ATMTOASC: Increase the string-size for the file name. Herring 110622
c         GRDTAB, GET_OTL_LIST: Fix bug due to change in format of OSO list files; fix the 
c             comments for the command line, incorrectly reflecting the old order.  King 120331
c    1.12 GRDTAB.f, MODEL.H: Modified number max entries maxatml (includes/model.h) that can be output
c         (366 days) and changed output 6-hr values to remove the tidal terms when a non-filtered  
c         input grid is used. Herring 120503
c    1.13 GET_MAP_GRID, GET_MAP_LIST: Fix bug in declaration of coefficients (integer but real in 
c          model.h). King 130116
c    1.14 ATL_ASC2BIN, ATLTOASC, GET_ATL_GRID, GRDTAB, INTERP_MONTHS, RD_ATL_GRID, Makefile, grdtab.h:  
c           Allow reading of a time-of-year dependent grid for atmospheric tidal loading.  King 130131/130204/120214
c         GET_MAP_GRID, GET_MAP_LIST: Fix (apparently innocuous) variable declaration of mapping 
c           coefficients.  King 130116
c         ATMFILT: Increase allowable filename length from 24 to 26 characters. King 120205
c         ATL_ASC2BIN: Fix format in reading name of model.  King 130212
c    1.15 GRDTAB: initialization of logical unit numbers in grdtab.h moved to block data at bottom of
c           file to satisfy the IBM AIX compiler. King 131011
c         ATL_ASC2BIN: Fix typo is call to report_stat. King 131011  
c         ATLTOASC: Add field width in format statement for logical variable to satisfy IBM AIX compiler. King 131011
c         GRDTAB: Fix mismatched variable types in call to intrinsic 'mod'. King 131011
c    1.16 UPDATE_VMFG, VMFTOASC, VMFASC2BIN:  Remove subroutine round6h and put it into gamit/lib. King 131025
c         RD_MAP_GRID, includes/grdtab.h: Save start, stop times in grdtab.h commons  for end-of-file checking King 131025
c    1.17 GET_OTL_LIST: Allow reading of a new OSO format, and change the logic to use the nearest
c           site rather than the first one within 10km.  Floyd 140209
c    1.18 GRDTAB, REP_S1S2COM (new), Makefile: Introduced new center of mass S1/S2 tide model based on 
c         analysis of  3-hr loading.u-strasbg.fr/ITRF/CM timeseries IB model.  The rep_S1S2CoM
c         subroutine replaces the CoM model used in the IERS series with the new one
c         developed here.  Use of the new model removes most the diurnal and semidiurnal
c         variations seen in the raw 6-hr time series.  To invoke the new model add a +
c         the name of the atl.grid file (i.e, use  atl.grid+ as name).  Added testing of
c         year and day for time range in map.grid files.  Herring 140326    
c    1.19 UPDATE_VMFG: Fix integer size mismatch in calling fix_y2k.  King 141129           
c         GRDTAB,: Change declarations to reflect change in model.h.  King 141206 
c         GET_ATML_GRID: Fix documentation for order and units of the ATML grid file. King 160311
c    1.20 TOC2DA, Makefile: New routine to convert MatLab produced 'toc' files to a GAMIT binary 
c           direct-access grid for ocean tidal loading.  Floyd/King 180217 
c         GET_OTL_GRID, RD_OTL_GRID: Read a version 2 OTL grid file created from a netCDF
c           file.  Floyd/King 180310
c    1.21 GET_MAP_LIST: Fix bug in initialzing the map_val array. S.H. Park/King 180827  
c    1.22 GRDTAB: Change 'maxnet' to 'maxsit' since multi-session no longer supported. King 190502
c    1.23 GET_OTL_LIST: Updated reading of TANG RADI lines for FES2014 format. Herring 200502
c    1.24 GET_OTL_LIST: Ensure longitude read from list file is also in positive range (0-360). Floyd 20201111


      return
      end

