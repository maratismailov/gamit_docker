c
c     Common block declarations for PLTSL. This common block is used
c     in conjunction with the PLOT common block.  This common contains
c     variables needed to decode the BCKFILE from the Kalman filter.
c
c     RESTRICTION:  This common block can not be currently used
c     with the OBS_@.FTNI commons used to define the KalObs files.
c     (Same variable names are used in both commons)
c
c                                  9:34 AM  TUE., 28  APR., 1987
 
*   data_type   - Indicates the type of data used in this solution.
*               - Bit Mapped (see OBS_VALUES.FTNI for details.)
*   ibepoch     - word count for start of epoch in BAKFILE
*   imar_val    - word count for start of markov values in BAKFILE
*   ipost_res   - word count for start of post fit residuals in BAKFILE
*   ibelev      - word count for start of elevation angles in BAKFILE
*   ibak_recs   - word count for start of record numbers from KalObs
*               - file
*   ibak_sou    - word count for start of source number in BAKFILE
*   ibak_unw    - word count for start of unwieght flag in BAKFILE
 
*   mar_elements(max_mar)   - the codes for the Markiv elements in
*               - the file.  These are generated from PLOT_TYPE and the
*               - site or soucre number in the left-hand 8 bits of
*               - the word. (Left_hand 8 bits is zero for site/source
*               - Independent markov variables)
*   num_mar     - Number of markov elements in current file.
*   num_site    - Number of sites in current file
*   num_souc    - Number of sources in current file
*   num_types   - Number of plot types available
*   run_time(6) - The runtime of the Kalman filter soultion
 
      integer*4 data_type, ibepoch, imar_val, ipost_res, ibelev,
     .    ibak_recs, ibak_sou, ibak_unw, mar_elements(max_mar),
     .    num_mar, num_site, num_souc, num_types, run_time(6)
 
*   conv_mar(max_mar_types) - The conversion factors from the internal
*               - units of the Kalman filter to the output units
*               - to be used on the plots.
 
      real*4 conv_mar(max_mar_types)
 
*   plot_types(max_types)   - The names of the various plot types which
*               - can be used (See PLTSL_BLOCKDATA.FTN)
*   site_names(max_site)    - Names of the sites for this file
*   souc_names(max_souc)    - Names of the sources for this file
 
 
      character*8 plot_types(max_types), site_names(max_site),
     .    souc_names(max_souc)
 
*   unit_label(max_types)   - The labels for each of the different
*               - units to be output
      character*6 unit_label(max_types)
 
*   data_file   - KalObs file name
 
 
      character*128 data_file
 
***********************************************************************

      common / pltsl_com / data_type, ibepoch, imar_val, ipost_res,
     .    ibelev, ibak_recs, ibak_sou, ibak_unw, mar_elements, num_mar,
     .    num_site, num_souc, num_types, run_time, 
     .    conv_mar, plot_types, site_names, souc_names, unit_label,
     .    data_file
 
*****************************************************************************
 
 
