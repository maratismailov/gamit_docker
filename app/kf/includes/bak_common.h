*@BAK_COMMON.FTNI
***********************************************************************
*                                                                     *
*     J.L. Davis 870413                                               *
*                                                                     *
*     Common block for handling back files                            *
*                                                                     *
***********************************************************************
 
*       bak_record(bak_record_size) - Array for reading BAKFIL
*   ,   bak_data_type               - Data type used in SOLVK run
*   ,   bak_dcb(bak_dcb_size)       - DCB for reading BAKFIL
*   ,   bak_ind_epoch               - Index for the epoch
*   ,   bak_ind_Markov              - Index for markov elements
*   ,   bak_ind_res                 - Index for residuals
*   ,   bak_ind_elev                - Index for elevations
*   ,   bak_ind_recs                - Index for records from Kalfil
*   ,   bak_ind_source              - Index ofr source numbers
*   ,   bak_ind_flags               - Index for flags
*   ,   bak_mar_elements(max_parn)  - Codes for Markov elements
*   ,   bak_num_epochs              - Number of epochs in back file
*   ,   bak_num_mar                 - Number of Markov elements
*                                   -   in back file
*   ,   bak_num_sites               - Number of sites in back file
*   ,   bak_num_sources             - Number of sources in back file
*   ,   bak_recl                    - Record length for bakfil
*   ,   bak_run_time(6)             - Run time of solution (YMDHMS)
 
      integer*4 bak_record(bak_record_size), bak_data_type,
     .    bak_dcb(bak_dcb_size), bak_ind_epoch, bak_ind_Markov,
     .    bak_ind_res, bak_ind_elev, bak_ind_recs, bak_ind_source,
     .    bak_ind_flags, bak_mar_elements(max_parn), bak_num_epochs,
     .    bak_num_mar, bak_num_sites, bak_num_sources, bak_recl,
     .    bak_run_time(6)
 
*       bak_site_names(max_sites)       - Names of sites read
*                                       -   from back file
*   ,   bak_source_names(max_sources)   - Names of sources read
*                                       -   from back file
*   ,   Mar_types(max_types)            - What types of parameters
*                                       -   are in BAKFIL?
 
      character*8 bak_site_names(max_sites),
     .    bak_source_names(max_sources), Mar_types(max_types)
 
*       bakfil_name                     - BAKFIL descriptor
 
      character*64 bakfil_name
 
***********************************************************************
*                                                                     *
*     Common blocks                                                   *
*                                                                     *
***********************************************************************
 
 
      common / bak_dcb_common / bak_dcb, bak_record
 
 
 
      common / bak_header / bak_data_type, bak_ind_epoch,
     .    bak_ind_Markov, bak_ind_res, bak_ind_elev, bak_ind_recs,
     .    bak_ind_source, bak_ind_flags, bak_mar_elements,
     .    bak_num_epochs, bak_num_mar, bak_num_sites, bak_num_sources,
     .    bak_recl, bak_run_time, bak_site_names, bak_source_names,
     .    Mar_types, bakfil_name
 
