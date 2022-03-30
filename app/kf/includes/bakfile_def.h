 
*
*     This file contains the definitions of the variables used
*     for decoding and reading the BAKFILE from the SOLVK Ver 2.0
*     Kalman filter.  BakFile is an epoch based file i.e., all
*     results about one epoch of data are stored in the one record
*     in BakFile.
*
*                                     09:35 PM SUN., 31 May., 1987
*
 
*   max_mar_types   - Maxiumum number of types of markov elements
 
      integer*4 max_mar_types
 
      parameter ( max_mar_types = 46 )
 
*------------------------------------------------------------------------
 
*   bak_recl    - The record length for the bakfile being
*               - read.
*   bdata_type  - Data types used in this solution, read from
*               - BakFile. (See Data_type in OBS_VALUES)
*   blen_KalObs - Length of the KalObs File name read from
*               - the BakFile
*   bmar_elements(max_parn) - The code for each of the markov
*               - elements.  This is split into two parts.  The
*               - first 8 bits give the Markov type.  The second
*               - 8 bits give the site or source number (if this
*               - is appropriate)
 
*   bnstep      - Number of types of observables used at each
*               - observation (eg. delay and rate bnstep=2)
 
*   bnbl        - Number of data per epoch which could be
*               - encounterd (=bnum_basel*bnum_dtypes)
*   bnum_basel  - Number of baselines for this solution
*   bnum_dtypes - Number of data types in this solution.
*   bnum_epochs - Number of epochs of data in this solution
*   bnum_headers- Number of header records in this BakFile.
*   bnum_mar    - Number of markov elements in this solution
 
*   bnum_sites  - Number of sites read from BakFile
*   bnum_sources    - Number of sources read from BakFile
 
*   brun_time(6)    - Runtime for Solvk solution read from
*               - BakFile
*   ib_epoch    - Index in bakarray for epoch (Real*8 JD)
*   ib_marval   - Index in bakarray for start of the markov
*               - elements.  Markov values are Real*4 and contain
*               - value and sigma.
*   ib_postres  - Index in bakkarray for the start of the post-fit
*               - residuals.  These are saved in a baseline order
*               - (computed from the site numbers).  The value
*               - and its sigma in real*4 are saved. If multiple
*               - data types are used these are saved next to each
*               - other.
*   ib_elev     - Index in bakarray for the start of the elevation
*               - angles for each site (rads. in site order)
*   ib_recs     - Index in bakkarray for the start of the record
*               - numbers from the KalObs file.  (Saved in
*               - baseline order, with a zero value if
*               - observation is not present.)
*   ib_source   - Index to the source number for this observation.
*   ib_unw      - Index to the unweight flags for each observation
*               - at this epoch.  Again in baseline order with
*               - multiple data types next to each other (I*2)
 
      integer*4 bak_recl, bdata_type, blen_KalObs,
     .    bmar_elements(max_parn), bnstep, bnbl, bnum_basel,
     .    bnum_dtypes, bnum_epochs, bnum_headers, bnum_mar, bnum_sites,
     .    bnum_sources, brun_time(6), ib_epoch, ib_marval, ib_postres,
     .    ib_elev, ib_recs, ib_source, ib_unw
 
*   bsite_names(max_sites)  - Site names read from BakFile
*   bsource_names(max_sources)  - Source names read from the
*               - BakFile.
*   bmar_types(max_mar_types)   - The list of markov types
*               - which are available for manipluation.  The
*               - order in this list sets the code for the
*               - type of markov element in the data records.
 
      character*8 bsite_names(max_sites), bsource_names(max_sources),
     .    bmar_types(max_mar_types)
 
*   bKalObs_file    - Name of the KalObs file used to generate the
*               - BakFile.
 
      character*64 bKalObs_file
 
*-------------------------------------------------------------------
*     Common declaration
 
*   bak_recl    - The record length for the bakfile being
*               - read.
*   bdata_type  - Data types used in this solution, read from
*               - BakFile. (See Data_type in OBS_VALUES)
*   blen_KalObs - Length of the KalObs File name read from
*               - the BakFile
*   bmar_elements - The code for each of the markov
*               - elements.  This is split into two parts.  The
*               - first 8 bits give the Markov type.  The second
*               - 8 bits give the site or source number (if this
*               - is appropriate)
*   bnstep      - Number of types of observables used at each
*               - observation (eg. delay and rate bnstep=2)
*   bnbl        - Number of data per epoch which could be
*               - encounterd (=bnum_basel*bnum_dtypes)
*   bnum_basel  - Number of baselines for this solution
*   bnum_dtypes - Number of data types in this solution.
*   bnum_epochs - Number of epochs of data in this solution
*   bnum_headers- Number of header records in this BakFile.
*   bnum_mar    - Number of markov elements in this solution
 
*   bnum_sites  - Number of sites read from BakFile
*   bnum_sources    - Number of sources read from BakFile
 
*   brun_time   - Runtime for Solvk solution read from
*               - BakFile
*   ib_epoch    - Index in bakarray for epoch (Real*8 JD)
*   ib_marval   - Index in bakarray for start of the markov
*               - elements.  Markov values are Real*4 and contain
*               - value and sigma.
*   ib_postres  - Index in bakkarray for the start of the post-fit
*               - residuals.  These are saved in a baseline order
*               - (computed from the site numbers).  The value
*               - and its sigma in real*4 are saved. If multiple
*               - data types are used these are saved next to each
*               - other.
*   ib_elev     - Index in bakarray for the start of the elevation
*               - angles for each site (rads. in site order)
*   ib_recs     - Index in bakkarray for the start of the record
*               - numbers from the KalObs file.  (Saved in
*               - baseline order, with a zero value if
*               - observation is not present.)
*   ib_source   - Index to the source number for this observation.
*   ib_unw      - Index to the unweight flags for each observation
*               - at this epoch.  Again in baseline order with
*               - multiple data types next to each other (I*2)
 
*   bsite_names - Site names read from BakFile
*   bsource_names  - Source names read from the
*               - BakFile.
*   bmar_types  - The list of markov types
*               - which are available for manipluation.  The
*               - order in this list sets the code for the
*               - type of markov element in the data records.
 
*   bKalObs_file    - Name of the KalObs file used to generate the
*               - BakFile.
 
      common / bak_def / bak_recl, bdata_type, blen_KalObs,
     .    bmar_elements, bnstep, bnbl, bnum_basel, bnum_dtypes,
     .    bnum_epochs, bnum_headers, bnum_mar, bnum_sites,
     .    bnum_sources, brun_time, ib_epoch, ib_marval, ib_postres,
     .    ib_elev, ib_recs, ib_source, ib_unw, bsite_names,
     .    bsource_names, bmar_types, bKalObs_file
 
