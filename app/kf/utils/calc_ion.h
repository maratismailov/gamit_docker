*@ calc_ion.h
************************************************************************
*                                                                      *
*     J.L. Davis                   2:12 PM  THU., 16  APR., 1987       *
*                                                                      *
*     EMA/VMA common block for CALC_ION.  As the KalObs files are      *
*     read, the pertinent information is stored in the variables       *
*     in this common.                                                  *
*                                                                      *
************************************************************************
 
*       max_obs                 - The maximum number of
*                               -   observations for any single
*                               -   KalObs file
 
      integer*4 max_obs
 
      parameter (max_obs = 50000)
 
***** Header information for "X-"  and "S-" files
 
*       calc_ion_cont(3)        - Flags for calculating
*                               -   group-phase-rate ionosphere
 
      logical calc_ion_cont(3)
 
*       sx_num_obs(2)            - Number of observations
 
      integer*4 sx_num_obs(2)
 
*       sx_site_names(max_sites,2)     - Site names
*   ,   sx_source_names(max_sources,2) - Source names
 
      character*8 sx_site_names(max_sites,2),
     .    sx_source_names(max_sources,2)
 
*       sx_KalObs(2)                  - KalObs file descriptor
 
      character*64 sx_KalObs(2)
 
 
      common / sx_header / calc_ion_cont, sx_num_obs, sx_site_names,
     .    sx_source_names, sx_KalObs
 
***** Record-by-record information from KalObs
 
*       rec_file_1(max_obs)          - rec_file_1(i) is the record
*                                    -   number in the first file
*                                    -   corresponding to the Ith
*                                    -   record in the second file
*   ,   sx_data_flag(max_obs,2)      - Data flag
*   ,   sx_FRNGE_code(max_obs,2)     - FRNGE quality codes
*   ,   sx_KalObs_rec(max_obs,2)     - KalObs record number
*   ,   sx_ion_flag(max_obs,2)       - Ion flag
*   ,   sx_site(2,max_obs,2)         - Site numbers
*   ,   sx_source(max_obs,2)         - Source number
 
      integer*4 rec_file_1(max_obs), sx_data_flag(max_obs,2),
     .    sx_FRNGE_code(max_obs,2), sx_KalObs_rec(max_obs,2),
     .    sx_ion_flag(max_obs,2), sx_site(2,max_obs,2),
     .    sx_source(max_obs,2)
 
*       sx_ion_sigma(3,max_obs)      - Ion cont. sigmas FREQ 1
*   ,   sx_sigma(3,max_obs,2)        - Observable sigmas
 
      real*4 sx_ion_sigma(3,max_obs), sx_sigma(3,max_obs,2)
 
*       sx_epoch(max_obs,2)          - Julian date
*   ,   sx_eff_frq(3,max_obs,2)      - Group-phase-rate effective
*                                    -  frequencies (MHz)
*   ,   sx_ion_cont(3,max_obs)       - Contribution FREQ #1
*   ,   sx_omc(3,max_obs,2)          - O - C
 
      real*8 sx_epoch(max_obs,2), sx_eff_frq(3,max_obs,2),
     .    sx_ion_cont(3,max_obs), sx_omc(3,max_obs,2)
 
      common / SX_records / rec_file_1, sx_data_flag, sx_FRNGE_code,
     .    sx_KalObs_rec, sx_ion_flag, sx_site, sx_source, sx_sigma,
     .    sx_ion_sigma, sx_omc, sx_epoch, sx_eff_frq, sx_ion_cont
 
