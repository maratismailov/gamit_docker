 
*     This include file contains the defintions for the main
*     markov part of the GLOBK_COMMON block.  Here we include
*     most of the information which will be needed by the various
*     part of the GLOBK program package.
*
*                                     09:46 PM TUE., 21 Jul., 1987
 
*   glb_markov  - First word of this part of the common
 
      integer*4 glb_markov
 
*   bak_opts    - Output options for the back solution
*   crt_opts    - Output options for user's terminal
*   grun_time(7)- Runtime of the last solution
*   guse_site(max_glb_site_wrds)    - Bit mapped array with bit
*               - set for each site to be used
*   guse_source(max_glb_sou_wrds)   - Bit mapped array with bit
*               - set for each site to be used
*   ltog_sites(max_glb_sites)   - Table to convert local site
*               - number to global site number.
*   ltog_sources(max_glb_sources)   - Table to convert local
*               - source number to global source number
*   ltog_svs(max_glb_svs)  - Table to convert local SV numbers to
*                 to global SV numbers.
*   parn_axo(2,max_glb_sites)       - Parameter numbers for
*               - axis offsets and rates.
*   parn_etd_coeff(2,max_etd_coeff,max_etd_sites) - parameter
*               - numbers for in and out of phase tidal
*               - coefficients.
*   parn_gamma  - Parameter number for gamma.
*   parn_nut_ang(8)- Parameter numbers for nutation in longitude
*               - and obliquity angles and seasonal model
*   parn_nut_coeff(2,max_nut_coeff) - Parameter numbers for
*               - in and out of phase nutation series coefficient
*               - corrections.
*   parn_ut1_coeff(2,max_ut1_coeff) - Parameter numbers for the
*                 ut1 coefficients
*   parn_xy_coeff(2,max_xy_coeff)   - Parameter numbers for the
*                 xy coefficients.
*   parn_rao    - Parameter number for RA origin.
*   parn_site(3,2,max_glb_sites)    - Parameter numbers for
*               - site XYZ and XYZ rate values.
*   parn_log(3,max_glb_sites)       - Parameter numbers for site log functions
*                 functions in NE and U.  (Actually computed in NEU frame.  Parameters
*                 appear after the site coordinates and velocities).
*   parn_source(2,2,max_glb_sources)    - Parameter numbers for
*               - the source Ra Dec and Ra Dec rate values.
*   parn_svs(max_svs_elem,max_glb_svs)     - Parameter numbers for the SVS
*   parn_tid(3,max_glb_sites)   - Tidal parameter numbers (h,l
*               - and lag)
*   parn_tran(3,2)    - Parameter numbers for three translation
*                   - parameters which handel shifts in origin
*                    (offset and rate)
*   parn_scale(2)   - Parameter numbers for scale and scale rate of
*                     change.
*   parn_wob(8)     - Parameter numbers for wobble x and y for
*               - offset, rate, seasonal offset and rate, and
*               - random walk and integrated random walk.
*   parn_ut1(6) - Parameter numbers for UT1-AT offset, rate
*               - and two different seasonal terms (usually
*               - annual and semiannual)
*   parn_eor_ut1(4) - Parameter numbers for the semi/diurn UT1 
*                 parameters
*   parn_eor_xy(4) - Parameter numbers for the semi/diurn xy
*                 parameters
*   parn_eor_etd(12,max_glb_sites) - Parameter numbers for the semi/diurn
*                 extended Earth tide parameters.

*   prt_opts    - Printer options for output
*   org_opts    - GLORG options for output.
 
*   sort_direction  - +1 mean sort values in increasing time
*               - order, -1 in decreasing time order.  The back
*               - solution will come out in the opposite direction
*               - (Passed through GLOBK runstring)

* MOD TAH 981020: New parns for multi-pmu
*   parn_mul_pmu(2,3,max_mul_pmu) -- Parameter numbers for the current
*                 set of multiday pmu values.
*   parn_rot(3,2) -- Parameter number so rotation parameters

* MOD TAH 991110: New estimated parameters for non-secular terms
*   parn_nonsec(2,max_nonsec) -- Two parameters per-nonsecular term
*   num_est_nons -- Number of estimated non-secular terms
*   param_est_nons(2,max_nonsec) -- Parameters of the estimated non-secular
*     terms
* MOD TAH 040703 -- New estimated parameters for atmospheric delays (average
*     zenith delays at sites)
*   parn_atm(max_glb_sites) -- Parameter numbers for the atmospheric delays
 
      integer*4 bak_opts, crt_opts, grun_time(7),
     .    guse_site(max_glb_site_wrds), guse_source(max_glb_sou_wrds),
     .    ltog_sites(max_glb_sites), ltog_sources(max_glb_sources),
     .    ltog_svs(max_glb_svs),
     .    parn_axo(2,max_glb_sites), 
     .    parn_etd_coeff(2,max_etd_coeff,max_etd_sites),
     .    parn_gamma, parn_nut_ang(8), parn_nut_coeff(2,max_nut_coeff),
     .    parn_ut1_coeff(2,max_ut1_coeff), 
     .    parn_xy_coeff(2,max_xy_coeff),
     .    parn_rao, parn_site(3,2,max_glb_sites), 
     .    parn_log(3,max_glb_sites),
     .    parn_source(2,2,max_glb_sources), 
     .    parn_svs(max_svs_elem,max_glb_svs),
     .    parn_tid(3,max_glb_sites),
     .    parn_tran(3,2), parn_scale(2), parn_wob(8), parn_ut1(6), 
     .    parn_eor_ut1(4), parn_eor_xy(6), 
     .    parn_eor_etd(12,max_glb_sites), prt_opts, org_opts,
     .    sort_direction, parn_mul_pmu(2,3,max_mul_pmu),
     .    parn_rot(3,2), parn_nonsec(2,max_nonsec), num_est_nons,
     .    param_est_nons(2,max_nonsec), parn_atm( max_glb_sites)
 
*   gsite_names(max_glb_sites)  - Names of the sites (in west
*               - longitude order)
*   gsource_names(max_glb_sources)  - Names of the sources (in
*               - increasing RA order)
*   gsvs_names(max_glb_svs)    - NAme of the SVS.
*   gsite_full(max_glb_sites)  - Full names of the sites (32 characters)
 
      character*8 gsite_names(max_glb_sites),
     .    gsource_names(max_glb_sources), 
     .    gsvs_names(max_glb_svs)

      character*32 gsite_full(max_glb_sites)
 
*   apr_table_file  - Name of file to output final site and source
*                   - positions.
*   glb_apr_file(max_apr_files)  - Names of files containing apriori site/source
*                   - positons, velocities and epochs
*   glb_svs_file    - File containing epoch, and SVS orbit
*                     information.
*   glb_bak_file    - Name of back solution file (also invokes
*                   - running of back solution)
*   glb_inp_file    - Name of the current global file being processed
*   glb_out_file    - Name of the combined output global file
*   glb_mar_file    - Name of the markov control file for GLOBK
*   glb_sol_file    - Name of the file containing the final
*                   - results and may be intermediate results if
*                   - a bak solution is run
*   nut_table_file  - Name of file for nutation angle corrections
*                   - (only for back solution)
*   pmu_table_file  - Name of file for PMU table (only for back
*                   - solution)
*   svs_mar_file    - Name of file containing time dependent markov
*                     statistics for satellite orbits.
 
*   list_file   - Name of the input ascii file with list of
*               - solutions to be combined
*   sort_file   - File with the list of input files sorted
*               - in ascending or descending time order. (Default
*               - GLBSRT::53, may be passed through globk
*               - runstring). This ia a type 2 file.
 
* SCM increased file name lengths from char*128 to char*256 - 090618 
      character*128 apr_table_file, glb_apr_file(max_apr_files), 
     .    glb_svs_file, glb_bak_file,
     .    glb_inp_file, glb_out_file, glb_mar_file, glb_sol_file,
     .    nut_table_file, pmu_table_file, svs_mar_file,
     .    list_file, sort_file
 
*   gdescription  - Description line for the solution being run.  It
*                 - is output as the first line of the back solution
*                 - file and the final solution.
 
      character*80 gdescription
 
*   last_glb_markov - Last word of this section of the common file
*   glb_markov_dummy(127)   - Padding at the end for when file is
*                   - read.
 
 
      integer*4 last_glb_markov, glb_markov_dummy(127)
 
*-----------------------------------------------------------
*     Common Declaration 
 
      common / glb_markov_block / glb_markov, bak_opts, crt_opts,
     .    grun_time, guse_site, guse_source, ltog_sites, ltog_sources,
     .    ltog_svs, 
     .    parn_axo, parn_etd_coeff, parn_gamma, parn_nut_ang,
     .    parn_nut_coeff, parn_ut1_coeff, parn_xy_coeff, 
     .    parn_rao, parn_site, parn_log, parn_source, 
     .    parn_svs, parn_tid, parn_tran, parn_scale,
     .    parn_wob, parn_ut1, parn_eor_ut1, parn_eor_xy,
     .    parn_eor_etd, prt_opts, org_opts, sort_direction, 
     .    parn_mul_pmu, parn_rot, parn_nonsec, num_est_nons,
     .    param_est_nons, parn_atm, 
     .    gsite_names, gsource_names, gsvs_names, gsite_full,
     .    apr_table_file, glb_apr_file, glb_svs_file, 
     .    glb_bak_file, glb_inp_file, glb_out_file, glb_mar_file,
     .    glb_sol_file, nut_table_file, pmu_table_file,
     .    svs_mar_file, list_file,
     .    sort_file, gdescription, last_glb_markov, glb_markov_dummy
 
