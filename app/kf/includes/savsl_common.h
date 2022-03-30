 
*     Include file for the SAVSL common with some definitions
*     needed for saving the global solution files.
*
*   max_apr_saves   - Maxium number of apriori values which
*                   - may need to be saved
      integer*4 max_apr_saves
 
      parameter ( max_apr_saves = 4*max_sites + 2*max_sources )
 
 
*   apr_codes(max_apr_saves)- List of codes for the apriori
*               - values saved in the global file
*   glb_codes(max_parn)     - List of the codes for the
*               - Global parameters stored in the global file.
*   min_parn                - Number of first parameter to be
*                           - saved
*   num_apr_codes           - Number of values in the apr_codes
*                           - array
*   num_glb_codes           - number of values in the glb_codes
*                           - array
*   pmu_code_loc(3)         - Location in the glb_codes array of
*                             of the pmu offset parameters.  Need if
*                             we are to suppress these values.
 
      integer*4 apr_codes(max_apr_saves), glb_codes(max_parn),
     .    min_parn, num_apr_codes, num_glb_codes, pmu_code_loc(3)
 
*------------------------------------------------------------------
*     The common declaration
* 
      common / glsave_com / apr_codes, glb_codes, min_parn,
     .    num_apr_codes, num_glb_codes, pmu_code_loc
 
