 
*     Include file for the SAVSL common with some definitions
*     needed for saving the global solution files.
*
*   max_apr_saves   - Maxium number of apriori values which
*                   - may need to be saved
      integer*4 max_apr_saves

* MOD TAH 090612: Inccreased size by max_svs_elem*max_glb_svs to account 
*     for satellites
 
      parameter ( max_apr_saves = 8*max_glb_sites + 4*max_glb_sources +
     .                            6*max_mul_pmu + 
     .                            max_svs_elem*max_glb_svs )
 
 
*   apr_codes(max_apr_saves)- List of codes for the apriori
*               - values saved in the global file
*   org_apr_codes(max_apr_saves) - Apriori codes using the
*                 original site numbers.  (need when some sites
*                 are not used but we want to save apriori value)
*   glb_codes(max_glb_parn)     - List of the codes for the
*               - Global parameters stored in the global file.
*   num_apr_codes           - Number of values in the apr_codes
*                           - array
*   num_glb_codes           - number of values in the glb_codes
*                           - array
*   ent_par_ep              - Number of parameter epochs needed
 
      integer*4 apr_codes(max_apr_saves), glb_codes(max_glb_parn),
     .    num_apr_codes, num_glb_codes, org_apr_codes(max_apr_saves),
     .    ent_par_ep
 
*   apr_codes               - List of codes for the apriori
*               - values saved in the global file
*   glb_codes               - List of the codes for the
*               - Global parameters stored in the global file.
*   num_apr_codes           - Number of values in the apr_codes
*                           - array
*   num_glb_codes           - number of values in the glb_codes
*                           - array
 
*   gecc_known(max_glb_site_wrds)  -- Bit mapped array that indicates
*     station height information is already known and should not be
*     updated.

      integer*4 gecc_known(max_glb_site_wrds)

*   mul_par_ep(max_glb_parn) -- Epochs parameter and apriori multi
*         day parameters.
*   apr_list(max_apr_saves)  -- List of apriori values

      real*8 mul_par_ep(max_glb_parn), apr_list(max_apr_saves) 

* MOD TAH 190621: Added back_type for when glsave run from glbak.
      logical back_type      ! Set true when 'BACK' type passed
 
 
      common / glsave_com / mul_par_ep, apr_list,
     .    apr_codes, glb_codes, num_apr_codes,
     .    num_glb_codes, org_apr_codes, ent_par_ep, 
     .    gecc_known, back_type
   
 
