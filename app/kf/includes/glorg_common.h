 
*     This is the include common block for the GLORG program
*     It is used locally in this program
 
*   num_glorg_commands  - Number of commands for GLORG
*   max_equates         - Maximum number of equates which can be
*                         set.
*   max_each_equate     - Maxmimum number of parameters for each
*                         equate (can be multiply invoked)
*   max_force           - Maximum number of parameters which can
*                         be forced to specific values

* MOD TAH 060511: Inceased file names from 64 to 256 characters.

* MOD TAH 040329: Added wild cards to plate and introduced new arguments to stab_min

* MOD TAH 040230: Added VXTRAN VYTRAN VTRANZ to force and equate commands to force Value
*   of X,Y,Z transaltion (needs apr_trans in globk command file), RXTRAN RYTRAN RTRANZ
*   will force/equate translation rate; VSCALE RSVALE forces/equates scale and scale 
*   rate (need apr_scale in globk).
*   Introduced output option SMAR (used in prt_opt, org_opt) to print markov process
*   noise on sites.
* MOD TAH 030513: Added explicit translation estimation for plate pole 
*       estimation as the default and introduced command to turn this
*       feature off (consistent with pre-
* MOD TAH 981105: Added feature to write out and solution file
*       so that it can be saved as a binary h-file.
* MOD TAH 980612: Added ASSIGN command to assign a site to plate
*       so that when residuals are computed, the new plate pole is
*       used.

* MOD TAH 980417: Added iteration of the stabalization calculations 
*       to edit and weight the sites.

* MOD TAH 970514: Added extra sigma checking for the pos_org and
*                 rate_org commands.
*   (1) hgt_var  command appened with two arguments which set the
*                height variation ratio allowed for position and
*                rate determination.

* MOD TAH 941227: Added three new commands:
*   (1) pos_org    -- List of parameters to be used in position
*                     constraint (xtran, ytran, ztran, xrot, yrot, zrot,
*                     scale)
*   (2) rate_org   -- List of parameter to be used in rate constraint
*   (3) hgt_var    -- Height variance pos and rate.

*   (4) eq_dist    -- Command to force equates between sites within
*                     specified distance.
 
      integer*4 num_glorg_commands, max_equates, max_each_equate,
     .          max_force

* MOD TAH 950109: Added plate command to estimate plate rotation
*                 vectors.

*   max_plates    - Maxiumum number of plates allowed

      integer*4 max_plates

      parameter ( num_glorg_commands =  32  )
* MOD TAH 060518: Increased numbers allowed
      parameter ( max_equates        = 8192 )
      parameter ( max_force          = 8192 )
      parameter ( max_each_equate    =  48  )
      parameter ( max_plates         =  64  )

*   val_force(max_force) - Values (adjustment values) to which
*                 parameters will be forced.
*   var_force(max_force) - Variance of the forced adjustments
*   cnd_hgt_var(2)       - Height variance to be used in constraint
*                 for position and rates.

*   eq_dist              - distance under which the specified
*                          parameter will be equated.
*   eq_var(max_equates)  - Variance to be assigned to each equate.
*                          (When the equate command is used, this
*                           is set to zero)
*   cond_var(7,2)        - Variances to be applied to the transaltion
*                          rotation and scale pos and velocity terms.
*                          (command COND_SIG)
*   use_ratio(2)         - Ratio of height variances to allow for use
*                          in fixing the origin of the system.
*                          (1.0 works well for loose solutions; 2-3
*                          is better for constrained systems i.e, CofM
*                          not free). Values for position and height.
* MOD TAH 980417: Added variables to allow iteration of orientation constraint
*   stab_site_var(max_glb_sites,2)  -- Variance assigned to site in the stabalization
*           for position and velocity.
*   stab_site_err(max_glb_sites)  -- Horizontal position error in each site.  If
*           too large then site is removed from list
*   stab_nsig                    -- Condition to remove a site from the list if
*           its residual is too large,
*   stab_rel                     -- Gives ratio of constant to site dependent
*           weighting, ie. 0 is constant weighting; 1 is site dependent only.
*   stab_min_rms(2)  -- Minimum RMS to use for testing outliers
*   stab_min_dh(2)   -- Minimun difference between min and median height contition.
*   stab_min_dne(2)  -- Mimimum difference between min and median NE sigmas for post
*                       stablization sigmas.

      real*8 val_force(max_force), var_force(max_force),
     .       cnd_hgt_var(2), eq_dist,
     .       eq_var(max_equates), cond_var(7,2), use_ratio(2) ,
     .       stab_site_var(max_glb_sites,2), 
     .       stab_site_err(max_glb_sites), stab_nsig, stab_rel,
     .       stab_min_rms(2), stab_min_dh(2), stab_min_dne(2)
 
*   use_sites(max_glb_site_wrds)    - Bit mapped word giving
*                       - the sites to be used
*   use_pos(max_glb_site_wrds)      - Bit mapped word giving
*                         the sites in position stabalization
*   use_rat(max_glb_site_wrds)      _ As above for rates.
*   cov_sites(max_glb_site_wrds)    - Bit mapped word giving
*                       - the sites to be used in computing the
*                         site covariance matrix
*   param_equates(max_each_equate, max_equates) - Parameter numbered to be
*                         equated.
*   num_equates - Total number of equates in this run
*   num_each_equate(max_equates) - number of parameters to be equated for
*                 each of the equates.
*   num_force   - Number of parameters to be forced
*   param_force(max_force) - Parameters to be forced

*   cnd_pos_bits, cnd_rat_bits -- Bits set for position and rate constraint

*   num_plates   - Number of plates to be estimated
*   num_plate_sites  - Number of sites to be used in the plate estimation
*   plate_number(max_glb_sites)  - Plate number assigned to plate (zero if
*                      site not allocated to a plate)
* MOD TAH 980611: Assign_number introduced.
*   assign_number(max_glb_sites) - Plate number assigned to a plate but 
*        these sites are not used in the computation of the plate pole.

* MOD TAH 980417: Added variables to allow iteration of orientation constraint
*   num_stab_iter  - Number of stabaliztion iterations to perform.  With each
*     iteration, station weights and editing of sites is performed
*   dcb(16)        - DCB buffer for reading covariance matrix file.

      integer*4 use_sites(max_glb_site_wrds),
     .          use_pos(max_glb_site_wrds), use_rat(max_glb_site_wrds),
     .          cov_sites(max_glb_site_wrds),
     .          param_equates(max_each_equate, max_equates),
     .          num_equates, num_each_equate(max_equates),
     .          num_force, param_force(max_force),
     .          cnd_pos_bits, cnd_rat_bits,
     .          num_plates, num_plate_sites, 
     .          plate_number(max_glb_sites),
     .          assign_number(max_glb_sites),
     .          num_stab_iter, dcb(16)

*   equate_loc - Logical which indicates that we should equate in
*                a local NUE frame and not in XYZ
*   first_eqf  - Logical which indicates that the eqautes and forces
*                should be done before the stabaliztion of the
*                corrdinate system
*   PlateTrans - Set to true to use translation estimation when
*                when estimating plate Euler poles.

      logical equate_loc, first_eqf, PlateTrans
 
*   glorg_commands(num_glorg_commands) - The GLORG commands
*                       - (see GLORG_BD.FTN)
*   org_types(7)        - List of origin constraint types
*   plate_names(max_plates)  - Names of plates
 
      character*8 glorg_commands(num_glorg_commands),
     .            org_types(7), plate_names(max_plates) 
 
*   glorg_apr_file      - Name of the glorg apr file (same format
*                       - as apriori file for GLOBK)
*   glorg_command_file  - Name of the glorg command file
*   glorg_plate_file    - Name of the file containing the plate
*                       - velocity file.
*   newoutfile             - Name of the output file being used in this run
 
      character*256 glorg_apr_file, glorg_command_file,
     .    glorg_plate_file, newoutfile

*   param_names  - list of parameter names or numerical
*                  values
      character*12 param_names(max_glb_parn)

 
*   parm_change(max_glb_parn)   - Changes to the parmeter values as
*                       - we update aprioris
 
 
      real*8 parm_change(max_glb_parn)

*   com_orgopt  -- Option for reading the glorg command file (normally
*     the same as comopt from the globk common)
* MOD TAH 190526: Increased size from 16 to 256 to be compatible with GLOBK
*     string comopt.   Also used in glbak to resolve frame in back solution.

      character*256 com_orgopt


*---------------------------------------------------------------------------
* COMMON DECLARATION
 
      common / glorg_common / val_force, var_force, parm_change,
     .    cnd_hgt_var, use_ratio, eq_dist, eq_var, cond_var, 
     .    stab_site_var, stab_site_err, stab_nsig, stab_rel,
     .    stab_min_rms, stab_min_dh, stab_min_dne, 
     .    use_sites, use_pos, use_rat, cov_sites, 
     .    param_equates, num_equates, num_each_equate,
     .    equate_loc, first_eqf, PlateTrans, num_force, param_force, 
     .    cnd_pos_bits, cnd_rat_bits, 
     .    num_plates, num_plate_sites, plate_number, assign_number,
     .    num_stab_iter, dcb, 
     .    glorg_commands, org_types, plate_names, 
     .    glorg_apr_file, glorg_command_file, glorg_plate_file,
     .    newoutfile, param_names, com_orgopt

      

*---------------------------------------------------------------------------
 
