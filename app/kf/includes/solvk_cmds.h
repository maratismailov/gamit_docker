 
*     Include file which contains the names of the commands used
*     by SOLVK. Commands contains the basic command set, the @_conts
*     values contain the names of the contributions which can
*     be applied.  (See also SOLVK_COMMANDS_BD.FTN, the blockdata
*     which initializes the commands)
 
*   commands(max_commands)  -  The array containing the commands
*               - for the Kalman filer.  (See kal_block_data for
*               - commands).
*   base_cont_types(max_base_conts)  - Contains the names of
*               - contributions which can be applied by applied
*               - baseline (the same contributions are applied
*               - for all baselines)
*   medium_cont_types(max_medium_conts)  - Contains the names of the
*               - atmospheric propagation delay quanities which
*               - can be applied to the theorectical delays.
*   site_cont_types(max_site_conts)  - contains the names of the
*               - contributions which can be applied by site.
*               - These include such things as cable calibration,
*               - feed rotation etc.
 
      character*8 commands(max_commands),
     .    base_cont_types(max_base_conts),
     .    medium_cont_types(max_medium_conts),
     .    site_cont_types(max_site_conts)
 
*   commands    -  The array containing the commands
*               - for the Kalman filer.  (See kal_block_data for
*               - commands).
*   base_cont_types - Contains the names of
*               - contributions which can be applied by applied
*               - baseline (the same contributions are applied
*               - for all baselines)
*   medium_cont_types  - Contains the names of the
*               - atmospheric propagation delay quanities which
*               - can be applied to the theorectical delays.
*   site_cont_types  - contains the names of the
*               - contributions which can be applied by site.
*               - These include such things as cable calibration,
*               - feed rotation etc.
 
      common / SOLVK_CMDS_BD / commands, base_cont_types,
     .    medium_cont_types, site_cont_types
 
