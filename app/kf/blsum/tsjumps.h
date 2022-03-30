*
*     Include file for tsjumps.f program
*
* PARAMETERS
* ----------

* max_ent  -- Maximum number of entries for a single time series
* max_jmp  -- Maximum number of breaks allowed in time series
* max_edt  -- Maximum number of edits allowed in time series

      integer*4 max_ent, max_jmp, max_edt, max_rn

      parameter ( max_ent = 36500 )  
      parameter ( max_jmp = 1024  )
      parameter ( max_edt = 1024  )
      parameter ( max_rn  = 1024  )

* COMMON BLOCK VARIABLES
* ----------------------

* Control parameters that user can specify
* ----------------------------------------
* scale_dres       -- Multiplier on chi**2 for checking changes in residuals
* scale_outlier    -- Multiplier on chi**2 to see if point is outlier (ie.,
*                     how close does a jump come back to its original value
* scale_sig        -- Multiplier on chi**2 for significance of jump.
* abs_dres         -- Absolute magntiude of dres
 
* abs_slope        -- Magnitude of largest slope before being set to
*                     to zero
* max_sig          -- Largest sigma of data to be checked for changes 

      real*8 scale_dres, scale_outlier, scale_sig, 
     .       abs_dres, abs_slope, max_sig


* ts_ep(max_ent)   -- Time series epochs (JD)
* ts_val(max_ent)  -- Time series values (units set by input)
* ts_sig(max_ent)  -- Time series sigmas on values


* rn_ep(3,max_ent) -- Start, Stop and refers to epochs from the
*                     org_in file.  Used to write the rename file.
*                     Order is:
*                     1 -- Commenced with
*                     2 -- Ended with
*                     3 -- Refers to epoch

* rn_jest(2,3,max_jmp)  -- Estimate and sigmas of jump estimates in NE and U

* norm_eq(max_jmp+2,max_jmp+2) -- Normal equations for estimating breaks
*           and rate.
* sol_vec(max_jmp+2)  -- Solution vector
* stats(4) -- Accumulation statistics
* dstat(4) -- Accumulation statistics for differences between adjacent
*             values
* pre_chi, pos_chi -- Prefit and postfit chi**2 per degree of freedom
* previous_chi     -- Chi**2/f from previous iteration.
* dif_chi  -- Chi**2 per degree of freedom for the difference values
*             after removing the slope from the data  

*rn_times -- Jumps file times 
*rn_dpos  -- Jumps file offsets  

      real*8 ts_ep(max_ent), ts_val(max_ent), ts_sig(max_ent), 
     .       rn_ep(3,max_ent),
     .       norm_eq(max_jmp+2,max_jmp+2), sol_vec(max_jmp+2),
     .       stats(4),  pre_chi, pos_chi, previous_chi,
     .       rn_jest(2,3,max_jmp),rn_times(2,max_rn),rn_dpos(3,max_rn)
              

* num_ent   -- Number of values in time series
* num_org   -- Number of entries for the epochs in the org file
* num_jmp   -- Number of breaks in time series  
* ts_jmp(max_jmp)  -- Epochs numbers of the jumps in the time series. 

* rn_jmp(max_jmp)  -- Epoch counter in the rn_ep array for the jump
* num_rn           -- Number of renames needed.
* rn_off    -- Offset number for renaming sites (default 00) 
* num_renames -- Number of jumps file rename entries

      integer*4 num_ent, num_org, num_jmp, ts_jmp(max_jmp),
     .          rn_jmp(max_jmp), num_rn, rn_off, num_renames 

* out_vals   -- Logical set true if we are writing a values files
* out_rn     -- Logical set true if we are writing a rename file 
* jmp_fl     -- Logical set true if we are reading a jumps file 
* in_org     -- Logical set true if we are reading org file
* app_rn     -- Logical set true if we are appending onto an exisiting
*               rename file (default)
* out_new    -- Output values only breaks are found
* debug      -- Set true for debug output
* rename_site  -- Logical set true to rename sites
* done_jmp -- Logical set true is jump read from rn_* array 

      logical out_vals, out_rn, in_org, app_rn, out_new, debug,
     .        rename_site, jmp_fl, done_jmp(3,max_rn) 

* values_in  -- Name of input values file
* values_out -- Name of output values file
* org_in     -- Name of input org/prt file from which values file was
*               generated (needed to make rename file)
* rename_file -- Output rename file.
* jumps_file  -- Input jumps file name.
* header(3)   -- Three header lines from each time series group 
*rn_types -- Jumps file types (NEU or XYZ)
*rn_codes  -- Jumps file rename site codes

      character*256 values_in, values_out, org_in, rename_file, 
     .              jumps_file
      character*128  header(3)

* site_name   -- Name of site
* comp_name   -- Name of component

      character*8 site_name, rn_codes(2,max_rn)
      character*4 comp_name
      character*3 rn_types(max_rn)

      common / tsj_com / scale_dres, scale_outlier, scale_sig, 
     .       abs_dres,  abs_slope, max_sig,
     .       ts_ep, ts_val, ts_sig, 
     .       rn_ep,  norm_eq, sol_vec,
     .       stats,  pre_chi, pos_chi, previous_chi, rn_jest, 
     .       rn_times, rn_dpos,
     .       num_ent, num_org, num_jmp, ts_jmp, rn_jmp, num_rn,
     .       rn_off, num_renames,
     .       out_vals, out_rn, in_org, app_rn, out_new, debug,
     .       rename_site, jmp_fl, done_jmp, 
     .       values_in, values_out, org_in, rename_file, jumps_file,
     .       header, site_name, rn_codes, comp_name, rn_types
 



