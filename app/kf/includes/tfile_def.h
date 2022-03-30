*
*     Definition of the t-file contents
*
* HEADER INFORMATION
*
* tf_head  -- Header line from the tfile (C*80)
* tf_comt(3) -- Three comment lines (C*80) each
* tf_start   -- Julian date of start (R*8)
* tf_ic      -- Julian date of IC integeration epoch (R*8)
* tf_end     -- Julian date of end (R*8)

* tf_ut1, tf_pole(2) -- UT1-UTC and pole positions saved in tfile
*               (If the tfile is inertial, these values are zero).
* tf_ics(max_svs_elem, max_glb_svs) -- The IC's at the IC epoch
*               (Upto 17 elements, and 32 satellites)
*
* tf_delt    --  Spacing of tehe tfile entries (seconds, R*8)
* tf_nep     --  Number of epochs in the tfile entries
* tf_nics    --  Number of IC constants (The exact types are deduced
*                from the second comment record; and the number
*                if read from the head record).
* tf_numsvs  --  Number of satellites in the tfile (I*4)
* tf_numpart --  Number of partial derivatives (I*4)
* tf_prns(max_glb_svs)  -- List of satellites in Tfile (PRN numbers).
* tf_numcrd  --  Number of coordinates given (3 if position only,
*                6 if position and velocity).

* tf_snames(max_glb_svs) -- Number of system name records by satellites
*                (C*16)

* tf_head_start -- Start of header
* tf_head_end   -- End of header
* tf_dummy      -- Dummy variable to keep real*8 on 8-byte boundaries.
* tf_head_len        -- Length of tfile head record.
* tf_ephem_len  -- Length of each ephereris record (computed).
* tf_codes(max_svs_elem) -- Code numbers for the parameters types
*                  in the tfile (see svel_to_code for list).

* tf_name       -- Name of t-file being read

      real*8 tf_start,  tf_ic, tf_end , tf_ut1, tf_pole(2), 
     .       tf_ics(max_svs_elem, max_glb_svs), tf_delt
      integer*4 tf_nep, tf_nics, tf_numsvs,  tf_numpart,
     .          tf_prns(max_glb_svs), tf_numcrd
      integer*4 tf_head_start, tf_head_end, tf_dummy, tf_head_len,
     .          tf_ephem_len,  tf_codes(max_svs_elem)
      character*80  tf_head, tf_comt(3)
      character*16  tf_snames(max_glb_svs)
      character*128 tf_name

      common / tf_head_rec / tf_head_start, tf_dummy,
     .       tf_start,  tf_ic, tf_end , tf_ut1, tf_pole, 
     .       tf_ics, tf_delt,
     .       tf_nep, tf_nics, tf_numsvs,  tf_numpart, tf_prns, 
     .       tf_numcrd, 
     .       tf_head, tf_comt,
     .       tf_snames, tf_name, 
     .       tf_head_len, tf_ephem_len,  tf_codes, 
     .       tf_head_end

*---------------------------------------------------------------------- 

* EPHEMERIS RECORD (one per epoch in the tfile)

* tf_epoch  -- Julian date of this entry
* tf_pos(6,max_glb_svs) -- Position and velocity of satellite at epoch.
* tf_parts(6,max_svs_elem,max_glb_svs) -- Partial derivatives of position
*              and velocity for each IC element for each satellite.
* tf_rec    -- TFile record number for current epoch. 

      real*8 tf_epoch,  tf_pos(6,max_glb_svs), 
     .       tf_parts(6,max_svs_elem,max_glb_svs)
      integer*4 tf_rec 

      common / tf_ephem / tf_epoch, tf_pos, tf_parts, tf_rec

*--------------------------------------------------------------------------------




