*
*     Include file for hfupd which makes hfile antenna information
*     and position estimate consistent with the standard igs.snx file. 
*
* GENERAL INFORMATION
* -------------------
* jd_zero  -- Value of JD when 0/0/0 is given as date

      real*8 jd_zero
      parameter ( jd_zero =  2415019.50d0 )

*
*-------------------------------------------------------------------------- 
* SINEX INFORMATION
* -----------------
* EXAMPLE OF CONTENTS OF IGS.SNX
*+SITE/ID
* albh  A 40129M003 P Victoria, Canada       236 30 45.4  48 23 23.3    32.0
*+SITE/RECEIVER
* albh  A ---- P 92:125:00000 92:128:85800 ROGUE SNR-8C         312   Meenix     
* albh  A ---- P 92:128:85800 92:351:72300 ROGUE SNR-8C         312   Meenix Upgr
*+SITE/ANTENNA
* albh  A ---- P 92:125:00000 94:104:74100 DORNE MARGOLIN B     91119
*+SITE/GPS_PHASE_CENTER
* 4000ST L1/L2 GEOD    -----  0.078  0.000 -0.003  0.074 -0.003 -0.001
* 4000ST L1/L2 GEOD    00080  0.078  0.000 -0.003  0.074 -0.003 -0.001
*+SITE/ECCENTRICITY
* albh  A ---- P 92:125:00000 92:146:00000 UNE   0.1290   0.0000   0.0000
* albh  A ---- P 92:146:00000 94:104:74100 UNE   0.1260   0.0000   0.0000
*
* Variables for saving site information

* max_snx_recs  --  Number of snx records to save (multiple records per site)
* max_snx_sites -- Max number of sites that can be handled
* max_snx_ent_per_site -- Max number of entries per sites (Separate entries for 
*                  recievers, antenna and eccentricity
* max_hfu_edits -- Maximum number of edits allowed
* max_hfu_edwrd -- Words need to set bit 1 to edit, bit 0 to restore
* max_snx_parn          -- Largest number of aprioris expected 
* max_rcvant    -- Maximum number of receiver or antenna names
 
      integer*4 max_snx_parn

      integer*4 max_snx_recs, max_snx_sites, max_snx_ent_per_site,
     .          max_hfu_edits, max_hfu_edwrd, max_rcvant

      parameter ( max_snx_recs  = 20000 )
      parameter ( max_snx_sites =  8192 )
      parameter ( max_snx_ent_per_site = 128 )
      parameter ( max_hfu_edits =  8192 )
      parameter ( max_hfu_edwrd = (max_hfu_edits-1)/32+1)
      parameter ( max_snx_parn  = 16384 )
      parameter ( max_rcvant    =  1024 )

* Site/ID information
* -------------------
* sx_site_code(max_snx_sites) -- 8-character code for each site
* sx_site_dome(max_snx_sites) -- Domes number for each site  ( 9 characters)
* sx_site_long(max_snx_sites) -- Long name for each site     (22 characters)
* sx_site_gpos(3,max_snx_sites) -- Geodetic position for each site (long degs
*                                latitudes degs, and height (m))

      real*8  sx_site_gpos(3,max_snx_sites)
      character*8  sx_site_code(max_snx_sites)
      character*12 sx_site_dome(max_snx_sites)
      character*24 sx_site_long(max_snx_sites) 

* sx_num_rec(max_snx_sites) -- Number of receiver entries for site
* sx_num_ant(max_snx_sites) -- Number of antenna entries for site
* sx_num_ecc(max_snx_sites) -- Number of eccetricity entries for site
* sx_rec_rec(max_snx_ent_per_site,max_snx_sites) -- Entry numbers for 
*                              each receiver entry
* sx_rec_ant(max_snx_ent_per_site,max_snx_sites) -- Entry numbers for 
*                              each antenna entry
* sx_rec_ecc(max_snx_ent_per_site,max_snx_sites) -- Entry numbers for 
*                              each eccentricity entry
* sx_ent_rec                -- Number of entries for recievers
* sx_ent_ant                -- Number of entries for anntennas
* sx_ent_ecc                -- Number of entries for eccentricity 

      integer*4 sx_num_rec(max_snx_sites), sx_num_ant(max_snx_sites), 
     .          sx_num_ecc(max_snx_sites), 
     .          sx_rec_rec(max_snx_ent_per_site,max_snx_sites),
     .          sx_rec_ant(max_snx_ent_per_site,max_snx_sites),
     .          sx_rec_ecc(max_snx_ent_per_site,max_snx_sites),
     .          sx_ent_rec, sx_ent_ant, sx_ent_ecc

* SITE/RECIEVER INFORMATION
* -------------------------
*+SITE/RECEIVER
* albh  A ---- P 92:125:00000 92:128:85800 ROGUE SNR-8C         312   Meenix     
* albh  A ---- P 92:128:85800 92:351:72300 ROGUE SNR-8C         312   Meenix Upg  

* sxrecv_st(max_snx_recs) -- Start time for reciever entry
* sxrecv_en(max_snx_recs) -- End time for reciever entry
* sxrecv_sty(max_snx_recs) -- Type of receiver
* sxrecv_sn(max_snx_recs)  -- Reciever serial number
* sxrecv_fw(max_snx_recs)  -- Reciever firmware

      real*8 sxrecv_st(max_snx_recs), sxrecv_en(max_snx_recs) 
      character*20 sxrecv_sty(max_snx_recs)
      character*5  sxrecv_sn(max_snx_recs)
      character*11 sxrecv_fw(max_snx_recs)

* SITE/ANTENNA INFORMATION
* ------------------------
*+SITE/ANTENNA
* albh  A ---- P 92:125:00000 94:104:74100 DORNE MARGOLIN B     91119
*+SITE/GPS_PHASE_CENTER
* 4000ST L1/L2 GEOD    -----  0.078  0.000 -0.003  0.074 -0.003 -0.001
* 4000ST L1/L2 GEOD    00080  0.078  0.000 -0.003  0.074 -0.003 -0.001

* These entries are generated from a combination of SITE/ANTENNA and
* SITE/GPS_PHASE_CENTER
* sxante_st(max_snx_recs) -- Start time for antenna entry
* sxante_en(max_snx_recs) -- End time for antenna entry
* sxante_sty(max_snx_recs) -- Antenna type
* sxante_sn(max_snx_recs)  -- Antenna Serial number
* sxL1a_ecc(3,max_snx_recs) -- Antenna L1 offsets (NEU: Note different
*                             order to SINEX file UNE)
* sxL2a_ecc(3,max_snx_recs) -- Antenna L2 offsets (NEU: Note different 
*                             order to SINEX file UNE) 

      real*8 sxante_st(max_snx_recs), sxante_en(max_snx_recs)
      real*4 sxL1a_ecc(3,max_snx_recs), sxL2a_ecc(3,max_snx_recs)
      character*20 sxante_sty(max_snx_recs)
      character*5  sxante_sn(max_snx_recs)

* SITE/ECCENTRICITY INFORMATION
* -----------------------------
*+SITE/ECCENTRICITY
* albh  A ---- P 92:125:00000 92:146:00000 UNE   0.1290   0.0000   0.0000
* albh  A ---- P 92:146:00000 94:104:74100 UNE   0.1260   0.0000   0.0000

* sxecce_st(max_snx_recs) -- Start time for eccentricity entry
* sxecce_en(max_snx_recs) -- End time for eccentricity entry 
* sxarp_ecc(3,max_snx_recs) -- Eccentricity values (NEU)

      real*8 sxecce_st(max_snx_recs), sxecce_en(max_snx_recs) 
      real*4 sxarp_ecc(3,max_snx_recs)
* MOD TAH 200205: Added sxantdaz for antenna azimuth
      real*4 sxantdaz(max_snx_recs) ! Antenna aligment from True N (deg). Value


 
* PROGRAM Variables
* -----------------

* num_snx_sites -- Number of sites in igs sinex file

      integer*4 num_snx_sites

* igs_snx_file -- Name of the IGS sinex file (default igs.snx)

      character*128 igs_snx_file


*-------------------------------------------------------------
* Common block declarations.
      common / sxinfo / sx_site_gpos,  sxrecv_st, sxrecv_en,
     .          sxante_st, sxante_en,  sxecce_st, sxecce_en,
     .          sxL1a_ecc, sxL2a_ecc,  sxarp_ecc, sxantdaz,
     .          sx_num_rec, sx_num_ant,  sx_num_ecc, 
     .          sx_rec_rec, sx_rec_ant,  sx_rec_ecc,
     .          sx_ent_rec, sx_ent_ant,  sx_ent_ecc,
     .          num_snx_sites, 
     .          igs_snx_file  

      common / sxchar / sx_site_code, sx_site_dome, sx_site_long, 
     .          sxrecv_sty, sxrecv_sn,  sxrecv_fw,
     .          sxante_sty, sxante_sn

****************************************************************

* VARIABLES NEEDED FOR UPDATING THE GLOBAL FILES
* ----------------------------------------------
 
* temp_r8(max_glb_parn) -- Temporary array for reading the apriori
*     and other real*8 arrays from the binary files

      real*8 temp_r8(max_snx_parn) 

      common / temp / temp_r8

***************************************************************
* VARIABLES NEEDED FOR PROGRAM CONTROL
* sx_num_upd   -- Number of stations that need updating
* dNEU_site(max_snx_sites) -- Site numbers to updated 

      integer*4 sx_num_upd, dNEU_site(max_snx_sites)

* dNEU_ecc(max_snx_sites) -- NEU position adjustments due to change
*     in antenna L1/L2 and ARP eccentricity

      real*8 dNEU_ecc(3,max_snx_sites)

* honly_opt  -- Options for just updating the header information
*               but not the solution.
*               Bit 1  -- Antenna information only
*               Bit 2  -- Mark pole tide as having been applied
* num_hfu_edits -- Number of edits listed by user
* hfu_edit_type(max_hfu_edwrd) -- Type of edit: bit 1 edit,
*                  bit 0 restore

      integer*4 honly_opt, num_hfu_edits,
     .          hfu_edit_type(max_hfu_edwrd)

* hfu_edit_start(max_hfu_edits) -- Start time for edit
* hfu_edit_end(max_hfu_edits)   -- End time for edit

      real*8  hfu_edit_start(max_hfu_edits),
     .        hfu_edit_end(max_hfu_edits)

*
* check_snx -- Set true if we are reading a sinex header file
* report_snx -- Set true is contents of sinex header file are
*              to be listed
* report_diff -- Set true if differences bewteen Hfiles and Sinex
*              header are to be reported
* upd_hf      -- Set true is hfiles are be updated
* edt_hf      -- Edit sites in hfiles
* app_ptide   -- Set true to update the pole tide (provided it
*                have not already been applied)
* honly_upd   -- Update headers only
* cwu_upd     -- Set true to Update CWU start and stop times to
*                the correct values (GIPSY bug fix: 091009)
* force_upd   -- Forces sub-range update when only one entry (case
*                when rename would split span but no data on initial
*                site)
* ps_ignore   -- Treat PS at end of name so that last two characters
*                of original name are copied).
* upd_needed  -- Set true when update needed (normally true but
*                for CWU case may be false if times are OK).
* report_rec  -- Report receiver differences (-rec option)

      logical check_snx,  report_snx, report_diff, upd_hf,
     .        edt_hf, app_ptide, honly_upd, cwu_upd, upd_needed,
     .        force_upd, ps_ignore, report_rec


* edit_file   -- Name of file with editing information (i.e., time
*                intervals over which specific sites should be removed
*                from h-file (or re-instated if previously removed)
* pmu_file    -- Optional polar motion file
 
      character*128 edit_file, pmu_file

* hfu_edit_names(max_hfu_edits) -- Names of sites and prns to
*     be edited
* hfu_edit_rname(max_hfu_edits) -- New name for a site (if blank
*     then site is edited)
* hfu_edit_code(max_hfu_edits)  -- Hfile code string that
*      must appear in hfile name. 

      character*8 hfu_edit_names(max_hfu_edits),
     .            hfu_edit_rname(max_hfu_edits)
      character*16  hfu_edit_code(max_hfu_edits)

      common / hfupd_ctrl /  dNEU_ecc, hfu_edit_start,
     .         hfu_edit_end, 

     .         sx_num_upd, dNEU_site,  honly_opt, num_hfu_edits,
     .         hfu_edit_type,

     .         check_snx,  report_snx, report_diff, 
     .         upd_hf, edt_hf, app_ptide, honly_upd, cwu_upd, 
     .         upd_needed, force_upd, ps_ignore, report_rec, 

     .         edit_file, pmu_file, hfu_edit_names,  hfu_edit_rname, 
     .         hfu_edit_code

*--------------------------------------------------------------------------

*  Declarations for saved  antenna and receiver information. These values
*  are read on the 1st call to routines  and then the saved values are used
*  in subsequent calls.

* Antenna and receiver codes.
* rcv_codes(max_rcvant)  -- Codes for receiver types
* rcv_names(max_rcvant)  -- Receiver full names
* ant_codes(max_rcvant)  -- Codes for antenna types
* ant_names(max_rcvant)  -- antenna full names

      character*6  rcv_codes(max_rcvant), ant_codes(max_rcvant)
      character*20 rcv_names(max_rcvant), ant_names(max_rcvant)

* arp_anthtcod(max_rcvant) -- Combined antname (16 chars) and 
*                 HTCODE  code (5-chars)
      character*21 arp_anthtcod(max_rcvant)

* num_rcvcode  -- Number of receiver codes in rcvant.dat
* num_antcode  -- Number of antenna codes in rcvant.dat
* num_arpcode  -- Number of antenna/ht codes found 
* num_offsave  -- Number of L1/L2 offsets saved.

      integer*4 num_rcvcode, num_antcode, num_arpcode
      integer*4 num_offSave 


* Antenna L1 and L2  phase center offsets
* L1_offSave,  L2_offSave -- L1/L2 phase center offsets
* arp_offsave -- ARP offsets: First is direct offset, second is
*     radius needed to convert sloped distance to vertical
* name_offsave(max_rcvant) -- Name of antenna and radome for
*     L1/L2 offsets.

      character*20 name_offsave(max_rcvant)

      real*8 L1_offSave(3, max_rcvant), L2_offSave(3,max_rcvant)
      real*8 arp_offsave(2,max_rcvant) 

      common / rcvant_i4 / num_rcvcode, num_antcode, num_offSave,
     .                     num_arpcode
      common / rcvant_ch / rcv_codes, ant_codes, rcv_names, ant_names,
     .                     arp_anthtcod, name_offsave

      common / rcvant_r8 / L1_offSave, L2_offSave, arp_offsave
