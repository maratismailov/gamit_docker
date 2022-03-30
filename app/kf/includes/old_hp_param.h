CTITLE OLD_HP_PARAM.H
*
*---------------------------------------------------------------
*     Include for the parameter values for the ol hp software.
*     These values currently have not changed in the new progam
*     but we keep them here as separate in case they do change
*
*----------------------------------------------------------------
 
*   omax_sites      - Maximum number of sites
*   omax_sources        - Maximum numbet of sources
*   owrd_site       - Number of 16 bit words needed for sites
*   owrd_source         - NUmber of 16 bit words needed for sources
 
*   omax_clk_brk        - Max number of clock breaks
*   omax_edit_types - Max number of edit types
*   omax_time_edits - MAx number of time edits
 
*   omax_obep       - Max number obervations per epoch
*   omax_clk_order  - Max order of clock polynomial (index starts
*                   - at zero).
 
*   omax_ephem_epoch    - Max number of epheremis epochs saved
*   omax_channels   - Max number of Mk3 channels
*   omax_frq_wvr        - Max number of WVR channels'
*   omax_glb_parn       - Maximumum number of global parameters
 
 
      integer*4 omax_sites, omax_sources, owrd_site , owrd_source ,
     .    omax_clk_brk, omax_edit_types, omax_time_edits, omax_obep,
     .    omax_clk_order, omax_ephem_epoch, omax_channels,
     .    omax_frq_wvr, omax_glb_parn
 
      parameter ( omax_sites   =  8 )
      parameter ( omax_sources = 32 )
      parameter ( owrd_site    = (omax_sites-1)/16 + 1)
      parameter ( owrd_source  = (omax_sources-1)/16 + 1)
 
      parameter ( omax_clk_brk     = 10 )
      parameter ( omax_edit_types  = 16 )
      parameter ( omax_time_edits  =  5 )
      parameter ( omax_obep = (omax_sites-1)*omax_sites )
      parameter ( omax_clk_order   =  2 )
      parameter ( omax_ephem_epoch =  3 )
      parameter ( omax_channels    = 14 )
      parameter ( omax_frq_wvr     =  3 )

      parameter ( omax_glb_parn    = 256 )
 
*-----------------------------------------------------------------
 
