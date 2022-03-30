 
*------------------------------------------------------------------
*                                                    OBS_NAMES.FTNI
*     The names part of the header records of the new Kalman filter
*     files.  Information about site and source names is contained
*     in this section.
*
*     T.Herring                   08:47 PM MON., 12 Jan., 1987
*------------------------------------------------------------------
 
*   names       - First word in the names block.  This values
*               - is equivalenced to MON_NAMES so that readf
*               - and writf can be used.
 
      integer*4 names
 
*   mon_names(max_sites)    - names of the monuments at each
*               - of the sites.  Used for mobile eccentricty
*               - monuments.
*   site_names(max_sites)   - the names of the sites
*   source_names(max_sources)   - Names of the radio sources
*               - used in this experiment
*   wvr_names(max_sites)    - Names of the WVR's used at each
*               - site.
 
      character*8 mon_names(max_sites), site_names(max_sites),
     .    source_names(max_sources), wvr_names(max_sites)
 
*   structure_root(max_sources) - roots of the files containing
*               - the maps used to correct for source structure.
 
      character*16 structure_root(max_sources)
 
*   alias_file  - Name of the file used to get the alias names
*               - for the sites and sources
 
      character*64 alias_file
 
*   spare_names - Spare space to be used if needed at some later time.
 
 
      integer*4 spare_names(128)
 
*   last_names_word - last word of the names block
*   dummy_names(127)    - 127 dummy words to ensure that we do
*               - overwrite values when the file is read.
 
      integer*4 last_names_word, dummy_names(127)
 
      equivalence ( names, mon_names )
 
*------------------------------------------------------------------
* Common declaration 
 
      common / names_block / mon_names, site_names, source_names,
     .    wvr_names, structure_root, alias_file, spare_names,
     .    last_names_word, dummy_names
 
*
