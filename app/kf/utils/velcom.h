 
*     VELCOM.H Include file containing the main common block
*     declarations
 
*     PARAMETER DECLARATIONS
 
*   max_files   - Maximum number of veloctity field files
*   max_sites   - Maximum number of sites allowed in total.
*               - (For sites to output they must appear in
*               -  the reference velocity field file.)
 
      integer*4 max_files, max_sites
 
      parameter ( max_files = 1024 )
      parameter ( max_sites = 4096 )
 
*     VARIABLE DECLARATIONS
 
*   num_files   - Number of files read (the first is the
*               - reference field)
*   tot_sites   - Total number of sites in reference field
*   num_sites(max_files)    - Number of sites in each of the fields
*   site(max_sites,max_files)   - Site numbers for each of the
*               - in the field being read.  (Only saved if
*               - match in the reference field)
 
      integer*4 num_files, tot_sites, num_sites(max_files),
     .    site(max_sites,max_files)
 
*   lat(max_sites), long(max_sites) - Latitude and longitude
*               - of the sites.
*   Nvel    - North velocity (mm/yr) for each
*               - site in each file
*   Evel    - East velocity (mm/yr)
*   Uvel    - Up velocity (mm/yr)
*   Nsig - North velocity sigma (mm/yr)
*   Esig - East  velocity sigma (mm/yr)
*   NErho   - Correlation between N and E
*                           - velocity
*   Usig - Up    velocity sigma (mm/yr)
 
      real*4 lat(max_sites), long(max_sites), Nvel(max_sites,max_files), 
     .    Evel(max_sites,max_files), Uvel(max_sites,max_files), 
     .    Nsig(max_sites,max_files), Esig(max_sites,max_files), 
     .    NErho(max_sites,max_files), Usig(max_sites,max_files)
 
*   sum_stat(8,max_sites)       - Summations statistics for
*                           - average field:
*                           - 1 - number of values
*                           - 2 - Sum dN
*                           - 3 - Sum dN**2
*                           - 4 - Sum dE
*                           - 5 - Sum dE**2
*                           - 6 - Sum dN*dE (for correlation)
*                           - 7 - Sum dU
*                           - 8 - Sum dU**2
 
      real*8 sum_stat(8,max_sites)
 
*   site_names(max_sites)       - Names of the sites
 
      character*8 site_names(max_sites)
 
*   file_names(max_files)       - Names of the input files with
*                           - first being the reference field.
*   average_file                - Name of file for average field.
*   extent                  - Name of the extent for the
*                           - difference files
*   header(max_files)           - First line for each velocity
*                           - field file.
 
      character*256 file_names(max_files), average_file, extent,
     .    header(max_files)
 
*----------------------------------------------------------------------
*     COMMON DECLARATION
 
      common / velcom_com / sum_stat, num_files, tot_sites, num_sites,
     .    site, lat, long, Nvel, Evel, Uvel, Nsig, Esig, NErho, Usig,
     .    site_names, file_names, average_file, extent, header
 
*-----------------------------------------------------------------------
 
 
