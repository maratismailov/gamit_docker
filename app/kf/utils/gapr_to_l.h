 
*     GAPR_TO_L include file
 
*   max_gtols       - Maxiumnum number of sites
 
      integer*4 max_gtols
 
      parameter ( max_gtols = 8192 )
 
* COMMON VARIABLES
 
*   lfile_epoch - Epoch to which the Lfile will referr
*               - (deciminal years)
*   lfile_jd    - Julian date of lfile time

      real*8 lfile_epoch, lfile_jd
 
*   num_full        - Number of sites with full site names
 
      integer*4 num_full
 
*   code_names(max_gtols)   - GAMIT codes for the sites (only
*                   - first four characters output)
 
      character*8 code_names(max_gtols)
 
*   full_names(max_gtols)   - Full names of the sites (only
*                   - 12 characters are used currently)
 
      character*32 full_names(max_gtols)
 
*   globk_apr_file  - Name of the globk apriori file
*   lfile           - Name of the Lfile (output)
*   full_name_file  - Name of the file with full site
*                   - names in it.
*   
 
      character*256 globk_apr_file, lfile, full_name_file
      character*128 site_comment(max_gtols)  ! Comments at end of
                     ! apriori lines.

* Added aprf_out  
      logical aprf_out   ! Set true when .apr file written instead
                     ! of old gamit format l-file
      logical ext_out    ! Set true to update EXTENDD entries and
                     ! out put new lines for renamed sites.

* Added  variables for extended terms
      integer*4 max_locsec  ! Maximum non-secular terms for one aite
      parameter ( max_locsec = 256 )
      integer*4 num_locsec  ! number non-secular terms for one aite
      integer*4 param_locsec(2,max_locsec) ! Site and type for local
      real*8 apr_val_locsec(8,max_locsec)  ! Paramters for local non-sec
 
*------------------------------------------------------------------
* COMMON DECLARATION
 
*   lfile_epoch - Epoch to which the Lfile will referr
*               - (deciminal years)
*   num_full        - Number of sites with full site names
 
*   code_names  - GAMIT codes for the sites (only
*                   - first four characters output)
 
*   full_names  - Full names of the sites (only
*                   - 12 characters are used currently)
 
*   globk_apr_file  - Name of the globk apriori file
*   lfile           - Name of the Lfile (output)
*   full_name_file  - Name of the file with full site
*                   - names in it.
      common / gapr_commR8 / lfile_epoch, lfile_jd, 
     .    apr_val_locsec
      common / gapr_commI4 / num_full, num_locsec,
     .     param_locsec, aprf_out, ext_out 
      common / gapr_commCH / code_names,
     .    full_names, globk_apr_file, lfile, full_name_file,
     .    site_comment
 
*--------------------------------------------------------------






