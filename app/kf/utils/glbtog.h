 
* DECALARATIONS FOR GLBTOG program
 
* PARAMETER STATEMENTS:
 
*   max_svs      - Maximum number of satellites in gfile
*   max_svs_elem - Maximum number of elements expected (should be
*                  the same as kalman_param.h)

      integer*4 max_svs, max_svs_elem
 
* MOD TAH 150825: Make consistent with GLOBK GNSS 
      parameter ( max_svs = 100 )
      parameter ( max_svs_elem   = 23 )

* VARIABLES
 
*   orb_svs(max_svs_elem,max_svs)  - Orbital elements of each satellite
*   orb_sig(max_svs_elem,max_svs)  - Sigmas on orbital elements.  These are
*                 m, mm/s and unitless.
*   exp_jd      - Julian date with fraction of day of experient
*   svs_jd     - Julian date with fraction of day of IC epoch
*   run_jd     - Julian date with fraction of day of run time
 
      real*8 orb_svs(max_svs_elem,max_svs), 
     .       orb_sig(max_svs_elem,max_svs), 
     .       exp_jd, svs_jd, run_jd
 
*   num_svs     - Number of satellites found in current IC epoch
*   prns(max_svs)   - PRN numbers of those satellites found
*   svns(max_svs)   - SV number of satellite (from svnav.dat Verion 2.0)
*   yr_last     - Last digit of year
*   gindx       - Index counter for number of gfiles written with
*               - the same root name
*   num_svs_elem - Number of orbital elements for the solution
 
      integer*4 num_svs, prns(max_svs), yr_last, gindx, 
     .          num_svs_elem, svns(max_svs) 

* MOD TAH 180401: New GNSS satellite naming schems
      logical new_gnss  ! Set true of new GNSS outout
 
*   curr_gfile      - Name of current gfile
*   prev_gfile      - Name of the last gfile written
 
      character*16 curr_gfile, prev_gfile

      character*8 btype(max_svs)   ! Block type (BLOCK removed) from svnav.dat
 
*   code_gfile      - Code four g-file names input by used
*                   - (could be 4 or 5 characters)
*   full_code       - Full code to be used for gfile (may have
*                   - last digit of year added.
*   gtime, gframe, gprec, gsrpmod - Models used in the GAMIT
*                     runs.
*   gnut, ggrav     - More gamit models (080103)
*   geradmod        - Earth-radiation model (Added 140327) 
*   gantradmod      - Antenna-radiation mode (Added 140327) 

      character*8 code_gfile, full_code, gtime, gframe, gprec, gsrpmod
      character*8 gnut, ggrav, geradmod, gantradmod
     
*   arels(max_svs_elem) - Orbital elements in gfile

      character*4 arels(max_svs_elem)
      
*   input_file      - Name of input file
*   comment_line    - Optional comment line for the gfile passed through
*                     runstring
      character*132 input_file, comment_line

      character*1 gnss_sys(max_svs)

* MOD TAH 200611: Added GNSS slection for output.  This option also
*     effects output g-file names with the 4th character replaced with
*     selected gnss
      character*1 sel_gnss   !  User selected GNSS
 
*
* COMMON BLOCK
* MOD TAH 180609: Added gnss_sys to common (previously defined)
 
      Common / glbtog_com / orb_svs, orb_sig, exp_jd, svs_jd, 
     .    run_jd, num_svs, prns, svns, yr_last, gindx,
     .    num_svs_elem, new_gnss, 
     .    curr_gfile, prev_gfile, code_gfile,
     .    full_code, gtime, gframe, gprec, gsrpmod, gnut, ggrav,
     .    geradmod, gantradmod, arels, input_file, btype, comment_line,
     .    gnss_sys, sel_gnss
 
