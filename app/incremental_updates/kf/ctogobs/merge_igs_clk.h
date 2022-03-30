*     Include file for the merge_igs_clk program
*
* PARAMETERS
      integer*4 max_ep   ! Maximum number of epochs
      integer*4 max_net  ! Maximum number of epochs
      integer*4 max_site ! Maximum number of sites
      integer*4 max_sat  ! Maximum number of satellites
*
C     parameter ( max_ep    = 5760 )
C     parameter ( max_net   =   16 )
C     parameter ( max_site  = 1024 )
* MOD TAH 200210: Increased from 32 to 35
* MOD TAH 200828: Increased from 35 to 100 to allow for full GNSS
      parameter ( max_sat   =  100 )
* MOD TAH 201211: Reduced size for cluster node limits
      parameter ( max_ep    = 2880 )
      parameter ( max_net   =   16 )
      parameter ( max_site  =  512 )
*
* Declarations.
*
* rcv_clk(max_site,max_net,max_ep) -- Receiver clocks
* rcv_sgc(max_site,max_net,max_ep) -- Sigmas of individual sites
* rcv_fin(max_site,max_ep)   -- Final merged receiver clocks
* rcv_sig(max_site,max_ep)   -- Sigmas for final clocks
* sat_clk(max_sat ,max_net,max_ep) -- Satellite clocks
* sat_sgc(max_sat ,max_net,max_ep) -- Sigmas for individual satellites
* sat_fin(max_sat, max_ep)   -- Final merged satellite clocks
* sat_sig(max_sat, max_ep)   -- Sigmas for final clocks
      real*8 rcv_clk(max_site,max_net,max_ep),
     .       rcv_sgc(max_site,max_net,max_ep)
      real*8 sat_clk(max_sat ,max_net,max_ep),
     .       sat_sgc(max_sat ,max_net,max_ep) 
      real*8 rcv_fin(max_site,max_ep), rcv_sig(max_site,max_ep)
      real*8 sat_fin(max_sat, max_ep), sat_sig(max_sat, max_ep)

* site_coords(3,max_site) ! XYZ of sites
      real*8 site_coords(3,max_site)

* start_ep -- MJD of start of clock values (initially set to zero)
* space_ep -- Spacing of values in days (initially zero): These
*             values are derived from the first two records of the
*             first clock file
      real*8 start_ep, space_ep

* rcv_fit(4,max_site,max_net) -- Fits to receivers (num, off, rate, rms)
* sat_fit(4,max_sat, max_net) -- Fits to satellites
      real*8 rcv_fit(4,max_site,max_net),
     .       sat_fit(4,max_sat, max_net) 

* net_mean(2,max_net) -- Mean offset and rate to applied to each 
*    net

      real*8 net_mean(2,max_net)


* rcv_flgc(max_site,max_net,max_ep) - Flags for input clocks
* rcv_flag(max_site,max_ep) - Flag for final clock
* sat_flgc(max_sat,max_net,max_ep) - Flags for input clocks satellites
* sat_flag(max_sat,max_ep) - Flag for final clock satellites

      integer*4 rcv_flgc(max_site,max_net,max_ep),
     .          rcv_flag(max_site,max_ep),
     .          sat_flgc(max_sat,max_net,max_ep),
     .          sat_flag(max_sat,max_ep)

* num_ep  ! Number of epochs
* num_net ! Number of networks
* num_site ! Number of sites
* num_ref  ! Number of reference sites
* num_sat  ! Number of satellites
      integer*4 num_ep, num_net, num_site, num_ref, num_sat

* debug  - Set true to get detailed output
      logical debug

* use_rcv(max_site), use_sat(max_sat) -- set true to use in the alignment
*     of the clocks
      logical use_rcv(max_site), use_sat(max_sat)

* fin_rcv(max_site), fin_sat(max_sat) -- Set true for final output
      logical fin_rcv(max_site), fin_sat(max_sat)

* ref_sites(max_site)  - Site numbers of the reference sites
* MOD TAH 200627: Converted prn to contain GNSS code based on conff function.
* prn_list(max_sat) - list of PRNS
      integer*4 ref_sites(max_site), prn_list(max_sat)

* best_ref          - Number of best reference clocks
      integer*4 best_ref

*  abs_mod          - Set true is -a option passed showing abs mode
      logical abs_mod
*
* clk_files(max_net) - Names of clock files
* out_file - names of output file
       character*128 clk_files(max_net), out_file

* site_codes(max_site) ! Codes for  the sites
       character*4 site_codes(max_site)
* site_names(max_site) ! Names of the sites
* autcln_ver   - version of autcln
       character*20 site_names(max_site), autcln_ver

* clk_frame_ver  ! Reference frame
       character*8 clk_frame_ver
* clk_phsmod ! Phase center model for clocks
       character*40 clk_phsmod
* clk_gamitver  ! Gamit version for clocks
       character*5 clk_gamitver

* MOD TAH 200623: Added user supplied refrence site list
       character*4 ref_codes(max_site)  ! 4-char codes input by user
       integer*4 num_usr_ref   ! Number of sites specified by user.

* MOD TAH 200723: Added clk_gnss for clock SYS (M for mixed; based on
*      satellure codes read
       character*1 clk_gnss   ! System G, R, E, C or M


*********************************************************************
* Common declaration

      common / mic_r8 / rcv_clk, rcv_sgc, sat_clk, sat_sgc
     .,                 rcv_fin, rcv_sig, sat_fin, sat_sig
     .,                 site_coords, start_ep, space_ep
     .,                 rcv_fit, sat_fit, net_mean

      common / mic_i4 / rcv_flgc, rcv_flag, sat_flgc, sat_flag
     .,                 num_ep, num_net, num_site, num_ref, num_sat
     .,                 ref_sites, prn_list, debug 
     .,                 use_rcv, use_sat, fin_rcv, fin_sat, best_ref 
     .,                 abs_mod, num_usr_ref

      common / mic_ch / clk_files, out_file, site_codes
     .,                 site_names, autcln_ver, ref_codes 
     .,                 clk_frame_ver, clk_phsmod, clk_gamitver
     .,                 clk_gnss

*********************************************************************





