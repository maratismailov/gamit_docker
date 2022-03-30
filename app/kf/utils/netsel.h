*     Include file for the netsel program
*
      integer*4  max_sites   ! Maximum of sites allowed
      integer*4  max_net_site     ! Maximum number of sites allowed
*                            ! allowed in each network
      integer*4  max_net     ! Maximum number of networks allowed

      parameter ( max_sites = 8192 )
!     parameter ( max_net_site = 80 )  ! This should be same as maxsit in gamit/includes/dimpar.h
      parameter ( max_net_site = 1024 )  ! Large increase for globk_US application to 
                                         ! split up velocity solution networks
      parameter ( max_net = 200 )  

      integer*4 nom_net_site ! Nominal number of sites per net
     .,         nom_tie_site ! Number of tie sites bewteen the networks
     .,          num_net_site ! Optimim number of sites per 
                              ! network 
     .,          num_site_ftp ! Number of sites in ftp log
     .,          num_site_vel ! number of sites in velocity file
     .,          network(4,max_sites)  ! Network to which a site
                              ! is assigned.  First entry is main
                              ! network, 2nd entry is rank in network
                              ! 3rd entry is the tie network number, 
                              ! 4th entry is any other network to which
                              ! this site is tied.
     .,          num_net      ! Actual number of networks 

      logical stinfav(max_sites) ! Set true if station.info available

      logical globk_US        ! Set true if being used for globk use_site list

      real*8 slon(max_sites),slat(max_sites) ! Long and lats of
*                             ! sites (deg)
     .,      flon(max_sites),flat(max_sites) ! Long and lats of
*                             ! sites (deg)
     .,      minlon, maxlon, minlat,maxlat ! Min and max long and
                              ! latitudes
     .,      maxuse_rw        ! Maximum horizontal RW to use (mm^2/yr).

      real*4 rw_stats(max_sites)    ! RW horizontal RW stats (mm^2/yr)
                              ! Used for weighting site selections.



      character*128 ftplog  ! Name of ftp log file
     .,             velfile ! Name of velocity coordinate files
     .,             sdfile  ! Name of site defaults files
     .,             getfile ! Name of file to issue commands to 
*                           ! get new coordinates
     .,             stinffile  ! Name of station.info file
     .,             rwfile  ! Name of randow walk sh_gen_stats file.

      character*8  snam(max_sites) ! Name of sites in velfile
     .,            fcod(max_sites) ! 4-char name from ftp log

      character*2  net_code        ! two character code for metwork (default is ne
                            ! and names will be neNN when NN is number

      character*128 ffile(max_sites) ! File file name for ftp site

****  COMMON DECLARARION
      common / net_i4 / nom_net_site, nom_tie_site, num_site_ftp, 
     .       num_site_vel, num_net_site, num_net, network, stinfav,
     .       globk_US, rw_stats

      common / net_r8 / slon,slat, flon,flat
     .,      minlon, maxlon, minlat,maxlat, maxuse_rw

      common / net_ch / ftplog, velfile, sdfile, getfile,
     .       snam, fcod, stinffile, net_code, rwfile, 
     .       ffile 
