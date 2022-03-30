*     Include for marge_apr_vel program
*

* PARAMETERS
      integer*4 max_site   ! Maximum number of sites allowed

      parameter ( max_site = 2048 )

* Common variables
      integer*4 num_vel_site       ! Number of velocity sites (unqiue)
     .,         num_apr_site       ! Number of sites in aprori
     .,         num_apr_nch        ! Number of sites changed based on name
     .,         num_apr_pch        ! Number of sites changed due to position
     .,         options            ! Options passed to program

      real*8 apr_pos(6, max_site)  ! Apriori site coordinates 
                                   ! and velocities
     .,      ep_pos(max_site)      ! Epoch of positions (JD)
     .,      vel_neu(3,max_site)   ! NEU velocity of sites from vel file
     .,      vel_ll (2,max_site)   ! Velocity file long and lat
     .,      dist_tol              ! Closeness of site to velocity site for
                                   ! velocity to be copied (m)

      logical zero_uvel            ! Set true if vertical velocity is to zerod

      character*8 apr_name(max_site)   ! Site names of apriori file
     .,           vel_name(max_site)   ! Site names from velocity file


      character*256 in_apr_file        ! Input apriori file name
     .,             out_apr_file       ! Output name
     .,             in_vel_file        ! Name of input velocity file 


      common / mav_r8 / apr_pos, ep_pos, vel_neu, vel_ll, dist_tol
      common / mav_i4 / num_vel_site, num_apr_site, num_apr_nch, 
     .                  num_apr_pch, options, zero_uvel
      common / mav_ch / apr_name, vel_name, in_apr_file, out_apr_file, 
     .                  in_vel_file

