 
      program velrot

      implicit none
 
*     Program for the comparison of site velocity estimated in
*     difference systems.  The runstring of the program is:
*
*     % velrot <sys1> <frame1> <sys2> <frame2> <outname> <out_frame> \
*              <fundamental sites> <height wght> <param_opt>
*
*     where <sys 1> is the name of the file containing the sites
*                 and their  velocity in system 1 (Standard vel file format)
*           <frame 1> is the frame for this first system, which should match
*                 the names in kf/gen_util/frame_to_fra.f.
*           <sys 2> is the second (comparison) system.  Sys 1 will be
*                 rotated and translated into this frame.  Differences
*                 will be given between sys1 and sys2 in this frame.
*           <frame 2> is the frame for this second system.  Same
*                 choices as above. 
*           <outname> is name of file for output of frame (overwritten)
*           <out_frame> is the frame for the output field (as above)
*           <fundamental sites> names of sites to be used for rotation
*                 and translation.  If not given then all sites are used. 
*                  ALL may be
*                 given as name and all common sites will be used.
*                 A '-' in front of name will stop it being used.
*           <height wght> Weight to be given to heights in the
*                 transformation determination.  If 1 then equal
*                 weight given.  If 0 then no weight is given
*                 (default is not to used heights)

      include 'velrot.h'
 
 
****  Read the runstring of the program
 
      write(*,100) velrot_ver
 100  format(/,' VELROT: Velocity field comparison and combination',
     .         ' Version ',a,/)
 
      call get_velrot_run
 
***** Read in system 1 coordinates
      call read_sys(sys1_file, sys1_coord, sys1_cov, sys1_names,
     .            num_sys1, max_sites)
 
*     Now rotate the system into the output frame
      frame_epoch = 2000.d0
      call frame_update(sys1_coord, sys1_cov, num_sys1, max_sites,
     .                sys1_frame, out_frame, frame_epoch)
 
*     Now read in system 2 coordinates
 
      call read_sys(sys2_file, sys2_coord, sys2_cov, sys2_names,
     .            num_sys2, max_sites)
 
*     Now rotate the system into the output frame
      call frame_update(sys2_coord, sys2_cov, num_sys2, max_sites,
     .                sys2_frame, out_frame, frame_epoch)
 
 
*     Now read in the sites to be used in the transformation
      call read_fund_sites
 
*     Now transform system 1 to system 2.
 
      call transframe
 
****  Finally write out the system 1 and system 2 frames transformed
*     into common system.
      
      call output_sum 
      call update_tran 

      if( av_dist.gt.0 ) call av_frame
 
      call output_frame
 
***** Thats all
      end
 
CTITLE GET_CORCOM_RUN
 
      subroutine get_velrot_run

      implicit none
 
*     Routine to read the runstring of the velrot program. (see
*     main for order)
 
      include 'velrot.h'
 
* LOCAL VARIABLES
 
*   len_run - Length of runstring parameter
*   rcpar   - Read runstring
*   trimlen - Length of string
*   ierr    - IOSTAT error.
*   indx    - Position in string
 
      integer*4 len_run, rcpar, trimlen, ierr
 
 
*   runstring   - Entry from runstring
 
      character*32 runstring
 
****  Start be getting the sys 1 file name
      eq_dist = 0
      cp_dist = 0
      param_opt = 'TR'

      len_run = rcpar(1, sys1_file )
      if( len_run.le.0 ) then
          call proper_runstring('velrot.hlp', 'velrot/sys1 file',-1)
      end if
 
****  Get frame for sys 1
      len_run = rcpar(2, sys1_frame)
      if( len_run.le.0 ) then
          call proper_runstring('velrot.hlp', 'velrot/sys1 frame',-1)
      end if
      call casefold(sys1_frame)
 
****  Now get the system 2 file name (strictly does not need to
*     be entered.)
      len_run = rcpar(3, sys2_file )
      if( len_run.le.0 ) then
 
*         Clear the name so that we know their is no second frame
          sys2_file = ' '
 
      end if
 
****  Get frame for sys 2
      len_run = rcpar(4, sys2_frame)
      if( len_run.le.0 ) then
*         Check to see if file given.  If no file given then
*         just to set none, otherwise give user help
          if( trimlen(sys2_file).eq.0 ) then
              sys2_frame = 'NONE'
          else
              call proper_runstring('velrot.hlp','velrot/sys2 frame',1)
          end if
      end if
      call casefold(sys2_frame)
 
*     Now get the name of the output file.  If none given then
*     write to screen.
      len_run = rcpar(5, out_file)
      if( len_run.le.0 ) out_file = '6'
 
*     Get the name of output frame.  (If not given, sys1 frame assumed)
*     Dec years is divided from frame by _ character
      len_run = rcpar(6, runstring)
      call casefold(runstring)
      if( len_run.eq.0 ) then
*        Nothing entered  so assume sys1 frame, and take 1990.0 as epoch
          out_frame = sys1_frame 
      else
          out_frame = runstring
      end if
      frame_epoch = 2000
      
 
****  Get name of file with fundamental sites
      len_run = rcpar(7, fund_file )
      if( len_run.le.0 ) fund_file = ' '
 
***** Get the weight of the heights
      len_run = rcpar(8, runstring )
      if( len_run.gt.0 ) then
          read(runstring,*,iostat=ierr ) height_weight
          call report_error('IOSTAT',ierr,'decod',runstring,
     .        0, 'velrot/Height Weight')
          if( ierr.ne.0 ) then
              call proper_runstring('velrot.hlp',
     .                        'velrot/Height Weight',1)
          end if
*         Don't let the weight get too small or we will have
*         numerical problems.
          if( height_weight.lt.1.d-6 ) height_weight = 1.d-6
      else
 
*         Set default weight 1 (not a value less than 1 downweights
*         the heights
          height_weight = 1.d0
      end if
 
****  Finally see what options have been passed for the type
*     transformation to be done.  Initially we assume a full
*     rotation and translation model.
      len_run = rcpar(9,param_opt)
      if( len_run.gt.0 ) then
          num_parn = 0
          call casefold(param_opt)
          if( index(param_opt,'T').gt.0 ) num_parn = 3
          if( index(param_opt,'R').gt.0 ) num_parn = 6
          if( index(param_opt,'S').gt.0 ) num_parn = 7
          if( index(param_opt,'L').gt.0 ) num_parn = 2
          if( num_parn.gt.7 ) then
              write(*,910) param_opt(1:len_run)
 910          format('**WARNING** Too many options in ',a,
     .              ' Setting options to TR')
              num_parn = 6
              param_opt = 'TR'
          endif
      else
          param_opt = 'TR'
          num_parn = 6
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_SYS
 
      subroutine read_sys(sys_file, sys_coord, sys_cov,
     .            sys_names, num_sys, max_sites)

      implicit none
 
*     This routine will read in a systems coordinates, velocities,
*     epochs and site names.  Any non-blank character is col 1 is
*     treated as a comment.  Any decoding errors will be igorned.
 
      include '../includes/const_param.h'

* PASSED VARIABLES
 
*   max_sites   - Max number of sites
*   num_sys     - Number of sites in this system
 
      integer*4 max_sites, num_sys
 
*   sys_coord(3,2,max_sites)    - Coordinates and velocities
*               - of sites (m and m/yr)
*   sys_cov(3,3,max_sites)    - velocity covariance matrix
 
      real*8 sys_coord(3,2,max_sites), sys_cov(3,3,max_sites)
 
*   sys_names(max_sites)    - Names of the sites
*   sys_file        - Name of file.
 
      character*(*) sys_names(max_sites), sys_file
 
* LOCAL VARIABLES
 
*   ierr, jerr  - IOSTAT errors
*   trimlen - Length of string.
*   i,j     - Loop counters
*   indx    - Position in string
 
      integer*4 ierr, jerr, trimlen, i, j, k, indx
 
*   values(12)   - Position and velocity read from file
*   pos_xyz(3)   - coordinates of site (height assumed zero)
*   vel_xyz(3)   - Velocity of site
*   geod_pos(3)   - Co-latitude, long and height
*   vel_neu(3)   - Velocity NEU
*   loc_coord(3) - Local coordinates returned from rot_geod
*   rot_mat(3,3) - Rotation matrix from NEU to XYZ
*   cov_neu(3,3), cov_xyz(3,3) -- Velocity covariance matrices.

 
      real*8 values(12), pos_xyz(3), vel_xyz(3), geod_pos(3),
     .       vel_neu(3), loc_coord(3), rot_mat(3,3)
      real*8 cov_neu(3,3)
 
*   tname   - Temp name read from file
 
      character*8 tname, cd
 
*   line    - Line read from file
 
      character*256 line
 
****  Try to open the file (see if name in non-blank)
      num_sys = 0
      if( trimlen(sys_file).gt.0 ) then
          open(100, file=sys_file, status='old', iostat=ierr)
          call report_error('IOSTAT',ierr,'open', sys_file,1,
     .                    'velrot/read_sys')
 
****      OK, FIle open now start reading
          i = 0
          do while( ierr.eq.0 .and. i.lt.max_sites)
              read(100,'(a)',iostat=ierr) line
              if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .            trimlen(line).gt.0 ) then
 
                  indx = 0

*                 Read all the values in the line upto the station
*                 name
                  do j = 1, 12
                     call read_line(line, indx, 'R8', jerr, values(j), 
     .                              cd )
                  end do

*                 Report if line is too short or if error occurred
                  if( jerr.ne.0 ) then
                      write(*,110) jerr, line(1:trimlen(line))
 110                  format('IOSTAT error ',I5,' reading ',a)
                  else

*                     If no error, continuing reading the line to get
*                     site name
                      call GetWord(line, tname, indx)

*                     Convert the geodetic coordinates to XYZ
                      geod_pos(1) = (90-values(2))*pi/180
                      geod_pos(2) = values(1)*pi/180
                      geod_pos(3) = 0.d0
                      call GEOD_to_XYZ(geod_pos, pos_xyz)

*                     Now rotate to velocites back into XYZ
                      vel_neu(1) = values(4)/1.d3
                      vel_neu(2) = values(3)/1.d3
                      vel_neu(3) = values(10)/1.d3
                      call rotate_geod(vel_neu,vel_xyz,'NEU','XYZ', 
     .                                 pos_xyz, loc_coord, rot_mat)
                      
*                     Make the NEU covariance matrix.
                      do j = 1,3
                         do k = 1,3
                            cov_neu(j,k) = 0.d0
                         end do
                      end do
                      cov_neu(1,1) = values(8)**2*1.d-6
                      cov_neu(2,2) = values(7)**2*1.d-6
                      cov_neu(1,2) = values(9)*values(7)*values(8)*1.d-6
                      cov_neu(2,1) = cov_neu(1,2)
                      cov_neu(3,3) = values(12)**2*1.d-6

*                     OK, Now save the values
                      i = i + 1
                      call check_lmax(i,'S')
                      sys_names(i) = tname
                      do j = 1,3
                         sys_coord(j,1,i) = pos_xyz(j)
                         sys_coord(j,2,i) = vel_xyz(j)
*                        Save the NEU version of the covariance matrix
                         do k = 1,3
                            sys_cov(k,j,i) = cov_neu(k,j)
                         end do
                      end do
*                             ! No decode error
                  end if
*                             ! Non-comment, and no error
              end if
*                             ! Looping over the input file
          end do
 
*****     Save the number of sites
          num_sys = i
          write(*,200) num_sys, sys_file(1:trimlen(sys_file))
 200      format(' There are ',i5,' sites in sys file ',a)
 
 
*                             ! Name actually given
      end if
 
****  Thats all
      return
      end
 
CTITLE FRAME_UPDATE
 
      subroutine frame_update(sys_coord, sys_epoch, num_sys,
     .            max_sites, sys_frame, out_frame, frame_epoch)
 
      implicit none

*     This routine will update the coordinate system frame from
*     the input one to the output frame.
 
      include '../includes/const_param.h' 
 
 
* PASSED VARIABLES
 
*   max_sites   - Max number of sites
*   num_sys     - Number of sites in this system
 
      integer*4 max_sites, num_sys
 
*   sys_coord(3,2,max_sites)    - Coordinates and velocities
*               - of sites (m and m/yr)
*   sys_epoch(max_sites)    - Epochs of the site positions
*   frame_epoch - Epoch for output frame
 
      real*8 sys_coord(3,2,max_sites), sys_epoch(max_sites),
     .    frame_epoch
 
*   sys_frame   - Frame for the input system
*   out_frame       - Frame for output system
 
      character*(*) sys_frame, out_frame
 
* LOCAL VARIABLES
 
*   i,j         - Loop counters
 
      integer*4 i,j
 
*   rot_vec(3)  - Rotation vector between the two frames.
*               - If the magnitude of this vetor is zero
*               - this routine simply returns
*   sang        - Sin of angle from cross product routine
*   xyz_vel(3)  - XYZ incremental velocity due frame
*               - difference.
 
      real*8 rot_vec(3), sang, xyz_vel(3)
 
****  First get the rotation vector bewteen the two frames
 
      call frame_to_frame(sys_frame, out_frame, rot_vec )
 
      write(*,100) sys_frame, out_frame,
     .             (rot_vec(j)*180.d6/pi,j=1,3)
 100  format(' Rotating from ',a10, ' to ',a10, ' using rotation',
     .       ' vector ',3F12.6,' degs/Myrs')
 
****  Check magntiudes.  (If very small then donot bother.
 
      if( abs(rot_vec(1)) + abs(rot_vec(2)) +
     .    abs(rot_vec(3)).lt.          1.d-10 ) RETURN
 
***** Loop over all of the sites and update position (if different
*     from frame_epoch) and velocity
 
      do i = 1, num_sys
          call cross_prod(rot_vec, sys_coord(1,1,i), xyz_vel, sang)
 
*         Update position and velocity
          do j = 1,3
              sys_coord(j,2,i) = sys_coord(j,2,i) + xyz_vel(j)
          end do
      end do
 
****  Thats all
      return
      end
 
 
CTITLE READ_FUND_SITES
 
      subroutine read_fund_sites

      implicit none
 
*     This routine will read the list of sites to be used in the
*     transformation between sys 1 and sys 2.  If no file name is
*     given then all possible sites will be used.
 
      include 'velrot.h'
 
* LOCAL VARIABLES
 
*   ierr, jerr  - IOSTAT errors
*   i,j         - Loop counter
*   ns1, ns2        - Matching sites from sys1 and sys2 lists
*   trimlen     - Length of string
 
      integer*4 ierr, i,j,k,l,n, ns1, ns2, trimlen, indx, jerr

*   link(2)     - Set true if sites are to be linked
      logical link(2)

*   dl          - site separtation (m)

      real*8 dl
 
*   tname(2)     - Name of fundamental site read from file.  May
*                 contain + or - at front to turn it off
*                 The two names are from sys 1 and sys 2.  If
*                 second name not given assumed to be same as first.

 
      character*16 tname(2), cd

*   line        - Line read from input file

      character*256 line
 
****  Try to open the file (see if name in non-blank)
      num_fund = 0 
      if( trimlen(fund_file).gt.0 ) then
          open(100, file=fund_file, status='old', iostat=ierr)
          call report_error('IOSTAT',ierr,'open', fund_file,1,
     .                    'velrot/read_fund')
 
****      OK, FIle open now start reading
          i = 0
          do while( ierr.eq.0 .and. i.lt.max_sites)
              read(100,'(a)', iostat=ierr) line
              if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .            trimlen(line).gt.0 ) then
 
                  indx = 0
                  call getword(line, tname(1), indx)
                  call casefold(tname(1))

****              See if distance option passed.  If it is passed
*                 then we used sites separated by upto this distance
                  if( tname(1)(1:7).eq.'EQ_DIST' ) then
                      call read_line(line, indx, 'R8', jerr, eq_dist, 
     .                              cd)
                      cp_dist = eq_dist

*                     Now loop over the sites seeing which one match
                      do j = 1, num_sys1
                         do k = 1, num_sys2
                            dl = 0
                            do l = 1,3
                               dl = dl + (sys1_coord(l,1,j)-
     .                                    sys2_coord(l,1,k))**2
                            end do
                            dl = sqrt(dl)
                            if( dl.lt.eq_dist .and.i.lt.max_links ) then
                                i = i + 1
                                call check_lmax(i,'L')
                                fund_links(1,i) = j
                                fund_links(2,i) = k
                            end if
                         end do
                      end do
                  else if( tname(1)(1:5).eq.'NAMES' ) then
*                     Assign all matching sites
                      i = 0
                      do j = 1, num_sys1
                          indx=0
                          call get_cmd(sys1_names(j), sys2_names,
     .                        -num_sys2,  ns2, indx)
                          if( ns2.gt.0 ) then
                              i = i + 1
                              call check_lmax(i,'L')
                              fund_links(1,i) = j
                              fund_links(2,i) = ns2
                          end if
                      end do
 
                  else if( tname(1)(1:7).eq.'CP_DIST' ) then 
                      call read_line(line, indx, 'R8', jerr, cp_dist, 
     .                              cd )
                  else if( tname(1)(1:7).eq.'AV_DIST' ) then
                      call read_line(line, indx, 'R8', jerr, av_dist, 
     .                              cd )
                  else

*                      See if second name passed  
                       tname(2) = ' '
                       call getword(line, tname(2), indx)
                       call casefold(tname(2))
                       if( tname(2)(1:1).eq.' ' ) tname(2) = tname(1)

*                      See if either name starts with a - sign (remove
*                      site from list)
                       link(1) = .true.
                       link(2) = .true.
                       if( tname(1)(1:1).eq.'-' ) then
                           link(1) = .false. 
                           cd = tname(1)
                           tname(1) = cd(2:)
                       else if ( tname(1)(1:1).eq.'+' ) then
                           cd = tname(1)
                           tname(1) = cd(2:)
                       end if
*                      Check the other sites
                       if( tname(2)(1:1).eq.'-' ) then
                           link(2) = .false. 
                           cd = tname(2)
                           tname(2) = cd(2:)
                       else if ( tname(2)(1:1).eq.'+' ) then
                           cd = tname(2)
                           tname(2) = cd(2:)
                       end if

*                      See if we can match sites
                       indx = 1
                       call get_cmd(tname(1), sys1_names,-num_sys1,
     .                            ns1, indx)
                       indx = 1
                       call get_cmd(tname(2), sys2_names,-num_sys2,
     .                            ns2, indx)

****                  See if we match sites sites
                      if( ns1.gt.0 .and. ns2.gt.0 .and. link(1)
     .                    .and. link(2) ) then
                           i = i + 1
                           call check_lmax(i,'L')
                           fund_links(1,i) = ns1
                           fund_links(2,i) = ns2
                      end if

                      if( ns1.gt.0 .and. .not.link(1) ) then
*                          Remove any link which has this site name
                           k = 0
                           do while ( k.lt.i )
                              k = k + 1
                              if( fund_links(1,k).eq.ns1 ) then
                                  do j = k,i-1
                                     fund_links(1,j) = fund_links(1,j+1)
                                     fund_links(2,j) = fund_links(2,j+1)
                                  end do
                                  i = i - 1
                                  k = k - 1
                              end if
                           end do
                      end if
                      if( ns2.gt.0 .and. .not.link(2) ) then
*                          Remove any link which has this site name
                           k = 0
                           do while ( k.lt.i )
                              k = k + 1
                              if( fund_links(2,k).eq.ns2 ) then
                                  do j = k,i-1
                                     fund_links(1,j) = fund_links(1,j+1)
                                     fund_links(2,j) = fund_links(2,j+1)
                                  end do
                                  i = i - 1
                                  k = k - 1
                              end if
                           end do
                      end if
                  end if
*                             ! Non-comment, and no error
              end if
*                             ! Looping over the input file
          end do
*                             ! Name actually given
      else 
 
*         Assign all matching sites
          i = 0
          do j = 1, num_sys1
              indx=0
              call get_cmd(sys1_names(j), sys2_names,-num_sys2,
     .            ns2, indx)
              if( ns2.gt.0 ) then
                  i = i + 1
                  call check_lmax(i,'L')
                  fund_links(1,i) = j
                  fund_links(2,i) = ns2
              end if
          end do
      end if

*     Clean up the fund links and remove any duplicates
      j = 0
      do while ( j .lt. i )
         j = j + 1
         n = 0
         do k = j+1, i
            if ( (fund_links(1,j).eq.fund_links(1,k) .and.
     .           fund_links(2,j).eq.fund_links(2,k)) .or.
     .           (fund_links(1,j).eq.fund_links(2,k) .and.
     .           fund_links(2,j).eq.fund_links(1,k)) ) then
                n = k
            end if
         end do
*        if n != 0 then we have a duplicate, so remove
         if( n.gt.0 ) then
             do k = n+1,i
                fund_links(1,k-1) = fund_links(1,k)
                fund_links(2,k-1) = fund_links(2,k)
             end do
             i = i - 1
         endif
      end do
 
***** Save the number of sites
      num_fund = i
      write(*,200) num_fund, fund_file(1:trimlen(fund_file))
 200  format(' There are ',i5,' matching sites in fundamental',
     .        ' file ',a)
 
      return
      end 

CTITLE CHECK_LMAX

      subroutine check_lmax(i,type)

      implicit none

*     Routine to make sure we do not exceed the number of fundamental
*     links or sites

      include 'velrot.h'

* i   -- Number of links so far
* type -- Type to check. L for links, S for sites

      integer*4 i
      character*(*) type

*     Check to see if above max and stop if need be
      if ( i.gt. max_links .and. type.eq.'L' ) then
         write(*,100) i
 100     format('**DISASTER** Too many link sites. ',i6,
     .          ' is maximum allowed.  Try reducing eq_dist')
         stop 'VELROT: Too many links, maximim allowed exceeded'
      end if

      if ( i.gt. max_sites .and. type.eq.'S' ) then
         write(*,120) i
 120     format('**DISASTER** Too many sites. ',i6,
     .          ' is maximum allowed.  Try using vel files with',
     .          ' fewer sites')
         stop 'VELROT: Too many sites, maximim allowed exceeded'
      end if

      return
      end
 
CTITLE TRANSFRAME
 
      subroutine transframe

      implicit none
 
*     Routine to compute the transformation parameters between
*     system 1 and system 2.  The loops over all of the fundamental
*     sites and computes the transformation from them.
 
      include '../includes/const_param.h'
      include 'velrot.h'
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   ierr    - invert_vis error
*   ipivot(7)   - Pivot array for inversion
*   ns1, ns2    - Sites numbers of the two being compared
 
 
      integer*4 i,j, ipivot(7), ns1, ns2
 
*   dprefit - Change between postfit and prefit sum of
*           - residuals squared
*   neu_part(3,7)   - Partials of NEU with respect to transformation
*           - parameters
*   xyz_part(3,7)   - Partials of parameters wrt XYZ coords.
*   rot_matrix(3,3) - Rotation between XYZ and NEU
*   dx(3), dn(3)    - Differences bwteen sys1 and sys2 in XYZ and
*           - NEU spaces
*   loc_coord(3)    - Local coordinates of site.
*   scale(7)        - Scale factors used by invert_vis.
 
 
      real*8 dprefit, neu_part(3,7), xyz_part(3,7), rot_matrix(3,3),
     .    dx(3), dn(3), loc_coord(3), scale(7), weights(3)
 
****  Initialize the system
 
      call clear_norm

****  Loop over the fundamental sites
      do i = 1, num_fund
 
*         get the transformation from XYZ to NEU for this site
          ns1 = fund_links(1,i)
          ns2 = fund_links(2,i)
 
*         Get velocity differences
          do j = 1,3
              dx(j) = sys2_coord(j,2,ns2) - sys1_coord(j,2,ns1) 
          end do
          call rotate_geod(dx, dn,  'XYZ', 'NEU', sys1_coord(1,1,ns1),
     .                    loc_coord, rot_matrix )
 
*****     Now generate the partials with respect to rotation,scale,and
*         translation.
          call get_parts(sys1_coord(1,1,ns1), rot_matrix, num_parn,
     .                neu_part, xyz_part, max_sites)

*         Get the combined covariance matrix of the velocities
          weights(1) = 1.d0/(sys1_cov(1,1,ns1)+sys2_cov(1,1,ns2))
          weights(2) = 1.d0/(sys1_cov(2,2,ns1)+sys2_cov(2,2,ns2))
          weights(3) = height_weight/
     .                      (sys1_cov(3,3,ns1)+sys2_cov(3,3,ns2)) 
   
*****     Now increment normal equations
          call increment_norm( dn, num_parn, neu_part, weights,
     .            norm_eq, bvec, sum_prefit, sum_weight, num_data )
 
      end do
 
***** Now solve the system
      do i = 1, num_parn
          trans_parm(i) = bvec(i)
      end do
       
      call invert_vis(norm_eq, trans_parm, scale, ipivot, num_parn,
     .                7, 1)
           
***** Compute the change in chi**2
      dprefit = 0
      do i = 1, num_parn
          do j = 1, num_parn
              dprefit = dprefit + bvec(i)*norm_eq(i,j)*bvec(j)
          end do
      end do
 
****  update the statistics
      sum_postfit = sum_prefit - dprefit
      if( num_data-num_parn.gt.0 ) then
          chi_fit = sqrt(sum_postfit/(num_data-num_parn))
          rms_fit = sqrt(num_data/sum_weight)*chi_fit
* MOD TAH 151203: If chi is <0.1 then reset to 1.0 
          if( chi_fit.lt.0.10 ) then
             write(*,220) chi_fit
 220         format('# Chi of fit less than 0.10 (',F6.4,').',
     .              ' Resetting to 1.0')
             chi_fit = 1.0d0
          endif
      else
          chi_fit = 1.d0
          rms_fit = 0.d0

*         OK: Not enough clear everything
          do i = 1, 7
             trans_parm(i) = 0.d0
             do j = 1,7
                norm_eq(i,j) = 0.d0
             end do
          end do
      end if
 
*     Now convert transformation parameters to output format
      trans_parm_out(1) = trans_parm(1) * 1000.d0
      trans_parm_out(2) = trans_parm(2) * 1000.d0
      trans_parm_out(3) = trans_parm(3) * 1000.d0
      trans_parm_out(4) = trans_parm(4) * rad_to_mas
      trans_parm_out(5) = trans_parm(5) * rad_to_mas
      trans_parm_out(6) = trans_parm(6) * rad_to_mas
      trans_parm_out(7) = trans_parm(7) * 1.d9
 
*     Now do sigmas
      trans_sigma(1) = sqrt(norm_eq(1,1)) * 1000.d0 * chi_fit
      trans_sigma(2) = sqrt(norm_eq(2,2)) * 1000.d0 * chi_fit
      trans_sigma(3) = sqrt(norm_eq(3,3)) * 1000.d0 * chi_fit
      trans_sigma(4) = sqrt(norm_eq(4,4)) * rad_to_mas * chi_fit
      trans_sigma(5) = sqrt(norm_eq(5,5)) * rad_to_mas * chi_fit
      trans_sigma(6) = sqrt(norm_eq(6,6)) * rad_to_mas * chi_fit
      trans_sigma(7) = sqrt(norm_eq(7,7)) * 1.d9 * chi_fit
 
*     Convert rms_fit to mm
      rms_fit = rms_fit * 1000.d0
 
****  Thats all
      return
      end 
 
CTITLE CLEAR_NORM
 
      subroutine clear_norm

      implicit none
 
*     This routine will initialize the estimation arrays.
 
      include 'velrot.h'
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
 
      integer*4 i,j
 
*     Initialize normal equations
      do i = 1, num_parn
          bvec(i) = 0.d0
          do j = 1, num_parn
              norm_eq(i,j) = 0.d0
          end do
      end do

****  See if some parameters are not be estimated
      if( index(param_opt,'T').eq. 0 .and. 
     .    index(param_opt,'L').eq. 0 ) then
          do i = 1,3
             norm_eq(i,i) = 1.d14
          end do
      end if

*     See if scale estimated
      if( index(param_opt,'S').eq.0 ) then
          norm_eq(7,7) = 1.d14
      end if      

 
*     Initialize the statistics variables
 
      sum_prefit = 0.d0
      sum_weight = 0.d0
      num_data = 0
 
 
*     Thats all
      return
      end
 
CTITLE GET_PARTS
 
      subroutine get_parts(coord, rot_matrix, num_parn,
     .                neu_part, xyz_part, max_sites)

      implicit none
 
*     Routine to compute the partials of NEU with respect to the
*     transformation parameters
 
 
* PASSED VARIABLES
 
*   max_sites       - MAx number of sites
*   num_parn            - Number of parameters to be estimated
*                   - IF 6 then no scale, if 7 then scale
 
      integer*4 max_sites, num_parn
 
*   coord(3)            - Coordinates of site
*   rot_matrix(3,3) - Rotation from XYZ to NEU
*   neu_part(3,7)       - Partials wrt NEU for this obervation
*   xyz_part(3,7)   - Partials wrt XYZ for this obs.
 
      real*8 coord(3), rot_matrix(3,3), neu_part(3,7), xyz_part(3,7)
 
* LOCAL VARIABLES
 
*   i,j             - Loop counters
 
 
      integer*4 i,j,k
 
****  First form the XYZ partials.  Clear first
      do i = 1,3
          do j = 1,7
              xyz_part(i,j) = 0.d0
          end do
      end do
 
****  Now put in the correct elements
*     Translation
      xyz_part(1,1) = 1.d0
      xyz_part(2,2) = 1.d0
      xyz_part(3,3) = 1.d0

*     Rotation partials
      xyz_part(2,4) = -coord(3)
      xyz_part(3,4) =  coord(2)
      xyz_part(1,5) =  coord(3)
      xyz_part(3,5) = -coord(1)
      xyz_part(1,6) = -coord(2)
      xyz_part(2,6) =  coord(1)
 
*     See if new scale partials
      if( num_parn.eq.7 ) then
          xyz_part(1,7) = coord(1)
          xyz_part(2,7) = coord(2)
          xyz_part(3,7) = coord(3)
      else
          xyz_part(1,7) = 0.d0
          xyz_part(2,7) = 0.d0
          xyz_part(3,7) = 0.d0
      end if
 
***** Now apply the rotation matrix to get NEU partials
      do i = 1,3
          do j = 1, num_parn
              neu_part(i,j) = 0.d0
              do k = 1,3
                  neu_part(i,j) = neu_part(i,j) +
     .                    rot_matrix(i,k)*xyz_part(k,j)
              end do
          end do
          do j = num_parn+1,7
             neu_part(i,j) = 0.d0
          end do
      end do
 
****  Thats all
      return
      end
 
CTITLE INCREMENT_NORM
 
      subroutine increment_norm( dn, num_parn, neu_part, wghs,
     .            norm_eq, bvec, sum_prefit, sum_weight, num_data )

      implicit none
 
*     Routine to increment normal equations and solution vector.
*     Prefit residuals are also accumulated for later computation
*     of the fit to the data.
 
 
* PASSED VARIABLES
 
*   num_parn        - Number of parameters being estimated
*   num_data        - Number of components summed into prefit
 
      integer*4 num_parn, num_data
 
*   dn(3)       - Difference in NE and U
*   wghs(3)   - Inverse sigma**2 on the components
*   neu_part(3,7)   - Partials of NEU wrt transformation
*               - parameters
*   norm_eq(7,7)    - Normal equations
*   bvec(7)     - Solution vector
*   sum_prefit  - Sum of prefit residuals
*   sum_weight  - Sum of weights
 
      real*8 dn(3), neu_part(3,7), norm_eq(7,7),
     .    bvec(7), sum_prefit, wghs(3), sum_weight
 
* LOCAL VARIABLES
 
*   i,j,k       - Loop counters
*   num_use     - Number of components to use
 
      integer*4 i,j,k, num_use
 
 
****  See what weights we should give (N and E always get 1)
* MOD TAH 040326: Changed tolerance on test to 1.d-5 so that
*     small height sigmas (compared to horizontal) will not
*     cause the number of components to count to get from 2 to 3.      
*     if( abs(wghs(3)/wghs(2)).lt.1.d-6 ) then
      if( abs(wghs(3)/wghs(2)).lt.1.d-5 ) then
          num_use = 2
      else
          num_use = 3
      end if
      

****  Now sum in the contributions from each component
      do i = 1, num_use
 
*         Now loop over parameters
          do j = 1,7
              bvec(j) = bvec(j) + neu_part(i,j)*dn(i)*wghs(i)
              do k = 1,7
                  norm_eq(j,k) = norm_eq(j,k) +
     .                neu_part(i,j)*wghs(i)*neu_part(i,k)
              end do
          end do
 
*         Sum the prefit statistics
          sum_prefit = sum_prefit + dn(i)**2*wghs(i)
          sum_weight = sum_weight + wghs(i)
          num_data   = num_data + 1
      end do
 
****  Thats all
      return
      end

CTITLE OUTPUT_SUM
 
      subroutine output_sum

      implicit none
 
*     This routine will output the results from the transformation.
      
      include '../includes/const_param.h'
      include 'velrot.h'
 
* LOCAL VARIABLES
 
*   i,j,k   - Loop counters
*   ierr    - IOSTAT error
*   trimlen - Length of string
*   run_time(7) - Time that program was run
 
      integer*4 i,j,k, ierr, trimlen, run_time(7), ns1, ns2
 
*   dx(3), dn(3)    - Changes in the position in XYZ and NEU
*   loc_coord(3)    - Local coordinates of sites
*   rot_matrix(3,3) - Rotation between XYZ and NEU
*   neu_part(3,7), xyz_part(3,7)    - Partials derivatives
*               - wrt to transformation parameters.
*   sectag      - Seconds part of time
 
 
      real*8 dx(3), dn(3), loc_coord(3), rot_matrix(3,3),
     .    neu_part(3,7), xyz_part(3,7), sectag, 
     .     cov_neu(3,3), temp_covar(7,7)

* ADDED Postfit sigma calculation
      real*8 summ(4)     ! Sum for weighted mean NEU  dp/sp^2
     .,      sumv(4)     ! Sum for variance (dp/sp)^2
     .,      sumw(4)     ! Sum of weights (1/sp)^2
      real*8 wmean(4), wrms(4), nrms(4)  ! Final stats; 4th 
                         ! entry is N+E 
 
*   unit_labs(7)        - Units for parameters
*   param_labs(7)   - Parameter names
 
      character*8 unit_labs(7), param_labs(7), comp_labs(4)

*   symbol - Character at end of name if separation < cp_dist

 
      data param_labs /   'X-Offset', 'Y-Offset', 'Z-Offset', 
     .                    'X-Rot', 'Y-Rot', 'Z-Rot','Scale' /
 
      data unit_labs /  '(mm/yr)', '(mm/yr)', '(mm/yr)', 
     .                  '(mas/yr)','(mas/yr)','(mas/yr)','(ppb/yr)' /

      data comp_labs /  'North','East','Up','Horz' /
 
 
****  Open up the output file
      out_unit = 200
      if( out_file(1:1).ne.'6' ) then
          open(out_unit, file = out_file, iostat=ierr, status='unknown')
          call report_error('IOSTAT',ierr,'open', out_file,1,
     .        'out_frame/open')
          if( ierr.ne.0 ) out_unit = 6
      else
          out_unit = 6
      end if
 
***** Write out the preamble file information.
      call systime( run_time, sectag)
 
      write(out_unit,100 ) (run_time(i),i=1,5), velrot_ver,
     .    sys1_frame, sys1_file(1:max(1,trimlen(sys1_file))),
     .    sys2_frame, sys2_file(1:max(1,trimlen(sys2_file))),
     .    fund_file(1:max(1,trimlen(fund_file))),
     .    out_frame, param_opt, eq_dist, cp_dist, av_dist,
     .    height_weight
 
 100  format('* VELROT Run on ',i4,2('/',i2),1x,i2,':',i2,
     .    ' Version ',a,/,
     .    '* SYSTEM 1 File    : ',a,1x,a,/,
     .    '* SYSTEM 2 File    : ',a,1x,a,/,
     .    '* FUNDAMENTAL File : ',11x,a,/,
     .    '* OUTPUT FRAME     : ',a,t32,' PARAM_OPT ',a,/,
     .    '* EQ_DIST          : ',6x,F12.1,' m,',
     .    ' CP_DIST : ',F12.1,' m*',/,
     .    '* AV_DIST          : ',6x,F12.1,' m',
     .    '* HEIGHT WEIGHT    : ',F10.6)
 
*     Now write out the transformation parameters
      if( num_data.gt.0 ) then
          write(out_unit,200)  num_data, num_fund, rms_fit,  chi_fit
200       format('* ',/, '* RMS fit for ',i8,' components',
     .        ' from ',i4,' stations was ', F10.2, 
     .        ' mm/yr, NRMS ',F8.2,/,
     .        '* Estimates of Transformation parameters are: ')

          do i = 1, num_parn
c             if( trans_sigma(i).gt.1.d-3 ) then
                 write(out_unit,220) i, param_labs(i),trans_parm_out(i),
     .                  trans_sigma(i), unit_labs(i)
220              format( '*',i5,1x,a,3x,F10.4,1x,F10.4,1x,a)
c             end if
          end do
      end if

*     Now print out the errors at the Fundamental stations
*     Clear stats
      summ = 0
      sumv = 0
      sumw = 0

 
      write(out_unit,300)
 300  format('*  Differences at the fundamental sites',/,
     .    '*   #  Name 1  Name Ref    dN (mm)    dE (mm)',
     .    '    dU (mm)    sN (mm)    sE (mm)    sU (mm)',
     .    '   sTN (mm)   sTE (mm)   sTU (mm)')
 
      do i = 1, num_fund
          ns1 = fund_links(1,i)
          ns2 = fund_links(2,i)
 
*         Get coordinate differences (before adjustment)
          do j = 1,3
              dx(j) = sys2_coord(j,2,ns2) - sys1_coord(j,2,ns1) 
          end do

****      Rotate into local frame to get rot_matrix for partials
          call rotate_geod(dx, dn,  'XYZ', 'NEU', sys1_coord(1,1,ns1),
     .                    loc_coord, rot_matrix )
 
*         Now get partials (in XYZ space )
          call get_parts(sys1_coord(1,1,ns1), rot_matrix, num_parn,
     .                neu_part, xyz_part, max_sites)
 
*         Update the adjustment
          do j = 1,3
              do k = 1, num_parn
                  dx(j) = dx(j) - xyz_part(j,k)*trans_parm(k)
              end do
          end do
 
*         Now rotate into local system
 
          call rotate_geod(dx, dn,  'XYZ', 'NEU', sys1_coord(1,1,ns1),
     .                    loc_coord, rot_matrix )

****      Now propagate the transformation uncertainity into the velocity
*         covariance matrix
          call var_comp(neu_part,norm_eq,cov_neu, temp_covar,3,7,1)
 
*         Now write the values
          write(out_unit,320) i, sys1_names(ns1),sys2_names(ns2),
     .            (dn(j)*1000.d0, j = 1,3),
     .            (sqrt(sys1_cov(j,j,ns1)+sys2_cov(j,j,ns2))*1000.d0,
     .                                                      j = 1,3),
     .            (sqrt(cov_neu(j,j))*1000.d0,j=1,3)
  320     format('A',i5,1x,a8,1x,a8,1x,9(F10.2,1x))

****      Sum up statistics
          do j = 1,3 
             summ(j) = summ(j) + 
     .                 dn(j)/(sys1_cov(j,j,ns1)+sys2_cov(j,j,ns2))
             sumv(j) = sumv(j) + 
     .                 dn(j)**2/(sys1_cov(j,j,ns1)+sys2_cov(j,j,ns2))
             sumw(j) = sumw(j) + 
     .                 1.d0/(sys1_cov(j,j,ns1)+sys2_cov(j,j,ns2))
          end do
 
      end do

***** Now finish statistics
      summ(4) = (summ(1)+summ(2))/2
      sumv(4) = (sumv(1)+sumv(2))/2
      sumw(4) = (sumw(1)+sumw(2))/2
      do j = 1,4 
         wmean(j) = summ(j)/sumw(j)*1000
         nrms(j) = sqrt(sumv(j)/num_fund)
         wrms(j) = sqrt((1.d0/sumw(j))*num_fund)*nrms(j)*1000
         write(*,420) comp_labs(j), num_fund, wmean(j), wrms(j), nrms(j)
         if( out_unit.ne.6 )  write(*,420) comp_labs(j), num_fund, 
     .                                        wmean(j), wrms(j), nrms(j)
 420     format('S Component ',a,' # ',i5,' WMean ',F6.2,' WRMS ',F6.2,
     .          ' mm/yr, NRMS ',F7.3)
      end do


****  Thats all
      return
      end 
 
CTITLE OUTPUT_FRAME
 
      subroutine output_frame

      implicit none
 
*     This routine will output the results from the transformation.
*     The coordinates from system 1 transformed to the system 2 frame
*     are first given, and the the system 2 values which were not in
*     system one are given.
      
      include '../includes/const_param.h'
      include 'velrot.h'
 
* LOCAL VARIABLES
 
*   i,j,k   - Loop counters
 
      integer*4 i,j, ns1
 
*   dx(3), dn(3)    - Changes in the position in XYZ and NEU
*   loc_coord(3)    - Local coordinates of sites
*   rot_matrix(3,3) - Rotation between XYZ and NEU
*   neu_part(3,7), xyz_part(3,7)    - Partials derivatives
*               - wrt to transformation parameters.
*   sectag      - Seconds part of time
 
 
      real*8  dn(3), loc_coord(3), rot_matrix(3,3),
     .    neu_part(3,7),  dsig(3),
     .    long, lat, rho, cov_neu(3,3), temp_covar(7,7)
 

*   symbol - Character at end of name if separation < cp_dist

      character*1 symbol
 
 
 
****  Now loop over all of the system 1 sites appying the
*     transformation.
      write(out_unit,400)
 400  format( /,'* SYSTEM 1 Velocities transformed to',
     .            ' SYSTEM 2 ',/,
     .          '*   Long.       Lat. ',8x,'E & N Rate ',4x,
     .          ' E & N Adj. ',4x,' E & N +-',2x,
     .          ' RHO ',6x,' H Rate   H adj.    +-',2x,'SITE',/,
     .          '* ',' (deg)      (deg) ',3x,3(7x,'(mm/yr)'),17x,
     .             '(mm/yr)' )


 
      do i = 1, num_sys1
 
****      Get the correction to the position due to the system
*         change. 
          call rotate_geod(sys1_coord(1,2,i), dn,  'XYZ', 'NEU', 
     .                     sys1_coord(1,1,i), loc_coord, rot_matrix )

*         Convert to mm/yr
          do j = 1,3
             dn(j) = dn(j)*1000.d0
          end do

****      Now propagate the transformation uncertainity into the velocity
*         covariance matrix
          call var_comp(neu_part,norm_eq,cov_neu, temp_covar,3,7,1)
           

*****     Convert back into velocity file format
          long = loc_coord(2)*180/pi
          lat  = (pi/2-loc_coord(1))*180/pi

          dsig(2) = sqrt(sys1_cov(2,2,i))*1000.d0
          dsig(1) = sqrt(sys1_cov(1,1,i))*1000.d0
          dsig(3) = sqrt(sys1_cov(3,3,i))*1000.d0 
          rho =  (sys1_cov(1,2,i))/(dsig(1)*dsig(2)/1.d6)

          call check_cp(i,0,cp_dist,symbol)

          write(out_unit,420) long, lat, dn(2), dn(1), dn(2), dn(1),
     .                        dsig(2), dsig(1), rho, dn(3), dn(3),
     .                        dsig(3), sys1_names(i), symbol
 420      format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,
     .               3(1x,f7.2), 1x,a8,a1)
 
      end do
 
***** Now output the system2 coordinates (comment out those which
*     are are already in the list.  NOTE: System 1 coordinates have
*     now been updated to system 2 frame.
 
      write(out_unit,500)
 500  format( /,'* SYSTEM 2 Velocities except those in',
     .            ' SYSTEM 1 ',/,
     .          '*   Long.       Lat. ',8x,'E & N Rate ',4x,
     .          ' E & N Adj. ',4x,' E & N +-',2x,
     .          ' RHO ',6x,' H Rate   H adj.    +-',2x,'SITE',/,
     .          '* ',' (deg)      (deg) ',3x,3(7x,'(mm/yr)'),17x,
     .             '(mm/yr)' )
 
      do i = 1, num_sys2
 
****      Get the correction to the position due to the system
*         change.  (If we want NEU we should call routine below
*         but we just want XYZ so don't bother.
 
*         Compute the difference in cooordinates if this site overlaps
          ns1 = 0
          do j = 1, num_sys1
              if( sys2_names(i).eq.sys1_names(j) ) then
                  ns1 = j
              end if
          end do
 
 
*         Now rotate into local system
          call rotate_geod(sys2_coord(1,2,i), dn,  'XYZ', 'NEU',
     .                sys2_coord(1,1,i),loc_coord, rot_matrix )
 
*          Convert to mm/yr
           do j = 1,3
             dn(j) = dn(j)*1000.d0
           end do
           

*****     Convert back into velocity file format
          long = loc_coord(2)*180/pi
          lat  = (pi/2-loc_coord(1))*180/pi

          dsig(2) = sqrt(sys2_cov(2,2,i))*1000.d0
          dsig(1) = sqrt(sys2_cov(1,1,i))*1000.d0
          dsig(3) = sqrt(sys2_cov(3,3,i))*1000.d0
          rho = (sys2_cov(1,2,i))/(dsig(1)*dsig(2)/1.d6)


          call check_cp(0,i,cp_dist,symbol)

          if( ns1.eq.0 ) then
             write(out_unit,520) ' ', long, lat, dn(2), dn(1), 
     .                        dn(2), dn(1),
     .                        dsig(2), dsig(1), rho, dn(3), dn(3),
     .                        dsig(3), sys2_names(i),symbol
! 520         format(a,f8.3,1x,f8.3,1x,6(1x,f7.2),1x,f6.3,2x,
!     .               3(1x,f7.2), 1x,a8,a1)
 520         format(a,f10.5,1x,F10.5,1x,6(1x,f7.2),1x,f6.3,2x,
     .               3(1x,f7.2), 1x,a8,a1)
          else
             write(out_unit,520) '-', long, lat, dn(2), dn(1), 
     .                        dn(2), dn(1),
     .                        dsig(2), dsig(1), rho, dn(3), dn(3),
     .                        dsig(3), sys2_names(i), symbol

          end if
      end do
 
****  Thats all
      return
      end 

CTITLE UPDATE_TRAN
 
      subroutine update_tran

      implicit none
 
*     This routine will update the velocities in each set of
*     coordinates so that they are in the correct frame.  It
*     will also modify the sigma to account for the uncertainity
*     in the transformations.
      
      include '../includes/const_param.h'
      include 'velrot.h'
 
* LOCAL VARIABLES
 
*   i,j,k   - Loop counters
*   ierr    - IOSTAT error
*   trimlen - Length of string
*   run_time(7) - Time that program was run
 
      integer*4 i,j,k 
*   dx(3), dn(3)    - Changes in the position in XYZ and NEU
*   loc_coord(3)    - Local coordinates of sites
*   rot_matrix(3,3) - Rotation between XYZ and NEU
*   neu_part(3,7), xyz_part(3,7)    - Partials derivatives
*               - wrt to transformation parameters.
 
 
      real*8 dx(3), dn(3), loc_coord(3), rot_matrix(3,3),
     .    neu_part(3,7), xyz_part(3,7),  dsig(3),
     .    long, lat, rho, cov_neu(3,3), temp_covar(7,7)
 
 
****  OK, adjust the values in sys 1 so that they are
*     in the same frame as sys 2.
 
      do i = 1, num_sys1
 
****      Get the correction to the position due to the system
*         change. 
          call rotate_geod(sys1_coord(1,2,i), dn,  'XYZ', 'NEU', 
     .                     sys1_coord(1,1,i), loc_coord, rot_matrix )
          call get_parts(sys1_coord(1,1,i), rot_matrix, num_parn,
     .                neu_part, xyz_part, max_sites)
 
*         Get adjustment to position
          do j = 1, 3
              dx(j) = 0.0d0
              do k = 1,num_parn
                  dx(j) = dx(j) + xyz_part(j,k)*trans_parm(k)
              end do
 
*             Now add to velocity
              sys1_coord(j,2,i) = sys1_coord(j,2,i) + dx(j)
          end do

          call rotate_geod(sys1_coord(1,2,i), dn,  'XYZ', 'NEU', 
     .                     sys1_coord(1,1,i), loc_coord, rot_matrix )

*         Convert to mm/yr
          do j = 1,3
             dn(j) = dn(j)*1000.d0
          end do

****      Now propagate the transformation uncertainity into the velocity
*         covariance matrix
          call var_comp(neu_part,norm_eq,cov_neu, temp_covar,3,7,1)
           

*****     Convert back into velocity file format
          long = loc_coord(2)*180/pi
          lat  = (pi/2-loc_coord(1))*180/pi

          dsig(2) = sqrt(sys1_cov(2,2,i)+cov_neu(1,1))*1000.d0
          dsig(1) = sqrt(sys1_cov(1,1,i)+cov_neu(2,2))*1000.d0
          dsig(3) = sqrt(sys1_cov(3,3,i)+cov_neu(3,3))*1000.d0 
          rho =  (sys1_cov(1,2,i)+cov_neu(1,2))/(dsig(1)*dsig(2)/1.d6)

****      Now save back in covariance matrix
          sys1_cov(1,1,i) = dsig(1)**2/1.d6
          sys1_cov(2,2,i) = dsig(2)**2/1.d6
          sys1_cov(3,3,i) = dsig(3)**2/1.d6
          sys1_cov(1,2,i) = rho*sqrt(sys1_cov(1,1,i)*sys1_cov(2,2,i))
      end do
 
 
****  Thats all
      return
      end


CTITLE CHECK_CP

      subroutine check_cp( ns1, ns2, dist, symbol)

      implicit none

*     Routine to see sites within in separation of dist

      include 'velrot.h'

* PASSED VARIABLES
* ns1, ns2 -- Site number in either system 1 or 2 (one should be
*             zero to say which system we want to check)
* dist     -- Distance that we must be less than to set symbol to *
* symbol   -- returned symbol (* if some site is < dist way, blank otherwise)

      integer*4 ns1, ns2
      real*8 dist
      character*(*) symbol

* LOCAL VARIABLE
* i, j  -- Loop counters
* dl    -- Length computation

      integer*4 i, j
      real*8 dl

****  Do each system separatley
      symbol = ' '
      if( ns1.gt.0 ) then
          do i = 1, num_sys2
             dl = 0.d0
             do j = 1,3
                dl = dl + (sys1_coord(j,1,ns1)-sys2_coord(j,1,i))**2
             end do
             if( sqrt(dl).lt.dist ) symbol = '*'
          end do
      else if( ns2.gt.0 ) then
          do i = 1, num_sys1
             dl = 0.d0
             do j = 1,3
                dl = dl + (sys2_coord(j,1,ns2)-sys1_coord(j,1,i))**2
             end do
             if( sqrt(dl).lt.dist ) symbol = '+'
          end do
      end if

****  Thats all
      return
      end


CTITLE AV_FRAME
 
      subroutine av_frame

      implicit none
 
*     This routine will average the velocities at points separated
*     by less than av_dist
      
      include '../includes/const_param.h'
      include 'velrot.h'
 
* LOCAL VARIABLES
 
*   i,j,k   - Loop counters
*   ierr    - IOSTAT error
*   trimlen - Length of string
*   run_time(7) - Time that program was run
 
      integer*4 i,j,k
 
*   dx(3), dn(3)    - Changes in the position in XYZ and NEU
*   loc_coord(3)    - Local coordinates of sites
*   rot_matrix(3,3) - Rotation between XYZ and NEU
*   neu_part(3,7), xyz_part(3,7)    - Partials derivatives
*               - wrt to transformation parameters.
*   sectag      - Seconds part of time
 
 
      real*8 dn(3), loc_coord(3), rot_matrix(3,3)

      real*8 sum_av(3), sum_var(3), dvar(3), blen
      real*8 sum_rho, rho
      integer*4 num_av

      logical kbit

****  Now loop over all sites seeing which ones match the
*     av_dist requirement and taking the averages. 
   
      do j = 1, max(num_sys1,num_sys2)/32+1
         av_sys1(j) = 0
         av_sys2(j) = 0
      end do 
 
      do i = 1, num_sys1
         num_av = 1
         call rotate_geod(sys1_coord(1,2,i), dn,  'XYZ', 'NEU', 
     .                     sys1_coord(1,1,i), loc_coord, rot_matrix )

         do j = 1, 3
            sum_av(j) = dn(j)/sys1_cov(j,j,i)
            sum_var(j) = 1.d0/sys1_cov(j,j,i)
         end do
         sum_rho = sys1_cov(1,2,i)/
     .             sqrt(sys1_cov(1,1,i)*sys1_cov(2,2,i))
         do j = 1, max(num_sys1,num_sys2)/32+1
            av_u1(j) = 0
            av_u2(j) = 0
         end do 
         call sbit(av_u1,i,1)

*        Loop over all remaining sites that have not
*        used and sum up for the average
         do k = i+1, num_sys1
            if( blen(sys1_coord(1,1,i),sys1_coord(1,1,k)).lt.av_dist 
     .          .and. .not.kbit(av_sys1,k) ) then
                num_av = num_av + 1
                call rotate_geod(sys1_coord(1,2,k), dn,  'XYZ', 'NEU', 
     .                     sys1_coord(1,1,k), loc_coord, rot_matrix )

               call sbit(av_u1,k,1)
               do j = 1,3
                   sum_av(j) = sum_av(j) + 
     .                         dn(j)/sys1_cov(j,j,k)
                   sum_var(j) = sum_var(j) + 1.d0/sys1_cov(j,j,k) 
                end do
                sum_rho = sum_rho + sys1_cov(1,2,k)/
     .             sqrt(sys1_cov(1,1,k)*sys1_cov(2,2,k))
            end if
         end do

*        Now repeat for sys 2 sites
         do k = 1, num_sys2
            if( blen(sys1_coord(1,1,i),sys2_coord(1,1,k)).lt.av_dist 
     .          .and. .not.kbit(av_sys2,k) ) then
                num_av = num_av + 1
                call rotate_geod(sys2_coord(1,2,k), dn,  'XYZ', 'NEU', 
     .                     sys2_coord(1,1,k), loc_coord, rot_matrix )
                call sbit(av_u2,k,1)
                do j = 1,3
                   sum_av(j) = sum_av(j) + 
     .                         dn(j)/sys2_cov(j,j,k)
                   sum_var(j) = sum_var(j) + 1.d0/sys2_cov(j,j,k) 
                end do
                sum_rho = sum_rho + sys2_cov(1,2,k)/
     .             sqrt(sys2_cov(1,1,k)*sys2_cov(2,2,k))
            end if
         end do

*****    OK, see if we have any thing to average
         if( num_av.gt.1 ) then
             do j = 1, 3
                dn(j) = sum_av(j)/sum_var(j)
                dvar(j) = 1.d0/sum_var(j)
             end do
             
             rho = sum_rho/num_av
                    
*            Map these back to XYZ and update the coordinates
             do k = 1, num_sys1
                if( kbit(av_u1,k) ) then
                    call rotate_geod(dn, sys1_coord(1,2,k),  
     .                     'NEU', 'XYZ', sys1_coord(1,1,k),
     .                     loc_coord, rot_matrix)
                    do j = 1, 3
                        sys1_cov(j,j,k) = dvar(j)
                    end do
                    sys1_cov(1,2,k) = rho*sqrt(dvar(1)*dvar(2))
                    call sbit(av_sys1,k,1)
c                   write(*,120) k, sys1_names(k), num_av, dn
 120                format('AV Sys 1 ',i5,1x,a8,1x,i5,1x,3f10.4)
                end if
             end do
*            Map these back to XYZ and update the coordinates
             do k = 1, num_sys2
                if( kbit(av_u2,k) ) then
                    call rotate_geod(dn, sys2_coord(1,2,k),  
     .                     'NEU', 'XYZ', sys2_coord(1,1,k),
     .                     loc_coord, rot_matrix)
                    do j = 1, 3
                        sys2_cov(j,j,k) = dvar(j)
                    end do
                    sys2_cov(1,2,k) = rho*sqrt(dvar(1)*dvar(2))
                    call sbit(av_sys2,k,1)
c                   write(*,140) k, sys2_names(k), num_av, dn
 140                format('AV Sys 2 ',i5,1x,a8,1x,i5,1x,3f10.4)
 
               end if
             end do
         end if
      end do   

*     Now loop over the sites in system 2.
      do i = 1, num_sys2
         if( .not.kbit(av_sys2,i) ) then 
            num_av = 1
            call rotate_geod(sys2_coord(1,2,i), dn,  'XYZ', 'NEU', 
     .                        sys2_coord(1,1,i), loc_coord, rot_matrix )

            do j = 1, 3
               sum_av(j) = dn(j)/sys2_cov(j,j,i)
               sum_var(j) = 1.d0/sys2_cov(j,j,i)
            end do
            sum_rho = sys2_cov(1,2,i)/
     .                sqrt(sys2_cov(1,1,i)*sys2_cov(2,2,i))

            do j = 1, num_sys2/32+1
               av_u2(j) = 0
            end do 
            call sbit(av_u2,i,1)

*           Loop over all remaining sites that have not
*           used and sum up for the average
            do k = i+1, num_sys2
               if( blen(sys2_coord(1,1,i),sys2_coord(1,1,k)).lt.av_dist 
     .             .and. .not.kbit(av_sys2,k) ) then
                   num_av = num_av + 1
                   call rotate_geod(sys2_coord(1,2,k), dn,  
     .                        'XYZ', 'NEU', 
     .                        sys2_coord(1,1,k), loc_coord, rot_matrix )
                   call sbit(av_u2,k,1)
                   do j = 1,3
                      sum_av(j) = sum_av(j) + 
     .                            dn(j)/sys2_cov(j,j,k)
                      sum_var(j) = sum_var(j) + 1.d0/sys2_cov(j,j,k) 
                   end do
                   sum_rho = sum_rho + sys2_cov(1,2,k)/
     .                sqrt(sys2_cov(1,1,k)*sys2_cov(2,2,k))
               end if
            end do

*****       OK, see if we have any thing to average
            if( num_av.gt.1 ) then
                do j = 1, 3
                   dn(j) = sum_av(j)/sum_var(j)
                   dvar(j) = 1.d0/sum_var(j)
                end do
                rho = sum_rho/num_av
                       
*               Map these back to XYZ and update the coordinates
                do k = 1, num_sys2
                   if( kbit(av_u2,k) ) then
                       call rotate_geod(dn, sys2_coord(1,2,k),  
     .                        'NEU', 'XYZ', sys2_coord(1,1,k),
     .                        loc_coord, rot_matrix)
                       do j = 1, 3
                           sys2_cov(j,j,k) = dvar(j)
                       end do
                       sys2_cov(1,2,k) = rho*sqrt(dvar(1)*dvar(2))
                       call sbit(av_sys2,k,1)
c                      write(*,220) k, sys2_names(k), num_av, dn
 220                   format('AV SYS 1 ',i5,1x,a8,1x,i5,1x,3f10.4)
                   end if
                end do
            end if
         end if
      end do 

 
****  Thats all
      return
      end
 
CTITLE BLEN

      real*8 function blen( x1, x2 )

      implicit none

*     Routine to return the length bewteen two sites
*
      real*8 x1(3), x2(3)

*     Straight calculation
      blen = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2 +
     .            (x1(3)-x2(3))**2 )

****  Thats all
      return
      end

