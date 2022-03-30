 
      program corcom

      implicit none 
 
*     Program for the comparison of site coordinates estimated by
*     difference systems.  The runstring of the program is:
*
*     % corcom <sys1> <frame1> <sys2> <frame2> <outname> <out_frame> \
*              <ties> <fundamental sites> <height wght> <scale>
*
*     where <sys 1> is the name of the file containing the sites
*                 and their coordinates and velocity in system 1
*           <frame 1> is the frame for this first system.  The name
*                 should match one of those in kf/gen_util/frame_to_fra.f
*           <sys 2> is the second (comparison) system.  Sys 1 will be
*                 rotated and translated into this frame.  Differences
*                 will be given between sys1 and sys2 in this frame.
*           <frame 2> is the frame for this second system.  Same
*                 choices as above.  
*           <outname> is name of file for output of frame (overwritten)
*           <out_frame> is the frame for the output field (as above)
*           <ties> if the file containing site ties.  The file is
*                 interpretted as using the second site name in the
*                 tie to generate the position of the first tie.  The
*                 ties are only applied to the sites in sys 1. (May
*                 be neglected in runstring---no ties will be used)
*           <fundamental sites> names of sites to be used for rotation
*                 and translation (names after ties have been applied)
*                 If not given then all sites are used.  ALL may be
*                 given as name and all common sites will be used.
*                 A '-' in front of name will stop it being used.
*           <height wght> Weight to be given to heights in the
*                 transformation determination.  If 1 then equal
*                 weight given.  If 0 then no weight is given
*                 (default is not to used heights)
*           <scale> indicates that scale should be estimated (Y will
*                 cause sacle to be estimated, unless height weight
*                 is zero, in which case scale can not be estimated.
*
* MOD TAH 920629: Added the EURA and NNR-NUVEL are reference frame choices
* MOD TAH 180718: Updated format to allow for fast motions of ice sites. 
*     Ice site velocities can exceed 100 m/yr.
*
      include 'corcom.h'
 
 
 
****  Read the runstring of the program
 
      write(*,100)
 100  format(/,' CORCOM: Coordinate system comparison program',/)
 
      call get_corcom_run
 
***** Read in system 1 coordinates
 
      call read_sys(sys1_file, sys1_coord, sys1_epoch, sys1_names,
     .            num_sys1, max_sites)
 
*     Now rotate the system into the output frame
      call frame_update(sys1_coord, sys1_epoch, num_sys1, max_sites,
     .                sys1_frame, out_frame, frame_epoch)
 
*     Now read in system 2 coordinates
 
      call read_sys(sys2_file, sys2_coord, sys2_epoch, sys2_names,
     .            num_sys2, max_sites)
 
*     Now rotate the system into the output frame
      call frame_update(sys2_coord, sys2_epoch, num_sys2, max_sites,
     .                sys2_frame, out_frame, frame_epoch)
 
*     Read in the ties
 
      call read_ties
 
***** Now map system 1 names to the system 2 names using the ties.
*     (NOTE: In general this will increase the number of sites in
*            system 1)
 
      call apply_ties
 
*     Now read in the sites to be used in the transformation
      call read_fund_sites
 
*     Now transform system 1 to system 2.
 
      call transframe
 
****  Finally write out the system 1 and system 2 frames transformed
*     into common system.
 
      call output_frame
 
***** Thats all
      end
 
CTITLE GET_CORCOM_RUN
 
      subroutine get_corcom_run

      implicit none 
 
*     Routine to read the runstring of the corcom program. (see
*     main for order)
 
      include 'corcom.h'
 
* LOCAL VARIABLES
 
*   len_run - Length of runstring parameter
*   rcpar   - Read runstring
*   trimlen - Length of string
*   ierr    - IOSTAT error.
*   indx    - Position in string
 
      integer*4 len_run, rcpar, trimlen, ierr, indx
 
*   decyrs  - Deciminal years for output frame.
 
      real*8 decyrs
 
*   runstring   - Entry from runstring
 
      character*32 runstring
 
****  Start be getting the sys 1 file name
      len_run = rcpar(1, sys1_file )
      if( len_run.le.0 ) then
          call proper_runstring('corcom.hlp', 'corcom/sys1 file',1)
      end if
 
****  Get frame for sys 1
      len_run = rcpar(2, sys1_frame)
      if( len_run.le.0 ) then
          call proper_runstring('corcom.hlp', 'corcom/sys1 frame',1)
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
              call proper_runstring('corcom.hlp','corcom/sys2 frame',1)
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
      if( len_run.gt.0 ) then
          indx = index(runstring,'_')
	  indx = indx + index(runstring(indx+1:),'_')
          if( indx.gt.1 ) then
              out_frame = runstring(1:indx-1)
              read(runstring(indx+1:),*,iostat=ierr) decyrs
              call report_error('IOSTAT',ierr,'decod',
     .                runstring(indx+1:),0,'Out frame epoch')
              if( ierr.ne.0 ) then
                  call proper_runstring('corcom.hlp',
     .                'corcom/out_frame epoch',1)
              end if
          else
              write(*,160) runstring(1:len_run)
 160          format(/,'CORCOM: Output frame not in correct format',
     .                ' Entry was ',a)
                  call proper_runstring('corcom.hlp',
     .                        'corcom/out_frame epoch',1)
          end if
 
*     Nothing entered  so assume sys1 frame, and take 1990.0 as epoch
      else
          out_frame = sys1_frame
          decyrs = 1990.0
      end if
      if( decyrs.lt.500 ) decyrs = decyrs + 1900.d0
      frame_epoch = decyrs
 
****  Get the name of the ties files
      len_run = rcpar(7, ties_file )
      if( len_run.le.0 ) ties_file = ' '
 
****  Get name of file with fundamental sites
      len_run = rcpar(8, fund_file )
      if( len_run.le.0 ) fund_file = ' '
 
***** Get the weight of the heights
      len_run = rcpar(9, runstring )
      if( len_run.gt.0 ) then
          read(runstring,*,iostat=ierr ) height_weight
          call report_error('IOSTAT',ierr,'decod',runstring,
     .        0, 'corcom/Height Weight')
          if( ierr.ne.0 ) then
              call proper_runstring('corcom.hlp',
     .                        'corcom/Height Weight',1)
          end if
      else
 
*         Set default weight 0
          height_weight = 0.d0
      end if
 
****  Finally see if scale to be estimated
      len_run = rcpar(10,runstring)
      if( len_run.gt.0 ) then
          if( runstring(1:1).eq.'y' .or.
     .        runstring(1:1).eq.'Y' ) then
              num_parn = 7
          elseif( runstring(1:1).eq.'r' .or.
     .            runstring(1:1).eq.'R' ) then
              num_parn = 3   ! Only do rotation
          else
              num_parn = 6
          end if
      else
 
*         Nothing given so use the default of no scale.
          num_parn = 6
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_SYS
 
      subroutine read_sys(sys_file, sys_coord, sys_epoch,
     .            sys_names, num_sys, max_sites)

      implicit none 
 
*     This routine will read in a systems coordinates, velocities,
*     epochs and site names.  Any non-blank character is col 1 is
*     treated as a comment.  Any decoding errors will be igorned.
 
 
* PASSED VARIABLES
 
*   max_sites   - Max number of sites
*   num_sys     - Number of sites in this system
 
      integer*4 max_sites, num_sys
 
*   sys_coord(3,2,max_sites)    - Coordinates and velocities
*               - of sites (m and m/yr)
*   sys_epoch(max_sites)    - Epochs of the site positions
 
      real*8 sys_coord(3,2,max_sites), sys_epoch(max_sites)
 
*   sys_names(max_sites)    - Names of the sites
*   sys_file        - Name of file.
 
      character*(*) sys_names(max_sites), sys_file
 
* LOCAL VARIABLES
 
*   ierr, jerr  - IOSTAT errors
*   trimlen - Length of string.
*   i,j     - Loop counters
*   indx    - Position in string
 
      integer*4 ierr, jerr, trimlen, i, indx
 
*   decyrs  - Deciminal years for site epoch
*   values(6)   - Position and velocity read from file
 
      real*8 decyrs, values(6)
 
*   tname   - Temp name read from file
 
      character*8 tname, cd
 
*   line    - Line read from file
 
      character*256 line
 
****  Try to open the file (see if name in non-blank)
      num_sys = 0
      if( trimlen(sys_file).gt.0 ) then
          open(100, file=sys_file, status='old', iostat=ierr)
          call report_error('IOSTAT',ierr,'open', sys_file,1,
     .                    'corcom/read_sys')
 
****      OK, FIle open now start reading
          i = 0
          do while( ierr.eq.0 .and. i.lt.max_sites)
              read(100,'(a)',iostat=ierr) line
              if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .            trimlen(line).gt.0 .and. 
     .            index(line,'EXTENDED').eq.0 ) then
 
                  indx = 0
                  call getword(line, tname, indx)
                  call casefold(tname)
                  call multiread(line, indx, 'R8', jerr, values, cd ,6)
                  call read_line(line, indx, 'R8', jerr, decyrs, cd )
C                 read(line(indx:),*,iostat=jerr) values, decyrs
                  call report_error('IOSTAT',jerr,'read', line,1,
     .                'read_sys/Coordinates')
 
*                 See error occurred, if not added to site.  See also if
*                 X and Y coords are less than 60 (if so this is a source
*                 position.  (means also that we can't work at the poles)
                  if( abs(values(1)).lt.25 .and.
     .                abs(values(2)).lt.61       ) jerr = 1

                  if( jerr.eq.0 ) then
                      i = i + 1
                      if( i.eq. max_sites ) then
                         write(*,150) max_sites
 150                     format('**WARNING** Maximum number of sites ',
     .                          I5,' exceeding.  Terminating input',/,
     .                      '            Update max_sites in corcom.h')
                      endif
                      sys_names(i) = tname
                      sys_coord(1,1,i) = values(1)
                      sys_coord(2,1,i) = values(2)
                      sys_coord(3,1,i) = values(3)
                      sys_coord(1,2,i) = values(4)
                      sys_coord(2,2,i) = values(5)
                      sys_coord(3,2,i) = values(6)
                      sys_epoch(i)     = decyrs
*                             ! No decode error
                  end if
*                             ! Non-comment, and no error
              end if
*                             ! Looping over the input file
          end do
          
*****     Save the number of sites
          num_sys = i
          write(*,200) num_sys, sys_file(1:trimlen(sys_file))
 200      format(' There are ',i4,' sites in sys file ',a)
 
 
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
     .    abs(rot_vec(3)).lt.          1.d-20 ) RETURN
 
***** Loop over all of the sites and update position (if different
*     from frame_epoch) and velocity
 
      do i = 1, num_sys
          call cross_prod(rot_vec, sys_coord(1,1,i), xyz_vel, sang)
 
*         Update position and velocity
          do j = 1,3
              sys_coord(j,1,i) = sys_coord(j,1,i) +
     .                        xyz_vel(j)*(sys_epoch(i)-frame_epoch)
              sys_coord(j,2,i) = sys_coord(j,2,i) + xyz_vel(j)
          end do
      end do
 
****  Thats all
      return
      end
 
CTITLE READ_TIES
 
      subroutine read_ties

      implicit none 
 
*     Routine to read in the ties information from the ties file
 
      include 'corcom.h'
 
* LOCAL VARIABLES
 
*   ierr, jerr  - IOSTAT errors
*   i           - Loop counter
*   trimlen     - Length of string
 
      integer*4 ierr, jerr, i, trimlen, indx
 
*   values(3)   - Value of tie
 
      real*8 values(3)
 
*   tname(2)        - Names of sites in tie file
*   ttype           - Type of tie: if not NEU than assumed to be XYZ
*   cd              - Dummy string
 
      character*8 tname(2), ttype, cd

*   line            - Line read from file

      character*256 line

 
****  Try to open the file (see if name in non-blank)
      num_ties = 0
      if( trimlen(ties_file).gt.0 ) then
          open(100, file=ties_file, status='old', iostat=ierr)
          call report_error('IOSTAT',ierr,'open', ties_file,1,
     .                    'corcom/read_ties')
 
****      OK, FIle open now start reading
          i = 0
          do while( ierr.eq.0 .and. i.lt.max_sites)
              read(100,'(a)', iostat=ierr) line
              if( ierr .eq.0 .and. line(1:1).eq.' ' .and.
     .            trimlen(line).gt.0 ) then
 
                  indx = 0
                  call getword(line, tname(1), indx)
                  call casefold(tname(1))
                  call getword(line, tname(2), indx)
                  call casefold(tname(2))
C                 read(line(indx:),*,iostat=jerr) values
                  call multiread(line,indx, 'R8', jerr, values, cd,3)
                  call report_error('IOSTAT',jerr,'read', line,1,
     .                'read_ties/Coordinates')

*                 See if coordinate type given.  If not NEU then
*                 assume XYZ
                  call getword(line, ttype, indx)
 
*                 See error occurred, if not added to site
                  if( jerr.eq.0 ) then
                      i = i + 1
                      ties_names(1,i) = tname(1)
                      ties_names(2,i) = tname(2)
                      call casefold(ttype)
                      if( ttype(1:4).ne.'NEU ' ) ttype = 'XYZ'
                      ties_type(i)    = ttype
                      ties_coord(1,i) = values(1)
                      ties_coord(2,i) = values(2)
                      ties_coord(3,i) = values(3)
*                             ! No decode error
                  end if
*                             ! Non-comment, and no error
              end if
*                             ! Looping over the input file
          end do
 
*****     Save the number of sites
          num_ties = i
          write(*,200) num_ties, ties_file(1:trimlen(ties_file))
 200      format(' There are ',i4,' sites in tie file ',a)
 
 
*                             ! Name actually given
      end if
 
****  Thats all
      return
      end
 
CTITLE APPLY_TIES
 
      subroutine apply_ties

      implicit none 
 
*     Routine to apply the ties read from the ties file.  These
*     are only applied to the system 1 sites, and are applied
*     recurcively mapping the site with the second name to the
*     first name in the ties file lists
 
      include 'corcom.h'
 
* LOCAL VARIABLES
 
*   i,j,k       - Loop counters
 
 
      integer*4 i,j,k,l

*   loc_coord(3)        - Local coordinates of the site
*   rot_matrix(3,3)     - Rotation matrix from rn_types system
*                       - to XYZ
*   dxyz(3)             - XYZ coordinate shift (m)
 
      real*8 loc_coord(3), rot_matrix(3,3), dxyz(3)
 
****  Start looping over the sites in system 1.  When a ties is
*     found, add the site to the bottom of the list.  (Latter the
*     ties will be checked to see if there are any ties for it.)
 
      i = 0
 
      do while ( i.lt. num_sys1 .and. num_sys1.lt.max_sites )
 
*         Increment site counter
          i = i + 1
 
*         Scan down list of ties to see if ant match.
          do j = 1, num_ties
              if( sys1_names(i).eq.ties_names(2,j) ) then
 
*                 Make sure that this is not a null tie (i.e.,
*                 names the same)
                  if( ties_names(1,j).ne.ties_names(2,j) ) then
*                     Add tie site to list
                      num_sys1 = num_sys1 + 1
                      k = num_sys1
                      sys1_names(k) = ties_names(1,j)
                      call rotate_geod(ties_coord(1,j), dxyz,
     .                     ties_type(j), 'XYZ', sys1_coord(1,1,i), 
     .                     loc_coord, rot_matrix)
 
                      do l = 1,3
                          sys1_coord(l,1,k) = sys1_coord(l,1,i) +
     .                                       dxyz(l)
                          sys1_coord(l,2,k) = sys1_coord(l,2,i)
                      end do
                      sys1_epoch(k) = sys1_epoch(i)
*                                 ! Tie sites have different name
                  end if
*                                 ! Sys1 site name match ties site
              end if
*                                 ! Looping over ties
          end do
*                                 ! Looping over sites in sys 1.
      end do
 
****  Report to user.
      write(*,200) num_sys1
 200  format(' After applying ties, there are ',i4,' sites in',
     .       ' system 1')
 
      return
      end
 
CTITLE READ_FUND_SITES
 
      subroutine read_fund_sites

      implicit none 
 
*     This routine will read the list of sites to be used in the
*     transformation between sys 1 and sys 2.  If no file name is
*     given then all possible sites will be used.
 
      include 'corcom.h'
 
* LOCAL VARIABLES
 
*   ierr, jerr  - IOSTAT errors
*   i,j         - Loop counter
*   ns1, ns2        - Matching sites from sys1 and sys2 lists
*   trimlen     - Length of string
 
      integer*4 ierr, i,j, k, ns1, ns2, trimlen, indx
 
*   tname       - Name of fundamental site read from file
 
      character*10 tname

*   line        - Line read from input file

      character*256 line

*   use_site 
      logical use_site  ! Set true if site to be used (- in front of
                        ! name causes non-use.
     .,       found     ! Set true if site found in fund_list
 
****  Try to open the file (see if name in non-blank)
      num_fund = 0
      if( trimlen(fund_file).gt.0 ) then
          open(100, file=fund_file, status='old', iostat=ierr)
          call report_error('IOSTAT',ierr,'open', fund_file,1,
     .                    'corcom/read_fund')
 
****      OK, FIle open now start reading
          i = 0
          do while( ierr.eq.0 .and. i.lt.max_sites)
              read(100,'(a)', iostat=ierr) line
              if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .            trimlen(line).gt.0 ) then
 
                  indx = 0
                  call getword(line, tname, indx)
                  call casefold(tname)
*                 See if we are to use site 
                  if( tname.eq.'ALL' ) then
*                     Assign all matching sites
                      i = 0
                      do j = 1, num_sys1
                          indx=0
                          call get_cmd(sys1_names(j), sys2_names,
     .                                  -num_sys2, ns2, indx)
                          if( ns2.gt.0 ) then
                              i = i + 1
                              fund_links(1,i) = j
                              fund_links(2,i) = ns2
                          end if
                      end do
                  else

*****                 Remove leading character
                      use_site = .true.
                      if( tname(1:1).eq.'-' ) then
                         use_site = .false.
                         tname = tname(2:9)
                      else if( tname(1:1).eq.'"' ) then
                         tname = tname(2:9)
                      endif
 
*                     See if we can match sites
                      indx = 1
                      call get_cmd(tname, sys1_names,-num_sys1,
     .                            ns1, indx)
                      indx = 1
                      call get_cmd(tname, sys2_names,-num_sys2,
     .                            ns2, indx)
 
****                  See if we match sites sites
                      if( ns1.gt.0 .and. ns2.gt.0 ) then
*                         See if already in list
                          found = .false.
                          do k = 1,i
                             if( fund_links(1,k).eq. ns1 .and.
     .                           fund_links(2,k).eq. ns2 ) then
                                 found = .true.
                                 exit
                             endif
                          end do
                          if( found ) then
*                            See what we should do (if use_site false
*                            then remove from list)
                             if( .not. use_site ) then
                                do j = k+1, i
                                   fund_links(1,j-1) = fund_links(1,j)
                                   fund_links(2,j-1) = fund_links(2,j)
                                enddo
                                i = i -1
                             endif
                          else
*                            Add site to list                                   
                             i = i + 1
                             fund_links(1,i) = ns1
                             fund_links(2,i) = ns2
                          endif
                      else
                          write(*,150) tname, ns1, ns2
 150                      format(' No match for ',a8,' in fundamental',
     .                        ' sites.  Returns from sys1 and sys2 are ',
     .                        2I4)
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
                  fund_links(1,i) = j
                  fund_links(2,i) = ns2
              end if
          end do
      end if
 
***** Save the number of sites
      num_fund = i
      write(*,200) num_fund, fund_file(1:trimlen(fund_file))
 200  format(' There are ',i4,' matching sites in fundamental',
     .        ' file ',a)
 
      return
      end
 
CTITLE TRANSFRAME
 
      subroutine transframe

      implicit none 
 
*     Routine to compute the transformation parameters between
*     system 1 and system 2.  The loops over all of the fundamental
*     sites and computes the transformation from them.
 
      include '../includes/const_param.h'
      include 'corcom.h'
 
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
     .    dx(3), dn(3), loc_coord(3), scale(7)
 
****  Initialize the system
 
      call clear_norm
 
****  Loop over the fundamental sites
      do i = 1, num_fund
 
*         get the transformation from XYZ to NEU for this site
          ns1 = fund_links(1,i)
          ns2 = fund_links(2,i)
 
*         Get coordinate differences
          do j = 1,3
              dx(j) = sys2_coord(j,1,ns2) - sys1_coord(j,1,ns1) +
     .            sys2_coord(j,2,ns2)*(sys1_epoch(ns1)-sys2_epoch(ns2))
          end do
          call rotate_geod(dx, dn,  'XYZ', 'NEU', sys1_coord(1,1,ns1),
     .                    loc_coord, rot_matrix )
 
*****     Now generate the partials with respect to rotation,scale,and
*         translation.
          call get_parts(sys1_coord(1,1,ns1), rot_matrix, num_parn,
     .                neu_part, xyz_part, max_sites)
 
*****     Now increment normal equations
          call increment_norm( dn, num_parn, neu_part, height_weight,
     .            norm_eq, bvec, sum_prefit, num_data )
 
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
          rms_fit = sqrt(sum_postfit/(num_data-num_parn))
      else
          rms_fit = 1.d0
      end if
 
*     Now convert transformation parameters to output format
      trans_parm_out(1) = trans_parm(1) * rad_to_mas
      trans_parm_out(2) = trans_parm(2) * rad_to_mas
      trans_parm_out(3) = trans_parm(3) * rad_to_mas
      trans_parm_out(4) = trans_parm(4) * 1000.d0
      trans_parm_out(5) = trans_parm(5) * 1000.d0
      trans_parm_out(6) = trans_parm(6) * 1000.d0
      trans_parm_out(7) = trans_parm(7) * 1.d9
 
*     Now do sigmas
      trans_sigma(1) = sqrt(norm_eq(1,1)) * rad_to_mas * rms_fit
      trans_sigma(2) = sqrt(norm_eq(2,2)) * rad_to_mas * rms_fit
      trans_sigma(3) = sqrt(norm_eq(3,3)) * rad_to_mas * rms_fit
      trans_sigma(4) = sqrt(norm_eq(4,4)) * 1000.d0 * rms_fit
      trans_sigma(5) = sqrt(norm_eq(5,5)) * 1000.d0 * rms_fit
      trans_sigma(6) = sqrt(norm_eq(6,6)) * 1000.d0 * rms_fit
      trans_sigma(7) = sqrt(norm_eq(7,7)) * 1.d9 * rms_fit
 
*     Convert rms_fit to mm
      rms_fit = rms_fit * 1000.d0
 
****  Thats all
      return
      end 
 
CTITLE CLEAR_NORM
 
      subroutine clear_norm

      implicit none 
 
*     This routine will initialize the estimation arrays.
 
      include 'corcom.h'
 
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
 
*     Initialize the statistics variables
 
      sum_prefit = 0.d0
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
      xyz_part(1,1) = -coord(3)
      xyz_part(3,1) =  coord(1)
      xyz_part(2,2) =  coord(3)
      xyz_part(3,2) = -coord(2)
      xyz_part(1,3) =  coord(2)
      xyz_part(2,3) = -coord(1)
 
*     Translation
      xyz_part(1,4) = 1.d0
      xyz_part(2,5) = 1.d0
      xyz_part(3,6) = 1.d0
 
*     See if new scale partials
      if( num_parn.eq.7 ) then
          xyz_part(1,7) = coord(1)
          xyz_part(2,7) = coord(2)
          xyz_part(3,7) = coord(3)
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
      end do
 
****  Thats all
      return
      end
 
CTITLE INCREMENT_NORM
 
      subroutine increment_norm( dn, num_parn, neu_part, height_weight,
     .            norm_eq, bvec, sum_prefit, num_data )

      implicit none 
 
*     Routine to increment normal equations and solution vector.
*     Prefit residuals are also accumulated for later computation
*     of the fit to the data.
 
 
* PASSED VARIABLES
 
*   num_parn        - Number of parameters being estimated
*   num_data        - Number of components summed into prefit
 
      integer*4 num_parn, num_data
 
*   dn(3)       - Difference in NE and U
*   height_weight   - Height weight (if set to zero heights
*               - are not used
*   neu_part(3,7)   - Partials of NEU wrt transformation
*               - parameters
*   norm_eq(7,7)    - Normal equations
*   bvec(7)     - Solution vector
*   sum_prefit  - Sum of prefit residuals
 
      real*8 dn(3), height_weight, neu_part(3,7), norm_eq(7,7),
     .    bvec(7), sum_prefit
 
* LOCAL VARIABLES
 
*   i,j,k       - Loop counters
*   num_use     - Number of components to use
 
      integer*4 i,j,k, num_use
 
*   wghs(3)     - Weights to be given to each component
 
 
      real*8 wghs(3)
 
****  See what weights we should give (N and E always get 1)
      wghs(1) = 1.d0
      wghs(2) = 1.d0
      wghs(3) = height_weight
 
      if( abs(wghs(3)).lt.1.d-6 ) then
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
          num_data   = num_data + 1
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
 
      include 'corcom.h'
 
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
     .    neu_part(3,7), xyz_part(3,7), sectag
 
*   unit_labs(7)        - Units for parameters
*   param_labs(7)   - Parameter names
 
      character*8 unit_labs(7), param_labs(7)
 
      data param_labs / 'X-Pole', 'Y-Pole', 'AT-UT1',
     .                'X-Offset', 'Y-Offset', 'Z-Offset', 'Scale' /
 
      data unit_labs / '(mas)','(mas)','(mas)',
     .                '(mm)', '(mm)', '(mm)', '(ppb)' /
 
 
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
 
      write(out_unit,100 ) (run_time(i),i=1,5),
     .    sys1_frame, sys1_file(1:max(1,trimlen(sys1_file))),
     .    sys2_frame, sys2_file(1:max(1,trimlen(sys2_file))),
     .    out_frame,
     .    ties_file(1:max(1,trimlen(ties_file))),
     .    fund_file(1:max(1,trimlen(fund_file))),
     .    out_frame, frame_epoch
 
 100  format('* CORCOM Run on ',i4,2('/',i2),1x,i2,':',i2,/,
     .    '* SYSTEM 1 File    : ',a,1x,a,/,
     .    '* SYSTEM 2 File    : ',a,1x,a,/,
     .    '+REFERENCE FRAME     ',a,/,
     .    '* TIES File        : ',9x,a,/,
     .    '* FUNDAMENTAL File : ',9x,a,/,
     .    '* ',/,
     .    '* OUTPUT FRAME     : ',a,t40,'REF. EPOCH ',f9.4  )
 
*     Now write out the transformation parameters
      if( num_data.gt.0 ) then
          write(out_unit,200)  num_data, num_fund, rms_fit,
     .        (i, param_labs(i),trans_parm_out(i),
     .        trans_sigma(i), unit_labs(i), i=1, num_parn)
 
200       format('*',/,
     .        '* RMS fit of transformation using ',i5,' components',
     .        ' from ',i4,' stations was ', F10.2, ' mm',/,
     .        '* Estimates of Transformation parameters are: ',
     .        7(/,'*',i4,1x,a,3x,F10.4,1x,F10.4,1x,a))
      end if
 
*     Now print out the errors at the Fundamental stations
 
      write(out_unit,300)
 300  format('*  Differences are the fundamental sites',/,
     .    '*   #    Name',4x,'dN (mm)',4x,'dE (mm)',4x,'dU (mm) ',
     .    '   dX (mm)',4x,'dY (mm)',4x,'dZ (mm)')
 
      do i = 1, num_fund
          ns1 = fund_links(1,i)
          ns2 = fund_links(2,i)
 
*         Get coordinate differences (before adjustment)
          do j = 1,3
              dx(j) = sys2_coord(j,1,ns2) - sys1_coord(j,1,ns1) +
     .            sys2_coord(j,2,ns2)*(sys1_epoch(ns1)-sys2_epoch(ns2))
          end do
 
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
 
*         Now write the values
          write(out_unit,320) i, sys1_names(ns1),
     .            (dn(j)*1000.d0, j = 1,3),
     .            (dx(j)*1000.d0, j = 1,3)
  320     format('*',i4,1x,a8,1x,6(F10.2,1x))
 
      end do
 
****  Now loop over all of the system 1 sites appying the
*     transformation.
 
      write(out_unit,400)
  400 format('*',/,'* SYSTEM 1 Coordinates in SYSTEM 2 Frame',/,
     .       '*  Name',11x,'X (m)', 9x,'Y (m)', 9x,'Z (m)',
     .       6x,'Xdot',6x,'Ydot',6x,'Zdot', 6x,'Epoch')
 
      do i = 1, num_sys1
 
****      Get the correction to the position due to the system
*         change.  (If we want NEU we should call routine below
*         but we just want XYZ so don't bother.
*         call rotate_geod(dx, dn,  'XYZ', 'NEU', sys1_coord(1,1,ns1),
*    .                    loc_coord, rot_matrix )
 
          call get_parts(sys1_coord(1,1,i), rot_matrix, num_parn,
     .                neu_part, xyz_part, max_sites)
 
*         Get adjustment to position
          do j = 1, 3
              dx(j) = 0.0d0
              do k = 1,num_parn
                  dx(j) = dx(j) + xyz_part(j,k)*trans_parm(k)
              end do
 
*             Now add to position
              sys1_coord(j,1,i) = sys1_coord(j,1,i) + dx(j)
          end do
 
*         Write out the results
          write(out_unit,420) sys1_names(i), (sys1_coord(j,1,i),j=1,3),
     .                     (sys1_coord(j,2,i),j=1,3), sys1_epoch(i)
* MOD TAH 180718: Updated format to 3(F10.5,1x)
 420      format(1x,a8,3(1x,F14.5),1x,3(F10.5,1x),F9.4)
      end do
 
***** Now output the system2 coordinates (comment out those which
*     are are already in the list.  NOTE: System 1 coordinates have
*     now been updated to system 2 frame.
 
      write(out_unit,500)
  500 format('*',/,'* SYSTEM 2 Coordinates except those in',
     .            ' SYSTEM 1 ',/,
     .       '*  Name',10x,'X (m)', 8x,'Y (m)', 8x,'Z (m)',
     .       5x,'Xdot',5x,'Ydot',5x,'Zdot',
     .       2x,'Epoch',2x,'dN (mm)',2x,'dE (mm)',2x,'dU (mm)')
 
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
 
*         Get coordinate differences (before adjustment)
          if( ns1.gt.0 ) then
              do j = 1,3
                  dx(j) = sys2_coord(j,1,i) - sys1_coord(j,1,ns1) -
     .                sys1_coord(j,2,ns1)*(sys2_epoch(i)-
     .                sys1_epoch(ns1))
              end do
 
*             Now rotate into local system
              call rotate_geod(dx, dn,  'XYZ', 'NEU',
     .                sys1_coord(1,1,ns1),loc_coord, rot_matrix )
 
*             write out the line
              write(out_unit,520) 'X',sys2_names(i),
     .                     (sys2_coord(j,1,i),j=1,3),
     .                     (sys2_coord(j,2,i),j=1,3), sys2_epoch(i),
     .                    (dn(j)*1000.d0,j=1,3)
* MOD TAH 180718: Updated format from 3(F8.5,1x) to 3(f10.5,1x)
 520          format(a1,a8,3(1x,F13.4),1x,3(F8.5,1x),F9.4,3(1x,F9.1))
          else
              write(out_unit,520) ' ',sys2_names(i),
     .                     (sys2_coord(j,1,i),j=1,3),
     .                     (sys2_coord(j,2,i),j=1,3), sys2_epoch(i)
 
          end if
      end do
 
****  Thats all
      return
      end
 
