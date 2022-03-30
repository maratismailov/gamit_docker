CTITLE READ_GEOC
 
      subroutine read_geoc( unit, line, np )
 
      implicit none

*     This routine will read the geodetic coordinates from the
*     hfile and save the scaling and parameter number information.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
*   np          - Running count on parameters
 
      integer*4 unit, np
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   indx        - Pointer in string
*   sn          - Site number
*   i           - Loop counter
*   nt          - Temporary parameter number read from file
 
      integer*4 ierr, indx, sn, i, nt, trimlen
 
*   apr, adj        - Apriori value and adjustment to parameter
 
      real*8 apr, adj
 
*   sr          - Site name read
 
      character*8 sr
 
*   est         - if * then estimated
 
 
      character*1 est
 
      read(line, 100, iostat=ierr) nt, est, sr, apr, adj
  100 format(i5,a1,a4,16x, f25.17, d30.17)
 
*     Find the site name
      indx = 1

      call get_cmd(sr, qsite_names, qnum_sites, sn,indx )
*                         ! Save values
      if( sn.gt.0 ) then
          llr(1,sn) = apr
 
*                                     ! Parameter estimated
          if( est.eq.'*' ) then
              np = np + 1
              qparn_sites(1,sn) = np
              qsol(np) = adj
              atoo(np) = nt
          end if
 
      else
          if( est.eq.'*' ) then
              write(*,150) line(1:trimlen(line))
  150         format(' ** ERROR *** decoding ',a,/,4x,'Could not find',
     .              ' site name',/,
     .               ' ****BINARY HFILE is Likely to be corrupt and ',
     .               ' not usuable' )
          end if

*         Must be no data on this site.  Continue and put the
*         site number as max_glb_sites.  Do nothing at the moment
C         sn = max_glb_sites
 
      end if
 
*     Now read the next two lines (with long and radius
      do i = 1, 2
          read(unit,'(a)', iostat=ierr ) line
          read(line, 100, iostat=ierr) nt, est, sr,
     .                    apr, adj
 
*         Save values
*                             ! Save values
          if( sn.gt.0 ) then
              llr(i+1,sn) = apr
              if( est.eq.'*' ) then
                  np = np + 1
                  qparn_sites(i+1,sn) = np
                  atoo(np) = nt
                  qsol(np) = adj
              end if
          end if
*                 ! Looping on longitude and radius
      end do
 
***** Now set up the scaling.  Convert lat, long, and radius to meter in
*     the N,E and U directions.
      if( sn.gt.0 ) then
          qscale(np-2) = llr(3,sn) * 1000.d0
          qscale(np-1) = llr(3,sn) * cos( llr(1,sn) )*1000.d0
          qscale(np)   = 1000.d0
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_ORBIT
 
      subroutine read_orbit ( unit, line, np, line_preread )
 
      implicit none

*     Routine to read the orbit block of parameters.  Line is passed
*     with the first entry (X position) and the remaining 8 lines are
*     read here.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
*   np          - Running count on parameters
 
      integer*4 unit, np

*   line_preread - Logical to indicate that Satellite X position
*                  line has been pre-read because we we searching
*                  for the radition parameters

      logical line_preread
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   vn          - Satellite vechile number
*   i           - Loop counter
*   nt          - Temporary parameter number read from file
*   svcode      - The svs element code number (see svel_to_code)
*   el_read     - Number of elements to read after the X satellite
*                 position
 
      integer*4 ierr, vn, i, nt, svcode, el_read  
 
*   apr, adj        - Apriori value and adjustment to parameter
*   scsv        - Scale from km to m and km/s to mm/s for the first
*                 six elements and 1 for the radiation parameters
*   scale_sv(max_svs_elem) - Scaling factors to use for sv orbits
 
      real*8 apr, adj, scsv, scale_sv(max_svs_elem)
 
*   est         - if * then estimated

      character*1 est

*   svp         - Character version of SV PRN number
 
      character*2 svp

*   svelem     - Name of the satellite element read
*   gkelem     - Dummy name for globk output
*   arelem     - Dummy name for arc output

      character*16 svelem
      character*14 gkelem
      character*4  arelem

      data scale_sv / 1000.d0, 1000.d0, 1000.d0,
     .                1000.d3, 1000.d3, 1000.d3,
     .                   1.d0,    1.d0,    1.d0,  
     .                   1.d0,    1.d0,    1.d0,  
     .                   1.d0,    1.d0,    1.d0,  
     .                   1.d0,    1.d0,    1.d0,
     .                   1.d0,    1.d0,    1.d0,
     .                   1.d0,    1.d0           /

 
***** Start decoding the line which we just read.  Then do the 8
*     remaining values

*     See if we know number of radiation parameters
      line_preread = .false.
      if( .not.rad_known ) then
          rad_num = 0
          el_read = 5
      else
*         Add 5 to number of radiation parameter to include
*         Y, Z position and XYZ velocities.
          el_read = rad_num + 5
      end if
 
      read(line, 200, iostat=ierr) nt, est,svelem,  svp, apr, adj
  200 format(i5,a1,a16,2x,a2, f25.16, d30.17)
 
*     Find the PRN number matched against the SVS name
      vn = -1
 
*                         ! change leading blank to zero
      if( svp(1:1).eq.' ' ) svp(1:1) = '0'
      do i = 1, qnum_svs
* MOD TAH 180401: Name format for GNSS changed: Now 6:7 not 5:6 
          if( svp.eq. qsvs_names(i)(6:7) ) then
              vn = i
          end if
      end do
 
*                         ! Save values
      if( vn.gt.0 ) then
          orb_est = .true.   ! Set to show orbits estimated,
          svs_pos(1,vn) = apr*scale_sv(1)
 
          if( est.eq.'*' ) then
              np = np + 1
              qparn_svs(1,vn) = np
              qsol(np) = adj
              atoo(np) = nt
              qscale(np) = scale_sv(1)
          end if
 
      else

*         No need to tell user any more.  (If no data then satellites
*         are not used and this message would be invoked.)
c         write(*,250) svp, line(1:trimlen(line))
c 250     format(' ** ERROR ** Could not find PRN ',a2,' from line:'
c    .            ,/,a)
*                                 ! Do nothing, but sure we check that
*                                 ! there is no satellite later. 
C         vn = max_glb_svs
 
      end if
 
*     Now read the next five lines which are the posiiton and velocity.
*     Then start decoding the radiation parameter models.
      do i = 1, el_read
          read(unit,'(a)', iostat=ierr ) line
          read(line, 200, iostat=ierr) nt, est,svelem,  svp, apr, adj
 
*         Save values
*                             ! Save values
          if( vn.gt.0 ) then

*             If i is less than of equal to 5, then these are  standard
*             elements other wize find out what it is
              if( i.le.5 ) then
*                 Rescale the apriori value
                  scsv = scale_sv(i+1)

                  svs_pos(i+1,vn) = apr * scsv
 
                  if( est.eq.'*' ) then
                      np = np + 1
                      qparn_svs(i+1,vn) = np
                      qsol(np) = adj
                      atoo(np) = nt
                      qscale(np) = scsv 
                  end if
              else
                  call svel_to_code(svelem, gkelem, arelem, 
     .                              svcode, 'NTOC')
                  if( svcode.ge.7 .and. svcode.le.max_svs_elem) then
*                     Rescale the apriori value
                      scsv = scale_sv(svcode)
                      svs_pos(svcode,vn) = apr * scsv
                      if( est.eq.'*' ) then
                          np = np + 1
                          qparn_svs(svcode,vn) = np
                          qsol(np) = adj
                          atoo(np) = nt
                          qscale(np) = scsv 
                      end if
                  end if
              end if
          end if
*                     ! Looping over orbital elements
      end do

****  Now see if we should search for the radiation parameter
*     model.
      if( .not.rad_known ) then

*         Start read
          line_preread = .false.
          do while ( .not.line_preread )
              read(unit,'(a)', iostat=ierr ) line
              read(line, 200, iostat=ierr) nt, est,svelem, svp, apr, adj

*             Get the parameter code
              call svel_to_code(svelem, gkelem, arelem, svcode, 'NTOC')
              if( svcode.ge.7 .and. svcode.le.max_svs_elem) then

*                 Sum number of radiation parameters
                  rad_num = rad_num + 1
*                 Rescale the apriori value
                  scsv = scale_sv(svcode)
                  svs_pos(svcode,vn) = apr * scsv
                  if( est.eq.'*' ) then
                      np = np + 1
                      qparn_svs(svcode,vn) = np
                      qsol(np) = adj
                      atoo(np) = nt
                     qscale(np) = scsv 
                  end if
              else
*                 We can not decode so must be next sattellite
                  line_preread = .true.
                  rad_known = .true.
                  write(*,300) rad_num
 300              format(i4,' Radiations parameter types found')
              end if
          end do
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_XYZ
 
      subroutine read_xyz( unit, line, np )
 
      implicit none

*     This routine will read the cartesian coordinates from the
*     hfile and save the scaling and parameter number information.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
*   np          - Running count on parameters
 
      integer*4 unit, np
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   indx        - Pointer in string
*   sn          - Site number
*   i           - Loop counter
*   nt          - Temporary parameter number read from file
 
      integer*4 ierr, indx, sn, i, nt
 
*   apr, adj        - Apriori value and adjustment to parameter
 
      real*8 apr, adj
 
*   sr          - Site name read
 
      character*8 sr
 
*   est         - if * then estimated
 
 
      character*1 est
 
      read(line, 100, iostat=ierr) nt, est, sr, apr, adj
  100 format(i5,a1,a4,18x, 2d30.17)
 
*     Find the site name
      indx = 1
      call get_cmd(sr, qsite_names, qnum_sites, sn,indx )
*                         ! Save values
      if( sn.gt.0 ) then
          site_pos(1,sn) = apr
 
*                                     ! Parameter estimated
          if( est.eq.'*' ) then
              np = np + 1
              qparn_sites(1,sn) = np
              atoo(np) = nt
              qsol(np) = adj
          end if
 
      end if
 
*     Now read the next two lines (with long and radius
      do i = 1, 2
          read(unit,'(a)', iostat=ierr ) line
          read(line, 100, iostat=ierr) nt, est, sr,
     .                    apr, adj
 
*         Save values
*                             ! Save values
          if( sn.gt.0 ) then
              site_pos(i+1,sn) = apr
              if( est.eq.'*' ) then
                  np = np + 1
                  qparn_sites(i+1,sn) = np
                  atoo(np) = nt
                  qsol(np) = adj
              end if
          end if
*                 ! Looping on longitude and radius
      end do
 
***** Now set up the scaling.  There is not change of units so set 1.00 
      if( sn.gt.0 ) then
          qscale(np-2) = 1.d0
          qscale(np-1) = 1.d0
          qscale(np)   = 1.d0 
      end if
 
****  Thats all
      return
      end
 
 
CTITLE READ_VEL
 
      subroutine read_vel( unit, line, np )
 
      implicit none

*     This routine will read the cartesian velocities 	from the
*     hfile and save the scaling and parameter number information.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
*   np          - Running count on parameters
 
      integer*4 unit, np
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   indx        - Pointer in string
*   sn          - Site number
*   i           - Loop counter
*   nt          - Temporary parameter number read from file
 
      integer*4 ierr, indx, sn, i, nt
 
*   apr, adj        - Apriori value and adjustment to parameter
 
      real*8 apr, adj
 
*   sr          - Site name read
 
      character*8 sr
 
*   est         - if * then estimated
 
 
      character*1 est
 
      read(line, 100, iostat=ierr) nt, est, sr, apr, adj
  100 format(i5,a1,a4,18x, 2d30.17)
 
*     Find the site name
      indx = 1
      call get_cmd(sr, qsite_names, qnum_sites, sn,indx )
*                         ! Save values
      if( sn.gt.0 ) then
          site_vel(1,sn) = apr
 
*                                     ! Parameter estimated
          if( est.eq.'*' ) then
              np = np + 1
              qparn_vel(1,sn) = np
              atoo(np) = nt
              qsol(np) = adj
          end if
 
      end if
 
*     Now read the next two lines (with long and radius
      do i = 1, 2
          read(unit,'(a)', iostat=ierr ) line
          read(line, 100, iostat=ierr) nt, est, sr,
     .                    apr, adj
 
*         Save values
*                             ! Save values
          if( sn.gt.0 ) then
              site_vel(i+1,sn) = apr
              if( est.eq.'*' ) then
                  np = np + 1
                  qparn_vel(i+1,sn) = np
                  atoo(np) = nt
                  qsol(np) = adj
              end if
          end if
*                 ! Looping on longitude and radius
      end do
 
***** Now set up the scaling.  There is not change of units so set 1.00 
      if( sn.gt.0 ) then
          qscale(np-2) = 1.d0
          qscale(np-1) = 1.d0
          qscale(np)   = 1.d0 
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_PMU
 
      subroutine read_pmu( unit, line, np )
 
      implicit none

*     This routine will read the polar motion and UT1 estimated
*     parameters in the hfile and 
*     save the scaling and parameter number information.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
*   np          - Running count on parameters
 
      integer*4 unit, np
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   type        - Type of polar motion/UT1:
*                 type 1 - X Pole
*                 type 2 - Y Pole
*                 type 3 - UT1-AT 
*   icmp         - Component (1= offset, 2=rate)
*   nt          - Temporary parameter number read from file
 
      integer*4 ierr, type, icmp, nt
 
*   apr, adj        - Apriori value and adjustment to parameter
 
      real*8 apr, adj
 
*   sr          - Parameter name read
 
      character*16 sr
 
*   est         - if * then estimated
 
      character*1 est
 
      read(line, 100, iostat=ierr) nt, est, sr, apr, adj
  100 format(i5,a1,a16, 4x, F25.16, D30.17)
 
*     See which component it is.
      type = 0
      icmp = 0
      if( index(sr,'X POLE').gt.0 ) type = 1
      if( index(sr,'Y POLE').gt.0 ) type = 2
      if( index(sr,'UT1-TAI').gt.0 ) type = 3

*     See if offset or rate
      icmp = 1
      if( index(sr,'RATE').gt.0 ) icmp = 2 
*                         ! Save values
      if( type.gt.0 .and. icmp.gt.0 ) then

*                                     ! Parameter estimated
          if( est.eq.'*' ) then
              np = np + 1
              qparn_pmu(icmp,type) = np
              atoo(np) = nt
              qsol(np) = adj

*****         Now set up the scaling.  
              if( (type.eq.1.or. type.eq.2) .and. icmp.gt.0 ) then
                  qscale(np)   = 1.d3 
              end if
              if( type.eq. 3 .and. icmp.gt.0 ) then
                  qscale(np)   = 15.d3 
              end if
          end if
      end if

*     Plus save the apriori values in the correct units.
      if( (type.eq.1.or. type.eq.2) .and. icmp.gt.0 ) then
          pmu_pos(icmp,type) = apr * 1.d3
      end if
      if( type.eq. 3 .and. icmp.gt.0 ) then
          pmu_pos(icmp,type) = apr * 15.d3
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_SVANT
 
      subroutine read_svant ( unit, line, np, line_preread )
 
      implicit none

*     Routine to read the satellite antenna offset parameters.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
*   np          - Running count on parameters
 
      integer*4 unit, np

*   line_preread - Logical to indicate that Satellite X position
*                  line has been pre-read because we we searching
*                  for the radition parameters

      logical line_preread
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   vn          - Satellite vechile number
*   i           - Loop counter
*   nt          - Temporary parameter number read from file
*   svcode      - The svs element code number (see svel_to_code)
*   el_read     - Number of elements to read after the X satellite
*   wn          - Test satellite number to see if changed
*                 position
 
      integer*4 ierr, vn, i, nt, svcode, el_read, wn  
 
*   apr, adj        - Apriori value and adjustment to parameter

      real*8 apr, adj                                        

*   est         - if * then estimated

      character*1 est

*   svp         - Character version of SV PRN number
*   swp         - Character version of SV PRS number
 
      character*2 svp, swp

*   svelem     - Name of the satellite element read
*   gkelem     - Dummy name for globk output
*   arelem     - Dummy name for arc output

      character*16 svelem
      character*14 gkelem
      character*4  arelem

 
 
***** Start decoding the line which we just read.  Then do the 8
*     remaining values

*     See if we know number of radiation parameters
      line_preread = .false.
      if( .not.svant_known ) then
          svant_num = 1
          el_read = 0
      else
*         Now get the number of elements minus 1 (since first is 
*         always read before hand).
          el_read = svant_num - 1
      end if
 
      read(line, 200, iostat=ierr) nt, est,svelem,  svp, apr, adj
  200 format(i5,a1,a16,2x,a2, f25.16, d30.17)
 
*     Find the PRN number matched against the SVS name
      vn = -1
 
*                         ! change leading blank to zero
      if( svp(1:1).eq.' ' ) svp(1:1) = '0'
      do i = 1, qnum_svs
          if( svp.eq. qsvs_names(i)(6:7) ) then
              vn = i
          end if
      end do
*                         ! Save values
      if( vn.gt.0 ) then
          call svel_to_code(svelem, gkelem, arelem, 
     .                      svcode, 'NTOC')
          svs_pos(svcode,vn) = apr
 
          if( est.eq.'*' ) then
              np = np + 1
              qparn_svs(svcode,vn) = np
              qsol(np) = adj
              atoo(np) = nt
              qscale(np) = 1.d0
          end if
 
      else

*         No need to tell user any more.  (If no data then satellites
*         are not used and this message would be invoked.)
c         write(*,250) svp, line(1:trimlen(line))
c 250     format(' ** ERROR ** Could not find PRN ',a2,' from line:'
c    .            ,/,a)
*                                 ! Do nothing, but sure we check that
*                                 ! there is no satellite later. 
C         vn = max_glb_svs
 
      end if
 
*     Now read the next five lines which are the posiiton and velocity.
*     Then start decoding the radiation parameter models.
      do i = 1, el_read
          read(unit,'(a)', iostat=ierr ) line
          read(line, 200, iostat=ierr) nt, est,svelem,  svp, apr, adj
 
*         Save values
*                             ! Save values
          if( vn.gt.0 ) then

*             Find out the type of element this is
              call svel_to_code(svelem, gkelem, arelem, 
     .                          svcode, 'NTOC')
              svs_pos(svcode,vn) = apr
              if( est.eq.'*' ) then
                   np = np + 1
                   qparn_svs(svcode,vn) = np
                   qsol(np) = adj
                   atoo(np) = nt
                   qscale(np) = 1.d0 
              end if
          end if
*                     ! Looping over orbital elements
      end do

****  Now see if we should search for the radiation parameter
*     model.
      if( .not.svant_known ) then

*         Start read
          line_preread = .false.
          do while ( .not.line_preread )
              read(unit,'(a)', iostat=ierr ) line
              read(line, 200, iostat=ierr) nt, est,svelem, swp, apr, adj 

*             Find the PRN Number (replace blank with zero)
              if( swp(1:1).eq.' ' ) swp(1:1) = '0'
              wn = -1
              do i = 1, qnum_svs
                  if( swp.eq. qsvs_names(i)(6:7) ) then
                     wn = i
                  endif
              enddo

*             Get the parameter code
              call svel_to_code(svelem, gkelem, arelem, svcode, 'NTOC')

*             See if satellite has changed.
* MOD TAH 020510: Check on character versions of PRN in case satellite not actually
*             used (in which case vn and wn will be zero).
              if( svp.eq.swp ) then

*                 Sum number of radiation parameters
                  svant_num = svant_num + 1
*                 Rescale the apriori value
                  if( vn.gt.0 ) svs_pos(svcode,vn) = apr
                  if( est.eq.'*' ) then
                      np = np + 1
                      qparn_svs(svcode,vn) = np
                      qsol(np) = adj
                      atoo(np) = nt
                      qscale(np) = 1.d0 
                  end if
              else
*                 We can not decode so must be next sattellite
                  line_preread = .true.
                  svant_known = .true.
                  write(*,300) svant_num
 300              format(i4,' Satellite Antenna offsets found')
              end if
          end do
      end if
 
****  Thats all
      return
      end

CTITLE READ_ATM
 
      subroutine read_atm( unit, line, np )
 
      implicit none

*     This routine will read the atmospheric delays from the
*     hfile and save the scaling and parameter number information.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
*   np          - Running count on parameters
 
      integer*4 unit, np
 
*   line        - Line read from file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   indx        - Pointer in string
*   sn          - Site number
*   nt          - Temporary parameter number read from file
 
      integer*4 ierr, indx, sn, nt
 
*   apr, adj        - Apriori value and adjustment to parameter
 
      real*8 apr, adj
 
*   sr          - Site name read
 
      character*8 sr
 
*   est         - if * then estimated
 
 
      character*1 est
 
      read(line, 100, iostat=ierr) nt, est, sr, apr, adj
  100 format(i5,a1,a4,16x, f25.17, d30.17)
 
*     Find the site name
      indx = 1
      call get_cmd(sr, qsite_names, qnum_sites, sn,indx )
*                         ! Save values
      if( sn.gt.0 ) then
          site_atm(sn) = apr
 
*                                     ! Parameter estimated
          if( est.eq.'*' ) then
              np = np + 1
              qparn_atm(sn) = np
              atoo(np) = nt
              qsol(np) = adj
              qscale(np) = 1.d0
          end if
 
      end if
 
****  Thats all
      return
      end














