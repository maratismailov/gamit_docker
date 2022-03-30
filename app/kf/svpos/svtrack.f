 
      program svtrack

      implicit none 
 
*     Program to compute the rise and set times of the GPS
*     satelites based on the ephemeris in the rinex navigiation
*     files
 
*     The runstring of the program is
*     % svtrack [navfile] site_lat site_long <'start yy mm dd hh \
*                     stop yy mm dd hh step (min)'>
*     where namvfile is nave file name
*         site_lat and site_long are the lat and long of the site
*                 (decimal degrees)
*         'start yy mm dd hh  stop yy mm dd hh step (min)' is the
*                 start and stop times (to nearest hour) and the
*                 step size in minutes.  These values are all enclosed
*                 in single quotes.  The defaults are ar for 1 day
*                 from the ephemeris epoch in 10 minute steps.
 
      include 'svtrack.h'
 
***** Get the runstring runstring
 
      call get_svsrun
 
*     Read in the ephemeris file
      call read_nav
 
*     Now prduce each of the types of data files we need
 
*     Inertial 3-d positions positions.  Name svs_[nav_file].i3d
      call svs_i3d
 
*     Terrestrial 3-d positions positions.  Name svs_[nav_file].t3d
 
       call svs_t3d
 
*     Ground track. Name svs_[nav_file].gtrck
 
       call svs_gtrck
 
*     Sky_map svs_[nav_file].sky
 
       call svs_skymap
 
*     Rise and set times svs_[nav_file].rs
 
       call svs_riseset
 
****  Thats all
      end
 
CTITLE GET_SVSRUN
 
      subroutine get_svsrun

      implicit none 
 
*     Routine to get the runstring.
 
      include '../includes/const_param.h'
      include 'svtrack.h'
 
* LOCAL VARIABLES
 
*  len_run     - Length of the runstring element
*  trimlen     - Length of string
* date_start(5), date_stop(5) - Start and stop dates
*               for generating results
*  i           - Loop counter
*  rcpar       - Gets runstring entry 
 
      integer*4 len_run,  date_start(5), date_stop(5),
     .          i, rcpar, ierr
 
      real*8 sectag
 
*  runstring   - Elements of runstring
 
 
      character*256 runstring
 
*     Get the first runstring parameter
      len_run = rcpar(1,nav_file)
      if( len_run.le.0 ) then
          call proper_runstring('svtrack.hlp','svtrack',1)
      end if
 
*     See if latitude and longitude passed
      len_run = rcpar(2,runstring)
      if( len_run.gt.0 ) then
          read(runstring,*,iostat=ierr) lat_deg
          call report_error('IOSTAT',ierr,'read',runstring,1,
     .                   'get_svsrun/latitude')
      else
          lat_deg = 45.d0
      end if
 
*     See if longitude passed
      len_run = rcpar(3,runstring)
      if( len_run.gt.0 ) then
          read(runstring,*,iostat=ierr) long_deg
          call report_error('IOSTAT',ierr,'read',runstring,1,
     .                   'get_svsrun/longitude')
      else
          long_deg = -70.d0
      end if
 
*     Convert the lat and long to XYZ
      site_loc(1) = pi/2 - lat_deg*pi/180.d0
      site_loc(2) = long_deg*pi/180.d0
      site_loc(3) = 0.d0
      call geod_to_xyz(site_loc, site_xyz)

*     Now get the transformation matrix from XYZ to NEU
      call xyz_to_geod(loc_rot, site_xyz, site_loc)
 
*     Now see if start and stop times passed
      len_run = rcpar(4,runstring)
      if( len_run.gt.0 ) then
          read(runstring,*,iostat=ierr) (date_start(i),i=1,4),
     .           (date_stop(i),i=1,4), step
          call report_error('IOSAT',ierr,'read',runstring,1,
     .                   'get_svsrun/dates and step')
          date_start(5) = 0
          sectag = 0.d0
          call ymdhms_to_jd(date_start, sectag, start)
          date_stop (5) = 0
          call ymdhms_to_jd(date_stop,sectag, stop)
 
*         Convert step from minutes to days
          step = step / 1440.d0
      else
 
*         Start the start and stop times to zero.  We will use
*         the ephemeris time later to get these values
          start = 0
          stop = 0
          step = 10.d0/1440.0
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_NAV
 
      subroutine read_nav

      implicit none 
 
*     Thuis routine will read the navigation file.  Only the first
*     occurence of a satellite ephemeris entry will be used.  (This
*     is set by the PRN number still being zero)
 
      include 'svtrack.h'
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   i           - Loop counter
*   rp          - Read Prn number
*   trimlen     - length of string
*   ll          - Length of line string
*   date(5)     - ymdhm
*   rinex_version - Version of rinex file
 
 
      integer*4 ierr, i, rp, trimlen, ll, date(5) 
      real*4 rinex_version
 
*   sectag      - Seconds tag in date
 
       real*8 sectag
 
*   still_header    - Indicates that we are sill in the header
*                   - part of file.
*   have            - Set true if we already have this satellite
 
       logical still_header, have
 
*   line            - line read from file
 
       character*256 line

*   cr - Carriage return (for handling dos files)
      character*1 cr

      cr = char(13)
 
****  Open the nav_file
 
      open(100, file=nav_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',nav_file,1,
     .            'svpos/read_nav')
 
*     Loop over the file reading the ephemeris entries.  First clear
*     all of the PRN numbers so we know when a PRN has been read
 
      do i = 1, max_sat
          prn(i) = 0
      end do
      rinex_version = 1
 
      still_header = .true.
      do while ( still_header )
          read(100,'(a)', iostat=ierr) line
          call sub_char(line,cr,' ')
          call casefold(line)
          ll = trimlen(line)
          if( ll.eq.0 .or. index(line,'END OF HEAD').gt.0 ) then
              still_header = .false.
          end if
          call report_error('IOSTAT',ierr,'read','NAV FILE HEADER',
     .                      1,'read_nav')
          if( index(line,'RINEX VERSION').gt.0 ) then
              read(line,*, iostat=ierr) rinex_version
              write(*,120) rinex_version
 120          format(' Rinex version ',F4.2,' Nav file found')
          end if
      end do
 
      num_sat = 0
 
*     Now start reading the entries
      do while ( ierr.eq.0 )
          read(100,'(a)',iostat=ierr) line
          call sub_char(line,cr,' ')
          if( trimlen(line).eq.0 ) ierr = -1
*                                     ! See which PRN
          if( ierr.eq.0 ) then
              read(line,*) rp
 
*             See if we already have this prn
              have = .false.
              do i = 1, num_sat
                  if( rp.eq.prn(i) ) then
                      have = .true.
                  end if
              end do
 
*                                         ! This is first
              if( .not.have ) then
*                                             ! on this prn so read
                  num_sat = num_sat + 1
                  rp = num_sat
                  read(line,200) prn(rp), date, sectag,
     .                        af0(rp),af1(rp), af2(rp)
 200              format(i2,5i3,f5.1,3d19.8)
 
                  call ymdhms_to_jd( date, sectag, toe_jd(rp))
 
*                 Read the rest of the entires
                  read(100,210) aode(rp), crs(rp), dn(rp), m0(rp)
                  read(100,210) cuc(rp), ecc(rp), cus(rp), art(rp)
                  read(100,210) toe(rp), cic(rp), om0(rp), cis(rp)
                  read(100,210) i0(rp) , crc(rp), w(rp)  , omd(rp)
                  read(100,210) idt(rp), cflg12(rp), weekno, pflg12(rp)
                  read(100,210) svacc, svhealth, tgd, aodc(rp)
 210              format(3x,4d19.8)
                  if( rinex_version.gt.1 ) read(100,'(a)') line
              else
*                 Skip 6 lines in file
                  do i = 1,6
                      read(100,'(a)') line
                  end do
                  if( rinex_version.gt.1 ) read(100,'(a)') line
              end if
          end if
      end do
 
      write(*,300) num_sat, nav_file(1:trimlen(nav_file))
 300  format('* ', i5,' satellites found in ',a)
 
      if( num_sat.eq.0 ) stop ' SVPOS: No satellites found'
 
      return
      end
 
CTITLE SVS_I3D
 
      subroutine svs_i3d

      implicit none 
 
*     Routine produce XYZ coordinates of GPS satellites in interial
*     space.  This makes file svs_[nav_file].i3d
 
      include 'svtrack.h'
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   i,j     - Loop counters
*   date(5) - ymdhm
 
      integer*4 ierr, i,j, date(5), trimlen, it
 
*   name    - Output file name
 
      character*256 name
 
*   t       - JD of calculation
*   sectag  - Seconds tag on date
 
 
      real*8 t, sectag
 
****  Generate the output file name and open the files
      call gen_name( nav_file, '.i3d', name)
 
      open(200, file = name, iostat=ierr,status='unknown')
      call report_error('IOSTAT',ierr,'open',name,1,'svs_i3d')

      write(*,100) name(1:trimlen(name))
 100  format('Generating ',a)
 
*     Now loop over all times and satellites
c     do t = start, stop, step
      do it = 0,nint((stop-start)/step)
          t = start + it*step
          do i = 1, num_sat
              call eph_to_xyz(t, i, 'I')
          end do
 
*         Write out the line into the file
          call jd_to_ymdhms( t, date, sectag)
          write(200,200) date, ((svs_xyz(j,i)/1.d6,j=1,3),i=1,num_sat)
 200      format(i5,4i3,32(3(F8.3,1x),2x))
      end do
 
      close(200)
      return
      end
 
CTITLE SVS_t3d
 
      subroutine svs_t3d

      implicit none 
 
*     Routine produce XYZ coordinates of GPS satellites in earthfixed 
*     frame.  This makes file svs_[nav_file].t3d
 
      include 'svtrack.h'
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   i,j     - Loop counters
*   date(5) - ymdhm
 
      integer*4 ierr, i,j, date(5), trimlen, it
 
*   name    - Output file name
 
      character*256 name
 
*   t       - JD of calculation
*   sectag  - Seconds tag on date
 
 
      real*8 t, sectag
 
****  Generate the output file name and open the files
      call gen_name( nav_file, '.t3d', name)
 
      open(200, file = name, iostat=ierr, status='unknown')
      call report_error('IOSTAT',ierr,'open',name,1,'svs_t3d')

      write(*,100) name(1:trimlen(name))
 100  format('Generating ',a)
 
*     Now loop over all times and satellites
c     do t = start, stop, step
      do it = 0,nint((stop-start)/step)
          t = start + it*step
          do i = 1, num_sat
              call eph_to_xyz(t, i, 'E')
          end do
 
*         Write out the line into the file
          call jd_to_ymdhms( t, date, sectag)
          write(200,200) date, ((svs_xyz(j,i)/1.d6,j=1,3),i=1,num_sat)
 200      format(i5,4i3,32(3(F8.3,1x),2x))
          write(*,220) date, ((svs_xyz(j,i),j=1,3),i=1,2)      
 220      format(i5,4i3,32(3(F15.3,1x),2x))
      end do

      close(200)
      return
      end

CTITLE SVS_gtrck
 
      subroutine svs_gtrck

      implicit none 
 
*     Routine produce spherical coord of GPS satellites in earthfixed 
*     frame.  This makes file svs_[nav_file].gtrck
      include '../includes/const_param.h' 
      include 'svtrack.h'
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   i,j     - Loop counters
*   date(5) - ymdhm
 
      integer*4 ierr, i, date(5), trimlen, it
 
*   name    - Output file name
 
      character*256 name
 
*   t       - JD of calculation
*   sectag  - Seconds tag on date
 
 
      real*8 t, sectag
 
****  Generate the output file name and open the files
      call gen_name( nav_file, '.gtrck', name)
 
      open(200, file = name, iostat=ierr, status= 'unknown')
      call report_error('IOSTAT',ierr,'open',name,1,'svs_gtrck')

      write(*,100) name(1:trimlen(name))
 100  format('Generating ',a)
 
*     Now loop over all times and satellites
c     do t = start, stop, step
      do it = 0,nint((stop-start)/step)
          t = start + it*step
          do i = 1, num_sat
              call eph_to_xyz(t, i, 'e')
          end do
 
*         Write out the line into the file
          call jd_to_ymdhms( t, date, sectag)
          write(200,200) date,
     .        (svs_loc(2,i)*180/pi, 90-svs_loc(1,i)*180/pi,
     .         svs_loc(3,i)/1.d6,i=1,num_sat)
 200      format(i5,4i3,32(3(F8.3,1x),2x))
      end do
 
      close(200)
      return
      end

CTITLE SVS_skymap
 
      subroutine svs_skymap

      implicit none 
 
*     Routine produce sky map of GPS satellites in earthfixed 
*     frame.  This makes file svs_[nav_file].skymap

      include '../includes/const_param.h' 
      include 'svtrack.h'
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   i,j     - Loop counters
*   date(5) - ymdhm
 
      integer*4 ierr, i, j, date(5), trimlen, it
 
*   name    - Output file name
 
      character*256 name
 
*   t       - JD of calculation
*   sectag  - Seconds tag on date
*   unit_xyz(3), unit_loc(3)  - Units vectors from site to satellite
*             in global and local frame
*   rang, hlen - Range to satellite and horizontal length of unit
*              vector
 
      real*8 t, sectag, unit_xyz(3), unit_loc(3), rang, hlen,
     .        svs_plt(3,max_sat)
 
****  Generate the output file name and open the files
      call gen_name( nav_file, '.skymap', name)
 
      open(200, file = name, iostat=ierr,status='unknown')
      call report_error('IOSTAT',ierr,'open',name,1,'svs_skymap')

      write(*,100) name(1:trimlen(name))
 100  format('Generating ',a)
 
*     Now loop over all times and satellites
c     do t = start, stop, step
      do it = 0,nint((stop-start)/step)
          t = start + it*step
          do i = 1, num_sat
              call eph_to_xyz(t, i, 'e')
*             Now convert the vector position of satellite into
*             topocentric frame so that we can AZ and El.
              rang = 0
              do j = 1,3
                 rang = rang + (svs_xyz(j,i)-site_xyz(j))**2
              end do
              rang = sqrt(rang)
              do j = 1,3
                 unit_xyz(j) = (svs_xyz(j,i)-site_xyz(j))/rang
              end do

*             now transform global XYZ into local NEU (loc_rot is
*             computed when the site coords are given in runstring.)
              do j = 1,3
                  call dvdot(unit_loc(j), loc_rot(j,1),3, unit_xyz,1,3)
              end do

*             Compute the horizontal length of the unit vector
              hlen = sqrt(unit_loc(1)**2+unit_loc(2)**2)
              svs_ang(1,i) = atan2(hlen, unit_loc(3))
              svs_ang(2,i) = atan2(unit_loc(2),unit_loc(1))
              svs_ang(3,i) = rang

*             Convert to plotable quanity.  North Up on the page.  This
*             way we can just plot the x y coordinates (i.e., mapping
*             to polar plot is done manually.)
              svs_plt(1,i) = min(svs_ang(1,i),pi/2)*sin(svs_ang(2,i))
              svs_plt(2,i) = min(svs_ang(1,i),pi/2)*cos(svs_ang(2,i))
              svs_plt(3,i) = svs_ang(3,i)
          end do
 
*         Write out the line into the file
          call jd_to_ymdhms( t, date, sectag)

c         write(200,200) date,
c    .        (svs_ang(2,i), min(svs_ang(1,i),pi/2)*180/pi,
c    .         svs_ang(3,i)/1.d6,i=1,num_sat)
          write(200,200) date, (svs_plt(1,i), svs_plt(2,i),
     .         svs_plt(3,i)/1.d6,i=1,num_sat)
 200      format(i5,4i3,32(3(F8.3,1x),2x))
      end do
 
      close(200)
      return
      end

CTITLE SVS_riseset
 
      subroutine svs_riseset

      implicit none 
 
*     Routine produce rise and set time of GPS satellites in earthfixed 
*     frame.  This makes file svs_[nav_file].riseset

      include '../includes/const_param.h' 
      include 'svtrack.h'
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   i,j     - Loop counters
*   date(5) - ymdhm
*   el_point(max_sat) - Point type to denote elevation
*             (Quantizied in bins of 5 degress)
*   el_out(max_sat)   - Y coordinate to be written
*             (Channel number if point is above the
*              horizon, -1 other wise)
 
      integer*4 ierr, i, j, date(5), trimlen, el_point(max_sat),
     .          el_out(max_sat), it
 
*   name    - Output file name
 
      character*256 name
 
*   t       - JD of calculation
*   sectag  - Seconds tag on date
*   unit_xyz(3), unit_loc(3)  - Units vectors from site to satellite
*             in global and local frame
*   rang, hlen - Range to satellite and horizontal length of unit
*              vector
 
      real*8 t, sectag, unit_xyz(3), unit_loc(3), rang, hlen
 
****  Generate the output file name and open the files
      call gen_name( nav_file, '.riseset', name)
 
      open(200, file = name, iostat=ierr, status='unknown')
      call report_error('IOSTAT',ierr,'open',name,1,'svs_riseset')

      write(*,100) name(1:trimlen(name))
 100  format('Generating ',a)
 
*     Now loop over all times and satellites
c     do t = start, stop, step
      do it = 0,nint((stop-start)/step)
          t = start + it*step
          do i = 1, num_sat
              call eph_to_xyz(t, i, 'e')
*             Now convert the vector position of satellite into
*             topocentric frame so that we can AZ and El.
              rang = 0
              do j = 1,3
                 rang = rang + (svs_xyz(j,i)-site_xyz(j))**2
              end do
              rang = sqrt(rang)
              do j = 1,3
                 unit_xyz(j) = (svs_xyz(j,i)-site_xyz(j))/rang
              end do

*             now transform global XYZ into local NEU (loc_rot is
*             computed when the site coords are given in runstring.)
              do j = 1,3
                  call dvdot(unit_loc(j), loc_rot(j,1),3, unit_xyz,1,3)
              end do

*             Compute the horizontal length of the unit vector
              hlen = sqrt(unit_loc(1)**2+unit_loc(2)**2)
              svs_ang(1,i) = atan2(hlen, unit_loc(3))
              svs_ang(2,i) = atan2(unit_loc(2),unit_loc(1))
              svs_ang(3,i) = rang

*             Convert the point to the type we need using the
*             zenith distance
              if( svs_ang(1,i).gt.pi/2 ) then
                  el_out(i) = -1
                  el_point(i) = 0
              else
                  el_out(i) = prn(i)
                  el_point(i) = int(((pi/2-svs_ang(1,i))*180/pi)/5)+1
              end if
          end do
 
*         Write out the line into the file
          call jd_to_ymdhms( t, date, sectag)

          write(200,200) date, (el_out(i), el_point(i), i=1,num_sat)
 200      format(i5,4i3,32(2i4,2x))
      end do
 
      close(200)

***** Now make the cplot comand file
      open(200,file='riseset.plt',status='unknown')

*     Write out the file command
      write(200,300) name(1:trimlen(name))
 300  format('* Cplot command file to plot rise and set times',/,
     .       '*',/
     .       ' head 0 0 ',/,
     .       ' file ',a,/,
     .       ' point -1',/,
     .       ' char 2.0',/,
     .       ' view 0.075 0.80 0.1 0.95') 

*     Now loop over all the satellites
      do i = 1, num_sat
         write(200,320) 6+(i-1)*2, 7+(i-1)*2
 320     format(' y_field 1 ',i3,' 0',/,
     .          ' p_field ',i3,' 0',/,
     .          ' read ',/,
     .          ' y_scale 0 33 ',/,
     .          ' draw')
      end do

*     Now do the axes
      write(200,340) name(1:trimlen(name)), lat_deg, long_deg
 340  format(' line 31',/,' char 2.0',/,
     .       ' label 0.01 33.2 1 0 "File ',a,' Lat ',
     .         F6.2,' Long ',f7.2,'"',/,
     .       ' ymn 1 2 "PRN Number"',/,
     .       ' xmn -1 1 "Date"',/,
     .       ' xmx -1 0',/,
     .       ' ymx 1 0')

****  Write out the legend
      write(200,360) 
 360  format(' view 0.81 0.99 0.1 0.95',/,
     .       ' char 2.0')
      write(200,370) (i/6.,i,i+5,char(i/5+ichar('A')),i=0,85,5)
 370  format(' Label 0 ',f4.1,' 1 0 "',i2,'-',i2,2x,a,'"')
      write(200,380) 
 380  format(' Label 0 15 1 0 "Range Symbol"',/,
     .       ' Label 0 16 1 0 "Elev."')

****  Thats all
      return
      end

CTITLE GEN_NAME
 
      subroutine gen_name(full_nav_file, ext, name)

      implicit none 
 
*     Rouitne to gerenrate file name form nav_file name and extent
 
 
*   full_nav_file       - nav file name
*   ext             - Extent for new name
*   name            - Generated name
 
      character*(*) full_nav_file, ext, name
 
*   start_root, end_root    - Start and end character numbers in
*                       - name root.
*   trimlen             - Length of string
*   i                   - loop counter
 
 
      integer*4 start_root, end_root, trimlen, i
 
*     Get length of full name
      end_root = trimlen(full_nav_file)
 
      i = end_root
      do while (i.ge.1 .and. full_nav_file(i:i).ne.'/' )
          i = i - 1
      end do
 
      start_root = i + 1
 
****  Construct the name
      name = 'svs_' // full_nav_file(start_root:end_root) // ext
 
****  Thats all
      return
      end
 
CTITLE EPH_TO_XYZ
 
      subroutine eph_to_xyz(t, i, sys)

      implicit none 
 
*     Routine to compute XYZ coordinates of satellite at time t for
*     sattellite number i.  SYS is 'I' or 'E' for inertial or Earth
*     fixed.
 
      include 'svtrack.h'
 
*   t       - Time for computation (day number)
 
      real*8 t
 
*   i       - Satellite number
 
      integer*4 i
 
*   sys     - System for results.
 
      character*(*) sys
 
* LOCAL VARIABLES
 
*   j           - Loop counter in Kepler's equation
 
      integer*4 j
 
*   gm      - GM
*   eom     - Earth rotation rate (rads/sec)
*   a       - Semimajor axis
*   n0      - Mean motion
*   tk      - Time of epoch from toe (seconds)
*   n       - COrrected mean motion
*   mk      - Mean anomaly
*   ek      - Eccentric anomaly
*   vk      - true anomaly
*   sinvk, cosvk    - Sin and cos of true anomaly
*   pk      - argument of latitude
*   duk, drk, dik   - Coorections to arg of lat, radius and
*           - inclinations
*   uk      - Argumenr of latitude
*   rk      - radius at time tk
*   ik      - Inclination at tk
*   xpk, ypk    - Inplane coordiantes
*   omk     - Longitude of the asecdning node
*   rot_mat(3,3)    - Rotation matrix from XYZ to NEU
 
 
 
      real*8 gm, eom, a, n0, tk, n, mk, ek, vk, sinvk, cosvk, pk,
     .    duk, drk, dik, uk, rk, ik, xpk, ypk, omk, rot_mat(3,3)
 
      gm = 3.986005d14
      eom = 7.2921151467d-5
 
      a = art(i)*art(i)
      n0 = sqrt(gm/a**3)
 
      tk = (t-toe_jd(i))*86400.0d0
      n = n0 + dn(i)
      mk = m0(i) + n*tk
 
****  Solve Keplers equation
      ek = mk
      do j = 1, 10
          ek = mk + ecc(i)*sin(ek)
      end do
 
***** Get the true anomaly
      sinvk = sqrt(1-ecc(i)**2)*sin(ek)/(1 - ecc(i)*cos(ek))
      cosvk = (cos(ek)-ecc(i))/(1-ecc(i)*cos(ek))
 
      vk = atan2(sinvk, cosvk)
 
*     Argument of latitude
      pk = vk + w(i)
 
*     Correction terms
      duk = cus(i)*sin(2*pk) +cuc(i)*cos(2*pk)
      drk = crs(i)*sin(2*pk) +crc(i)*cos(2*pk)
      dik = cis(i)*sin(2*pk) +cic(i)*cos(2*pk)
 
      uk = pk + duk
      rk = a*(1-ecc(i)*cos(ek)) + drk
      ik = i0(i) + dik + idt(i)*tk
 
*     Get inplane coordinates
      xpk = rk*cos(uk)
      ypk = rk*sin(uk)
 
*     Compute the longitude of the ascending node
      omk = om0(i) + omd(i)*tk
 
*     If we are in Earth fixed frame account for rotation of Earth
      if( sys(1:1).eq.'E' .or. sys(1:1).eq.'e') then
          omk = omk - eom*(tk+toe(i))
      end if
 
*     Get three_d coordinates
      svs_xyz(1,i) = xpk*cos(omk) - ypk*sin(omk)*cos(ik)
      svs_xyz(2,i) = xpk*sin(omk) + ypk*cos(omk)*cos(ik)
      svs_xyz(3,i) = ypk*sin(ik)
 
*     Now convert to latitude and longitude
      call xyz_to_geod( rot_mat, svs_xyz(1,i), svs_loc(1,i))
 
****  Thats all
      return
      end
 
 
 
