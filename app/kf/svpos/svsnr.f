 
      program svsnr

      implicit none 
 
*     Program to output SNR by time, satellite and elevation angle/
 
*     The runstring of the program is
*     % svsnr [navfile] [data file]
*     where nav_file is nave file name
*         data_file is the name of the Rinex data file.
 
      include 'svsnr.h'
 
* Main pogram variables

      integer*4 i

*   rad     - Radius of site coordinates (used to check if resonably
*             valid).

      real*8 rad
 
*   eof     - Indicates end of file.
 
      logical eof
 
***** Get the runstring runstring
 
      call get_svsrun
 
*     Read in the ephemeris file
      call read_nav
 
*     Now loop over the data and get an estimate of the site position.
 
      call read_data_head
 
      write(*,100) 
 100  format('* SVSNR Output:',/,
     .       '* Code:  0 - Pcode, 1 - Xcorrelation',/,
     .       '*    Date           Sec.   SecofD  PRN      Az (deg)',
     .       '  El (deg)  L1    L2    Code')

      eof = .false.
      i = 0

*     Check that we have coordinates for the site.
      rad = sqrt(site_xyz(1)**2 + site_xyz(2)**2 + site_xyz(3)**2)
*     If the radius is less than 6,000 km then stop.
      if( rad.lt. 6.d6 ) then
          write(*,200) site_xyz
 200      format(/,'***DISASTER*** Site coordinates are not valid',/,
     .        15x, 'X = ',F12.2,' Y = ',F12.2,' Z = ',F12.2,' m.',/,
     .        15x, 'Use runstring option to give valid coordinates')
          stop 'SVSNR: Invalid site coordinates'
      endif

      call xyz_to_geod(loc_rot, site_xyz, site_loc)
      do while ( .not.eof)
          call read_range(eof)
          i = i + 1 
          if( .not.eof) call write_out
      end do
 
      end

CTITLE write_out
 
      subroutine write_out

      implicit none 
 
*     Routine to increment the solution for the new data.
 
      include  '../includes/const_param.h'
      include 'svsnr.h'

* sec10 - minutes and seconds in tenths of seconds units.

      integer*4 i,j, num_av
      
      real*8 pc(max_sat), average ,  svclk(max_sat)

      real*8  clock_epoch
 
****  compute the epheremis position at the measurement
*     time
 
      do i = 1, num_sat
 
          call eph_to_xyz( data_epoch-0.066666/86400.d0, i, 'E')
          svclk(i) = (af0(i)+
     .                af1(i)*(data_epoch-toe_jd(i))*86400.d0)*
     .               vel_light

      end do
 
      do i = 1, num_chan
 
          omc_OK(i) = .false.
          do j = 1, num_sat
              if( chan(i).eq.prn(j) ) then
 
*                 Compute a rough range
                  p1c(i) = sqrt( (site_xyz(1)-svs_xyz(1,j))**2+
     .                           (site_xyz(2)-svs_xyz(2,j))**2+
     .                           (site_xyz(3)-svs_xyz(3,j))**2)
 
*                 accumulate the average clock offset
                  omc(i) =  p1o(i)-p1c(i)+svclk(j)
                  omc_OK(i) = .true.
               end if
           end do
      end do 

*     Now get the average clock offset
      average = 0.d0
      num_av  = 0
      do i = 1, num_chan
         if( omc_OK(i) ) then
             average = average + omc(i)
             num_av  = num_av + 1
         end if
      end do
 
      if( num_av.gt.0 ) average = average / num_av
C     average = 0
 
****  Now do the final range computation

      do i = 1, num_chan
          do j = 1, num_sat
              if( chan(i).eq.prn(j) ) then
 
*                 Compute the range accounting for the propagation
*                 delay.
                  clock_epoch = data_epoch -
     .                          (average+p1c(i))/vel_light/86400.d0
                  svclk(i) = (af0(i)+
     .                        af1(i)*(clock_epoch -
     .                                toe_jd(i))*86400.d0)*
     .                        vel_light
                  call eph_to_xyz( data_epoch-
     .                            (average+p1c(i))/vel_light/86400.d0,
     .                             j, 'E')
                  p1c(i) = sqrt( (site_xyz(1)-svs_xyz(1,j))**2+
     .                           (site_xyz(2)-svs_xyz(2,j))**2+
     .                           (site_xyz(3)-svs_xyz(3,j))**2)
                  if( p2o(i).gt.0 ) then
                      pc(i) = (p1o(i)*(77.d0/60.d0)**2-p2o(i))/
     .                     ((77.d0/60.d0)**2-1.d0)
                  else
                      pc(i) = p1o(i)
                  end if

****              compute O-C
                  omc(i) = pc(i) -p1c(i)+svclk(j)         
              end if
          end do
      end do

      call svs_skymap(data_epoch)

****  Thats all
      return
      end

CTITLE READ_RANGE
 
      subroutine read_range(eof)

      implicit none 
 
*     This routine will read the next group of ranges from the rinex file
*     eof if returned true if we run out of data.
 
      include '../includes/const_param.h'
      include 'svsnr.h'
 
*   date(5)     - Ymdhm of observation
*   flags(5)    - Flags read from file.
*   ierr        - IOSTAT Error
*   i,j         - Loop counter
*   id          - Dummy entry for power fail flag
 
      integer*4 date(5), flags(5), ierr, i,j, id, trimlen
 
*   sectag      - Seconds tag
*   vals(5)     - Range values read
 
      real*8 sectag, vals(5)
 
*   eof     - Indicates end of file.
 
      logical eof
 
      character*256 line
 
****  Read in the next line
      read(101,'(a)', iostat=ierr) line
      if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
          eof = .true.
          RETURN
      end if
 
c Changed to a formatted read to avoid problems with rinex files that specify
c "G" for GPS satellites in the date/time/svs line of the rinex file. Simon 970723
* MOD TAH 971013: Corrected obvious problem with using 12(a,i2) with no argument
*     read.
      read(line,100,iostat=ierr) date, sectag, id, num_chan,
     .            (chan(i),i=1,num_chan)
 100  format(5i3,f11.7,i3,i3,12(1x,i2))
      call ymdhms_to_jd( date, sectag, data_epoch)
 
*     Now loop over the data records.
      do i = 1, num_chan
          read(101,120, iostat=ierr) (vals(j), flags(j), j =1,5)
 120      format(10(f14.3,1x,i1))
 
*         Now assign the phase and range measurements
          code_type(i) = 1
          do j = 1,5
              if((data_types(j).eq.'C1' .or.
     .            data_types(j).eq.'P1').and. vals(j).ne.0) 
     .                       P1o(i) = vals(j)
              if( data_types(j).eq.'P1' .and. vals(j).ne.0)
     .                       code_type(i) = 0
              if( data_types(j).eq.'L1' ) L1o(i) = vals(j)/
     .                    (2*77*10.23d6)*vel_light
              if( data_types(j).eq.'L2' ) L2o(i) = vals(j)/
     .                    (2*60*10.23d6)*vel_light
              if( data_types(j).eq.'P2' ) P2o(i) = vals(j)
              if( data_types(j).eq.'C2' ) P2o(i) = vals(j)
              if( data_types(j).eq.'S1' ) S1(i) = vals(j)
              if( data_types(j).eq.'S2' ) S2(i) = vals(j)
          end do
          if( num_data_types.gt.5 ) then
              read(101,120, iostat=ierr) (vals(j), flags(j), 
     .                   j =1,num_data_types-5)
              do j = 1,num_data_types-5
                  if( data_types(j+5).eq.'S1' ) S1(i) = vals(j)
                  if( data_types(j+5).eq.'S2' ) S2(i) = vals(j)
              end do
          end if
          
*         See if we need to convert log SNR
          if( rcv_type.eq.1 ) then
              S1(i) = exp(S1(i)/25.d0)/100.d0
              S2(i) = exp(S2(i)/25.d0)/100.d0
          else if ( rcv_type.eq.2 ) then
c              S1(i) = exp(S1(i)/5.d0)/100.d0
c              S2(i) = exp(S2(i)/5.d0)/100.d0
* MOD TAH 040315: Mod'd uZ to be db-Hz (power to voltage)
	      S1(i) = 10.d0**(S1(i)/20.d0)
	      S2(i) = 10.d0**(S2(i)/20.d0)
          else
* MOD TAH 040315: Mod'd all others to db-Hz (power to voltage)
	      S1(i) = 10.d0**(S1(i)/20.d0)
	      S2(i) = 10.d0**(S2(i)/20.d0)
          end if          
      end do

      if( ierr.ne.0 ) eof = .true.
 
****  Thats all
      return
      end
 
 
CTITLE READ_DATA_HEAD
 
      subroutine read_data_head

      implicit none 
 
*     Routine to read the header from the data files.
*     Returns will be:
*     Approximate site position
*     station name
 
      include 'svsnr.h'
 
* Local variables
 
*   ierr    - IOSTAT error
*  trimlen  - Length of string
*   indx    - Position of substring in a string
 
      integer*4 ierr, trimlen, indx, i
      integer*4 blks   ! Numnber of data type blocks based on num_data_type
     .,         it     ! Counter reading through blocks
     .,         str, rem   ! Start element (9 per line) and
                       ! remaining number on line.
 
*   eoh     - End of header flags
*   rinex   - True if rinex file
 
      logical eoh, rinex
 
*   line    - Line read from file
 
 
      character*256 line
      character*1 cr
 
****  Open the data file
      cr = char(13) 
      open(101,file=data_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',data_file, 1,
     .        'read_data_file/svsnr')
 
*     Start looping over the header
      eoh = .false.
      rinex = .false.
      rcv_type = 0
      
      do while ( .not.eoh )
          read(101,'(a)', iostat=ierr) line
          call sub_char(line,cr,' ')
          call casefold(line)
          call report_error('IOSTAT',ierr,'read', data_file,1,
     .                'End of hedader not found')
          if( trimlen(line).lt.10) eoh = .true.
          if( index(line,'END OF HEAD').gt.0 ) eoh = .true.
 
*         See if we find Rinex file
          indx = index(line,'OBSERVATION DATA')
          if( indx.gt.0 ) rinex = .true.
 
*         See if marker name
          indx = index(line,'MARKER NAME')
          if( indx.gt.0 ) site_name = line(1:8)
 
*         See if position
          indx = index(line,'APPROX POSITION')
          if( indx.gt.0 .and. site_xyz(1).eq.0 ) then
              read(line,*,iostat=ierr) site_xyz
          end if
          
*         See if Log SNR values
          indx = index(line,'REC # / TYPE / VERS')
          if( indx.gt.0 ) then
              if( index(line,'XII').gt.0 .or.
     .            index(line,' Z-12').gt.0      ) then
                  rcv_type = 1
                  write(*,'(a)') '* Log SNR receiver found'
              else if( index(line,'UZ-12').gt.0  ) then
                  rcv_type = 2
                  write(*,'(a)') '* Log uZ SNR receiver found'
              end if 
          end if
 
*         See if data types
          indx = index(line,'TYPES OF OBS')
           if( indx .gt.59 ) then
              read(line,'(i6,9a6)', iostat=ierr) num_data_types, 
     .               (data_types(i), i=1,min(9,num_data_types))
* MOD TAH 130330: Allow more data types
              if( num_data_types.gt.max_data_types) then
                  call report_stat('FATAL','SVPOS','Read TYPES OF OBS',
     .               data_file,'Too many types of data',
     .               max_data_types)
              end if
              if( num_data_types.gt.9 ) then
* MOD TAH 151020: Loop over the remaining blocks needed
                  blks = (num_data_types-1)/9+1  ! Includes count for values
                                                 ! already read
                  do it = 2,blks
                     rem = min(9,num_data_types-(it-1)*9)
                     str = it*9
                     read(101,'(a)', iostat=ierr) line
                     read(line,'(6x,9a6)', iostat=ierr)  
     .                  (data_types(str+i), i=1,rem)
                  enddo
              endif
          end if
      end do
 
*     See if we found rinex file
      if( rinex ) then
          write(*,150) data_file(1:trimlen(data_file)),
     .        site_name, site_xyz
 150      format('* For RINEX data file ',a,/,
     .        '* Site ',a,' Aprrox. position ',3F15.3)
      else
          write(*,170) data_file(1:trimlen(data_file))
 170      format(a,' Does not appear to be RINEX file')
          stop 'svsnr: Wrong type data file'
      end if
 
****  Thats all
      return
      end
 
 
 
CTITLE GET_SVSRUN
 
      subroutine get_svsrun

      implicit none 
 
*     Routine to get the cunstring.
 
      include '../includes/const_param.h'
      include 'svsnr.h'
 
* LOCAL VARIABLES
 
*  len_run     - Length of the runstring element
*  trimlen     - Length of string
* date_start(5), date_stop(5) - Start and stop dates
*               for generating results
*  i           - Loop counter
*  rcpar       - Gets runstring entry
 
 
      integer*4 len_run, i, rcpar
 
*  runstring   - Elements of runstring
 
      character*256 runstring
 
*     Get the first runstring parameter
      len_run = rcpar(1,nav_file)
      if( len_run.le.0 ) then
          call proper_runstring('svsnr.hlp','svsnr/nav file',1)
      end if
 
*     See if latitude and longitude passed
      len_run = rcpar(2,data_file)
      if( len_run.le.0 ) then
          call proper_runstring('svsnr.hlp','svsnr/data file',1)
      end if

*     Now see if positions passed 
      do i = 1,3
         len_run = rcpar(2+i,runstring)
         if(  len_run.gt.0 ) then
             read(runstring,*) site_xyz(i)
         else
             site_xyz(i) = 0.d0
         end if
      end do

****  Thats all
      return
      end
 
CTITLE READ_NAV
 
      subroutine read_nav

      implicit none 
 
*     Thuis routine will read the navigation file.  Only the first
*     occurence of a satellite ephemeris entry will be used.  (This
*     is set by the PRN number still being zero)
 
      include 'svsnr.h'
 
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
     .            'svsnr/read_nav')
 
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
 120          format('* Rinex version ',F4.2,' Nav file found')
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
 
      if( num_sat.eq.0 ) stop ' svsnr: No satellites found'
 
      return
      end
CTITLE EPH_TO_XYZ
 
      subroutine eph_to_xyz(t, i, sys)

      implicit none 
 
*     Routine to compute XYZ coordinates of satellite at time t for
*     sattellite number i.  SYS is 'I' or 'E' for inertial or Earth
*     fixed.
 
      include 'svsnr.h'
 
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
 
CTITLE SVS_skymap
 
      subroutine svs_skymap(t)

      implicit none 
 
*     Routine produce sky map of GPS satellites in earthfixed 
*     frame.  This makes file svs_[nav_file].skymap

      include '../includes/const_param.h' 
      include 'svsnr.h'
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   i,j     - Loop counters
*   date(5) - ymdhm
 
      integer*4  i, j, date(5), k
 
*   t       - JD of calculation
*   sectag  - Seconds tag on date
*   unit_xyz(3), unit_loc(3)  - Units vectors from site to satellite
*             in global and local frame
*   rang, hlen - Range to satellite and horizontal length of unit
*              vector
 
      real*8 t, sectag, unit_xyz(3), unit_loc(3), rang, hlen,
     .       sec_of_day

*     Write out the line into the file
      call jd_to_ymdhms( t, date, sectag)
      sec_of_day = date(4)*3600 + date(5)*60 + sectag

      do k = 1, num_chan
          do i = 1, num_sat
              if( chan(k).eq.prn(i) ) then

*                 Now convert the vector position of satellite into
*                 topocentric frame so that we can AZ and El.
                  rang = 0
                  do j = 1,3
                     rang = rang + (svs_xyz(j,i)-site_xyz(j))**2
                  end do
                  rang = sqrt(rang)
                  do j = 1,3
                     unit_xyz(j) = (svs_xyz(j,i)-site_xyz(j))/rang
                  end do

*                 now transform global XYZ into local NEU (loc_rot is
*                 computed when the site coords are given in runstring.)
                  do j = 1,3
                      call dvdot(unit_loc(j), loc_rot(j,1),3, unit_xyz,
     .                           1,3)
                  end do

*                 Compute the horizontal length of the unit vector
                  hlen = sqrt(unit_loc(1)**2+unit_loc(2)**2)
                  svs_ang(2,i) = atan2(unit_loc(3),hlen)*180/pi
                  svs_ang(1,i) = atan2(unit_loc(2),unit_loc(1))*180/pi
                  if(svs_ang(1,i).lt.0 ) svs_ang(1,i) = svs_ang(1,i)+360
                  svs_ang(3,i) = rang

*                 Convert to plotable quanity.  North Up on the page.  This
*                 way we can just plot the x y coordinates (i.e., mapping
                  write(*,200)  date, sectag, sec_of_day, prn(i), 
     .                  svs_ang(1,i), svs_ang(2,i), s1(k), s2(k),
     .                  code_type(k)
 200              format(i5,4i3,F6.1,f10.1,1x,'PRN ',i2.2,1x, 
     .                   2f10.3,1x,2f6.1,3x,i2)
              end if
          end do
      end do
      return
      end

