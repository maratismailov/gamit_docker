CTITLE READ_NAV
 
      subroutine read_nav
 
      implicit none

*     Thuis routine will read the navigation file.  Only the first
*     occurence of a satellite ephemeris entry will be used.  (This
*     is set by the PRN number still being zero) (from Rinex)
 
      include 'track_com.h'
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   i           - Loop counter
*   rp          - Read Prn number
*   trimlen     - length of string
*   ll          - Length of line string
*   date(5)     - ymdhm
*   rinex_version - Version of rinex file
 
 
      integer*4 ierr, i, rp, trimlen, ll, date(5), rinex_version,
     .          j,k, jerr, iun, it
 
      integer*4 max_fnd_prn ! Largest PRN seen.  This set the num_sat
               ! we need since slots are saved by prn number

*   sectag      - Seconds tag in date
 
       real*8 sectag, start_sp3, t, nav_mjd, r8weekno
 
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
*     Initialize the start and end times
      start_nav = 1.d30
      end_nav   = 0.d0
      iun = 102
      open(iun, file=nav_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',nav_file,1,
     .            'svpos/read_nav')
 
*     Loop over the file reading the ephemeris entries.  First clear
*     all of the PRN numbers so we know when a PRN has been read
 
      do i = 1, max_sat
          prn_sp3(i) = 0
      end do
      rinex_version = 1
 
      still_header = .true.
      do while ( still_header )
          read(iun,'(a)', iostat=ierr) line
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
          end if
      end do
 
      num_sat = 0
 
*     Now start reading the entries
      do while ( ierr.eq.0 )
          read(iun,'(a)',iostat=ierr) line 
          call sub_char(line,cr,' ')
          if( trimlen(line).lt.70 ) then
*             This is a teqc file: skip this line
* MOD TAH 000819: Check longer line since there can be an
*             extra value on this line.
              read(iun,'(a)',iostat=ierr) line
              call sub_char(line,cr,' ')
          end if    
          if( trimlen(line).eq.0 ) ierr = -1
*                                     ! See which PRN
          if( ierr.eq.0 ) then
              read(line,*,iostat=jerr) rp
 
*             See if we already have this prn
              have = .false.
              do i = 1, num_sat
                  if( rp.eq.prn_sp3(i) ) then
                      have = .true.
                  end if
              end do

*             Get the time range in the nav file
              read(line,200) j, date, sectag
 
              call ymdhms_to_mjd( date, sectag, nav_mjd)
              if( nav_mjd.lt.start_nav ) start_nav = nav_mjd
              if( nav_mjd.gt.end_nav   ) end_nav = nav_mjd
 
*                                         ! This is first
              if( .not.have ) then
*                                             ! on this prn so read
                  num_sat = num_sat + 1
                  rp = num_sat
                  read(line,200,iostat=jerr) prn_sp3(rp), date, sectag,
     .                        af0(rp),af1(rp), af2(rp)
 200              format(i2,5i3,f5.1,3d19.8)
 
                  call ymdhms_to_mjd( date, sectag, toe_jd(rp))
                  if( toe_jd(rp).lt.start_nav ) start_nav = toe_jd(rp)
                  if( toe_jd(rp).gt.end_nav   ) end_nav = toe_jd(rp)

 
*                 Read the rest of the entires
                  read(iun,210) aode(rp), crs(rp), dn(rp), m0(rp)
                  read(iun,210) cuc(rp), ecc(rp), cus(rp), art(rp)
                  read(iun,210) toe(rp), cic(rp), om0(rp), cis(rp)
                  read(iun,210) i0(rp) , crc(rp), w(rp)  , omd(rp)
                  read(iun,210) idt(rp), cflg12(rp),r8weekno,pflg12(rp)
                  weekno = nint(r8weekno)
                  read(iun,210) svacc, svhealth, tgd, aodc(rp)
 210              format(3x,4d19.8)
                  if( rinex_version.gt.1 ) read(iun,'(a)') line
              else
*                 Skip 6 lines in file
                  do i = 1,6
                      read(iun,'(a)') line
                  end do
                  if( rinex_version.gt.1 ) read(iun,'(a)') line
              end if
          end if
      end do
 
      write(*,300) num_sat, nav_file(1:trimlen(nav_file))
 300  format('* ', i5,' satellites found in ',a)
      write(*,310) start_nav, end_nav
      if( end_nav-start_nav.lt.0.25 ) then
          write(*,305) end_nav, start_nav
 305      format('Short duration nav file: Times ',2F20.3,
     .           ' increasing duration by 0.25 days')
          start_nav = start_nav - 0.125
          end_nav = end_nav + 0.125
      end if
 310  format('Entries start at MJD ',f20.3,' and end at JD ',f20.3) 
      start_sp3 = (int(start_nav*96.d0)-15)/96.d0
      sp3_refmjd = int(start_sp3)
 
      max_fnd_prn = 0
      do i = 1, num_sat
         j = 0 
C        do t = start_sp3, end_nav+15.d0/24.d0, 15.d0/1440.d0
         do it = 0,nint((end_nav+15.d0/24.d0-start_sp3)/(15.d0/1440.d0))
            t = start_sp3 + it*15.d0/1440.d0
            j = j + 1
            call  neph_to_xyz(t, prn_sp3(i), 'E')
            do k = 1, 3
               sp3_xyz(k,prn_sp3(i),j) = svs_xyz(k,prn_sp3(i)) 
            end do
            sp3_clk(prn_sp3(i),j) = af0(i) + 
     .                           af1(i)*(t-toe_jd(i))*86400.d0 
            sp3_time(j) = t
            sp3_sec(j)  = (t-sp3_refmjd)*86400.d0
         end do
         num_sp3 = j
         if ( prn_sp3(i).gt.max_fnd_prn ) max_fnd_prn = prn_sp3(i)
C        do j = 1, num_sp3
C           write(*,420) prn_sp3(i), (sp3_xyz(k,prn_sp3(i),j), k = 1,3),
C    .                   sp3_clk(prn_sp3(i),j), sp3_time(j), sp3_sec(j)
C420        format('P ',i2.2,3F14.3,f14.6,f20.6,f12.4)
C        end do          
      end do 


*     Save num_sat as largest PRN
      write(*,520 ) num_sat, max_fnd_prn
 520  format('Found ',i4,' satellites with max_prn ',i4)
      num_sat = max_fnd_prn     
 
      if( num_sat.eq.0 ) stop ' SVPOS: No satellites found'
      close(iun)
      
      return
      end
 
