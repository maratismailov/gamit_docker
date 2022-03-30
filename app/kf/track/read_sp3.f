CTITLE READ_SP3
 
      subroutine read_sp3

      implicit none
 
*     Thuis routine will read sp3 epehemeris file. All ephemeris
*     entries are read
 
      include '../includes/const_param.h'
      include 'track_com.h'
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   i           - Loop counter
*   rp          - Read Prn number
*   trimlen     - length of string
*   ll          - Length of line string
*   date(5)     - ymdhm
*   num_nosp3clk - Number of missing SP3 clock values 
 
      integer*4 ierr, i, j,k,  trimlen,  date(5), pn, num_nosp3clk
      integer*4 iun   ! Unit number for SP3 file
 
*   sectag      - Seconds tag in date
      real*8 sectag
 
*   still_header    - Indicates that we are sill in the header
*                   - part of file.
*   have            - Set true if we already have this satellite
      logical still_header
 
*   line            - line read from file
      character*256 line
 
*   cr - Carriage return (for handling dos files)
      character*1 cr

* MOD TAH 060302: Updated for SP3 file versions
      character*1 sp3_ver
 
      cr = char(13)
 
****  Open the sp3_file
      iun = 102
      open(iun, file=sp3_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',sp3_file,1,
     .            'svsp3/read_sp3')
 
*     Loop over the file reading the ephemeris entries.  First clear
*     all of the PRN numbers so we know when a PRN has been read

      do i = 1, max_sat
          prns(i) = 0
      end do

****  Clear the SP3 XYZ entries.  This was we can later tell if we read a
*     value or not
      do i = 1, max_eph_eps
         do j = 1, max_prn
            do k = 1,3
               sp3_xyz(k,j,i) = 0.d0
            end do
         end do
      end do
 
      still_header = .true.

* MOD TAH 060302: Get SP3 version
      read(iun,'(a)', iostat=ierr) line
      sp3_ver = line(2:2)
      if( sp3_ver .eq. 'c' ) then
         write(*,'(a)') 'SP3c version file'
      endif 
 
      do while ( still_header )
          read(iun,'(a)', iostat=ierr) line
          if( line(1:1).eq.'*' ) still_header = .false.
          if( ierr.ne.0 ) then
              write(*,120) sp3_file(1:trimlen(sp3_file))
 120          format(' Error reading ',a,' before end of ',
     .               'header found')
              stop 'svsp3: Reading sp3 file'
           end if
      end do
 
      num_sat = 0
      num_sp3 = 1
      do i = 1, max_sat
         sp3_clk(i,1) = 0.d0
      end do


* MOD TAH 060302: Start at character 4 instead of 2      
      read(line(4:),*) date, sectag
      call ymdhms_to_mjd( date, sectag, sp3_time(num_sp3))
      sp3_refmjd = int(sp3_time(num_sp3)+1.d-6)
      sp3_sec(num_sp3) = date(4)*3600.0d0 + date(5)*60.d0 +
     .                   sectag
       
*     Now start reading the entries
      num_nosp3clk = 0
      do while ( ierr.eq.0 )
 
*         Read down the elements of the satellites
          read(iun,'(a)',iostat=ierr) line
          if( line(1:3).eq.'EOF') then
C              Do not set error so that we can read through concatinated
C              files.
C              ierr = -1
* MOD TAH 060302: Only get GPS satellites
          else if (line(1:2).eq.'P ' .or.line(1:2).eq.'PG' ) then

* MOD TAH 060302: Start at character 3 instead of 2       
             read(line(3:),*) pn,
     .              sp3_xyz(1,pn,num_sp3),
     .              sp3_xyz(2,pn,num_sp3),
     .              sp3_xyz(3,pn,num_sp3),
     .              sp3_clk(pn,num_sp3)
              sp3_xyz(1,pn,num_sp3) = sp3_xyz(1,pn,num_sp3)*1000.d0
              sp3_xyz(2,pn,num_sp3) = sp3_xyz(2,pn,num_sp3)*1000.d0
              sp3_xyz(3,pn,num_sp3) = sp3_xyz(3,pn,num_sp3)*1000.d0

*             Check to make sure we have valid clock
              if( sp3_clk(pn,num_sp3).eq.999999.999999d0 ) then
                  num_nosp3clk = num_nosp3clk + 1
                  if( 1.eq.2 )
     .            write(*,250) pn, num_sp3, 
     .                         sp3_clk(pn,num_sp3-1)*1.d6
 250              format('**WARNING** No clock estimate for PRN',
     .                   I2.2,' at SP3 Epoch ',i5,' Using ',F10.3,
     .                   ' usec')
                   
                  sp3_clk(pn,num_sp3) = sp3_clk(pn,num_sp3-1)*1.d6
              end if
*             Convert clock to seconds
              sp3_clk(pn,num_sp3)   = sp3_clk(pn,num_sp3)*1.d-6
              prns(pn) = pn
 
              
              num_sat = max(num_sat,pn)
          else if (line(1:1).eq.'*' ) then
              num_sp3 = num_sp3 + 1
              read(line(2:),*) date, sectag
              call ymdhms_to_mjd( date, sectag, sp3_time(num_sp3))
              sp3_sec(num_sp3) = date(4)*3600.0d0 + date(5)*60.d0 +
     .                   sectag + 
     .                   int(sp3_time(num_sp3)-sp3_refmjd)*86400.d0

*             Clear next set of values
              do i = 1, max_sat
                 sp3_clk(i,num_sp3) = 0.d0
              end do   
          end if
      end do
      
       
      write(*,300) num_sat, num_sp3, sp3_file(1:trimlen(sp3_file)),
     .             num_nosp3clk
  300 format('* ', i5,' satellites, at ',i4,' epochs found in ',a,/,
     .       '* ', i5,' entries had no clock values')
 
      if( num_sat.eq.0 ) stop 'READ_SP3: No satellites found'
 
      close(iun)
 
      call get_svs_blk
      return
      end

CTITLE GET_SVS_BLK

      subroutine get_svs_blk

      implicit none

*     Routine to get the block numbers for the PRNs in the SP3 file.  
*     tries to read svnav.dat.  If it can't find it uses, values valid
*     in July 2003.

      include 'track_com.h'

* LOCAL VARIABLES

      integer*4 i, ierr, jerr, iun
      integer*4 prn, svn, blk, date(5), chan
      integer*4 prn_blk_072003(32)
      integer*4 trimlen
      real*4 svnav_ver  ! Version of svnav.dat (0.0 or 2.0 for GNNS)

      real*8 sec, lmjd

      character*256 line
      character*256 home_dir, svnav
      character*22 blk_name


*     Block numbers: 1 Blk 1; 2 Blk II, 3 Blk IIA, 4 Blk IIR
      data prn_blk_072003 /  3, 2, 3, 3, 3, 3, 3, 3, 3, 3,
     .                      4, 1, 4, 4, 2, 4, 2, 4, 2, 4,
     .                      4, 2, 3, 3, 3, 3, 3, 4, 3, 3,
     .                      3, 3 / 

****  Try to open local copy of svnav.dat
      svnav = 'svnav.dat'
      iun = 103
      open(iun,file=svnav, iostat=ierr, status='old')
      if ( ierr.ne.0 ) then   ! Try gg tables
          call getenv('HOME',home_dir)
          svnav = home_dir(1:max(1,trimlen(home_dir))) //
     .                   '/EDDY_LOCAL/gg/tables/svnav.dat'
          open(iun,file=svnav, iostat=ierr, status='old')
      end if
      sec = 0.0d0
      svnav_ver = 0.0
      if( ierr.eq.0 ) then   ! Open OK, so read file
*         MOD TAH 160809: See if new version of svnav.day
          read(iun,'(a)',iostat=ierr ) line
          if( line(1:9).eq.'svnav.dat' ) then
              read(line(19:25),*,iostat=jerr) svnav_ver 
              write(*,110) svnav_ver, jerr
 110          format('Found Ver ',F6.1,' svnav.dat. Jerr ',I4)
          endif
          if( svnav_ver.eq. 0.0 ) then
             do while ( ierr.eq.0 )
                read(iun,'(a)',iostat=ierr ) line
                read(line,120, iostat=jerr) prn, blk, date
 120            format(i2,6x,i1,28x,i4,1x,4i2)
                if( ierr.eq.0 .and. jerr.eq.0 .and.
     .              prn.gt.0 .and. prn.le.32 ) then
                    call ymdhms_to_mjd(date,sec, lmjd)
                    if( lmjd.le.sp3_refmjd ) then
                        prn_blk(prn) = blk
                    end if
                end if
             end do
          else
*            Read GNSS version which contains GPS plus other systems. 
             do while ( ierr.eq.0 )
                read(iun,'(a)',iostat=ierr ) line
*               Currently read only G line
                if( line(2:2).eq.'G' ) then
                   date(2) = 1   ! Date is given as DOY
                   read(line,140, iostat=jerr) svn, prn, chan, blk_name, 
     .                                         date(1),date(3:5)
 140               format(5x,I3,1x,I3,1x,I3,2x,a20,31x,i4,1x,i3,1x,i2,
     .                    1x,i2) 
!G     1   4   0  BLOCK I                 453800.     U      0.1999  1978  53  0  0  1985 199  0  0 
                   if( ierr.eq.0 .and. jerr.eq.0 ) then
                      call ymdhms_to_mjd(date,sec, lmjd)
                      call name_to_blk(+1, blk_name, blk)
                      if( lmjd.le.sp3_refmjd ) then
                          call name_to_blk(+1, blk_name, blk)
                          prn_blk(prn) = blk
                       end if
                   endif
                endif
             enddo
          endif  
      else     ! No svnav.dat file, so assign values
          write(*,220) 
 220      format('No SVNAV.DAT file found. Using July 2003 Blk numbers')
          do i = 1, num_sat
             prn_blk(i) = prn_blk_072003(i)
          end do
      end if

***** Report the values
      write(*,320) (i,prn_blk(i),i=1,num_sat)
 320  format(('Blk Numbers: ',8(:,1x,'PRN',i2.2,1x,i2)))

****  Thats all
      close(iun)
      return
      end
          

 
CTITLE CLEAN_SP3_CLK

      subroutine clean_sp3_clk
      
      implicit none

*     Routine to clean the 99999's from the SP3 clocks by linear
*     interpolation

      include 'track_com.h'
      
* LOCAL VARIABLES

* nb, nf  - Indices of god clock before and after bad clock
* i,j     - Loop over satelllites and epochs

      integer*4 nb, nf, i, j
      
* dt      - TIme difference from back point
* dt_sp3, dclk_sp3 - Interpolation values
      
      real*8    dt, dt_sp3, dclk_sp3
      
* bood, fgood - Indicates good clock found in back and forward
*               directions
      
      logical bgood, fgood
      

***** Loop over all the satellites accross all the sp3 epochs
      do i = 1, num_sat
          do j = 1, num_sp3
             if( sp3_clk(i,j).gt.99999.d-6 ) then
             
*                This clock epoch needs interpolation
                 nb = j
                 bgood = .false.
                 do while ( .not.bgood .and. nb.gt.1 )
                     nb = nb - 1
                     if( sp3_clk(i,nb).lt.99999.d-6 ) then
                         bgood = .true.
                     end if
                 end do
                 if( .not.bgood ) then
                     write(*,110) i,j
 110                 format('WARNING: Cant find first good clock',
     .                      ' for PRN ',i2.2,' at epoch ',i4)
                 end if
                 
****             now do foward search
                 nf = j
                 fgood = .false.
                 do while ( .not.fgood .and. nf.lt.num_sp3 )
                     nf = nf + 1
                     if( sp3_clk(i,nf).lt.99999.d-6 ) then
                         fgood = .true.
                     end if
                 end do
                 if( .not.fgood ) then
                     write(*,120) i,j
 120                 format('WARNING: Cant find last good clock',
     .                      ' for PRN ',i2.2,' at epoch ',i4)
                 end if
                                    
****             Now do interpolation
                 if( bgood .and. fgood ) then
                     dt_sp3 = sp3_time(nf) - sp3_time(nb)
                     dclk_sp3 = sp3_clk(i,nf) - sp3_clk(i,nb)
                     dt = sp3_time(j) - sp3_time(nb)
                     sp3_clk(i,j) =  sp3_clk(i,nb) +  
     .                              (dclk_sp3/dt_sp3)*dt
                     write(*,200) i, j,  sp3_clk(i,j)*1.d6
 200                 format(' Interpolating clock for PRN    ',i2.2,
     .                      ' at epoch ',i4,' Value ',f10.3,' usec')
                 else if ( bgood ) then
                     sp3_clk(i,j) =  sp3_clk(i,nb)
                     write(*,210) i, j,  sp3_clk(i,j)*1.d6
 210                 format(' Back estimate of clock for PRN ',i2.2,
     .                      ' at epoch ',i4,' Value ',f10.3,' usec')
                 else if ( fgood ) then
                     sp3_clk(i,j) =  sp3_clk(i,nf)
                     write(*,220) i, j,  sp3_clk(i,j)*1.d6
 220                 format(' Forw estimate of clock for PRN ',i2.2,
     .                      ' at epoch ',i4,' Value ',f10.3,' usec')
                 else
                     sp3_clk(i,j) = 0
                 end if
              end if
          end do
      end do

      
****  Thats all
      return
      end 
 
     
                      
