      subroutine update_sp3

*     Routine to update the SP3 file entries if needed at each new time 
*     that a measurement becomes available.

      implicit none


      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* LOCAL
      logical update   ! Set true if update needed

****  See if we need to update
      if( num_sp3.eq.0 ) then   ! No files at all yet, so get entries.
*        Read the data we need before the current time
         call get_newsp3('start')
      end if

*     See if need to extent current entries.  Since ultrarapids
*     are upated every 6 hours with a latency of 3 hours, see if the 
*     end of the SP3 entries suggests we can update. When there is less
*     21-hrs of data left, we can upate.  (epoch+24 - 3hrs latency)
*      print *,'SP3 Times ',sp3_time(num_sp3),RT_MJD_obs(sblk), sblk
      if ( sp3_time(num_sp3) - RT_MJD_obs(sblk).lt. 21/24.d0 ) then
*        If update is needed then get new files
         call get_newsp3('end')
      end if

****  Thats all
      return
      end

      subroutine get_newsp3(option)

*     Routine to read new SP3 files as needed.  Two options are used:
*     start -- initial SP3 file.  This option simply reads the first file
*     end   -- This option tries to read the next file (due to timing it
*              may not be available and either an older file or the existing
*              entries will be used.

      implicit none

      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED
      character*(*) option   ! Option: start or end

* LOCAL
      real*8 obs_mjd  ! Observation time used to generate file names

      integer*4  nf   ! File number found (1--final, 2--rapid, 2>Ultrarapid)
     .,      ierr     ! IOSTST error
     .,      uns      ! Unit number for SP3 file.
     .,      j        ! Loop counter
     .,      trimlen  ! String length

      character*256 sp3_names(10)  ! Names of possible SP3 files to use
                      ! in order of perference (final, rapid and upto
                      ! ultra-rapids with degress if latency (in real-time
                      ! first four files are in the future).
      logical found   ! Set true when able to find SP3 file


****  See which option so we know what file name to generate
C      if( option(1:1).eq.'s' ) then

*        Start option:  Get time 120-minutes before currunt and see
*        which SP3 file will have that name. (14-point interpolator)
         obs_mjd = RT_MJD_obs(sblk)
*        Generate file names that could contain this epoch of data
         call gen_sp3name(obs_mjd, sp3_names, option)

*        Now find which sp3 file is available and use the best
*        one.
         found = .false.
         uns = 80
         nf = 0
         do while ( .not.found .and. nf.lt.10 )
            nf = nf + 1
            open(uns,file=sp3_names(nf),status='old',iostat=ierr)
            if( ierr.eq.0 ) then  ! File found
               found = .true.
            end if
         enddo

*        See if found
         if( .not.found ) then
*           We have a fatal error: No SP3 orbit file can be found
            write(*,120) obs_mjd
     .,          (j,sp3_names(j)(1:trimlen(sp3_names(j))),j=1,10)
 120        format('FATAL No SP3 files are MJD ',F15.6,/,
     .             'Choices were: ',/,
     .              10(i3,1x,a,/))
            call report_stat('FATAL','TrackRT','get_newsp3',
     .           ' ','No SP3 files found',0)
         endif

         if( sp3_file.ne.sp3_names(nf) ) then
             write(*,140) sp3_names(nf)(1:trimlen(sp3_names(nf)))
 140         format('Reading SP3 file ',a)
             call report_stat('STATUS','TrackRT','get_newsp3',
     .            sp3_names(nf),'opening',nf)
             call read_sp3RT(uns, nf)
             sp3_file = sp3_names(nf)
         else
             close(uns)
         endif

         
C      else
*        We already have sp3 orbit information and we 

C      endif

      return
      end 

      subroutine eval_endt( update )


      implicit none


      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED 
      logical update    ! Set true if the need to update SP3 files.


      end

      subroutine gen_sp3name(obs_mjd, sp3_names, option) 

      implicit none


      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED
      real*8 obs_mjd  ! Time of observation
      character*(*) sp3_names(10)  ! Possible file names (final rapid and
                       ! 8 ultra rapids) listed in preference order
      character*(*) option  ! Either 'start' or 'end' depending on if 
                       ! this is first sp3 file or extending the end

* LOCAL
      real*8 early_mjd ! Earliest MJD needed (2-hrs before current time)
     .,      ultra_mjd ! MJD for nominal Ultra rapid file
     .,      test_mjd, mjd       ! Generic MJD

      integer*4 imjd, jmjd   ! Integer MJD
     .,      week, dow, hr  ! GPS week, doy of week and hr (0,6,12,18)
     .,      ultra_hrs ! Nominal hr for ultra file for this time (0,6,12,18)
     .,      trimlen   ! Length of string
     .,      j         ! Loop counter
     .,      lsp3d     ! Length of SP3 file name.

****  Get the week and day
      if( option(1:1).eq.'s' ) then
         test_mjd = obs_mjd - 2/24.d0  ! Go back 2-hour
      else  
         test_mjd = obs_mjd + 2/24.d0  ! Go forward 2-hour
      endif
  
      imjd = int(test_mjd)
      week = (imjd - 44244)/7
      dow = imjd - (week*7+44244)
     
      if( trimlen(sp3_dir).eq.0 ) sp3_dir = './'
      lsp3d = trimlen(sp3_dir)
      if( sp3_dir(lsp3d:lsp3d).ne.'/' ) then
          sp3_dir(lsp3d+1:) = '/'
          lsp3d = lsp3d + 1
      endif
      
      write(sp3_names(1),120) sp3_dir(1:lsp3d), sp3_root(1:2), 
     .                        sp3_root(3:3), week, dow
 120  format(a,a2,a1,I4.4,I1,'.sp3')
      write(sp3_names(2),120) sp3_dir(1:lsp3d), sp3_root(1:2), 'r', 
     .                        week, dow

****  Now generate the names for the ultra-rapid files.  These file arrive 
*     every 6 hrs at 3:00, 9:00, 15:00, 21:00 GPST with predictions starting
*     at 0:00, 6:00, 12:00 and 18:00 hrs.  Normally we go back 3-hrs and the
*     latest 6-hr block should be available.
      ultra_mjd = obs_mjd
      imjd = int(ultra_mjd)
      ultra_hrs = int((ultra_mjd -imjd)*4.d0)*6
*     Advance MJD by 24 hrs to generate name of file where this time 
*     will be in the observed segment (only available when not in realtime).
      imjd = imjd + 1
      do j = 1,8 
          mjd = imjd + ultra_hrs/24.d0 - (j-1)*0.25d0
          jmjd = int(mjd)
          week = (jmjd - 44244)/7
          dow = jmjd - (week*7+44244)
          hr = nint((mjd - (week*7+44244) - dow)*24.d0) 
          write(sp3_names(j+2),160) sp3_dir(1:lsp3d),sp3_root(1:2), 'u', 
     .          week, dow, hr
 160      format(a,a2,a1,I4.4,I1,'_',i2.2,'.sp3')
      enddo

***** Thats all
C     write(*,220) obs_mjd
C    .,      (j,sp3_names(j)(1:trimlen(sp3_names(j))),j=1,10)
C220  format('SP3 Time ',F15.6,/,10(i3,1x,a,/))

      return
      end

CTITLE READ_SP3RT
 
      subroutine read_sp3RT( uns, nf )

      implicit none
 
*     This routine will read sp3 epehemeris file connected to 
*     unit uns.  The file number nf indicates if ultrarapid or 
*     rapid/final to be read.  (Ultras are 48-hrs long and replace
*     the whole ephemeris (with check on older values), where as
*     rapid/final move the eariler version down. (48-hrs of SP3
*     are present in both case). 
 
      include '../includes/const_param.h'
      include 'trackRT.h'

* PASSED
      integer*4 uns  ! Unit number of SP3 file
     .,         nf   ! File number, 1 and 2 are final/rapid; >2 ultra rapid
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   i           - Loop counter
*   rp          - Read Prn number
*   trimlen     - length of string
*   ll          - Length of line string
*   date(5)     - ymdhm
*   num_nosp3clk - Number of missing SP3 clock values 
 
      integer*4 ierr, i, j,k,  trimlen,  date(5), pn, num_nosp3clk
      integer*4 sp3_sig_code(3)  ! 2**(code) mm sigma on orbits
      integer*4 isp3   ! MJD offset from ref MJD
      integer*4 off_ind  ! Offset in SP3 files, used to compute
                       ! RMS difference when file updated.
      integer*4 old_ref  ! Old reference time when SP3 updated.
      integer*4 dsec     ! Seconds offset to account for different
                         ! day in SP3 file.
      integer*4 num      ! Number of epochs in SP3 file (not used)
 
      real*8 sectag     ! Seconds tag in date
     .,      sp3_start  ! Start MJD of SP3 file
 
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

***** Read first line of file to get file type and the
*     start time of this file.
      read(uns,'(a)', iostat=ierr) line
      sp3_ver = line(2:2)
C     if( sp3_ver .eq. 'c' ) then
C        write(*,'(a)') 'SP3c version file'
C     endif 
*     Get the start time of the file
      read(line,100,iostat=ierr) date, sectag, num
 100  format(3x,I5,4(1x,I2),1x,f11.8,1x,I7)
      call ymdhms_to_mjd(date, sectag, sp3_start)

****  Before we start: If r/f then move the end of the current
*     file back 24-hrs; 
      if( nf.le.2 .and. num_sp3.gt.96 ) then
*         Move values back by 24-hrs (96 15-minute epochs)
          old_ref = sp3_refmjd 
          sp3_refmjd = nint(sp3_time(97))
          dsec = (sp3_refmjd-old_ref)*86400  ! Seconds change
          do i = 97, num_sp3
             sp3_time(i-96) = sp3_time(i)
             sp3_sec(i-96) = sp3_sec(i) - dsec
             do j = 1, max_sat
                sp3_clk(j,i-96) = sp3_clk(j,i)
                do k = 1,3
                   sp3_xyz(k,j,i-96) = sp3_xyz(k,j,i)
                   sp3_sig(k,j,i-96) = sp3_sig(k,j,i)
                end do
             end do
          end do
          num_sp3 = num_sp3 - 96  ! Move counter back
       end if
*      See we need an offset index for the ultra-rapid comparison
       off_ind = 0
       if( nf.gt.2 .and. num_sp3.gt.0 ) then
*         See if the can find start time in current entries
          do i = 1, num_sp3
             if( sp3_time(i).eq.sp3_start ) then
                 off_ind = i
             end if
          end do
          num_sp3 = 0   ! Reset back to zero
       end if

***   Now start reading file 
      still_header = .true.

      do while ( still_header )
          read(uns,'(a)', iostat=ierr) line
          if( line(1:1).eq.'*' ) still_header = .false.
          if( ierr.ne.0 ) then
              write(*,120) sp3_file(1:trimlen(sp3_file))
 120          format(' Error reading ',a,' before end of ',
     .               'header found')
              stop 'svsp3: Reading sp3 file'
           end if
      end do

* MOD TAH 060302: Start at character 4 instead of 2
      num_sp3 = num_sp3 + 1   ! First value      
      read(line(4:),*) date, sectag
      call ymdhms_to_mjd( date, sectag, sp3_time(num_sp3))
      if( nf.gt.2 .or. num_sp3.eq.1 ) then  ! Set reference time
          sp3_refmjd = int(sp3_time(num_sp3)+1.d-6)
      end if
      isp3 = int(sp3_time(num_sp3))-sp3_refmjd
      sp3_sec(num_sp3) = date(4)*3600.0d0 + date(5)*60.d0 +
     .                   sectag + isp3*86400
      do i = 1, max_sat
         sp3_clk(i,num_sp3) = 0.d0
         do j = 1, 3
            sp3_sig(j,i,num_sp3) = 1.e6   ! Large error
         enddo
      end do
      
*     Now start reading the entries
      num_nosp3clk = 0
      do while ( ierr.eq.0 )
 
*         Read down the elements of the satellites
          read(uns,'(a)',iostat=ierr) line
          if( line(1:3).eq.'EOF') then
C              Do not set error so that we can read through concatinated
C              files.
C              ierr = -1
* MOD TAH 060302: Only get GPS satellites
          else if (line(1:2).eq.'P ' .or.line(1:2).eq.'PG' ) then

* MOD TAH 060302: Start at character 3 instead of 2   
             read(line,220) pn,
     .             (sp3_xyz(j,pn,num_sp3),j=1,3),
     .              sp3_clk(pn,num_sp3),
     .             (sp3_sig_code(j),j=1,3)
 220         format(2x,I2.2,4f14.6,3(1x,i2),1x,i3)
*            Convert to meeters for position and sigma codes
             do j = 1,3
                sp3_xyz(j,pn,num_sp3) = sp3_xyz(j,pn,num_sp3)*1000.d0
                if( sp3_sig_code(j).ne.0 ) then
                   sp3_sig(j,pn,num_sp3) = 2**sp3_sig_code(j)/1000.d0
                else  ! Set error 1-meter when not known
                   sp3_sig(j,pn,num_sp3) = 100.d0
                end if
             end do
 
*             Check to make sure we have valid clock
* CODE Replaced with call to clean_sp3_clk
             if( sp3_clk(pn,num_sp3).eq.999999.999999d0 ) then
                 num_nosp3clk = num_nosp3clk + 1
C                 if( 1.eq.2 )
C    .            write(*,250) pn, num_sp3, 
C    .                         sp3_clk(pn,num_sp3-1)*1.d6
C250              format('**WARNING** No clock estimate for PRN',
C    .                   I2.2,' at SP3 Epoch ',i5,' Using ',F10.3,
C    .                   ' usec')
C                 if( num_sp3.gt.1 ) then 
C                    sp3_clk(pn,num_sp3) = sp3_clk(pn,num_sp3-1)*1.d6
C                 endif 
              end if

*             Convert clock to seconds
              sp3_clk(pn,num_sp3)   = sp3_clk(pn,num_sp3)*1.d-6
              prns(pn) = pn
              num_sat = max(num_sat,pn)

****          See if we can compare value with old values (for
*             ulra-rapid orbits)
              if( off_ind.gt.0 ) then
C                  call check_sp3(pn, num_sp3, off_ind, stats)
              endif
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
                 do j = 1, 3
                    sp3_sig(j,i,num_sp3) = 1.e6   ! Large error
                 enddo
             end do   
          end if
      end do
      
       
      if( debug(9).gt.0 ) 
     .write(*,300) num_sat, num_sp3, sp3_file(1:trimlen(sp3_file)),
     .             num_nosp3clk
  300 format('* ', i5,' max PRN with ',i4,' epochs found in ',a,/,
     .       '* ', i5,' entries had no clock values')
 
      if( num_sat.eq.0 ) stop 'READ_SP3: No satellites found'
 
      close(uns)

      call clean_sp3_clkRT

      call get_svs_blkRT

      close(uns)

      return
      end

CTITLE GET_SVS_BLKRT

      subroutine get_svs_blkRT

      implicit none

*     Routine to get the block numbers for the PRNs in the SP3 file.  
*     tries to read svnav.dat.  If it can't find it uses, values valid
*     in July 2003.

      include 'trackRT.h'

* LOCAL VARIABLES

      integer*4 i, ierr, jerr, unt
      integer*4 prn, blk, date(5)
      integer*4 prn_blk_072003(32)
      integer*4 trimlen

      real*8 sec, lmjd

      character*256 line
      character*256 home_dir, svnav


*     Block numbers: 1 Blk 1; 2 Blk II, 3 Blk IIA, 4 Blk IIR
      data prn_blk_072003 /  3, 2, 3, 3, 3, 3, 3, 3, 3, 3,
     .                      4, 1, 4, 4, 2, 4, 2, 4, 2, 4,
     .                      4, 2, 3, 3, 3, 3, 3, 4, 3, 3,
     .                      3, 3 / 

****  Try to open local copy of svnav.dat
      svnav = 'svnav.dat'
      unt = 82
      open(unt,file=svnav, iostat=ierr, status='old')
      if ( ierr.ne.0 ) then   ! Try gg tables
          call getenv('HOME',home_dir)
          svnav = home_dir(1:max(1,trimlen(home_dir))) //
     .                   '/gg/tables/svnav.dat'
          open(unt,file=svnav, iostat=ierr, status='old')
      end if
      sec = 0.0d0
      if( ierr.eq.0 ) then   ! Open OK, so read file
          do while ( ierr.eq.0 )
             read(unt,'(a)',iostat=ierr ) line
             read(line,120, iostat=jerr) prn, blk, date
 120         format(i2,6x,i1,28x,i4,1x,4i2)
             if( ierr.eq.0 .and. jerr.eq.0 .and.
     .           prn.gt.0 .and. prn.le.32 ) then
                 call ymdhms_to_mjd(date,sec, lmjd)
                 if( lmjd.le.sp3_refmjd ) then
                     prn_blk(prn) = blk
                 end if
             end if
          end do
      else     ! No svnav.dat file, so assign values
          write(*,220) 
 220      format('No SVNAV.DAT file found. Using July 2003 Blk numbers')
          do i = 1, num_sat
             prn_blk(i) = prn_blk_072003(i)
          end do
      end if

***** Report the values
      if( debug(9).gt.0 ) 
     .write(*,320) (i,prn_blk(i),i=1,num_sat)
 320  format(('Blk Numbers: ',8(:,1x,'PRN',i2.2,1x,i2)))

****  Thats all
      close(unt)
      return
      end
          

 
CTITLE CLEAN_SP3_CLKRT

      subroutine clean_sp3_clkRT
      
      implicit none

*     Routine to clean the 99999's from the SP3 clocks by linear
*     interpolation

      include 'trackRT.h'
      
* LOCAL VARIABLES

* nb, nf  - Indices of god clock before and after bad clock
* i,j     - Loop over satelllites and epochs

      integer*4 nb, nf, i, j,k 
      
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
                 if( .not.bgood .and. debug(1).gt.0 ) then
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
                 if( .not.fgood .and. debug(1).gt.0 ) then
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
                     if( debug(9).gt.0 )
     .               write(*,200) i, j,  sp3_clk(i,j)*1.d6
 200                 format(' Interpolating clock for PRN    ',i2.2,
     .                      ' at epoch ',i4,' Value ',f10.3,' usec')
                 else if ( bgood ) then
                     sp3_clk(i,j) =  sp3_clk(i,nb)
                     if( debug(9).gt.0 )
     .               write(*,210) i, j,  sp3_clk(i,j)*1.d6
 210                 format(' Back estimate of clock for PRN ',i2.2,
     .                      ' at epoch ',i4,' Value ',f10.3,' usec')
                 else if ( fgood ) then
                     sp3_clk(i,j) =  sp3_clk(i,nf)
                     if( debug(9).gt.0 )
     .               write(*,220) i, j,  sp3_clk(i,j)*1.d6
 220                 format(' Forward Estimate for       PRN ',i2.2,
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
 
     
                      
    




