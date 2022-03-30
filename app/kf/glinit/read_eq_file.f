CTITLE READ_EQ_FILE
 
      subroutine read_eq_file

      implicit none
 
*     This routine will check to see of the earthquake list
*     file have been given and then go read the file.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/globk_cmds.h'
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   trimlen - Length of string
*   iel     - Command number from list of earthquake commands
*   indx    - Position in string that we have decoded to.
*   eol     - End of line pointer for removing ! an d# comments
 
 
 
      integer*4 ierr, trimlen, iel, indx, jndx, eol, i
 
*   buffer  - Line read from input file
 
 
      character*256 buffer
 
****  See if the file name has been given
      num_eq = 0
      num_renames = 0

      do i = 1, num_eqfiles 
* MOD TAH 160324: See if file name is RESET
          if( eq_inp_file(i)(1:6).eq.'RESET ' ) then
             eq_reset = .true.
          else
             open(106,file=eq_inp_file(i), status='old', iostat=ierr)
             call report_error('IOSTAT',ierr,'open',eq_inp_file(i),0,
     .                       'read_eq_file/glinit')
             if( ierr.ne.0 ) eq_inp_file(i) = ' '
 
***          Now start reading the file and extracting commands
 
             do while ( ierr.eq.0 )
                 read(106,'(a)',iostat=ierr) buffer

*                Clear any comments
                 eol = index( buffer,'!' )
                 if( eol.eq.0 ) eol = index(buffer,'#')
                 if( eol.gt.0 ) then
                     buffer(eol:) = ' '
                 end if

                 if( ierr.eq.0 .and. buffer(1:1).eq.' ' .and.
     .               trimlen(buffer).gt.0 )  then
 
*                    OK, Valid line for command.  Process:
                     indx = 0
                     call get_cmd(buffer, eq_commands, max_eq_cmd, iel,
     .                           indx)
* MOD TAH 021009: Remove any comments after !
                     jndx = index(buffer,'!')
                     if( jndx.gt.1 ) then
                         buffer = buffer(1:jndx-1)
                     end if 

*                    Process the command
                     call proc_eq_cmd(iel, indx, buffer )
 
                 end if
             end do
             if( ierr.eq.-1 ) close(106)
          end if
      end do
 
****  Close the input file
      return
      end
 
CTITLE PROC_EQ_CMD
 
      subroutine proc_eq_cmd(iel, indx, buffer )

      implicit none
 
*     Routine to process the earthquake commands
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/const_param.h'
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   trimlen - Length of string
*   iel     - Command number from list of earthquake commands
*   indx    - Position in string that we have decoded to.
 
*   jel     - Earthquake number
*   jndx    - Index when looking up code
*   date(5) - YMDHM of earthquake
*   jerr    - Error to test if hfile name passed in rename
*   indx_curr - Current value of indx incase no hfile name passed
 
      integer*4 ierr, trimlen, iel, indx, jel, jndx, date(5), j,
     .          jerr, indx_curr
 
*   geod(3) - Geodetic postion of earthquake (co-lat, long,
*           - ellipsiodal height.
*   sectag  - Seconds tag
*   lat, long       - Decimal degrees for epicenter
*   rad_km, dep_km  - Radius of influence and depth in
*           - km (converted to meters for internal use)
 
      real*8 geod(3), sectag, lat, long, rad_km, dep_km
 
*   buffer  - Line read from input file
 
      character*(*) buffer
 
*   test_code   - TEst code for Earth quake to see if
*           - duplicate or new event.
*	cdum    - Dummy string for readline
*   rn_test     - Test string to see if hfile name passed in rename
 
      character*8 test_code, cdum
      character*16 rn_test
 
****  See if there was an error
      if( iel.le.0 ) then
          call report_error('Get_command',ierr,'decod',buffer,0,
     .            'proc_eq_cmd')
          RETURN
      end if
 
****  Process the command:
      goto (  100,  200,  300,  400,  500,  600, 700, 800, 900 ) iel
 
*     EQ_DEF: Define the parameters of a specifiic Earth quake.
*         Command must be given before the other commands
*         referring to this earthquake
 100  continue
 
*         Get the Earthquake code
          call getword(buffer, test_code, indx)
          call casefold(test_code)
 
*         See if we have used this code before
          jndx = 1
          call get_cmd(test_code, eq_codes, num_eq, jel, jndx)
          if( jel.gt.0 ) then
              write(*,120) buffer(1:trimlen(buffer))
 120          format('**WARNING** Earthquake defined by: ',a,/,
     .               14x,'has been previously defined.  Using new',
     .                ' definition')
          else if ( jel.lt.-1 ) then
              write(*,140) buffer(1:trimlen(buffer))
 140          format('**WARNING** Earthquake defined by: ',a,/,
     .                14x,'is not UNIQUE code. Event being ignored')
          else
 
*             Add this event
              num_eq = num_eq + 1
              jel = num_eq
              if( num_eq.gt.max_eq ) then
                  write(*,160) max_eq, buffer(1:trimlen(buffer))
 160              format('**WARNING** Maxiumum number of earthquakes',
     .                 '(',i4,') Exceeded.  Ignoring ',/,a)
                  jel = 0
              end if
          end if
 
****      Now process the command if jel is greater than 0 indicating
*         that we have a valid earthquake
          if( jel.gt.0 ) then

              eq_rename(jel) = .false.
              eq_co_applied(jel) = .false.

              eq_codes(jel) = test_code
              read(buffer(indx:), *, iostat=ierr) lat, long, rad_km,
     .                dep_km, date
              geod(1) = pi/2 - lat*pi/180
              geod(2) = long*pi/180
              geod(3) = -dep_km*1000.d0
              call GEOD_to_XYZ( geod, eq_pos(1,jel) )
 
              eq_rad(jel) = rad_km*1000.d0
              eq_depth(jel) = dep_km*1000.d0
 
              sectag = 0.d0
              call ymdhms_to_jd( date, sectag, eq_epoch(jel))

*             Initialize the values associated with the earthquake
              do j = 1, 6
                 eq_apr_coseismic(j,jel) = 0.d0
                 eq_mar_post(j,jel) = 0.d0
                 eq_mar_pre(j,jel) = 0.d0
                 eq_log_sig(j,jel) = 0.d0
              end do

          end if
      RETURN
 
****  EQ_RENAM: Rename option for Earthquake.
 200  continue
 
*         Get the Earthquake code
          call getword(buffer, test_code, indx)
          call casefold(test_code)
          jndx = 1
          call get_cmd(test_code, eq_codes, num_eq, jel, jndx)
          if( jel.gt.0 ) then
              eq_rename(jel) = .true.
          else
              write(*,210) buffer(1:trimlen(buffer))
 210      format('**WARNING** Earthquake code not found. Command:',
     .                a)
          end if
          RETURN
 
****  EQ_COSEI: Set the coseimic displacement apriori sigmas
 300  continue
 
*         Get the Earthquake code
          call getword(buffer, test_code, indx)
          call casefold(test_code)
          jndx = 1
          call get_cmd(test_code, eq_codes, num_eq, jel, jndx)
          if( jel.gt.0 ) then
 
*             Now read the rest of the line.  sigmas in NEU (m) and
*             spatial sigmas (NEU at 1 earthquake depth).
              read(buffer(indx:),*,iostat=ierr) 
     .                (eq_apr_coseismic(j,jel), j = 1,6)
              call report_error('IOSTAT',ierr,'read',buffer,0,
     .                'PROC_EQ_CMD/COSEISMIC')
          else
              write(*,210) buffer(1:trimlen(buffer))
          end if
          RETURN
 
****  EQ_POST : Set the postseismic markov process noise in
*                 m**2/yr in NEU . Read in mm**2/day and coverted
*                 here. 
 400  continue
 
*         Get the Earthquake code
          call getword(buffer, test_code, indx)
          call casefold(test_code)
          jndx = 1
          call get_cmd(test_code, eq_codes, num_eq, jel, jndx)
          if( jel.gt.0 ) then
 
*             Now read the rest of the line. Markov process values
*             (fixed and spatially dependent)
              read(buffer(indx:),*,iostat=ierr) eq_dur(2,jel),
     .                       (eq_mar_post(j,jel), j = 1,6)
              call report_error('IOSTAT',ierr,'read',buffer,0,
     .                'PROC_EQ_CMD/POSTSEISMIC')
*             Convert from mm**2/day to m**2/yr
              do j =1,6
                 eq_mar_post(j,jel) = eq_mar_post(j,jel)*365.25d-6
              end do
          else
              write(*,210) buffer(1:trimlen(buffer))
          end if
          RETURN
 
****  EQ_PRE  : Set the preseismic markov process noise in
*                 m**2/yr in NEU
 500  continue
 
*         Get the Earthquake code
          call getword(buffer, test_code, indx)
          call casefold(test_code)
          jndx = 1
          call get_cmd(test_code, eq_codes, num_eq, jel, jndx)
          if( jel.gt.0 ) then
 
*             Now read the rest of the line. Markov process values
*             (fixed and spatially dependent)
              read(buffer(indx:),*,iostat=ierr) eq_dur(1,jel),
     .                     (eq_mar_pre(j,jel), j = 1,6)
              call report_error('IOSTAT',ierr,'read',buffer,0,
     .                'PROC_EQ_CMD/PRESEISMIC')
*             Convert from mm**2/day to m**2/yr
              do j =1,6
                 eq_mar_pre(j,jel) = eq_mar_pre(j,jel)*365.25d-6
              end do
          else
              write(*,210) buffer(1:trimlen(buffer))
          end if
          RETURN
 
***** RENAME : General feature for renaming sites
 600  CONTINUE
 
*         Make sure we do not have too many renames already.
          if( num_renames+1.gt.max_rn ) then
              write(*,610) max_rn, (buffer(1:trimlen(buffer)))
 610          format('**WARNING** Too many renames specified.',
     .            ' Maximum allowed ',i5,/,
     .            ' Ignoring ',a)
              RETURN
          end if
 
          num_renames = num_renames + 1
          call getword(buffer, rn_codes(1,num_renames),indx)
          call getword(buffer, rn_codes(2,num_renames),indx)
          call casefold(rn_codes(1,num_renames))
          call casefold(rn_codes(2,num_renames))

* MOD TAH 971112: Check to see if a name restriction has been included
*         in the line
          indx_curr = indx
          call getword(buffer, rn_test, indx )
*         See if number or string
          call check_num( rn_test, jerr )
          if( jerr.ne.0 ) then
              rn_hfiles(num_renames) = rn_test
          else
              rn_hfiles(num_renames) = ' '
              indx = indx_curr
          end if

 
*         Now start reading the oprional parts of the command.
*         Start time for the rename (ymdhms)
          do j = 1,5
              call read_line(buffer,indx,'I4',ierr,date(j),cdum)
              if( ierr.ne.0 ) date(j) = 1
          end do
*         If there is an error, default the start to 1900
          if( ierr.ne.0 ) date(1) = 1900
 
*         Process the starting data
          sectag = 0.0d0
          call ymdhms_to_jd(date, sectag, rn_times(1,num_renames))
 
****      Get the end time:
          do j = 1,5
              call read_line(buffer,indx,'I4',ierr,date(j),cdum)
              if( ierr.ne.0 ) date(j) = 1
          end do
*         If there is an error, default the end to 2100
          if( ierr.ne.0 ) date(1) = 2100
 
*         Process the starting data
          sectag = 0.0d0
          call ymdhms_to_jd(date, sectag, rn_times(2,num_renames))
 
****      See if there is a displacement associated with this
*         rename:
          do j = 1,3
              call read_line(buffer,indx,'R8',ierr,
     .            rn_dpos(j,num_renames),cdum)
              if( ierr.ne.0 ) rn_dpos(j,num_renames) = 0.d0
          end do
 
****      Now see what frame this is (XYZ/NEU)
          call read_line(buffer,indx,'CH',ierr,sectag,
     .            rn_types(num_renames))
          call casefold(rn_types(num_renames))
          if( ierr.ne.0 ) rn_types(num_renames) = 'XYZ'
          RETURN

***** EQ_LOG : Read the log estimation information
 700  CONTINUE

*         Get the Earthquake code
          call getword(buffer, test_code, indx)
          call casefold(test_code)
          jndx = 1
          call get_cmd(test_code, eq_codes, num_eq, jel, jndx)
          if( jel.gt.0 ) then
 
*             Now read the rest of the line. Markov process values
*             (fixed and spatially dependent)
              read(buffer(indx:),*,iostat=ierr) eq_log_tau(jel),
     .                     (eq_log_sig(j,jel), j = 1,6)
              call report_error('IOSTAT',ierr,'read',buffer,0,
     .                'PROC_EQ_CMD/EQ_LOG')
          else
              write(*,210) buffer(1:trimlen(buffer))
          end if
          RETURN

****  BREAK : Read the time and optional name and hfile code for
*         break in the time series
 800  CONTINUE

*         Start decoding the line and make this look like a regular
*         rename command except that the new name contains a ? character
*         that is updated in new routine sort_unf_ren
          if( num_renames+1.gt.max_rn ) then
              write(*,610) max_rn, (buffer(1:trimlen(buffer)))
              RETURN
          end if

*         Get site name
          num_renames = num_renames + 1
          call getword(buffer, rn_codes(1,num_renames),indx)
          call casefold(rn_codes(1,num_renames))

*         Now set the rn_codes(2,num_renames) to denote that
*         this is a break command
          rn_codes(2,num_renames) = rn_codes(1,num_renames)(1:4) //
     .          '_?PS' 

*         See if h-file restriction has been added
          indx_curr = indx
          call getword(buffer, rn_test, indx )
*         See if number or string
          call check_num( rn_test, jerr )
          if( jerr.ne.0 ) then
              rn_hfiles(num_renames) = rn_test
          else
              rn_hfiles(num_renames) = ' '
              indx = indx_curr
          end if

*         Now get the date
          call multiread(buffer,indx,'I4',ierr,date,cdum,5)
*         Process the starting data
          sectag = 0.0d0
          call ymdhms_to_jd(date, sectag, rn_times(1,num_renames))
*         set an open ended end-date
          date(1) = 2100
          do j = 2,5
             date(j) = 1
          end do
          call ymdhms_to_jd(date, sectag, rn_times(2,num_renames))

****      Finally see if open re-name string has been passed
          call getword(buffer, rn_test, indx )
          call casefold(rn_test)
          jel = trimlen(rn_test)
          if( jel.gt.0 ) then
              rn_codes(2,num_renames)(5:5+jel) = '_' // rn_test(1:jel)
          endif

****      Set the offsets to zero for this rename
          do j = 1,3
              rn_dpos(j,num_renames) = 0.d0
          end do
          rn_types(num_renames) = 'NEU'
          RETURN

****  RESET: Causes all name extents to be reset to _GPS
 900  CONTINUE
          eq_reset = .true.
          RETURN

 
***** Thats all
      END

CTITLE SORT_UNF_RN

      subroutine sort_unf_rn

      implicit none

****  Routine to time sort, name sort and unify the rename entries.
*     In this routine we adjust any open ended entries.

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* LOCAL VARIABLES
      integer*4 i,j, n, sel, icc, jend
      character*1 cc
      character*4 ccode 
      logical kbit, fnd
  

      if ( num_renames.le.1 ) RETURN

****  Now sort by site name (based on 4-char code)
      do i = 1, num_renames-1
         do j = i+1,num_renames
            if( rn_codes(1,j)(1:4).lt.rn_codes(1,i)(1:4) ) then
                call swap_rn(i,j)
            endif
         end do
      end do

****  First start by time sorting all of the entries for sites
*     that have the names.

****  Run through the list getting groups of the same site
      sel = 1  ! Start element
      ccode = rn_codes(1,sel)(1:4)
      n = 1

      do while ( n.lt. num_renames )
          n = n + 1
*         Time sort with new name or end of list
          if( rn_codes(1,n)(1:4).ne. ccode(1:4) .or. 
     .        n .eq. num_renames ) then
* MOD TAH 161227: Moved increment of n code up so that sort
*             will check all entries.

*****         If we reached the end of the list increment n
*             so that loop below will get last entry (normally
*             n is the first entry of the next station).
              if ( n.eq.num_renames ) n=n+1

*             Name has changed, time sort current group
              do i = sel, n-2
                 do j = i+1,n-1
                   if( rn_times(1,j).lt. rn_times(1,i)  ) then
                       call swap_rn(i,j)
                   endif
                end do
              end do


*****         Now that we have a time sort list unify the names
*             and open ended times
              cc = '0'  ! Current code: Work here to Z
              do i = sel, n-1   ! Loop over entries at this site
*                 See if we have ? char
                  if ( rn_codes(2,i)(6:6).eq.'?' ) then
*                     We need a new name here.  Scan to see what we
*                     can use
                      fnd = .false.
                      do while ( .not. fnd )
                         icc = ichar(cc)+1
                         if( icc.eq.58 ) icc = 65
                         if( icc.gt.90 ) then
                             call report_stat('FATAL','GLOBK',
     .                            'sort_unf_rn',rn_codes(1,i),
     .                            'Too many BREAKS at site',0)
                         endif
                         cc = char(icc)
                         fnd = .true.
                         do j = sel,n-1
                            if( rn_codes(2,j)(6:6).eq.cc ) fnd = .false.
                         end do
                      end do
*                     Assign the name
                      rn_codes(2,i)(6:6) = cc(1:1)
                  end if
*                 Now unify the times
                  if( rn_times(2,i).gt. 2488000.0d0 ) then  ! Past 2099
*                     Find the next enntry that is not an XPS,XCL
                      j = i
                      jend = 0
                      do while ( j.lt.n-1 )
                         j = j + 1
* MOD TAH 140319: Verify that the site code has not changed while 
*                        looking for excluded data. 
                         if( rn_codes(2,j)(6:6).ne.'X' .and.
     .                       rn_codes(1,j)(1:4).eq.ccode(1:4) ) then
                             jend = j
                             j = n
                         endif
                      end do
                      if( jend.gt.0 ) rn_times(2,i) =  rn_times(1,jend)
                  endif
              end do  

*             Save current (only save if not at end of list)
              if( n.le. num_renames ) then
                 ccode = rn_codes(1,n)(1:4)
                 sel = n
              endif 
          end if
      end do

****  Thats all
      return
      end

CTITLE SWAP_RN

      subroutine swap_rn(i,j)

      implicit none

*     Routine to swap two rename entries.
      
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

* PASSED
      integer*4 i, j  ! Entries to be swapped 


* LOCAL VARIABLES

      character*8 tm_codes(2), tm_types
      character*16 tm_hfiles
      real*8 tm_times(2), tm_dpos(3)

      integer*4 k

****  Save the ith entry
      tm_codes(1) = rn_codes(1,i)
      tm_codes(2) = rn_codes(2,i)
      tm_hfiles   = rn_hfiles(i)
      tm_types    = rn_types(i)
      do k = 1,2
         tm_times(k) = rn_times(k,i)
      enddo
      do k = 1,3
         tm_dpos(k) = rn_dpos(k,i)
      enddo

****  Replace i with jth entry
      rn_codes(1,i) = rn_codes(1,j)
      rn_codes(2,i) = rn_codes(2,j)
      rn_hfiles(i)  = rn_hfiles(j)
      rn_types(i)   = rn_types(j)
      do k = 1,2
         rn_times(k,i) = rn_times(k,j)
      enddo
      do k = 1,3
         rn_dpos(k,i) = rn_dpos(k,j)
      enddo

****  Now save the tmp entries to j entry
      rn_codes(1,j) = tm_codes(1)
      rn_codes(2,j) = tm_codes(2)
      rn_hfiles(j)  = tm_hfiles
      rn_types(j)   = tm_types
      do k = 1,2
         rn_times(k,j) = tm_times(k)
      enddo
      do k = 1,3
         rn_dpos(k,j) = tm_dpos(k)
      enddo

****  Thats all
      return
      end

