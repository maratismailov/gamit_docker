
      program excite_pmu
 
*     Program to compare tables of polar motion UT1 series.  The
*     runstring is:
*     % excite_pmu <output excite>  <in primary> <UT1_def> <SD def>
*
*     where  <output excite> is the name of the excitation file
*            <in primary> is the name of the file for the input
*                 series.  Must be uniformly spaced.
*            <UT1_def> is the definition of UT1 in the input series
*                 (R-denoted regularized)
*            <SD def> indicates that series still have short period
*                 terms (tidal) that should be removed.

      include 'excite_pmu.h'
 
*   date_start(5), date_stop(5) - Start and stop calender
*               - dates (generic)
*   trimlen     - Length of string
 
      integer*4 date_start(5), date_stop(5), trimlen
 
*   sectag      - Seconds tag of JD
 
 
      real*8 sectag
 
****  Start program.  Start decoding the runstring upto the end of
*     of primary file entries.  The rest will be decoded later.
 
      write(*,100)
 100  format(/' EXCITE_PMU: Compare Polar Motion/UT1 excitations',
     .        ' functions from oberved data',/)
 
      call get_cp_runstring
 
*     Now read the primary data set.  Tell user what is going on.
      write(*,120) in_pr(1:trimlen(in_pr))
 120  format(' Starting to read primary data set in ',a)
 
      call read_primary
 
*     Tell user about primary date set.
      call jd_to_ymdhms(pr_start+0.0001d0, date_start, sectag)
      call jd_to_ymdhms(pr_stop+0.0001d0,  date_stop,  sectag)
      write(*,200) num_pr, date_start, date_stop
 200  format(' Primary data set has ',i4,' entries',/,
     .       ' Start date ',i4,2('/',i2),1x,i2,':',i2, 5x,
     .       ' Stop date  ',i4,2('/',i2),1x,i2,':',i2)

      call gen_exite
 
****  Thats all
      end
 
CTITLE GET_CP_RUNSTRING
 
      subroutine get_cp_runstring
 
*     This routine will read the first part of the runstring up to
*     the output primary file name.  The rest of the run string will
*     be read later
 
      include 'excite_pmu.h'
 
* LOCAL VARIABLES
 
*   rcpar   - Runstring reader
*   len_run - Length of runstring
 
 
      integer*4 rcpar, len_run
 
*   runstring   - Temporary storage of runstring
 
      character*40 runstring
 
****  Get the name of the output smoothed file
      len_run = rcpar(1, out_ex)
      if( len_run.le.0 ) then
          call proper_runstring('excite_pmu.hlp','excite_pmu',1)
      end if

*     See if we should generate the smooth values of just use the
*     primary data set.
      runstring = out_ex
 
****  Get the primary input series
      len_run = rcpar(2, in_pr )
      if( len_run.le.0 ) then
          call proper_runstring('excite_pmu.hlp','compare_pmu',1)
      end if
 
*     Get the UT1 definition for primary input
      len_run = rcpar(3, UT1_pr)
      if( len_run.le.0 ) then
          call proper_runstring('excite_pmu.hlp','compare_pmu',1)
      end if
      call casefold(UT1_pr)
 
*     Get the SD definition for primary input
      len_run = rcpar(4, SD_pr)
      if( len_run.le.0 ) then
          SD_pr = 'N'
      end if
      call casefold(SD_pr)
****  Thats all
      return
      end
 
CTITLE READ_PRIMARY
 
      subroutine read_primary
 
*     Routine to read the primary UT1 series.  If UT1 is not
*     regularized, it is regularized here.  We leave the series
*     as either UT1-UTC or UT1-AT (if UT1-UTC need to worry about
*     the jumps, when we smooth the data.)
 
      include 'excite_pmu.h'
 
* LOCAL VARIABLES
 
*   ierr, jerr  - IOSTAT error
*   date(5)     - Date and hr:min of the value read
 
      integer*4 ierr, jerr, date(5)
 
*   sectag      - Seconds tag for jdate
*   dut1        - Short period UT1 contribution.
*   dx, dy      - Tidal contributions to x and y
 
      real*8 sectag, dut1, dx, dy
 
*   line    - Line read from input
 
      character*128 line
 
****  Open the input file
      open(100, file=in_pr, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',in_pr,1, 'read_primary')
 
****  Start reading the file.
      num_pr = 0
      sectag = 0
      do while ( ierr.eq.0 )
          read(100,'(a)',iostat=ierr) line
          write(*,'(a)') line(1:80)
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
*             Read the line
              num_pr = num_pr + 1
              read(line,*,iostat=jerr) date, pr_xp(num_pr),
     .            pr_xsig(num_pr), pr_yp(num_pr), pr_ysig(num_pr),
     .            pr_ut(num_pr)  , pr_usig(num_pr)
              call report_error('IOSTAT',jerr,'read',line,
     .            0, 'read_primary')
              write(*,*) date, num_pr, pr_ut(num_pr)
              if( jerr.eq.0 ) then
                  call ymdhms_to_jd( date, sectag, pr_jd(num_pr))
 
*                 See if we need to remove tidal UT1
                  if( ut1_pr(1:1).ne.'R' ) then
                      call short_period_ut1(pr_jd(num_pr), dut1)
                      pr_ut(num_pr) = pr_ut(num_pr) - dut1
                  end if
                  if( sd_pr(1:1).eq.'Y' ) then
                      call sd_comp(pr_jd(num_pr), dx, dy, dut1)
                      pr_xp(num_pr) = pr_xp(num_pr) - dx
                      pr_yp(num_pr) = pr_yp(num_pr) - dy
                      pr_ut(num_pr) = pr_ut(num_pr) - dut1
                  end if
              else
*                 Skip this entry
                  num_pr = num_pr - 1
              end if
*                     ! File read OK and blank in col 1.
          end if
*                     ! Looping over the input file.
      end do
 
****  Save the start and stop times
      pr_start = nint(pr_jd(1)     -0.5d0) + 0.5d0
      pr_stop  = nint(pr_jd(num_pr)-0.5d0) + 0.5d0
      write(*,*) pr_start, pr_stop, pr_jd(1), pr_jd(num_pr)
     
 
****  Thats all
      close(100)
      return
      end 
 
CTITLE gen_exite
 
      subroutine gen_exite
 
*     This routine will generate the exciation values of the
*     series given the primary data set
 
      include 'excite_pmu.h'
 
* LOCAL VARIABLES
 
*   start_pt        - First point to be used in primary data
*               - for smoothing
*   num_pt      - Number of points to use
*   ns          - Smoothed point number
*   date_start(5), date_stop(5) - Start and stop calender
*               - dates (generic)
*   iout        - Output unit number
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   i           - Loop counter
*   inf(3)      - Information table for Lagrangian interpolation
 
      integer*4 start_pt, num_pt, ns, date_start(5), date_stop(5),
     .    iout, ierr, trimlen, i

      real*8 inf(3)
 
*   sectag      - Seconds tag for JD
*   dut1        - Tide contribution to tide.
*   ut_tab(max_tab) - UT1 table to tab out discontinuties in
*               - UT1-UTC
*   px_tab(max_tab), py_tab(max_tab) - Tabular values for removing
*                 trend
*   jd          - Date of the smoothed value
 
      real*8  xp(2),yp(2),ut(2), sectag, dut1, 
     .    ut_tab(max_tab), jd,
     .    px_tab(max_tab), py_tab(max_tab), dt
      integer*4 ijd
 
***** Start looping over the times to be computed with the
*     smooth series
 
      start_pt = 1
      ns = 0
      dt = (pr_jd(2)-pr_jd(1))*2
 
C     do jd = pr_start, pr_stop, dt/2 
      do ijd = 0, nint((pr_stop-pr_start)/(dt/2))
          jd = pr_start + ijd * (dt/2)   
 
*****     Increment smoothed point number
          ns = ns + 1
 
*****     Find the range of data which stradles the data.  The limits
*         are at least 2 points on either sides or no more than
*         10 days
 
          call get_sample( jd, dt, start_pt, num_pt )
 
*         Now make a continuous UT1 series
 
          call make_ut1_tab( jd, pr_jd(start_pt), num_pt,
     .                    pr_ut(start_pt), ut_tab )
          do i = 1, num_pt
             px_tab(i) = pr_xp(start_pt+i-1)
             py_tab(i) = pr_yp(start_pt+i-1)
          end do
          if ( num_pt.lt.4 ) then
              write(*,101) jd, num_pt, pr_jd(start_pt)
 101          format('JD, NUM, START ',f12.2,1x,i4,f12.2)
          end if

 
*         Use Largrangain interpolation.
          inf(1) = pr_jd(start_pt)
          inf(2) = dt/2.d0
          inf(3) = num_pt
          call Lagrange_intp(inf,px_tab, jd, xp, 1)
          call Lagrange_intp(inf,py_tab, jd, yp, 1)
          call Lagrange_intp(inf,ut_tab, jd, ut, 1)

*         Now compute the excitations 	
          sm_xp(ns) = xp(1)
          sm_yp(ns) = yp(1)
          sm_ut(ns) = ut(1)	
          
          sm_xsig(ns) = xp(1) + yp(2)*433.d0/6.283185d0
          sm_ysig(ns) =-yp(1) + xp(2)*433.d0/6.283185d0
          sm_usig(ns) = ut(2)

          sm_jd(ns)   = jd
 
*         end loop getting smooth data
      end do
 
      num_sm = ns
 
****  Now write the smoothed values to file.
      iout = 200
      call open_lu(iout, out_ex, ierr, 'append' )
      call report_error('IOSTAT',ierr,'open', out_ex,0, 'gen_exite')
      if( ierr.ne.0 ) iout = 6
 
***** Write some header information.
      call jd_to_ymdhms(pr_start, date_start, sectag)
      call jd_to_ymdhms(pr_stop,  date_stop,  sectag)
      write(iout,100) in_pr(1:trimlen(in_pr)), 
     .               date_start, date_stop 
 100  format('* PMU Excitation from ',a,/,
     .       '* Start date ',i4,2('/',i2),1x,i2,':',i2, 5x,
     .       ' Stop date  ',i4,2('/',i2),1x,i2,':',i2,
     .       '*   Date            Xp       Xex       Yp      Yex   ',
     .       '     UT1           Lod',/,
     .       '*                    "        "        "        "   ',
     .       '      s             s' )
 
      do i = 1, ns
          call jd_to_ymdhms( sm_jd(i)+0.0001d0, date_stop, sectag)
          if( ut1_sm(1:1).ne.'R' ) then
              call short_period_ut1( sm_jd(i), dut1 )
              ut(1) = sm_ut(i) + dut1
          else
              ut(1) = sm_ut(i)
          end if
 
*         Now write out values
          write(iout, 200) date_stop, sm_xp(i), sm_xsig(i),
     .            sm_yp(i), sm_ysig(i), ut(1), sm_usig(i)
  200     format(1x,i4,4i3,1x,2(f9.5,1x,f8.5,1x),1x,f10.6,1x,f9.6)
      end do
 
****  Thats all.
      return
      end
 
CTITLE GET_SAMPLE
 
      subroutine get_sample( jd, dt, start_pt, num_pt )
 
*     Routine to get the range of data to be smoothed.  We start at
*     start_pt and search in both directions.
 
      include 'excite_pmu.h'
 
*   start_pt        - Place to start searching, and to return
*               - value in.
*   num_pt      - number of points.
 
      integer*4 start_pt, num_pt
 
*   jd          - JD of point we want to interpolate
*   dt          - Time spacing to seach over
 
      real*8 jd, dt
 
* LOCAL VARIABLES
 
*   i,j         - Looping values finding range
 
 
      integer*4 i,j
 
****  Start at start_pt.  We know that we do not need to go back
*     from this point.
 
      i = start_pt
      do while( jd-pr_jd(i).gt. dt    .and. i.lt.num_pr )
          i = i + 1
      end do
 
*     OK, we have found the first point within 10 days.  Now find
*     the last point
      j = i
      do while ( pr_jd(j)-jd.le. dt   .and. j.lt.num_pr )
          j = j + 1
      end do
 
****  Now see is we have enough points
      if( j-i.lt.4 ) then
          if( j.le.num_pr-2) then
              j = j + 2
          else
              j = num_pr
          end if
 
          if( i.gt.2 ) then
              i = i - 2
          else
              i = 1
          end if
      end if
 
****  Now save values and return
      start_pt = i
      num_pt = j - i + 1
 
****  Thats all
      return
      end
 
CTITLE MAKE_UT1_TAB
 
      subroutine make_ut1_tab( jd, pr_jd, num_pt, pr_ut, ut_tab )
 
*     Routine to make a continuous table of UT1 (relative to the
*     point nearest the value at jd.
 
*   num_pt  - Number of points being considered
 
      integer*4 num_pt
 
*   jd      - JD of evalution point
*   pr_jd(num_pt)   - List of JD's in primary series
*   pr_ut(num_pt)   - Primary series values of UT
*   ut_tab(num_pt)  - Continuous UT1 series.
 
      real*8 jd, pr_jd(num_pt), pr_ut(num_pt), ut_tab(num_pt)
 
* LOCAL VARIABLES
 
*   i       - Loop counters
*   ic      - index of pr_jd closest to jd.
 
      integer*4 i, ic
 
*   djd     - Differenece in jd
*   secs    - number of seconds difference between ut at
*           - each time and ut neasrest jd.
 
      real*8 djd, secs
 
****  First find closest point to jd
      djd = 1.d20
      ic  = 1
      do i = 1, num_pt
          if( abs(jd-pr_jd(i)).lt.djd ) then
              djd = abs(jd-pr_jd(i))
              ic = i
          end if
      end do
 
****  Now make the table so that all values are closest to the
*     UT1 value at this time.
 
      do i = 1, num_pt
          secs = nint( pr_ut(i)-pr_ut(ic) )
          ut_tab(i) = pr_ut(i) - secs
      end do
 
****  Thats all
      return
      end
 

