 
      program compare_pmu
 
*     Program to compare tables of polar motion UT1 series.  The
*     runstring is:
*     % compare_pmu <output smooth> <spacing> <smoothing> <UT1_def> \
*                   <in primary>   <UT1_def> <out primary> \
*                   <in secondary> <UT1_def> <out secondary> \
*                        ""           ""         ""
*     where <output smooth> is the names of the output smoothed filed
*                 (which can be used in update_pmu) If the name
*                 USEPR is entered then the primary series will be
*                 used as the smoothed series.  Spacing, smmothing,
*                 and UT1_def will then have no effect.
*           <spacing> is the spacing of the smoothed entries (days)
*           <smoothing> is a multiplier to the FWHM computed to pass
*                 through the primary data. A value >1 will tend to
*                 smooth the primary series.
*           <UT1_def> is the definition for UT1 or UT1R in the output
*                 (UT1R will be output if UT1_def=R, otherwise UT1
*                  will be output)
*           <in primary> is the input primary series which will be
*                  smoothed (IRIS formar)
*           <UT1_def> is UT1 definition in primary series.
*           <out primary> is name of the output file for the
*                    differences (primary-smoothed) (mas for output)
*         Repeating group of elements
*         | <in secondary> is the input secondary series which will be
*         |        smoothed (IRIS format)
*         | <UT1_def> is UT1 definition in primary series.
*         | <out secondary> is name of the output file for the
*         |         differences (secondary-smoothed) (mas for output)
*         The last group can be repeated for as many secondary files
*         as desired.
*
 
      include 'compare_pmu.h'
 
*   date_start(5), date_stop(5) - Start and stop calender
*               - dates (generic)
*   trimlen     - Length of string
*   nr          - Runstring element at start of next secondary
*               - file group (actually -1 from this value)
*   iout        - Unit numbers for output
 
      integer*4 date_start(5), date_stop(5), trimlen, nr, iout
 
*   finished        - Indicates that we exhausted the supply of
*               - runstring parameters for secondary files.
 
      logical finished
 
*   sectag      - Seconds tag of JD
 
 
      real*8 sectag
 
****  Start program.  Start decoding the runstring upto the end of
*     of primary file entries.  The rest will be decoded later.
 
      write(*,100)
 100  format(/' COMPARE_PMU: compare tables of polar motion UT1',
     .        ' series',/)
 
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
 
****  Now generate the smoothed series.  If the output file name is
*     USEPR then we have simply copied the smoothed over.  Note in this
*     the smoothed series may not be regularly spaced (so don't try to
*     take advantage of this)
 
      if( gen_sm ) then
          call gen_smooth
      else
          call copy_smooth
      end if
 
***** We have now generated and written the smoothed solution.  Now form
*     the difference between smoothed and primary if they are different
 
      if( gen_sm ) then
          call diff_pr( iout )
          call out_stats(6, prx_stats, 'Primary X Pole','mas')
          call out_stats(6, pry_stats, 'Primary Y Pole','mas')
          call out_stats(6, pru_stats, 'Primary UT1-AT','mts')
          call out_stats(iout, prx_stats, 'Primary X Pole','mas')
          call out_stats(iout, pry_stats, 'Primary Y Pole','mas')
          call out_stats(iout, pru_stats, 'Primary UT1-AT','mts')
      end if
 
***** Now we finish looping over the secondary files.
      finished = .false.
      nr = 7
      do while ( .not.finished )
 
          call get_sec_run( nr, finished )
          if( .not.finished ) then
              call diff_sec(iout)
              call write_sec(iout)
              call out_stats(6, secx_stats, 'Secondary X Pole','mas')
              call out_stats(6, secy_stats, 'Secondary Y Pole','mas')
              call out_stats(6, secu_stats, 'Secondary UT1-AT','mts')
              call out_stats(iout, secx_stats, 'Secondary X Pole','mas')
              call out_stats(iout, secy_stats, 'Secondary Y Pole','mas')
              call out_stats(iout, secu_stats, 'Secondary UT1-AT','mts')
          end if
      end do
 
****  Thats all
      end
 
CTITLE GET_CP_RUNSTRING
 
      subroutine get_cp_runstring
 
*     This routine will read the first part of the runstring up to
*     the output primary file name.  The rest of the run string will
*     be read later
 
      include 'compare_pmu.h'
 
* LOCAL VARIABLES
 
*   rcpar   - Runstring reader
*   ierr    - IOSTAT error
*   len_run - Length of runstring
 
 
      integer*4 rcpar, ierr, len_run
 
*   runstring   - Temporary storage of runstring
 
      character*40 runstring
 
****  Get the name of the output smoothed file
      len_run = rcpar(1, out_sm)
      if( len_run.le.0 ) then
          call proper_runstring('compare_pmu.hlp','compare_pmu',1)
      end if

*     See if we should generate the smooth values of just use the
*     primary data set.
      runstring = out_sm
      call casefold(runstring)
      if( runstring(1:5).eq.'USEPR' ) then
          gen_sm = .false.
      else
          gen_sm = .true.
      end if 

 
*     Get spacing (days)
      len_run = rcpar(2,runstring)
      if( len_run.le.0 ) then
          call proper_runstring('compare_pmu.hlp','compare_pmu',1)
      else
          read(runstring,*, iostat=ierr) sm_step
          call report_error('IOSTAT',ierr,'decod',runstring,
     .            1, 'get_cp_runstring')
 
      end if
*     Get smoothing multiplier (>1 will smooth data)
      len_run = rcpar(3,runstring)
      if( len_run.le.0 ) then
          call proper_runstring('compare_pmu.hlp','compare_pmu',1)
      else
          read(runstring,*, iostat=ierr) sm_fwhm_scale
          call report_error('IOSTAT',ierr,'decod',runstring,
     .            1, 'get_cp_runstring')
 
      end if
 
*     Get the UT1 definition for smooth output
      len_run = rcpar(4, UT1_sm)
      if( len_run.le.0 ) then
          call proper_runstring('compare_pmu.hlp','compare_pmu',1)
      end if
      call casefold(UT1_sm)
 
****  Get the primary input series
      len_run = rcpar(5, in_pr )
      if( len_run.le.0 ) then
          call proper_runstring('compare_pmu.hlp','compare_pmu',1)
      end if
 
*     Get the UT1 definition for primary input
      len_run = rcpar(6, UT1_pr)
      if( len_run.le.0 ) then
          call proper_runstring('compare_pmu.hlp','compare_pmu',1)
      end if
      call casefold(UT1_pr)
 
****  Get the primary output residuals series file
      len_run = rcpar(7, out_pr )
      if( len_run.le.0 ) then
          call proper_runstring('compare_pmu.hlp','compare_pmu',1)
      end if
 
****  Thats all
      return
      end
 
CTITLE READ_PRIMARY
 
      subroutine read_primary
 
*     Routine to read the primary UT1 series.  If UT1 is not
*     regularized, it is regularized here.  We leave the series
*     as either UT1-UTC or UT1-AT (if UT1-UTC need to worry about
*     the jumps, when we smooth the data.)
 
      include 'compare_pmu.h'
 
* LOCAL VARIABLES
 
*   ierr, jerr  - IOSTAT error
*   date(5)     - Date and hr:min of the value read
 
      integer*4 ierr, jerr, date(5)
 
*   sectag      - Seconds tag for jdate
*   dut1        - Short period UT1 contribution.
 
      real*8 sectag, dut1, dx, dy
 
*   line    - Line read from input
 
      character*128 line
 
****  Open the input file
      open(100, file=in_pr, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',in_pr,1, 'read_primary')
 
****  Start reading the file.
      num_pr = 0
      do while ( ierr.eq.0 )
          read(100,'(a)',iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
*             Read the line
              num_pr = num_pr + 1
              read(line,*,iostat=jerr) date, pr_xp(num_pr),
     .            pr_xsig(num_pr), pr_yp(num_pr), pr_ysig(num_pr),
     .            pr_ut(num_pr)  , pr_usig(num_pr)
              call report_error('IOSTAT',jerr,'read',line,
     .            0, 'read_primary')
              if( jerr.eq.0 ) then
                  call ymdhms_to_jd( date, sectag, pr_jd(num_pr))
 
*                 See if we need to remove tidal UT1
                  if( ut1_pr(1:1).ne.'R' ) then
                      call short_period_ut1(pr_jd(num_pr), dut1)
                      pr_ut(num_pr) = pr_ut(num_pr) - dut1
                  end if
                  if( ut1_pr(2:2).eq.'S' ) then
                      call sd_comp(pr_jd(num_pr), dx,dy,dut1)
                      pr_xp(num_pr) = pr_xp(num_pr) - dx/1000.d0
                      pr_yp(num_pr) = pr_yp(num_pr) - dy/1000.d0
                      pr_ut(num_pr) = pr_ut(num_pr) - dut1/1000.d0
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
 
****  Thats all
      close(100)
      return
      end 
 
CTITLE GEN_SMOOTH
 
      subroutine gen_smooth
 
*     This routine will generate the smoothed values of the
*     series given the primary data set
 
      include 'compare_pmu.h'
 
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
 
      integer*4 start_pt, num_pt, ns, date_start(5), date_stop(5),
     .    iout, ierr, trimlen, i
 
*   fwhmx, fwhmy, fwhmu - Full widths at half max for
*               - x, y pole and UT1. (days)
*   trendx(2), trendy(2), trendu(2) - Linear fit through
*               - the pr data (weighted by sigma and distance
*               - from central epoch)
*   xp,yp,ut        - Smoothed values for x,y pole and ut1
*   xps,yps,uts - Sigmas for smoothed values.
*   sectag      - Seconds tag for JD
*   dut1        - Tide contribution to tide.
*   ut_tab(max_tab) - UT1 table to tab out discontinuties in
*               - UT1-UTC
*   px_tab(max_tab), py_tab(max_tab) - Tabular values for removing
*                 trend
*   jd          - Date of the smoothed value
 
      real*8 fwhmx, fwhmy, fwhmu, trendx(2), trendy(2), trendu(2),
     .    xp,yp,ut, xps,yps,uts, sectag, dut1, ut_tab(max_tab), jd,
     .    px_tab(max_tab), py_tab(max_tab)
      integer*4 ijd
 
***** Start looping over the times to be computed with the
*     smooth series
 
      start_pt = 1
      ns = 0
 
C     do jd = pr_start, pr_stop+3*sm_step, sm_step
      do ijd = 0, nint(((pr_stop+3*sm_step)-pr_start)/sm_step)
          jd = pr_start + ijd*sm_step
 
*****     Increment smoothed point number
          ns = ns + 1
 
*****     Find the range of data which stradles the data.  The limits
*         are at least 2 points on either sides or no more than
*         10 days
 
          call get_sample( jd, start_pt, num_pt )
 
*         Now make a continuous UT1 series
 
          call make_ut1_tab( jd, pr_jd(start_pt), num_pt,
     .                    pr_ut(start_pt), ut_tab )
          do i = 1, num_pt
             px_tab(i) = pr_xp(start_pt+i-1)
             py_tab(i) = pr_yp(start_pt+i-1)
          end do

 
*****     Compute the FWHM of guassian filter to use for each of the
*         data types
          call wt_filter_fwhm(pr_jd(start_pt),pr_xsig(start_pt),
     .                      num_pt,fwhmx)
          call wt_filter_fwhm(pr_jd(start_pt),pr_ysig(start_pt),
     .                      num_pt,fwhmy)
          call wt_filter_fwhm(pr_jd(start_pt),pr_usig(start_pt),
     .                      num_pt,fwhmu)
 
*****     Modify the fwhm's by the scaling factors
          fwhmx = fwhmx*sm_fwhm_scale
          fwhmy = fwhmy*sm_fwhm_scale
          fwhmu = fwhmu*sm_fwhm_scale/2.d0
 
*****     Now remove the trend about the time to be considerd
          call remove_trend( px_tab, pr_jd(start_pt),
     .            pr_xsig(start_pt), num_pt, jd, jd, trendx)
          call remove_trend( py_tab, pr_jd(start_pt),
     .            pr_ysig(start_pt), num_pt, jd, jd, trendy)
          call remove_trend( ut_tab(1)      , pr_jd(start_pt),
     .            pr_usig(start_pt), num_pt, jd, jd, trendu)
 
*****     Now get smoothed value with Gaussian
          call wt_gauss_filter(pr_jd(start_pt), px_tab,
     .            pr_xsig(start_pt), num_pt, jd, xp, xps, fwhmx)
          call wt_gauss_filter(pr_jd(start_pt), py_tab,
     .            pr_ysig(start_pt), num_pt, jd, yp, yps, fwhmy)
          call wt_gauss_filter(pr_jd(start_pt), ut_tab,
     .            pr_usig(start_pt), num_pt, jd, ut, uts, fwhmu)
 
****      Now add back in trend (Only the offset need by added
*         since we referred linear fit to epcoh of smoothed point)
          sm_xp(ns) = xp + trendx(1)
          sm_yp(ns) = yp + trendy(1)
          sm_ut(ns) = ut + trendu(1)
*         Copy the sigmas and save date
          sm_xsig(ns) = xps
          sm_ysig(ns) = yps
          sm_usig(ns) = uts
          sm_jd(ns)   = jd

          sm_xfwhm(ns) = fwhmx
          sm_yfwhm(ns) = fwhmy
          sm_ufwhm(ns) = fwhmu
 
*         end loop getting smooth data
      end do
 
      num_sm = ns
 
****  Now write the smoothed values to file.
      iout = 200
      call open_lu(iout, out_sm, ierr, 'append' )
      call report_error('IOSTAT',ierr,'open', out_sm,0, 'gen_smooth')
      if( ierr.ne.0 ) iout = 6
 
***** Write some header information.
      call jd_to_ymdhms(pr_start, date_start, sectag)
      call jd_to_ymdhms(pr_stop,  date_stop,  sectag)
      write(iout,100) in_pr(1:trimlen(in_pr)), sm_step,
     .               date_start, date_stop, sm_fwhm_scale
 100  format('* Smoothed PMU from ',a,' with ',F4.2,' day spacing',/,
     .       '* Start date ',i4,2('/',i2),1x,i2,':',i2, 5x,
     .       ' Stop date  ',i4,2('/',i2),1x,i2,':',i2,/,
     .       '* Smoothing factor applied is ',f9.2)
 
      do i = 1, ns
          call jd_to_ymdhms( sm_jd(i)+0.0001d0, date_stop, sectag)
          if( ut1_sm(1:1).ne.'R' ) then
              call short_period_ut1( sm_jd(i), dut1 )
              ut = sm_ut(i) + dut1
          else
              ut = sm_ut(i)
          end if
 
*         Now write out values
          write(iout, 200) date_stop, sm_xp(i), sm_xsig(i),
     .            sm_yp(i), sm_ysig(i), ut, sm_usig(i),
     .            sm_xfwhm(i), sm_yfwhm(i), sm_ufwhm(i)
  200     format(1x,i4,4i3,1x,2(f9.5,1x,f8.5,1x),1x,f10.6,1x,f9.6,
     .           3f5.1)
      end do
 
****  Thats all.
      return
      end
 
CTITLE COPY_SMOOTH
 
      subroutine copy_smooth
 
*     Routine to copy the primary arrays to the smoothed arrays
 
      include 'compare_pmu.h'
 
*   i       - Loop counter
 
      integer*4 i
 
****  Simply copy all of the values over.
      do i = 1, num_pr
          sm_jd(i) = pr_jd(i)
          sm_xp(i) = pr_xp(i)
          sm_yp(i) = pr_yp(i)
          sm_ut(i) = pr_ut(i)
 
          sm_xsig(i) = pr_xsig(i)
          sm_ysig(i) = pr_ysig(i)
          sm_usig(i) = pr_usig(i)
      end do
 
      num_sm = num_pr
 
****  Thats all
      return
      end
 
 
CTITLE GET_SAMPLE
 
      subroutine get_sample( jd, start_pt, num_pt )
 
*     Routine to get the range of data to be smoothed.  We start at
*     start_pt and search in both directions.
 
      include 'compare_pmu.h'
 
*   start_pt        - Place to start searching, and to return
*               - value in.
*   num_pt      - number of points.
 
      integer*4 start_pt, num_pt
 
*   jd          - JD of point we want to interpolate
 
      real*8 jd
 
* LOCAL VARIABLES
 
*   i,j         - Looping values finding range
 
 
      integer*4 i,j
 
****  Start at start_pt.  We know that we do not need to go back
*     from this point.
 
      i = start_pt
      do while( jd-pr_jd(i).gt.sm_fwhm_scale*10.d0 .and. i.lt.num_pr )
          i = i + 1
      end do
 
*     OK, we have found the first point within 10 days.  Now find
*     the last point
      j = i
      do while ( pr_jd(j)-jd.le.sm_fwhm_scale*10.d0 .and. j.lt.num_pr )
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
 
CTITLE DIFF_PR
 
      subroutine diff_pr( iout )
 
*     Routine to difference the smoothed tables from the primary
*     data set and to accumulate the statistics of difference.
 
      include 'compare_pmu.h'
 
* LOCAL VARIABLES
 
*   next_pt     - next point to be used in interpolating
*               - the smoothed series
*   date_start(5), date_stop(5) - Start and stop calender
*               - dates (generic)
*   iout        - Output unit number
*   i           - Loop counter over primary data set
*   ierr        - IOSTAT error
*   trimlen     - Length of string
 
      integer*4 next_pt, date_start(5), date_stop(5), iout, i, ierr,
     .    trimlen
 
*   dxp,dyp,dut - differences for x,y pole and ut1 (mas,mts)
*   sxp, syp, sut   - Sigmas for differences (mas, mts)
*   dt          - Time difference between primary and
*               - first bounding smoothed data point
*   ut1(2)      - Two values of ut1 for interpolation.  (This is
*                 to account for possible leap seconds accross
*                 the boundary.
*   sectag      - Seconds tag for JD
 
      real*8 dxp,dyp,dut, sxp, syp, sut, dt, sectag, ut1(2)
 
****  Clear statistics
 
      call clear_stats( prx_stats )
      call clear_stats( pry_stats )
      call clear_stats( pru_stats )
 
****  Create the output difference file
      iout = 200
      call open_lu(iout, out_pr, ierr, 'append' )
      call report_error('IOSTAT',ierr,'open', out_pr, 0, 'diff_pr')
      if( ierr.ne.0 ) iout = 6
 
***** Write some header information.
      call jd_to_ymdhms(pr_start, date_start, sectag)
      call jd_to_ymdhms(pr_stop,  date_stop,  sectag)
      write(iout,100) in_pr(1:trimlen(in_pr)), sm_step,
     .               date_start, date_stop, sm_fwhm_scale
 100  format('* Residuals in ',a,' with ',F4.2,' day spacing',/,
     .       '* Start date ',i4,2('/',i2),1x,i2,':',i2, 5x,
     .       ' Stop date  ',i4,2('/',i2),1x,i2,':',i2,/,
     .       '* Smoothing factor applied is ',f6.3)
 
 
****  Loop over the primary data set.
      next_pt = 2
      do i = 1, num_pr
 
*         Find the first point in the smmothed data after the
*         primary JD
 
          do while ( pr_jd(i).gt.sm_jd(next_pt) .and.
     .              next_pt.lt.num_sm )
              next_pt = next_pt + 1
          end do
 
****      Now interpolate linearly betwwen the smoothed values
          dt = pr_jd(i) - sm_jd(next_pt-1)
 
          dxp = pr_xp(i) - ((sm_xp(next_pt)-sm_xp(next_pt-1))/
     .                    (sm_jd(next_pt)-sm_jd(next_pt-1)) )*dt
     .                   - sm_xp(next_pt-1)
          dyp = pr_yp(i) - ((sm_yp(next_pt)-sm_yp(next_pt-1))/
     .                    (sm_jd(next_pt)-sm_jd(next_pt-1)) )*dt
     .                   - sm_yp(next_pt-1)

*         dut = pr_ut(i) - ((sm_ut(next_pt)-sm_ut(next_pt-1))/
*    .                    (sm_jd(next_pt)-sm_jd(next_pt-1)) )*dt
*    .                   - sm_ut(next_pt-1)
          ut1(1) = sm_ut(next_pt-1)
          ut1(2) = sm_ut(next_pt) - 
     .                 nint(sm_ut(next_pt)-sm_ut(next_pt-1))

          dut  = pr_ut(i)  - ((ut1(2)-ut1(1))/
     .                      (sm_jd(next_pt)-sm_jd(next_pt-1)) )*dt
     .                     - ut1(1) 
 
****      Now convert to milliarc sec and millitime secs
          dxp = dxp * 1000.d0
          dyp = dyp * 1000.d0
*                                         ! Account for leapsecond
          dut = mod(dut,1.d0)* 1000.d0
 
*         Do sigmas
          sxp = pr_xsig(i) * 1000.d0
          syp = pr_ysig(i) * 1000.d0
          sut = pr_usig(i) * 1000.d0
 
*****     Now write values.
          call jd_to_ymdhms( pr_jd(i)+0.0001d0, date_stop, sectag)
          write(iout, 200) date_stop, dxp, sxp, dyp, syp, dut, sut
  200     format(1x,i4,4i3,1x,2(f9.3,1x,f8.3,1x),1x,f10.3,1x,f9.3)
 
*****     Accumulate statistics
 
          call accum_stats(dxp,sxp, prx_stats )
          call accum_stats(dyp,syp, pry_stats )
          call accum_stats(dut,sut, pru_stats )
      end do
 
****  Thats all
      return
      end
 
CTITLE CLEAR_STATS
 
      subroutine clear_stats ( stats )
 
*     routine initializes the statistics variables.
 
*   stats(4)        - Stats to be clears
      real*8 stats(4)
 
*    i          - Loop counter
 
      integer*4 i
 
****  Just loop to clear
      do i = 1,4
          stats(i) = 0
      end do
 
****  Thats all
      return
      end
 
CTITLE ACCUM_STATS
 
      subroutine accum_stats (dp,sp, stats )
 
*     Routine to accumulate statistics for pole positions
 
*   dp      - Differences
*   sp      - Sigma
*   stats(4)    - Accumulation array
 
      real*8 dp, sp, stats(4)
 
* LOCAL
 
*   wgh     - Weight (1/sigma**2)
 
      real*8 wgh
 
****  Compute weight for statistics
 
      if( sp.gt.0 ) then
          wgh = 1.d0/sp**2
      else
          RETURN
      end if
 
      stats(1) = stats(1) + dp*wgh
      stats(2) = stats(2) + 1.d0*wgh
      stats(3) = stats(3) + dp**2*wgh
      stats(4) = stats(4) + 1.d0
 
****  Thats all
      return
      end
 
CTITLE OUT_STATS
 
      subroutine out_stats ( unit, stats, title, units)
 
*     Routine to finish statistics and output values.
 
*   unit        - Unit number for output
 
      integer*4 unit
 
*   stats(4)        - The accumulated stats
 
      real*8 stats(4)
 
*   title       - title for stats output.
*   units       - Units to be written in the output
 
      character*(*) title, units
 
* LOCAL
 
*   wmean       - Weighted mean difference
*   nrms_mean   - NRMS scatter about mean
*   wrms_mean   - WRMS scatter about mean
 
 
      real*8 wmean, nrms_mean, wrms_mean
 
      if( stats(4).gt.1 ) then
c
          wmean = stats(1)/stats(2)
          nrms_mean = sqrt( (stats(3)-stats(1)*wmean)/(stats(4)-1) )
          wrms_mean = sqrt( stats(4)/stats(2) )*nrms_mean
c
*              ! only one data so we can not compute nrms scatter
      else
c
          wmean     = 0.d0
          nrms_mean = 1.d0
          wrms_mean = sqrt( stats(4)/stats(2) )*nrms_mean
c
      end if
 
*     Now write out the resulys
      write(unit,100) title, nint(stats(4)), wmean, units,
     .        wrms_mean, units, nrms_mean
 100  format('* For ',a,' statistics are from ',i4,' data',/,
     .        '* Weighted Mean ',f10.3,1x,a,' Wrms ',f10.3,1x,a,
     .        ' Nrms ',f8.3)
 
****  Thats all
      return
      end
 
CTITLE GET_SEC_RUN
 
      subroutine get_sec_run( nr, finished )
 
*     Routine to decode the secondary runstring to file in, out, and
*     UT definition
 
      include 'compare_pmu.h'
 
*   nr          - Current number in runstring
 
      integer*4 nr
 
*   finished        - Set true when we run out of runstring parameters
 
      logical finished
 
* LOCAL
 
*   rcpar   - Runstring reader
*   ierr    - IOSTAT error
*   len_run - Length of runstring
 
      integer*4 rcpar, len_run
 
****  Get the secondary input series
      nr = nr + 1
      len_run = rcpar(nr, in_sec )
      if( len_run.le.0 ) then
          finished = .true.
          RETURN
      end if
 
*     Get the UT1 definition for secondary input
      nr = nr + 1
      len_run = rcpar(nr, UT1_sec)
      if( len_run.le.0 ) then
          finished = .true.
          RETURN
      end if
      call casefold(UT1_sec)
 
****  Get the secondary output residuals series file
      nr = nr + 1
      len_run = rcpar(nr, out_sec )
      if( len_run.le.0 ) then
          finished = .true.
          RETURN
      end if
 
****  Thats all
      return
      end
 
CTITLE DIFF_SEC
 
      subroutine diff_sec(iout)
 
*     Routine to compute the differences from the secondary inputs.
*
      include 'compare_pmu.h'
 
* LOCAL VARIABLES
 
*   next_pt     - next point to be used in interpolating
*               - the smoothed series
*   date(5)     - Generic date array
*   iout        - Output unit number
*   i           - Loop counter over primary data set
*   ierr        - IOSTAT error
*   trimlen     - Length of string
 
      integer*4 next_pt, iout, ierr,  date(5), jerr
 
*   dxp,dyp,dut - differences for x,y pole and ut1 (mas,mts)
*   sxp, syp, sut   - Sigmas for differences (mas, mts)
*   dt          - Time difference between primary and
*               - first bounding smoothed data point
*   sectag      - Seconds tag for JD
*   dut1        _ Tide Ut1 contribtion

*   inf         - Polar/ut1 information for Lagrange interpolation.
*                 (jd, spacing and number)
*   wx(2), wy(2), wu(2) - Lagrangian interpolated values and rates
 
      real*8 dxp,dyp,dut, sxp, syp, sut, sectag, dut1, inf(3),
     .       wx(2), wy(2), wu(2), dx, dy

*   line        - Line read from input

      character*128 line
 
****  Clear statistics
 
      call clear_stats( secx_stats )
      call clear_stats( secy_stats )
      call clear_stats( secu_stats )

****  Set up for Lagrangian interpolation

      inf(1) = sm_jd(1)
      inf(2) = sm_jd(2)-sm_jd(1)
      inf(3) = num_sm
 
****  OPen the input file
 
      open(100,file=in_sec, status='old', iostat=ierr )
      call report_error('IOSTAT',ierr,'open',in_sec,0, 'diff_sec')
      if( ierr.ne.0 ) RETURN
 
 
****  Create the output difference file
      iout = 200

****  Loop over the secimary data set.
      next_pt = 2
      do while ( ierr.eq.0 )
 
          read(100,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
*             Decode the line
              read(line,*,iostat=jerr) date, sec_xp, sec_xsig, sec_yp,
     .                        sec_ysig, sec_ut, sec_usig
              call report_error('IOSTAT',jerr,'decod', line,
     .            0, 'read_secimary')
              if( jerr.ne.0 ) write(*,'(a)') line(1:70)
              if( jerr.eq.0 ) then
                  call ymdhms_to_jd( date, sectag, sec_jd)
 
*                 See if we need to remove tidal UT1
                  if( ut1_sec(1:1).ne.'R' ) then
                      call short_period_ut1(sec_jd, dut1)
                      sec_ut = sec_ut - dut1
                  end if
                  if( ut1_sec(2:2).eq.'S' ) then
                      call sd_comp(sec_jd, dx,dy,dut1)
                      sec_xp = sec_xp - dx/1000.d0
                      sec_yp = sec_yp - dy/1000.d0
                      sec_ut = sec_ut - dut1/1000.d0
                  end if
              end if
 
*             Find the first point in the smmothed data after the
*             secimary JD
 
****          Use Lagrangian interpolation
              call lagrange_intp(inf,sm_xp, sec_jd, wx, 1)
              call lagrange_intp(inf,sm_yp, sec_jd, wy, 1)
              call lagrange_intp(inf,sm_ut, sec_jd, wu, 1)

              dxp  = sec_xp - wx(1)
              dyp  = sec_yp - wy(1)
              dut  = sec_ut - wu(1)
 
****          Now convert to milliarc sec and millitime secs
              dxp = dxp * 1000.d0
              dyp = dyp * 1000.d0
*                                             ! Account for leapsecond
              dut = (dut - nint(dut)) * 1000.d0
 
*             Do sigmas
              sxp = sec_xsig * 1000.d0
              syp = sec_ysig * 1000.d0
              sut = sec_usig * 1000.d0
 
*****         Now accumulate statistics
              if( jerr.eq.0 ) then

*****             Accumulate statistics
 
                  call accum_stats(dxp,sxp, secx_stats )
                  call accum_stats(dyp,syp, secy_stats )
                  call accum_stats(dut,sut, secu_stats )
              end if
          end if
 
      end do
 
****  Thats all
      return
      end

CTITLE REMOVE_TREND

      subroutine remove_trend( values, dates, sigmas, num, start, mid,
     .                         trend )

*     Routine to remove a linear trend from values.  The slope estimate
*     is weighted by the sigmas and by the time from mid(_epoch).

*   i,j         - Loop counters
*   num         - Number of values in values

      integer*4 i, num

*   a(3), b(2)  - Nomral equations and solution vector
*   dates(1)    - Dates of the values
*   det         - Determinate of normal equations
*   dt          - Time difference between start and date (days)
*   mid         - Mid epoch of VLBI data
*   sigmas(1)   - Sigmas of the values
*   start       - Epoch to which trend will be referred
*   trend(2)    - Offset and rate
*   values(1)   - Values to be fitted
*   wgh         - Weight given to each point

      real*8 a(3), b(2), dates(1), det, dt, mid, sigmas(1), start,
     .    trend(2), values(1), wgh

***** First check to see if we have enough data

*                             ! NOT ENOUGH
      if( num.lt.2 ) then
          stop ' UPDATE_PMU Abort: Not enough data in tables'
      end if

*     Clear normal equations
      do i = 1,3
          a(i) = 0.d0
      end do
      do i = 1,2
          b(i) = 0.d0
      end do

****  Start incrementing normal equations

      do i = 1, num
          dt = dates(i) - start
          wgh = 1.d0/ (sigmas(i)**2 * (abs(dates(i)-mid)+.1d0)**2 )

          a(1) = a(1) + wgh
          a(2) = a(2) + wgh*dt
          a(3) = a(3) + wgh*dt*dt

          b(1) = b(1) + wgh*values(i)
          b(2) = b(2) + wgh*dt*values(i)
      end do

*     Invert and get the trend
      det = a(1)*a(3) - a(2)*a(2)

*     Chech the det
*                                     ! Not a good solution
      if( det.lt.1.d-7*a(1) ) then
          stop ' COMPARE_PMU aborted:  Trend removal singular'
      end if

      trend(1) = (a(3)*b(1) - a(2)*b(2) )/det
      trend(2) = (a(1)*b(2) - a(2)*b(1) )/det

****  Now remove trend from data
      do i = 1,num
          values(i) = values(i) - trend(1) - trend(2)*(dates(i)-start)
      end do

****  Thats all
      return
      end

CTITLE WRITE_SEC
 
      subroutine write_sec(iout)
 
*     Routine to compute the differences from the secondary inputs.
*     and write out the differencces with the mean removed.
*
      include 'compare_pmu.h'
 
* LOCAL VARIABLES
 
*   next_pt     - next point to be used in interpolating
*               - the smoothed series
*   date(5)     - Generic date array
*   iout        - Output unit number
*   i           - Loop counter over primary data set
*   ierr        - IOSTAT error
*   trimlen     - Length of string
 
      integer*4 next_pt, iout, ierr, trimlen, date(5), jerr
 
*   dxp,dyp,dut - differences for x,y pole and ut1 (mas,mts)
*   sxp, syp, sut   - Sigmas for differences (mas, mts)
*   dt          - Time difference between primary and
*               - first bounding smoothed data point
*   sectag      - Seconds tag for JD
*   dut1        _ Tide Ut1 contribtion

*   inf         - Polar/ut1 information for Lagrange interpolation.
*                 (jd, spacing and number)
*   wx(2), wy(2), wu(2) - Lagrangian interpolated values and rates
*   wmx, wmy, wmu  - Weighted means for XY and UT1
 
      real*8 dxp,dyp,dut, sxp, syp, sut,  sectag, dut1, inf(3),
     .       wx(2), wy(2), wu(2), wmx, wmy, wmu, dx, dy

*   line        - Line read from input

      character*128 line

***** Get the wweighted means
      if( secx_stats(2).ne.0 ) then
          wmx = secx_stats(1)/secx_stats(2)
          wmy = secy_stats(1)/secy_stats(2)
          wmu = secu_stats(1)/secu_stats(2)
      else
          wmx = 0.d0 
          wmy = 0.d0 
          wmu = 0.d0 
      end if

****  Set up for Lagrangian interpolation

      inf(1) = sm_jd(1)
      inf(2) = sm_jd(2)-sm_jd(1)
      inf(3) = num_sm
 
****  OPen the input file
      rewind(100) 
c     open(100,file=in_sec, status='old', iostat=ierr )
c     call report_error('IOSTAT',ierr,'open',in_sec,0, 'write_sec')
c     if( ierr.ne.0 ) RETURN
 
 
****  Create the output difference file
      iout = 200
      call open_lu(iout, out_sec, ierr, 'append' )
      call report_error('IOSTAT',ierr,'open', out_sec,0, 'write_sec')
      if( ierr.ne.0 ) iout = 6
 
***** Write some header information.
      write(iout,100) in_sec(1:trimlen(in_sec)), sm_step,
     .              in_pr(1:trimlen(in_pr)), sm_fwhm_scale
      write( 6  ,100) in_sec(1:trimlen(in_sec)), sm_step,
     .              in_pr(1:trimlen(in_pr)), sm_fwhm_scale
 100  format('* Residuals in ',a,' with ',F4.2,' day spacing',/,
     .       '* Smoothed series computed from ',a,/,
     .       '* Smoothing factor applied is ',f6.3)
 
****  Loop over the secimary data set.
      next_pt = 2
      do while ( ierr.eq.0 )
 
          read(100,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
*             Decode the line
              read(line,*,iostat=jerr) date, sec_xp, sec_xsig, sec_yp,
     .                        sec_ysig, sec_ut, sec_usig
              call report_error('IOSTAT',jerr,'decod', line,
     .            0, 'read_secimary')
              if( jerr.ne.0 ) write(*,'(a)') line(1:70)
              if( jerr.eq.0 ) then
                  call ymdhms_to_jd( date, sectag, sec_jd)
 
*                 See if we need to remove tidal UT1
                  if( ut1_sec(1:1).ne.'R' ) then
                      call short_period_ut1(sec_jd, dut1)
                      sec_ut = sec_ut - dut1
                  end if
                  if( ut1_sec(2:2).eq.'S' ) then
                      call sd_comp(sec_jd, dx,dy,dut1)
                      sec_xp = sec_xp - dx/1000.d0
                      sec_yp = sec_yp - dy/1000.d0
                      sec_ut = sec_ut - dut1/1000.d0
                  end if
              end if
 
*             Find the first point in the smmothed data after the
*             secimary JD
 
****          Use Lagrangian interpolation
              call lagrange_intp(inf,sm_xp, sec_jd, wx, 1)
              call lagrange_intp(inf,sm_yp, sec_jd, wy, 1)
              call lagrange_intp(inf,sm_ut, sec_jd, wu, 1)

              dxp  = sec_xp - wx(1)
              dyp  = sec_yp - wy(1)
              dut  = sec_ut - wu(1)
 
****          Now convert to milliarc sec and millitime secs
              dxp = dxp * 1000.d0 - wmx
              dyp = dyp * 1000.d0 - wmy 
*                                             ! Account for leapsecond
              dut = (dut - nint(dut)) * 1000.d0 - wmu
 
*             Do sigmas
              sxp = sec_xsig * 1000.d0
              syp = sec_ysig * 1000.d0
              sut = sec_usig * 1000.d0
 
*****         Now write values.
              if( jerr.eq.0 ) then
                  call jd_to_ymdhms( sec_jd+0.0001d0, date, sectag)
                  write(iout, 200) date, dxp, sxp, dyp, syp, dut, sut
c                 write( 6  , 200) date, dxp, sxp, dyp, syp, dut, sut
  200             format(1x,i4,4i3,1x,2(f9.3,1x,f8.3,1x),1x,
     .                   f10.3,1x,f9.3)
              end if
          end if
 
      end do
 
****  Thats all
      return
      end

