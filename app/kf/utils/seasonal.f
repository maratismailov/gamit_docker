 
      program seasonal

      implicit none 
 
*     This program will read a series of summap files and produce
*     a seasonal model for temperature (in this model).  This model
*     is then fitted to a latitude and height function.
 
      include '../includes/const_param.h'
      include 'seasonal.h'
 
* DECLARATIONS for MAIN Program
 
*   rcpar               - Reads the runstring
*   ierr                - IOSTAT error on file open
*   i,j,k               - Loop counters
*   ir                  - Runstring counter
 
      integer*4 rcpar, ierr, i,ir, trimlen
 
*   infile                  - Name of current input file
 
 
 
      character*132 infile
 
***** Loop over the runstring and read all of the data
 
      write(*,100)
 100  format(/' Program SEASONAL:',/,
     .        ' Runstring: seasonal <list of input files>',/)
 
      ir = 0
      call init_nobs( max_rsites, num_robs )
 
      do while (rcpar(ir+1,infile) .gt.0 )
          ir = ir + 1
          write(*,150) ir, infile(1:trimlen(infile))
 150      format(' Processing file ',i3,1x,a)
 
          open(100, file = infile, iostat=ierr, status='old')
          call report_error('IOSTAT',ierr,'open', infile, 0,
     .                    'SEASONAL')
          if( ierr.eq.0 ) then
              call read_data( 100 )
          end if
      end do
 
***** Now summarise the data contents
      call sumdat
 
***** Now loop over site and fit the seasonal terms to data
 
      write(*,200)
 200  format(/' Seasonal model fit to each site',/,
     .        ' -------------------------------',/,
     .    ' Site name ','Total Used ', 10x,' Projected (C)',10x,
     .    ' Tbias (C)', 10x, 'RMS Tp (C)',1x,'RMS Tb(C)',
     .    2x,'Latitude',2x,'Height (km)',/,
     .    20x,' Mean     Cos(month)    Sin(month)',
     .     5x,' Mean     Cos(month)    Sin(month)',
     .    20x,' (deg)',6x,'(km)')
 
 
      do i = 1, max_rsites
 
          if( num_robs(i).gt.0 ) then
 
              call fit_seas(i)
 
              call out_seas(i)
 
          end if
 
      end do
 
****  Now we have all the site dependendent values, fit to a latitude
*     height model
 
      call fit_gmean
      call fit_gseas 
 
****  Write out the global model
 
      write(*,300)
 300  format(/,' Global Model',/,
     .        ' ------------',/,
     .        ' Term ',15x,' Mean ',5x, ' Cos/Sin(lat) ',2x,' ht (km)',
     .        2x,' RMS ')
 
      call out_global
 
****  Thats all
      end
 
CTITLE READ_DATA
 
      subroutine read_data( unit )
 
      implicit none 

*     This routine will read the input files.  If a station is found
*     which is not the list of known stations it is ignored.  If there
*     is a read error on the line, it is also ignored.
 
      include 'seasonal.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number for reading the file
 
      integer*4 unit
 
* LOCAL VARIABLES
 
*   ierr, jerr  - IOSTAT errors
*   is          - Station number (from list)
*   jn          - Temporary storage of obsevration number
*   indx        - Position in string
*   trimlen     - Length of line
*   date(5)     - Date ymdmh
*   iscode      - Integer version of code
 
      integer*4 ierr, jerr, is, jn, indx, trimlen, date(5), iscode
 
*   ah,bh,ch, rmsh  - ABC and RMS for Hydrostatic delay
*   aw,bw,cw, rmsw  - ABC and RMS for Wet delay
*   th, tw          - Hydrostatic and Wet zenith delays (mm)
*   es              - Surface partial pressure of water vapor
*   ts, tp, dt      - Surface temp, projected temp and
*                   - difference (C)
*   beta, ht            - Lapse rate and height of trop.
 
      real*4 ah,bh,ch, rmsh, aw,bw,cw, rmsw, th, tw, es, ts, tp, dt,
     .    beta, ht
 
*   sectag          - Seconds tag
 
      real*8 sectag
 
*   scode           - Code for the station
 
      character*8 scode
 
*   line        - Line read from file
 
      character*256 line
 
****  Start reading the file
 
      ierr = 0
      sectag = 0.d0
 
      do while ( ierr.eq.0 )
          read(unit,'(a)',iostat=ierr) line
          if( trimlen(line).gt.0 .and. ierr.eq.0 ) then
              read(line,*,iostat=jerr) date, iscode, ah,bh,ch, rmsh,
     .                aw,bw,cw, rmsw, th, tw, es, ts, tp, dt,
     .                beta, ht
 
*             If read OK, continue processing
              if( jerr.eq.0 ) then

                  write(scode,'(i5)') iscode
 
*                 See if we can find the station
                  indx = 1
                  call get_cmd( scode, site_codes, max_rsites,
     .                    is, indx )
 
*                 See if we found station, and not too many obs
                  if( is.gt.0 .and.
     .                num_robs(is)+1 .lt. max_robs ) then

                      num_robs(is) = num_robs(is) + 1
                      jn = num_robs(is)
                      temps(jn,is) = tp
                      tbias(jn,is) = dt
                      call ymdhms_to_jd( date, sectag, epochs(jn,is))
                      edit_flag(jn,is) = 0
 
*                     A quick edit check
                      if( tp.gt.40.0 .or. tp.lt. -40 ) then
                          edit_flag(jn,is) = 1
                      end if
                  end if
              end if
          end if
      end do
 
****  We have now finished
      close(unit)
      return
      end
 
CTITLE INIT_NOBS
 
      subroutine init_nobs( max_rsites, num_robs )

      implicit none 
 
*     Routine to initialize the number of observations
 
*   max_rsites  - MAx number of sites
*   num_robs(max_rsites)    - Number of observations
 
 
      integer*4 max_rsites, num_robs(max_rsites)
 
*   i           - Loop counter
 
      integer*4 i
 
      do i = 1, max_rsites
          num_robs(i) = 0
      end do
 
****  Thats all
      return
      end
 
CTITLE SUMDAT
 
      subroutine sumdat

      implicit none 
 
*     Routine to summarize the contents of the data files read
 
      include 'seasonal.h'
 
 
*   i,j     - Loop counters
 
 
      integer*4 i,j
 
***** Loop over all possible stations and see what is available
      j = 0
      do i = 1, max_rsites
          if( num_robs(i).gt.0 ) then
              j = j + 1
              write(*,100) j, site_names(i), site_codes(i), num_robs(i)
  100         format(2x,i3,'. ',a8,1x,a8, I5)
          end if
      end do
 
      write(*,150) j
  150 format(' There are ',i3,' Sites in the data files')
 
****  Thats all
      return
      end
 
 
CTITLE FIT_SEAS
 
      subroutine  fit_seas( is )

      implicit none 
 
****  Routine to fit a model of the form
*
*     T = Tm + Tc * cos(month) + Ts * sin(month)
*
*     to the temperature data
 
 
      include '../includes/const_param.h'
      include 'seasonal.h'
 
* Passed variables
 
*   is      - Site number
 
      integer*4 is
 
* Local variables
 
*   i,j,k   - Loop counters
*   iter    - Number of iterations
*   nt      - Number of data for T and Bias
*   ne      - Nunber of new edits.
*   ipvt(3) - Pivot elements
 
      integer*4 i,j, iter, nt, ipvt(3), ne
 
*   normt                   - Normal equations
*   bvec(3,2)               - Bvector for temps and bias
*   sol(3,2)                    - Solution results for T and B
*   y(2)                    - Observations
*   dt                      - Time difference from DJ2000 in
*                           - radians
*   rms(2)                  - RMS scatters for each.
*   partials(3)             - Partial derivatives
*   scale(3)                    - Scale factors for inversion
*   drms                    - Change in rms
*   dsol                    - Solution estimate
 
      real*8 norm(3,3), bvec(3,2), sol(3,2), y(2), dt, rms(2),
     .     partials(3), scale(3), drms, dsol
 
*   converged   - Indicates that the values have converged after
*               - editing.
 
      logical converged
 
****  Loop, editing data until the results are completed
 
      converged = .false.
      iter = 0
      nt = 0
      do while ( .not.converged .and. iter.lt.10 )
          iter = iter + 1
 
*         Loop over the data
          call clear_fit( norm, bvec, rms, 3, 2)
          nt = 0
          do i = 1, num_robs(is)
 
*             Convert the time to days since Jan 1, 2000
              dt = epochs(i,is) - DJ2000
              dt = mod(dt/365.25d0,1.d0)*2*pi
              partials(1) = 1.d0
              partials(2) = cos(dt)
              partials(3) = sin(dt)
 
              y(1) = temps(i,is)
              y(2) = tbias(i,is)
 
*             Increment the normal equations is data OK.
              if( edit_flag(i,is).eq.0 ) then
                  nt = nt + 1
                  call incr_norm( norm, bvec, partials, y, rms, 3,2)
              end if
          end do
 
****      Now solve the system of equations.  Save B vector so that we
*         can compute the adjustment to the rms
          do i = 1,3
              do j = 1,2
                  sol(i,j) = bvec(i,j)
              end do
          end do
 
          call invert_vis( norm, sol, scale, ipvt, 3, 3, 2 )
 
*****     Now get the adjustment to the rms
          do i = 1,2
              drms = 0.d0
              do j = 1,3
                  drms = drms + sol(j,i)*bvec(j,i)
              end do
              rms(i) = sqrt((rms(i) - drms)/nt)
          end do
 
*****     Now loop and edit the data.
          ne = 0
          do i = 1, num_robs(is)
 
*             Convert the time to days since Jan 1, 2000
              dt = epochs(i,is) - DJ2000
              dt = mod(dt/365.25d0,1.d0)*2*pi
              partials(1) = 1.d0
              partials(2) = cos(dt)
              partials(3) = sin(dt)
 
              y(1) = temps(i,is)
              y(2) = tbias(i,is)
 
              do j = 1,2
 
                  dsol = partials(1)*sol(1,j) + partials(2)*sol(2,j)
     .                + partials(3)*sol(3,j)
 
                  if( abs(temps(i,is)-dsol).le.3*rms(j) ) then
                      edit_flag(i,is) = 0
                  end if
 
                  if( abs(y(j)-dsol).gt.3*rms(j) ) then
                      ne = ne + 1
                      edit_flag(i,is) = 1
                  end if
              end do
          end do
 
*****     See if we did not edit
          if( ne.eq.0 ) converged = .true.
      end do
 
***** Copy the final results
      do j = 1,3
          seas_temp(j,is) = sol(j,1)
          seas_bias(j,is) = sol(j,2)
      end do
 
      num_used(is) = nt
      rms_temp(is) = rms(1)
      rms_bias(is) = rms(2)
 
***** Thats all
      return
      end
 
CTITLE CLEAR_FIT
 
      subroutine clear_fit ( a,b, rms, np, ns)

      implicit none 
 
*     Routine to clear the normal equations and setup the rms
 
*   np, ns  - Number of parameters and number of solutions
 
      integer*4 np, ns
 
*   a(np,np), b(np,ns)  - Normal equations and solution vector
*   rms(ns)             - Rms of each solution.
 
      real*8 a(np,np), b(np,ns), rms(ns)
 
*   i,j                 - Loop counters
 
      integer*4 i,j
 
****  Clear everything
 
      do i = 1, np
          do j = 1, ns
              b(i,j) = 0.d0
          end do
          do j = 1, np
              a(i,j) = 0.d0
          end do
      end do
 
      do i = 1, ns
          rms(i) = 0.d0
      end do
 
***** Thats all
      return
      end
 
CTITLE INCR_NORM
 
      subroutine incr_norm( a,b, p, y, rms, np, ns )

      implicit none 
 
*     Routine to increment normal equations and solution vector.
 
*   np, ns  - Number of parameters and number of solutions
 
      integer*4 np, ns
 
*   a(np,np), b(np,ns)  - Normal equations and solution vector
*   rms(ns)             - Rms of each solution.
*   p(np)               - Partial derivatives
*   y(ns)               - Data
 
      real*8 a(np,np), b(np,ns), rms(ns), p(np), y(ns)
 
*   i,j                 - Loop counters
 
      integer*4 i,j
 
*                                 T       T
****  Loop over partials forming P P and P y
 
      do i = 1, np
          do j = 1, ns
              b(i,j) = b(i,j) + p(i)*y(j)
          end do
          do j = 1, np
              a(i,j) = a(i,j) + p(i)*p(j)
          end do
      end do
 
      do i = 1, ns
          rms(i) = rms(i) + y(i)**2
      end do
 
***** Thats all
      return
      end
 
CTITLE OUT_SEAS
 
      subroutine out_seas( is )

      implicit none 
 
*     Routine to output the site dependent seasonal results
 
      include 'seasonal.h'
 
*   is      - Site number to output
*   j       - Loop counter
 
 
      integer*4 is, j
 
****  Just write out the lines
 
      write(*,100) site_names(is), num_robs(is), num_used(is),
     .    (seas_temp(j,is), j=1,3), (seas_bias(j,is), j=1,3),
     .    rms_temp(is), rms_bias(is), site_lat(is), site_ht(is)
  100 format(1x,a8,1x,2i5,1x,6(F8.2,1x),1x,2(f8.2,1x), F5.1,1x,
     .       f5.3)
 
***** Thats all
      return
      end
 
CTITLE fit_gmean
 
      subroutine fit_gmean
 
      implicit none 

 
*     This routine will fit a global model to the seasonal coefficients
*     for the mean values.
*
*     The model form is (for mean temperature)
*
*     Ti = Tm + Ts * cos (latitude) + Th * h

*     Where the fit is to each of the model components.
 
      include '../includes/const_param.h'
      include 'seasonal.h'
 
 
* Local variables
 
*   i,j,k   - Loop counters
*   iter    - Number of iterations
*   is      - Site number
*   nt      - Number of stations
*   ipvt(3) - Pivot elements
 
      integer*4 i,j,is, nt, ipvt(3)
 
*   norm(3,3)               - Normal equations
*   bvec(3,2)               - Bvector for temps and bias
*   sol(3,2)                    - Solution results for T and B
*   y(2)                    - Observations
*   dt                      - Time difference from DJ2000 in
*                           - radians
*   rms(2)                  - RMS scatters for each.
*   partials(3)             - Partial derivatives
*   scale(3)                    - Scale factors for inversion
*   drms                    - Change in rms
 
 
      real*8 norm(3,3), bvec(3,2), sol(3,2), y(2), rms(2), 
     .    partials(3), scale(3), drms
 
*     Set up and Loop over the stations
 
      call clear_fit( norm, bvec, rms, 3, 2)
      nt = 0
      do is = 1, max_rsites
 
*         See if we have data
          if( num_used(is).gt.0 ) then
              nt = nt + 1
              partials(1) = 1.d0
              partials(2) = cos(site_lat(is)*pi/180.d0)
              partials(3) = site_ht(is)
 
              y(1) = seas_temp(1,is)
              y(2) = seas_bias(1,is)
 
 
              call incr_norm( norm, bvec, partials, y, rms, 3,2)
          end if
      end do
 
****  Now solve the system of equations.  Save B vector so that we
*     can compute the adjustment to the rms
      do i = 1,3
          do j = 1,2
              sol(i,j) = bvec(i,j)
          end do
      end do
 
****  See if we have enough stations
      if( nt.lt.3 ) then
          write(*,*) ' Only ', nt, ' stations '
          stop ' STOP: Not enough stations to fit a global model'
      end if
 
      call invert_vis( norm, sol, scale, ipvt, 3, 3, 2 )
 
***** Now get the adjustment to the rms
      do i = 1,2
          drms = 0.d0
          do j = 1,3
              drms = drms + sol(j,i)*bvec(j,i)
          end do
          if( nt.gt.3 ) then
              rms(i) = sqrt((rms(i) - drms)/(nt-3))
          else
              rms(i) = sqrt((rms(i) - drms))
          end if
 
      end do
 
***** Save the results
 
      do i = 1,3
          temp_mean(i)  = sol(i,1)
          tbias_mean(i) = sol(i,2)
      end do
 
      do i = 1,2
          rms_global((i-1)*3 + 1) = rms(i)
      end do
 
****  Print out residuals to global model
      write(*,100)
  100 format(/' RESIDUALS TO GLOBAL MEAN FIT',/,
     .        ' ----------------------------',/,
     .    ' Site name ', 10x,' Projected (C)',10x, ' Tbias (C)',/,
     .    20x,' Mean ',
     .     5x,' Mean ')
 
      do is = 1, max_rsites
 
*         See if we have data
          if( num_used(is).gt.0 ) then
              partials(1) = 1.d0
              partials(2) = cos(site_lat(is)*pi/180.d0)
              partials(3) = site_ht(is)
 
              y(1) = seas_temp(1,is)
              y(2) = seas_bias(1,is)
 
              do i = 1,2
                  y(i) = y(i) - sol(1,i)*partials(1)
     .                       - sol(2,i)*partials(2)
     .                       - sol(3,i)*partials(3)
              end do
 
              write(*,200) site_names(is), (y(i), i = 1,2)
 200          format(1x,a8,1x,6(F8.2,1x))
          end if
      end do
 
 
***** Thats all
      return
      end
 
CTITLE OUT_GLOBAL
 
      subroutine out_global

      implicit none 
 
*     Routine to output the global model results
 
      include 'seasonal.h'
 
 
*   i       - Loop counter
 
      integer*4 i
 
***** Now write each of the types
 
      write(*,100) 'Tp_mean', (temp_mean(i),  i=1,3), rms_global(1)
      write(*,100) 'Tp_cosm', (temp_cosm(i),  i=1,3), rms_global(2)
      write(*,100) 'Tp_sinm', (temp_sinm(i),  i=1,3), rms_global(3)
      write(*,100) 'dT_mean', (tbias_mean(i), i=1,3), rms_global(4)
      write(*,100) 'dT_cosm', (tbias_cosm(i), i=1,3), rms_global(5)
      write(*,100) 'dT_sinm', (tbias_sinm(i), i=1,3), rms_global(6)
 
  100 format(1x,a8,3(F8.3,1x),1x,F8.2)
 
***** Thats all
      return
      end
 
CTITLE SEAS_BD
 
      block data seas_bd

      implicit none 
 
*     Block data with site names, codes, locations
 
      include 'seasonal.h'
 
      data site_names / 'ALBANY  ',
     .                'CHATHAM ',
     .                'FAIRBNKS',
     .                'NOME    ',
     .                'PT.MUGU ',
     .                'PORTLAND',
     .                'SAN DIEG',
     .                'SAN NIC ',
     .                'WEST PB ',
     .                'ALBUQUER',
     .                'DESRT RK',
     .                'EDW AFB ',
     .                'EL PASO ',
     .                'VANDNBRG'  /
 
      data site_lat /   42.5,    41.7,    64.8,    64.5,    34.2,
     .                  45.6,    32.8,    33.3,    26.8,    35.0,
     .                  36.6,    34.9,    31.5,    34.8  /
 
      data site_ht  /   0.089,   0.023,   0.137,   0.014,   0.003,
     .                  0.008,   0.010,   0.173,   0.077,   1.619,
     .                  1.007,   0.706,   1.198,   0.116  /
 
      data site_codes / '72518', '74494', '70261', '70200', '72391',
     .                '72606', '72290', '72291', '72203', '72365',
     .                '72387', '72381', '72270', '72393'  /
 
      end
 
 
CTITLE fit_gseas
 
      subroutine fit_gseas

      implicit none 
 
 
*     This routine will fit a global model to the seasonal coefficients
*
*     and for the month dependent arguments.
*     Ti = Tm + Ts * sin (latitude) + Th * h
*
*     Where the fit is to each of the model components.
 
      include '../includes/const_param.h'
      include 'seasonal.h'
 
 
* Local variables
 
*   i,j,k   - Loop counters
*   iter    - Number of iterations
*   is      - Site number
*   nt      - Number of stations
*   ipvt(3) - Pivot elements
 
      integer*4 i,j,is, nt, ipvt(3)
 
*   norm(3,3)               - Normal equations
*   bvec(3,4)               - Bvector for temps and bias
*   sol(3,4)                    - Solution results for T and B
*   y(4)                    - Observations
*   dt                      - Time difference from DJ2000 in
*                           - radians
*   rms(4)                  - RMS scatters for each.
*   partials(3)             - Partial derivatives
*   scale(3)                    - Scale factors for inversion
*   drms                    - Change in rms
 
 
      real*8 norm(3,3), bvec(3,4), sol(3,4), y(4), rms(4), 
     .    partials(3), scale(3), drms
 
*     Set up and Loop over the stations
 
      call clear_fit( norm, bvec, rms, 3, 4)
      norm(1,1) = 1.d0
      nt = 0
      do is = 1, max_rsites
 
*         See if we have data
          if( num_used(is).gt.0 ) then
              nt = nt + 1
              partials(1) = 0.d0
              partials(2) = sin(site_lat(is)*pi/180.d0)
              partials(3) = site_ht(is)*sin(site_lat(is)*pi/180.d0)
 
              y(1) = seas_temp(2,is)
              y(2) = seas_temp(3,is)
 
              y(3) = seas_bias(2,is)
              y(4) = seas_bias(3,is)
 
 
              call incr_norm( norm, bvec, partials, y, rms, 3,4)
          end if
      end do
 
****  Now solve the system of equations.  Save B vector so that we
*     can compute the adjustment to the rms
      do i = 1,3
          do j = 1,4
              sol(i,j) = bvec(i,j)
          end do
      end do
 
****  See if we have enough stations
      if( nt.lt.3 ) then
          write(*,*) ' Only ', nt, ' stations '
          stop ' STOP: Not enough stations to fit a global model'
      end if
 
      call invert_vis( norm, sol, scale, ipvt, 3, 3, 4 )
 
***** Now get the adjustment to the rms
      do i = 1,4
          drms = 0.d0
          do j = 1,3
              drms = drms + sol(j,i)*bvec(j,i)
          end do
          if( nt.gt.3 ) then
              rms(i) = sqrt((rms(i) - drms)/(nt-3))
          else
              rms(i) = sqrt((rms(i) - drms))
          end if
 
      end do
 
***** Save the results
 
      do i = 1,3
          temp_cosm(i)  = sol(i,1)
          temp_sinm(i)  = sol(i,2)
 
          tbias_cosm(i) = sol(i,3)
          tbias_sinm(i) = sol(i,4)
      end do
 
      rms_global(2) = rms(1)
      rms_global(3) = rms(2)
      rms_global(5) = rms(3)
      rms_global(6) = rms(4)
 
****  Print out residuals to global model
      write(*,100)
  100 format(/' RESIDUALS TO GLOBAL SEASONAL FIT',/,
     .        ' -------- -----------------------',/,
     .    ' Site name ', 10x,' Projected (C)',10x, ' Tbias (C)',/,
     .    20x,'  Cos(month)    Sin(month)',
     .     5x,'  Cos(month)    Sin(month)')
 
      do is = 1, max_rsites
 
*         See if we have data
          if( num_used(is).gt.0 ) then
              partials(1) = 0.d0
              partials(2) = sin(site_lat(is)*pi/180.d0)
              partials(3) = site_ht(is)*sin(site_lat(is)*pi/180.d0)
 
              y(1) = seas_temp(2,is)
              y(2) = seas_temp(3,is)
 
              y(3) = seas_bias(2,is)
              y(4) = seas_bias(3,is)
 
              do i = 1,4
                  y(i) = y(i) - sol(1,i)*partials(1)
     .                       - sol(2,i)*partials(2)
     .                       - sol(3,i)*partials(3)
              end do
 
              write(*,200) site_names(is), (y(i), i = 1,4)
 200          format(1x,a8,1x,6(F8.2,1x))
          end if
      end do
 
 
***** Thats all
      return
      end
 
