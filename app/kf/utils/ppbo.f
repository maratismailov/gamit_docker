      Program ppb
 
*     Program to compute the average and site dependent average
*     values for the horizontal and vertical rms scatters from
*     baseline length repeatibilty results (obtained from blsum)
*
*     Runstring:
*     % ppb <sum_file> <min_number> <max_len>
*
*     where <sum_file> is the name of a blsum summary file, and
*           <min_number> is the minimum number of determinations to
*                 to consider in the computation.
*           <max_len> maximum baseline length to be considered.
 
****  Include the main declarations
 
      include 'ppb.h'
 
****  Utility type declarations
 
*   rcpar       - Runstring values
*   decimaltoint    - Get interval value from string
*   len_run     - Length of string
*   ierr,jerr       - IOSTAT error
*   iel         - Site number return from list
*   indx        - Pointer in string
*   i,j,k           - Loop counters
*   iter        - iteration counter
 
 
      integer*4 rcpar, decimaltoint, len_run, ierr,jerr, iel, indx,
     .    i, trimlen, iter
 
*   values(11)  - The eleven values needed from the summary
*               - file lines (last entry is duration)
 
 
      real*8 values(11)
 
*   header      - Set true while header is being read
*   converged   - Indicates converged for all sites
 
 
      logical header, converged
 
*   rname       - Name of site read from line
 
 
      character*8 rname
 
*   line    - String read from runstring
 
 
 
      character*256 line
 
***** Decode the runstring
 
      len_run = rcpar(1, sum_file)
      if( len_run.le.0 ) then
          call proper_runstring('ppb.hlp', 'ppb', 1)
      else
*         OPen summaary file
          open(100, file = sum_file, iostat=ierr, status='old')
          call report_error('IOSTAT',ierr,'open', sum_file, 1,
     .                    'ppb/main')
      end if
 
****  Now see if min_num passed
 
      len_run = rcpar(2, line )
      if( len_run.gt.0 ) then
          min_num = decimaltoint( line, ierr)
          if( ierr.ne.0 ) then
              call report_error('decimaltoint',ierr,'decod', line,
     .                        1, 'ppb/main')
          end if
      else
         min_num = 3
      end if
 
***** See if max length passed
      len_run = rcpar(3, line )
      if( len_run.gt.0 ) then
          read(line,*,iostat=ierr) max_len
          max_len = max_len * 1000.d0
          if( ierr.ne.0 ) then
              call report_error('IOSTAT',ierr,'decod', line,
     .                        1, 'ppb/main')
          end if
      else
          max_len = 13000000.d0
      end if
 
***** Now readin the summary file, get the site names and statistics.
*     Get the start of the file from (mm/yr) line
 
      write(*,150) sum_file(1:trimlen(sum_file))
 150  format('* PPB: Reading summary file ',a)
 
      header = .true.
      do while ( header )
          read(100,'(a)', iostat=ierr) line
          indx = index( line, '(mm/yr)' )
          if( indx.gt.0 ) header = .false.
*                                     ! Error reading file
          if( ierr.ne.0 ) then
              call report_error('IOSTAT',ierr,'read', sum_file,0,
     .                        'ppb/main')
              write(*,*) ' *** Check that this is a summary file'
              stop ' PPB: End of header not found in summary file'
          end if
      end do
 
***** Now start reading the repeatability data.
      peak_soln = 1
      num_sites = 0
      num_ent = 0
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr) line
 
*         Get the first site name
          if( ierr.eq.0 ) then
 
*****         Now decode the rest of lines
              indx = 18
              read(line(18:),*,iostat=jerr) (values(i),i=1,11)
              call report_error('IOSTAT',jerr,'read',line(18:),0,
     .                        'ppb/main')
              if( jerr.eq.0 ) then
 
*****             See if the data satisfies the conditions set.
                  if( values(2).ge. min_num .and.
     .                values(3).le. max_len ) then
                      num_ent = num_ent + 1
                      rname = line(1:8)
                      indx = 1
                      call get_cmd(rname, site_names, num_sites,
     .                     iel, indx)
                      if( iel.le.0 ) then
                          num_sites = num_sites + 1
                          site_names(num_sites) = rname
                          iel = num_sites
                      end if
                      site(1,num_ent) = iel
 
*                     Now get second site
                      rname = line(10:17)
                      indx = 1
                      call get_cmd(rname, site_names, num_sites,
     .                     iel, indx)
                      if( iel.le.0 ) then
                          num_sites = num_sites + 1
                          site_names(num_sites) = rname
                          iel = num_sites
                      end if
                      site(2,num_ent) = iel
                      num_soln(num_ent) = nint(values(1))
                      peak_soln = max(peak_soln,num_soln(num_ent))
                      length(num_ent) = values(3)
 
*                     Use rms about slope of duration > 1 year
*                     and number of values greater than 3
                      if( values(11).gt. 1.d0 .and.
     .                    values(2).gt. 3 ) then
                          wrms(num_ent) = values(9)
                          chi(num_ent) = values(10)
                      else
                          wrms(num_ent) = values(5)
                          chi(num_ent) = values(6)
                      end if
                      num_est(num_ent) = nint(values(2))
*                             ! Does not meet condition, skip
                  end if
*                             ! Decode OK.
              end if
*                             ! EOF not reached on file
          end if
*                             ! Looping reading files
      end do
 
***** Give a short summary
      write(*,200) num_sites, num_ent, peak_soln
 200  format('* There are ',i3,' sites, ',i6,' entries, and ',i3,
     .    ' peak number of solutions in input')
 
****  Now loop over the peak number of solutions and estimate the
*     paramters of the fit
 
      do i = 1, peak_soln
 
          converged = .false.
          iter = 0

          call init_apr
 
          do while ( .not.converged .and. iter.lt.30  )
              call clear_norm(norm_av, bvec_av, 1, norm_st, bvec_st,
     .                       num_sites )
              call get_ppb_est(i, converged)
              iter = iter + 1
c         call out_ppb_est(i,iter)
c     call out_repeat
          end do
 
          call out_ppb_est(i,iter)
 
      end do
 
***** Now write out the repeatabiliies
 
      call out_repeat
 
****  Thats all
      end
 
CTITLE CLEAR_NORM
 
      subroutine clear_norm( norm_av, bvec_av, nav,
     .                       norm_st, bvec_st, nst )
 
*     Routine to clear the normal equations for the estimation of
*     the horizonatal and vertical site rms

* nav, nst - number of sites in average (1) and number of sites
*            in all.

      integer*4 nav, nst

* norm_av(2*nav,2*nav), bvec_av(2*nav) - Normal eqna and bvector of
*            average
* norm_st(2*nst,2*nst), bvec_st(2*nst) - Normal eqna and bvector of   
*            the site dependent averages.

      real*8 norm_av(2*nav,2*nav), bvec_av(2*nav), 
     .       norm_st(2*nst,2*nst), bvec_st(2*nst)

 
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
 
 
      integer*4 i,j
 
      do i = 1,2*nav
          bvec_av(i) = 0.d0
          do j = 1,2*nav
              norm_av(i,j) = 0.d0
          end do
      end do
 
      do i = 1, 2*nst
          bvec_st(i) = 0.d0
          do j = 1, 2*nst
              norm_st(i,j) = 0.d0
          end do
          norm_st(i,i) = 0.01d0 
      end do
 
***** Thats all
      return
      end
 
CTITLE GET_PPB_EST
 
      subroutine get_ppb_est(is, converged)
 
*     Routine to accumulate and solve the normal equations for
*     the average and site dependent horizontal and up noise
 
      include 'ppb.h'
 
* PASSED VARIABLES
 
*   is      - Solution number
*   converged - Indicates solution converegd
 
      integer*4 is
 
 
      logical converged
 
* LOCAL VARIABLES
 
*   i       - Loop counter
*   pivot(2*max_sites)  - Pivot elements.
 
 
      integer*4 i, pivot(2*max_sites), j
 
*   dbdh, dbdu  - derivatives of length with repect to horz and
*           - up changes (dimensionless)
*   scale(2*max_sites)  - Scaling for the inversion
*     cc(2*max_sites)  - One column of the covariance matrix
*                 (Needed for fixing values to positive numbers)
*     gn(2*max_sites)  - Gain array for fixing parameters
 
 
      real*8 dbdh, dbdu, scale(2*max_sites), cc(2*max_sites),
     .    gn(2*max_sites)
 
****  Loop over all of the data that we have
 
      do i = 1, num_ent
 
*         See if this is correct solution
*                                         ! Yes it is, so increent
          if( num_soln(i).eq.is) then
              call bl_part( length(i), dbdh, dbdu, Re )
 
****          Now increment the average normal equatiosn
              call inc_norm( 1, 1, wrms(i), num_est(i), dbdh, dbdu,
     .                bvec_av, norm_av, hs_apr_av(is), us_apr_av(is),
     .                hs_apr_av(is), us_apr_av(is), 1)
 
*****         Now increment site dependend values
              call inc_norm( site(1,i), site(2,i), wrms(i), num_est(i),
     .                dbdh, dbdu, bvec_st, norm_st,
     .                hs_apr_st(is,site(1,i)), us_apr_st(is,site(1,i)),
     .                hs_apr_st(is,site(2,i)), us_apr_st(is,site(2,i)),
     .                num_sites)
 
          end if
      end do
 
 
****  Now solve the system of equations.  Do the Average values first.
 
      call invert_vis( norm_av, bvec_av, scale, pivot, 2, 2,1)
 
c     hs_av(is) = sign(1.d0,bvec_av(1))*sqrt(abs(bvec_av(1)))
c     us_av(is) = sign(1.d0,bvec_av(2))*sqrt(abs(bvec_av(2)))
 
      converged = .true.
***** Get updated estimated and 
*     try to a stabilize the system if we have negative values
      hs_av(is) = hs_apr_av(is) + bvec_av(1)
      call fix_neg( norm_av, bvec_av, 2, hs_av(is), 1, 5.0d0,cc,gn )
      us_av(is) = us_apr_av(is) + bvec_av(2)
      call fix_neg( norm_av, bvec_av, 2, us_av(is), 2,15.0d0,cc,gn )
 
      if( abs(bvec_av(1)).gt. 0.1 ) converged = .false.
      if( abs(bvec_av(2)).gt. 0.1 ) converged = .false.
 
      hs_apr_av(is) = hs_av(is)
      us_apr_av(is) = us_av(is)
 
***** Now do the site dependent values
      call invert_vis( norm_st, bvec_st, scale, pivot, 2*num_sites,
     .                2*num_sites,1 )
 
      do i = 1, num_sites
c         hs_st(is,i) = sign(1.d0,bvec_st(2*i-1))*
c    .                       sqrt(abs(bvec_st(2*i-1)))
c         us_st(is,i) = sign(1.d0,bvec_st(2*i))*sqrt(abs(bvec_st(2*i)))

*         Get updated estimate and see if we need to fix negatives
          hs_st(is,i) = hs_apr_st(is,i) + bvec_st(2*i-1)
          call fix_neg( norm_st, bvec_st, 2*num_sites, hs_st(is,i),
     .                  2*i-1, hs_av(is),cc,gn )
          us_st(is,i) = us_apr_st(is,i) + bvec_st(2*i)
          call fix_neg( norm_st, bvec_st, 2*num_sites, us_st(is,i),
     .                  2*i,   us_av(is),cc,gn )
 
      end do

*     Do four iterations to make sure we get all negative
      do j = 1, 4 

c     write(*,'(''Second pass'',i4)') j
      do i = 1, num_sites
c         hs_st(is,i) = sign(1.d0,bvec_st(2*i-1))*
c    .                       sqrt(abs(bvec_st(2*i-1)))
c         us_st(is,i) = sign(1.d0,bvec_st(2*i))*sqrt(abs(bvec_st(2*i)))

*         Get updated estimate and see if we need to fix negatives
          hs_st(is,i) = hs_apr_st(is,i) + bvec_st(2*i-1)
          call fix_neg( norm_st, bvec_st, 2*num_sites, hs_st(is,i),
     .                  2*i-1, hs_av(is),cc,gn )
          us_st(is,i) = us_apr_st(is,i) + bvec_st(2*i)
          call fix_neg( norm_st, bvec_st, 2*num_sites, us_st(is,i),
     .                  2*i,   us_av(is),cc,gn )
 
      end do
      end do
c     write(*,'(''Third  pass'')')
      do i = 1, num_sites
c         hs_st(is,i) = sign(1.d0,bvec_st(2*i-1))*
c    .                       sqrt(abs(bvec_st(2*i-1)))
c         us_st(is,i) = sign(1.d0,bvec_st(2*i))*sqrt(abs(bvec_st(2*i)))

*         Get updated estimate and see if we need to fix negatives
          hs_st(is,i) = hs_apr_st(is,i) + bvec_st(2*i-1)
          call fix_neg( norm_st, bvec_st, 2*num_sites, hs_st(is,i),
     .                  2*i-1, hs_av(is),cc,gn )
          us_st(is,i) = us_apr_st(is,i) + bvec_st(2*i)
          call fix_neg( norm_st, bvec_st, 2*num_sites, us_st(is,i),
     .                  2*i,   us_av(is),cc,gn )
 
*         See if we have converged.
          if( abs(bvec_st(2*i-1)).gt.0.1 ) converged = .false.
          if( abs(bvec_st(2*i  )).gt.0.1 ) converged = .false.
          hs_apr_st(is,i) = hs_st(is,i)
          us_apr_st(is,i) = us_st(is,i)
      end do
 
***** Thats all.  Later we may want to add something for negative values
      return
      end

CTITLE INIT_APR

      subroutine init_apr

*     Rouitne to set the initial apriori values for the horzontal
*     and up sigmas

      include 'ppb.h'

* LOVAL VARIABLES

      integer*4 i, j

***** Set the average values
      do j = 1, max_soln
         hs_apr_av(j) =  5.d0
         us_apr_av(j) = 15.d0
      end do

      do i = 1, num_sites
         do j = 1, max_soln
            hs_apr_st(j,i) = hs_apr_av(j)
            us_apr_st(j,i) = us_apr_av(j)
         end do
      end do

***** Thats all
      return
      end
 
CTITLE FIX_NEG
 
      subroutine fix_neg( cov_parm, sol_parm, num_parm, force_value,
     .                nc, default,cov_col, equ_gn )
 
*     This routine to check the value so the force_value (parameter
*     in solution) and see if it is negative.  If it is then then
*     value will be forced to 'default'
 
*   num_parm    - Dimension of the solution system
*   nc          - parmater number to be checked
 
      integer*4 num_parm, nc
 
*   cov_parm(num_parm,num_parm) - Covariance matrix
*   sol_parm(num_parm)          - Solution vector.
*   force_value                 - Final estimate of the parameter
*                               - to be checked.
*   default                     - Default value to use fr the
*                               - for force_value if it is negative.
*   cov_col(num_parm)           - One column of covariance matrix
*   equ_gn(num_parm)            - Gain for KF implementation
 
      real*8 cov_parm(num_parm,num_parm), sol_parm(num_parm),
     .    force_value, default, cov_col(num_parm), equ_gn(num_parm)
 
* LOCAL VARIABLES
 
*   ipivot      - Part of Kalman gain calcultion.
 
      integer*4 ipivot
 
*   sol_col     - Value used to parameter to be fixed.
*   dchi        - Change is chi**2 when forces
*   avat        - Part of Kalmn gain calcultion
*   scale       - PArt of Kalmn gain calcultion
*   new_bv       - New value of solution vector such that
*               - apriori + new_bv = default.  apriori is
*               - implicitly computed from force_value -
*               - currect bvec value
 
 
      real*8 sol_col, dchi, avat, scale, new_bv
 
***** See if the parameter of concern is negative
 
      if( force_value.gt.0 ) RETURN
 
*     Value is negative. See how mucg much we need to change it
*     by.
c     write(*,100) nc, force_value, sol_parm(nc), default
c100  format(' Fixing P# ',i3,' from ',f8.3,'(dest ',f8.3,
c    .       ') to ',f8.3,$)
 
      new_bv = default - force_value + sol_parm(nc)
 
*     Now force value this value
 
      call force_parms( cov_parm, sol_parm, num_parm,
     .    cov_col, sol_col, nc, 1, new_bv, 0.d0, dchi, avat, equ_gn,
     .    ipivot, scale )
 
*     Update the value of the value of the forced parameter
      force_value = default
c     write(*,120) sol_parm(nc), sol_parm(nc+1)
c120  format(' After ',2f8.3)
 
***** Thats all
      return
      end
 
CTITLE  INC_NORM
 
      subroutine inc_norm( s1,s2, lwrms, num, dbdh, dbdu,
     .                bvec, norm, hs_apr1, us_apr1,
     .                hs_apr2, us_apr2, ns)
 
*     Subroutine to increment the normal equations
 
      include 'ppb.h'
 
* PASSED VARIABLES
*   s1,s2       ! The two sites numbers (if s1=s2=1 then
*               ! average values will be determined (single
*               ! site))
*   ns          ! Number of sites we are dimensioned for
*   num     ! Numebr of values in repeatability (used for
*               ! weight)
 
 
      integer*4 s1,s2, ns, num
 
*   dbdh, dbdu  - Baseline partials for horizontal and up
*   bvec(2*ns), norm(2*ns,2*ns) - B vector and normal
*               - equations
*   lwrms        - WRMS scatter
*   hs_apr1, us_apr1, hs_apr2, us_apr2 - Apriori values of the
*                  fits to the model.
 
 
      real*8 dbdh, dbdu, bvec(2*ns), norm(2*ns,2*ns), lwrms,
     .    hs_apr1, us_apr1, hs_apr2, us_apr2
 
* LOCAL VARIABLES
 
*   iels(4)     - Elements to be incremented
*   i,j         - Loop counter
 
 
      integer*2 iels(4), i,j
 
*   wgh         - Weight in least sqaures estimate
*   OmC         - wrms - apriorui value
*   parts(4)        - partials for incrmenting normals
*   denom       - Denominator for computing partials
 
      real*8 wgh, OmC, parts(4), denom
 
****  Get the weiught and "data"
      wgh = num
      OmC = lwrms - sqrt((hs_apr1**2+hs_apr2**2)*dbdh**2 +
     .                   (us_apr1**2+us_apr2**2)*dbdu**2  )
      denom = sqrt((hs_apr1**2+hs_apr2**2)*dbdh**2 +  
     .             (us_apr1**2+us_apr2**2)*dbdu**2  ) 
 
      parts(1) = dbdh**2*hs_apr1/denom
      parts(2) = dbdu**2*us_apr1/denom 
      iels(1) = s1*2 - 1
      iels(2) = s1*2
 
      parts(3) = dbdh**2*hs_apr2/denom
      parts(4) = dbdu**2*us_apr2/denom
      iels(3)  = s2*2 - 1
      iels(4)  = s2*2
 
****  Now incremnt the bvector and normal equations
      do i = 1,4
          bvec(iels(i)) = bvec(iels(i)) + OmC*wgh*parts(i)
          do j = 1, 4
              norm(iels(i),iels(j)) = norm(iels(i),iels(j)) +
     .            parts(i)*wgh*parts(j)
          end do
      end do
 
****  Thats all
      return
      end
 
CTITLE BL_PART
 
      subroutine bl_part ( locl, dbdh, dbdu, Re )
 
 
*     Routine to compute the partial of length with reprect to the
*     horizontal and up changes
 
*   locl  - Baseline length (m)
*   dbdh, dbdu  - derivatives of length with repect to horz and
*           - up changes (dimensionless)
*   Re      - Radius of Earth (m)
 
 
      real*8 locl, dbdh, dbdu, Re
 
      dbdu = locl / (2.d0*Re)
      dbdh = sqrt(1.d0-dbdu**2)
 
***** Thats all
      return
      end
 
CTITLE OUT_PPB_EST
 
      subroutine out_ppb_est(is, iter)
 
*     Routitine to write out the estimates of the horizontal and
*     up errors (the Up errors are also expressed as ppb number)
 
 
 
      include 'ppb.h'
 
* PASSED VARIABLES
 
*   is      - Solution number
*   iter    - NUmber of iterations for solution
 
 
      integer*4 is, iter
 
* LOCAL VARIABLES
 
*   i   - Loop counter
 
      integer*4 i

*   r2   - Sqrt(2) for getting baseline results

      real*8 r2

      data r2 / 1.414213563237d0 /
 
***** First output average values
 
      write(*,100) is, iter
 100  format('* For solution number ',i3,' after ',i3,' iterations',/,
     .       '* ',20x,'SINGLE SITE VALUES ',10x,
     .       '       BASELINE VALUES ',/,
     .       '*   Site        Horizontal         Up       PPB ',
     .       '    Horizontal         Up       PPB ',/,
     .       '*                 (mm)            (mm)          ',
     .       '       (mm)            (mm)')
      write(*,150) 'AVERAGE', hs_av(is), us_av(is),
     .            us_av(is)/(2*Re)* 1.d6, r2*hs_av(is),
     .            r2*us_av(is), r2*us_av(is)/(2*Re)* 1.d6
 150  format('*',2x,a8,3x,F10.2,3x,f10.2,3x,f8.3, 3x,
     .          F10.2,3x,f10.2,3x,f8.3)
 
 
****  Loop over sites
      do i = 1, num_sites
          write(*,150) site_names(i), hs_st(is,i), us_st(is,i),
     .            us_st(is,i)/(2*Re)* 1.d6
      end do
 
****  Thats all
      return
      end
 
CTITLE OUT_REPEAT
 
      subroutine out_repeat
 
*     This rouitne will write out the data set along the estimates
*     of expected repeatability from the average results and the
*     site dependent results.
 
      include 'ppb.h'
 
* LOCAL VARIABLES
 
*   i       - Loop counter
*   ns      - Solution number (local for convenience)
 
 
      integer*4 i, ns
 
*   dbdh, dbdu      - Partials of baseline length wrt
*           - horizontal and up.
*   wrms_av - Expected Wrms from average
*   wrms_st - expected wmrs from site dependent results.
 
 
      real*8 dbdh, dbdu, wrms_av, wrms_st
 
****  Loop over all of the entries
      write(*,100)
 100  format('*  Length    wrms   nrms   #  wrms (av) '
     .        ' wrms (site) Soln. Baseline',/,
     .       '*   (km)     (mm)               (mm)    ',
     .        '   (mm)')
 
      do i = 1, num_ent
 
*         Get the average and site dependent wrms expected
 
          call bl_part( length(i), dbdh, dbdu, Re )
          ns = num_soln(i)
 
          wrms_av = sqrt(2*(hs_av(ns)*dbdh)**2 +
     .                   2*(us_av(ns)*dbdu)**2 )
          wrms_st = sqrt((hs_st(ns,site(1,i))*dbdh)**2 +
     .                    (hs_st(ns,site(2,i))*dbdh)**2 +
     .                    (us_st(ns,site(1,i))*dbdu)**2 +
     .                    (us_st(ns,site(2,i))*dbdu)**2 )
 
*****     Now write out the resulys
          write(*,150) length(i)/1000.d0, wrms(i), chi(i), num_est(i),
     .                wrms_av, wrms_st, num_soln(i),
     .                site_names(site(1,i)), site_names(site(2,i))
 150      format( F10.3, F7.2,1x,f6.3,1x,i4,1x,2(F7.2,1x), I3,2x,
     .            a8,1x,a8)
      end do
 
***** Thats all
      return
      end
 
