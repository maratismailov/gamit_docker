 
CTITLE READ_PMU_INP     
 
      subroutine read_pmu_inp

      implicit none 
 
*     Routine to read the Polar motion UT1 values from a file.  The
*     routine assumes that the tabular points are unformly spaced
*     although the spacing may change through the file provided that
*     we do not try to interpolate accross the boundary between spcing
*     changes.
 
*     The also assumes that the UT1 values have been regularized before
*     being used here.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
* LOCAL VARIABLES:
 
*   ierr, jerr  - IOSTAT error
*   trimlen     - Length of string
*   date(5)     - Date of Pole position estimate
*   ntab        - Number of tabular points passed to
*               - Lagrange_inpt routine
 
      integer*4 ierr, jerr, trimlen, date(5), ntab, i, j, k
 
*   sectag      - Seconds tag for JD comps
*   jdpole      - Julian date of pole position/UT1 in file
*   prev_jd     - Julian date of previous value (used to get
*               - spacing
*   spacing     - Spaceing between table values (days)
*   values(6)   - X pole +- Ypole +- and UT1-TAI +- values
*   pmu_vals(3,2)   - Interpolated XY and UT1 values (mas) and
*               - their derivatives (mas/day)
*   inf(3)      - Information array for interpolation.
*   pts(3,max_mul_pmu+20)   - Tabular points to be used (XY and UT1)
*   lps         - Estimate of leapsec size.

      real*8 sectag, jdpole, prev_jd, spacing, values(6),
     .    pmu_vals(3,2), inf(3), pts(3,max_mul_pmu+20), start, stop, lps 

      real*8 gep_utc  ! Epoch converted to UT1 by truncating to nearest
     .                ! minute

*   Leap_sec    - Leap second from tables.
      integer*4 leap_sec

*   line        - Line read from input file
 
 
      character*256 line
 
****  Initialize
      prev_jd = 0.d0
      sectag  = 0.d0
      spacing = 0.d0
      ntab    = 0 

      if( num_mul_pmu.gt.0 ) then
          start = min(gepoch_expt,gmul_pmu_ep(1))
          stop  = max(gepoch_expt,gmul_pmu_ep(num_mul_pmu))

*         See if we need to read new values
C         if( .not.kbit(mul_pmu_opt,32) ) RETURN
      else
          start = gepoch_expt
          stop  = gepoch_expt
      end if
       
 
****  If we applying pmu open file and get values.
!     print *,'DEBUG: Start GWOB ',gwob_apr(1:4), gut1_apr(1:2)
!     print *,'DEBUG: Start CWOB ',cwob_apr(1:4), cut1_apr(1:2)
!     print *,'DEBUG:    mul_pmu ',apr_val_mul_pmu(:,:,1) 


      if( trimlen(pmu_inp_file).gt.0 ) then

          open(105, file=pmu_inp_file, status='old', iostat = ierr)
          call report_error('IOSTAT',ierr,'open',pmu_inp_file,0,
     .                    'GLB_UPR_APR')
 
*         Do nothing if there is an error.  Setting name blank will
*         turn off this feature.
          if( ierr.ne.0 ) then
              pmu_inp_file = ' '
          end if
 
*         Start reading the file
          do while ( ierr.eq.0 )
              read(105,'(a)', iostat=ierr) line
              call report_error('IOSTAT',ierr,'read',pmu_inp_file,0,
     .                    'GLB_UPR_APR')
              if( ierr.ne.0 ) then
                  write(*,200)
  200             format('*** WARNING *** Polar Motion/UT1 values',
     .                   ' will not be updated')
              end if
 
              if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
*****             See if we are in correct time range.
                  read(line,*,iostat=jerr) date, values
 
                  call ymdhms_to_jd( date, sectag, jdpole)
 
*****             Check to get spacing of tables
                  if( prev_jd.ne.0.d0 ) then
                      spacing = jdpole - prev_jd
                  end if
                  prev_jd = jdpole
 
*****             See if we are close to correct position.
                  if( jdpole.ge.start-2*spacing .and.
     .                jdpole.le.stop +2*spacing ) then
 
*****                 Time in correct range so say tabular points
*                     (Skip over sigmas in values)
                      ntab = ntab + 1 
                      if( ntab.gt.max_mul_pmu+20 ) then
                          print *,'MAXPMU: ',jdpole, start, stop, 
     .                              spacing
                          call report_stat('FATAL','GLFOR',
     .                         'read_pmu_inp',' ',
     .                        'Too many PMU tabular points: Max ',
     .                        max_mul_pmu+20)
                      end if
                      pts(1,ntab) = values(1)*1000.d0
                      pts(2,ntab) = values(3)*1000.d0

*                     Check to see UT1-UTC or UT1-TAI (convert to
*                     UT1-TAI if it is the other. Convert from time
*                     seconds to mas.
C                     dut = nint(cut1_apr(1)/15.d3-values(5))
* MOD TAH 990222: Check to see if UT1-UTC; if then add in leap seconds
                      if( abs(values(5)).lt.2.d0 ) then
                         call get_leapsec(jdpole, 1.d5, leap_sec)
                      else
                         leap_sec = 0.d0
                      endif 
                      values(5) = (values(5) + leap_sec)*15.d3
 
*                     Check value
                      pts(3,ntab) = values(5)
                      if( ntab.eq.1 ) then
*                         Save time of first point
                          inf(1) = jdpole
                          inf(2) = spacing
                      end if
 
                  end if
 
*****             Now see if we have found all of the values
                  if( jdpole.ge.stop +2*spacing ) then
*                                 ! Force exit
                      ierr = -1

****                  Interpolate the values and save them
                      if( ntab.ge.4 ) then 
                          inf(3) = ntab
                          gep_utc = nint((gepoch_expt-2450000)*1440)/
     .                        1440.d0+2450000
!                         call Lagrange_intp(inf, pts, gepoch_expt,
!    .                        pmu_vals, 3)
                          call Lagrange_intp(inf, pts, gep_utc,
     .                        pmu_vals, 3)
!      print *,'PMU ',inf, pts(:,1:4), gepoch_expt, pmu_vals, gep_utc
*****                     See if UT1-AT versus UT1-UTC problem
                          lps = nint((pmu_vals(3,1)-
     .                                cut1_apr(1))/15.d3)*15.d3
 
*****                     Save them
*                         Check values to make sure OK.

                          if(abs(pmu_vals(3,1)-cut1_apr(1)-lps).gt.
     .                                                100.d0 ) then 

*                             If the value 0.d0 then most likely
*                             not used, so don't warn user
                              if( cut1_apr(1).ne.0.d0 )
     .                        write(*,250) pmu_vals(3,1)/15.d0,
     .                                     cut1_apr(1)/15.d0,
     .                                     lps/15.d3
 250                          format(' *** ERROR **** in PMU tables:',
     .                            ' Interpolated UT1-TAI is ',f10.2,
     .                            ' Apriori value is ',f10.2,' LPS ',
     .                            F8.1,' NOT USED')
                          else

* DEBUG:
                              if( glb_com_file(1:8).eq.'test.com' ) then
                                  write(*,*) 'DEBUG: Apriori EOP ',
     .                                    gwob_apr, ' UT1 ',gut1_apr
                                  write(*,*) 'DEBUG: New EOP ',lps, 
     .                                    pmu_vals
                              end if
                          end if    
                          gwob_apr(1) = pmu_vals(1,1)
                          gwob_apr(2) = pmu_vals(2,1)
                          gwob_apr(3) = pmu_vals(1,2)
                          gwob_apr(4) = pmu_vals(2,2)
                          gut1_apr(1) = pmu_vals(3,1)
                          gut1_apr(2) = pmu_vals(3,2)

c                          end if

****                      Now interpolate the multi-pmu epochs
                          do i = 1, num_mul_pmu
                              call Lagrange_intp(inf, pts, 
     .                              gmul_pmu_ep(i),  pmu_vals, 3)
 
*****                         Save them
*****                         See if UT1-AT versus UT1-UTC problem
                              lps = nint((pmu_vals(3,1)-
     .                            apr_val_mul_pmu(1,3,i))/15.d3)*15.d3
*                             Check values to make sure OK.
                              if(abs(pmu_vals(3,1)-
     .                           apr_val_mul_pmu(1,3,i)-lps).gt.300.
     .                           .and. apr_val_mul_pmu(1,3,i)
     .                            .ne. -9999.d0 )then
                                  write(*,250) pmu_vals(3,1)/15.d0,
     .                                  apr_val_mul_pmu(1,3,i)/15.d0,
     .                                  lps/15.d3
                              else 
                                  do j = 1,3
                                     do k = 1,2
                                       apr_val_mul_pmu(k,j,i) = 
     .                                             pmu_vals(j,k)
                                     end do
                                  end do
                              end if
                          end do 
  
                      else

                          write(*,270) pmu_inp_file(1:
     .                                 trimlen(pmu_inp_file))
  270                     format('READ_PMU_INP: Experiment not in',
     .                           ' time range of tables in ', a)

* ADD DEBUG TAH 131223: See what statis
                          if( glb_com_file(1:8).eq.'test.com' ) then
                              write(*,*) 'DEBUG: Apriori EOP ',
     .                                gwob_apr, ' UT1 ',gut1_apr
                              do i = 1, num_mul_pmu
                                 write(*,*) 'DEBUG: MultiPMU ',i,
     .                                ' apr_val_mul_pmu(1) ',
     .                                apr_val_mul_pmu(1,:,i),
     .                                ' apr_val_mul_pmu(2) ',
     .                                apr_val_mul_pmu(2,:,i)
                              end do

                           end if
                      end if
                  end if
*                         ! Error zero and line not a comment
              end if
*                         ! Looping over lines in file
          end do
*                         ! Close the input file.
          close(105)
*                         ! If we reading the polar motion/ut1 file
      end if
!     print *,'DEBUG: End   GWOB/UT1 ',gwob_apr(1:4), gut1_apr(1:2)
!     print *,'DEBUG: End   CWOB/UT1 ',cwob_apr(1:4), cut1_apr(1:2)

***** Thats all
      return
      end
 
