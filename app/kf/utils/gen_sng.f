      Program gen_sng

      implicit none 
 
*     This program will read a pait of files with the quadrature
*     components of the diurnal and semi-diurnal UT1 and Pole position
*     variations and output the instantaneous pole position and UT1
*     values.
 
*   rcpar       - Reads the runstring
*   len_run     - Length of runstring
*   ierr        - IOSTAT error
*   date(5)     - Date ymdhms
*   trimlen     - Length of string
*   i           - Loop variable
 
      integer*4 rcpar, len_run, ierr, jerr, date(5), trimlen, i
 
*   ut1_tot(4) - Diurnal and Semidiurnal totals
*   ut1_adj(4)  - Diurnal and semidiurnal adjustments (mas)
*   ut1_sig(4)  - Sigmas on the adjustments
*   xy_tot(6)  - Diurnal and Semidiurnal totals for Pole
*               - position
*   xy_adj(6)   - Diurnal and semidiurnal adjustments (mas)
*   xy_sig(6)   - Sigmas on the adjustments
 
*   sectag      - Seconds tag
*   xtot, ytot, utot    - Total positions of pole and UT1
*   xady, yadj, uadj    - Adjuments to the apriori values
*   xsig, ysig, usig    - Sigmas of the adjustments (mas/mts)
*   gst         - Greenwich Sidereal time (rads)
*   cd,sd, cs, ss       - Cos and Sin of diurnal and semidiurna;
*               - arguments
 
*   jd          - General Julian date
*   jd_mid      - JD of the values read from file
 
*   step        - Step size for values in days
*   duration    - NUmber of days to generate each value for.
 
 
      real*8 ut1_tot(4), ut1_adj(4), ut1_sig(4), xy_tot(6),
     .    xy_adj(6), xy_sig(6), sectag, xtot, ytot, utot,
     .    xadj, yadj, uadj, cd,sd, cs, ss, jd, jd_mid, step, duration,
     .    gst, xsig, ysig, usig

      integer*4 ijd, njd  ! Counter for JD loop 
 
*   in_ut1      - Input UT1 file
*   in_xy       - Input XY file
*   runstring   - RUnstring for getting step and duration.
*   line        - Line read from input file
 
      character*256 in_ut1, in_xy, runstring, line
 
***** Start decoding the runstring, Get input UT1 file
      len_run = rcpar(1, in_ut1)
      if( len_run.le.0 ) then
          call proper_runstring( 'gen_sng.hlp','gen_sng',1)
      end if
 
      open(100, file=in_ut1, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',in_ut1,0, 'gen_sng')
      if( ierr.ne.0 ) then
          call proper_runstring( 'gen_sng.hlp','gen_sng',1)
      end if
 
****  Get the xy file
      len_run = rcpar(2, in_xy)
      if( len_run.le.0 ) then
          call proper_runstring( 'gen_sng.hlp','gen_sng',1)
      end if
 
      open(101, file=in_xy, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',in_xy,0, 'gen_sng')
      if( ierr.ne.0 ) then
          call proper_runstring( 'gen_sng.hlp','gen_sng',1)
      end if
 
****  Now get the spacing.  If none passed then use 0.1 days
      len_run = rcpar(3, runstring)
      if( len_run.gt.0 ) then
          read(runstring,*,iostat=ierr) step
          call report_error('IOSTAT',ierr,'decod',runstring,1,
     .                    'gen_sng/step size')
      else
          step = 0.1d0
      end if
 
*     Get the optional duration
      len_run = rcpar(4, runstring)
      if( len_run.gt.0 ) then
          read(runstring,*,iostat=ierr) duration
          call report_error('IOSTAT',ierr,'decod',runstring,1,
     .                    'gen_sng/duration')
      else
          duration = 1.d0
      end if
 
 
***** Now start reading the input files.  These are assumed to
*     synchronized so that each can be read in parallel.
 
      write(*,100 ) in_ut1(1:trimlen(in_ut1)), in_xy(1:trimlen(in_xy))
 100  format('* Short period UT1/polar motion genereted from inputs:',/,
     .       '* ',a,' and ',a,/,
     .       '*      Date       xtot (mas) +-  ytot (mas) +- ',
     .       '  Utot (mts) +-',
     .       '   xadj (mas)   yadj (mas)   Uadj (mts)')
 
      do while ( ierr.eq.0 )
          line = '*'
          do while ( (line(1:1).ne.' ' .or. trimlen(line).eq.0 )
     .               .and. ierr.eq.0  )
              read(100,'(a)', iostat=ierr ) line
          end do
 
          read(line,*,iostat=jerr) date, (ut1_tot(i), ut1_adj(i),
     .                                ut1_sig(i), i=1,4)
 
*         If no error continue
          line = '*'
          do while ( (line(1:1).ne.' ' .or. trimlen(line).eq.0 )
     .               .and. ierr.eq.0  )
              read(101,'(a)', iostat=ierr ) line
          end do
 
          read(line,*,iostat=jerr) date, (xy_tot(i), xy_adj(i),
     .                                xy_sig(i), i=1,6)
 
 
*****     OK, If no error then generate the values
          if( ierr.eq.0 .and. jerr.eq.0 ) then
              sectag = 0.d0
              call ymdhms_to_jd( date, sectag, jd_mid)
              njd = nint(duration/step) 
!             do jd = jd_mid-duration/2, jd_mid+duration/2, step
              do ijd = 0, njd 
                  jd = jd_mid-duration/2+ijd*step

                  call gst_jd( jd, gst)
 
*****             Now compute the values
                  cd = cos(gst)
                  sd = sin(gst)
                  cs = cos(2*gst)
                  ss = sin(2*gst)
 
*                 Now do the sums.  UT1 first
                  utot = ut1_tot(1)*cd + ut1_tot(2)*sd +
     .                    ut1_tot(3)*cs + ut1_tot(4)*ss
                  uadj = ut1_adj(1)*cd + ut1_adj(2)*sd +
     .                    ut1_adj(3)*cs + ut1_adj(4)*ss
                  usig = sqrt((ut1_sig(1)*cd)**2 + 
     .                        (ut1_sig(2)*sd)**2 +
     .                        (ut1_sig(3)*cs)**2 +
     .                        (ut1_sig(4)*ss)**2  )
 
*****             Now do the pole positions
                  xtot = -xy_tot(1)*cd + xy_tot(2)*sd
     .                - xy_tot(3)*cs - xy_tot(4)*ss
     .                - xy_tot(5)*cs + xy_tot(6)*ss
                  xadj = -xy_adj(1)*cd + xy_adj(2)*sd
     .                - xy_adj(3)*cs - xy_adj(4)*ss
     .                - xy_adj(5)*cs + xy_adj(6)*ss
                  xsig  = sqrt((xy_sig(1)*cd)**2 +
     .                         (xy_sig(2)*sd)**2 +
     .                         (xy_sig(3)*cs)**2 +
     .                         (xy_sig(4)*ss)**2 +
     .                         (xy_sig(5)*cs)**2 +
     .                         (xy_sig(6)*ss)**2   )  
 
                  ytot = xy_tot(1)*sd + xy_tot(2)*cd
     .                - xy_tot(3)*ss + xy_tot(4)*cs
     .                + xy_tot(5)*ss + xy_tot(6)*cs
                  yadj = xy_adj(1)*sd + xy_adj(2)*cd
     .                - xy_adj(3)*ss + xy_adj(4)*cs
     .                + xy_adj(5)*ss + xy_adj(6)*cs
                  ysig  = sqrt((xy_sig(1)*sd)**2 +
     .                         (xy_sig(2)*cd)**2 +
     .                         (xy_sig(3)*ss)**2 +
     .                         (xy_sig(4)*cs)**2 +
     .                         (xy_sig(5)*ss)**2 +
     .                         (xy_sig(6)*cs)**2   )  
 
 
*****             Now write out the results
                  call jd_to_ymdhms( jd, date, sectag)
                  write(*,200) date, xtot, xsig, ytot, ysig, 
     .                         utot/15.d0, usig/15.d0,  
     .                         xadj, yadj, uadj/15.d0
 200              format(I5,4i3,4F8.3, 2F9.4, 2f8.3,f9.4)
              end do
          end if
      end do
 
***** Thats all
      close(100)
      close(101)
 
      end
 
 
 
 
 
 
 
