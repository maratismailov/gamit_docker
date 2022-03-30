 
      program gamit_to_iris
 
*     This program will take as input the two files used by GAMIT for
*     it polar motion and uT1 tables (nominally pole. and ut1.) and
*     output a file in standard IRIS format.  (For use in SOLVK)
*
 
      integer*4 max_iers
 
*                                     ! Max number of entries allowed
      parameter ( max_iers = 10000 )
 
*   ierr    - IOSTAT errors
*   rcpar   - Read runstring
*   date(5) - Date in ymdhms
*   len_inpm, len_inut, len_outiris - Lengths of the pm and
*           - ut file names
*   npm     - Number of pm values
*   nut     - Number of ut values
*   i,j     - Loop counters
*   iout    - Output unit number
*   gen_time(7) - Time at which table generated
*   imjd        - Integer mjd read from files
*   ivalues(12) - Values read from tables.
*   trimlen     - Length of used portion of string
*   spacing     - Spacign between values (days) read from format line
 
      integer*4 ierr, rcpar, date(5), len_inpm, len_inut, len_outiris,
     .    npm, nut, i,j, iout, gen_time(7), imjd, ivalues(12), trimlen,
     .    spacing
 
*   found   - Indicates that we have found UT and PM value if in
*           - search mode.
 
 
      logical found
 
*   sectag  - Seconds tag for jd to date conversion
*   mjd_pm(max_iers)    - Julian dates for the pole positions
*   mjd_ut(max_iers)    - Julian dates for UT1-AT
*   pm(2,max_iers)  - X and Y pole positions
*   ut(max_iers)        - UT1-AT values
*   pole_scael, ut1_scale - Scale factors to asec and tsec.
 
      real*8 sectag, mjd_pm(max_iers), mjd_ut(max_iers),
     .    pm(2,max_iers), ut(max_iers), pole_scale, ut1_scale
 
*   line            - Line read from files
*   in_pm           - Name of in pm file
*   in_ut           - Name of in UT1-AT file
*   out_iris            - Output IRIS format file
*   format_pm, format_ut    - Formats for the pm and ut files
*   ut_header, pm_header   - Header lines for ut and pole tables.
 
 
      character*128 in_pm, in_ut, out_iris, format_pm, format_ut,
     .    ut_header, pm_header
 
****  Start decoding the runstring.  Get PM file
      len_inpm = rcpar(1,in_pm)
      if( len_inpm.le.0 ) then
          call proper_runstring('gamit_to_iris.hlp','gamit_to_iris',1)
      end if
 
      open(100, file = in_pm, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',in_pm,1, 'gamit_to_iris')
 
*     Get UT file
      len_inut = rcpar(2,in_ut)
      if( len_inut.le.0 ) then
          call proper_runstring('gamit_to_iris.hlp','gamit_to_iris',1)
      end if
 
      open(101, file = in_ut, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',in_ut,1,'gamit_to_iris')
 
*     Get the output file name
      len_outiris = rcpar(3,out_iris)
      if( len_outiris.le.0 ) then
          call proper_runstring('gamit_to_iris.hlp','gamit_to_iris',1)
      end if
 
      iout = 200
      call open_lu(iout, out_iris, ierr, 'unknown')
      call report_error('IOSTAT',ierr,'open',out_iris,1,'gamit_to_iris')
 
****  Now loop over the PM values
 
      read(100,'(a)', iostat=ierr) pm_header
      read(100,'(a35,18x,I3,6x,F12.6)', iostat=ierr) format_pm, spacing,
     .                                      pole_scale
 
      npm = 0
 
      do while ( ierr.eq.0 )
 
          read(100,format_pm,iostat=ierr) imjd, (ivalues(i),i=1,12)
 
****      If no error then assign to tables
          if( ierr.eq.0 ) then
              do i = 1,6
                  mjd_pm(npm+i) = imjd + spacing*(i-1)
                  pm(1,npm+i) = ivalues(2*(i-1)+1)
                  pm(2,npm+i) = ivalues(2*(i-1)+2)
              end do
              npm = npm + 6
          end if
      end do
 
****  Now loop over the UT values
 
      read(101,'(a)', iostat=ierr) ut_header
      read(101,'(a30,32x,f12.6)', iostat=ierr) format_ut, ut1_scale
 
      nut = 0
 
      do while ( ierr.eq.0 )
 
          read(101,format_ut,iostat=ierr) imjd, (ivalues(i),i=1,6)
 
****      If no error then assign to tables
          if( ierr.eq.0 ) then
              do i = 1,6
                  mjd_ut(nut+i) = imjd + spacing*(i-1)
*                 Change sign to make UT1-AT.
                  ut(nut+i) = -ivalues(i)
              end do
              nut = nut + 6
          end if
      end do
 
****  Tell use what is happening
 
      write(*,200) npm, in_pm(1:len_inpm), nut, in_ut(1:len_inut)
 200  format(/' gamit_to_iris: Conversion program',/,
     .        i4,' XY Pole doublets read from ',a,/,
     .        i4,' UT1-AT values    read from ',a)
 
****  Write the header for the output
 
      call systime( gen_time, sectag)
      gen_time(6) = nint(sectag)
 
      write(iout,250) in_pm(1:len_inpm), in_ut(1:len_inut),
     .            gen_time, pole_scale, ut1_scale,
     .            pm_header(1:trimlen(pm_header)),
     .            ut_header(1:trimlen(ut_header))
 250  format('* IRIS Format table form ',a,' and ',a,/,
     .    '* Generated at ',I4,2('/',i2),1x,i2,':',i2,1x,
     .    i2,'.',i2,' Scales: ',2d9.2,/,
     .    '* ',a,/,'* ',a,/,
     .    '*    Date ',t20,'X-Pole (arcsec) +-',
     .                 t40,'Y-Pole (arcsec) +-',
     .                 t60,'UTC-AT (timesec)  +-'  )
 
 
****  Now merge the two tables and write the output.
      j = 0
      do i = 1, npm
 
*         Find the correct ut entry
          j = j + 1
*                                                     ! Enter search
          if( abs(mjd_pm(i)-mjd_ut(j)).gt.0.001d0 ) then
*                                                 ! mode
              j = 0
              write(*,300) i,j
  300         format(' Search mode for ',i4,'th pm value.',
     .               ' At ',i4,'th UT value...',$)
              j = 0
              found = .false.
              do while ( .not.found .and. j.lt.nut)
                  j = j + 1
                  if( abs(mjd_pm(i)-mjd_ut(j)).gt.0.001d0 ) then
                      found = .true.
                      write(*,320) j
  320                 format('Found at ',i4)
                  end if
              end do
          end if
 
****      Now convert to date and write
          call jd_to_ymdhms( mjd_pm(i)+0.0001d0, date, sectag)
          write(iout,400) date, pm(1,i)*pole_scale,
     .                          pm(2,i)*pole_scale,
     .                          ut(j)*ut1_scale
 400      format(1x,i4,4i3,1x,2(F9.6,' 0.0010 '),1x,
     .            f12.6,1x,' 0.00010 ')
      end do
 
****  Thats all
      end
 
