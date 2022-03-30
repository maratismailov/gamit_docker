CTITLE REPORT_EQ
 
      subroutine report_eq(iout, opt)

      implicit none 
 
*     Routine to report the values of the Earthquake parameters
*     used.
* MOD TAH 110528: Added option to report all renames.  (Mainly for debug
*     testing of code).
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/globk_common.h'
 
* PASSED variables
 
*   iout        - Output unit
 
      integer*4 iout
      character*(*) opt  ! Option: ALL will output all renames 
 
* LOCAL variables
 
*   i,j     - Loop counters
*   date_start(5), date_end(5)  - Start and stop times of
*           - renames
 
      integer*4 i,j, date_start(5), date_end(5)
      integer*4 num_rnout, num_eqout  ! Number of renames and 
                                      ! earthquakes reported
 
*   loc_coord(3)        - Colat, long and height
*   rot_mat(3,3)        - Rotatopn from xyz to neu
 
      real*8 sectag, loc_coord(3), rot_mat(3,3)

*   kbit  -- Test bit
      logical kbit
 
*   ren         - YES or NO to rename option
 
 
      character*4 ren
 
***   Start with the renames:
      num_rnout = 0 
      num_eqout = 0

****  See if global eq_reset used
      if( eq_reset ) then
          write(iout,80)
 80       format(/,'EQ_FILE STATUS: All name extents set to _GPS')
      endif

****  Count how many used
      do i = 1, num_renames

* MOD TAH 980514: See if rename really used.
          if( kbit(rn_used,i) .or. index(opt,'ALL').gt.0 ) then
             num_rnout = num_rnout + 1
         endif
      end do
      if( num_renames.gt.0 ) then
          write(iout, 100) num_renames, num_rnout, opt
 100      format(/' There were ',i6,' site renames input ',i6,
     .            ' were used. ',
     .            ' Reporting ',a,'. Renames used are: ',/,
     .        '   #   Orig     New          Specific ',
     .        ' Period from   -----> ',
     .        '     To              Position change (m)       Type')
      end if
 
      do i = 1, num_renames

* MOD TAH 980514: See if rename really used.
          if( kbit(rn_used,i) .or. index(opt,'ALL').gt.0 ) then
             num_rnout = num_rnout + 1
             call jd_to_ymdhms(rn_times(1,i), date_start, sectag)
             call jd_to_ymdhms(rn_times(2,i), date_end, sectag)
             write(iout,120) i, rn_codes(1,i), rn_codes(2,i),
     .               rn_hfiles(i),
     .               date_start, date_end, (rn_dpos(j,i), j=1,3),
     .               rn_types(i)
120          format(i6,1x,a8,'->',a8,1x,a16,i4,'/',i2,'/',i2,1x,i2,':',
     .           i2,'  ',i4,'/',i2,'/',i2,1x,i2,':',i2,2x,
     .           3(F10.4,1x),a3)
          endif
      end do
 
****  Now do Earthquakes:
      do i = 1,num_eq
          if( kbit(eq_used,i) .or. index(opt,'ALL').gt.0 ) then
              num_eqout = num_eqout + 1
          end if
      end do

      if( num_eq.gt.0 ) then
          write(iout,200) num_eq, num_eqout
 200      format(/,' There were ',i3,' earthquakes input; ',i3,
     .        ' were used. Earthquakes used are:',/,
     .        '  #  CODE     Lat (deg)  Long (deg) ',
     .        'Radius (km) Depth (km)    Date       Rename?')
      end if
 
      do i = 1,num_eq
          if( kbit(eq_used,i) .or. index(opt,'ALL').gt.0 ) then
             call XYZ_to_GEOD(rot_mat, eq_pos(1,i), loc_coord)
             if ( eq_rename(i) ) then
                 ren = 'YES'
             else
                 ren = 'NO'
             end if
             call jd_to_ymdhms(eq_epoch(i), date_start,sectag)
             write(iout,220) i, eq_codes(i), 
     .           (pi/2-loc_coord(1))*180./pi,
     .           loc_coord(2)*180./pi, eq_rad(i)/1000.0, 
     .           eq_depth(i)/1000.,
     .           date_start, ren
 220         format(i3,3x,a2,3x,F10.4,1x,f10.4,3x,f10.4,1x,f10.4,1x,
     .               i4,'/',i2,'/',i2,1x,i2,':',i2,2x,a3)
          end if
      end do
 
*     Now give coseimic sigmas
      if( num_eq.gt.0 ) then
          write(iout,250)
 250      format(/,' COSEISMIC characteristics',/,
     .        ' #  CODE              Static sigma        ',
     .        '       Spatial Sigma (Depth/Dist)^2',/,
     .        '               North      East       Height (m) ',
     .        '   North     East       Height (m)')
          do i = 1, num_eq
            if( kbit(eq_used,i) ) 
     .      write(iout,260) i,eq_codes(i), 
     .                      (eq_apr_coseismic(j,i),j=1,6)
 260        format((i3,3x,a3,1x,3(f10.5,1x),2x,3(f10.5,1x)))
          end do
          write(iout,270)
 270      format(/,' PRE-SEISMIC characteristics',/,
     .         ' #  CODE    Dur              Static Process   ',
     .         '        Spatial Process (Depth/Dist)^2',/,
     .         '          (days)    North      East       Height  ',
     .         '      North     East       Height',/,
     .         '                             (mm^2/day)           ',
     .         '            (mm^2/day)')
          do i = 1, num_eq
            if( kbit(eq_used,i) ) 
     .      write(iout,280) i,eq_codes(i), eq_dur(1,i),
     .                    (eq_mar_pre(j,i)/365.25d-6,j=1,6)
          enddo
 280      format((i3,3x,a3,1x,F5.1,1x, 3(f10.5,1x),2x,3(f10.5,1x)))
          write(iout,290)
 290      format(/,' POST-SEISMIC characteristics',/,
     .         ' #  CODE    Dur              Static Process   ',
     .         '        Spatial Process (Depth/Dist)^2',/,
     .         '          (days)    North      East       Height  ',
     .         '      North     East       Height',/,
     .         '                             (mm^2/day)           ',
     .         '            (mm^2/day)')
          do i = 1, num_eq
            if( kbit(eq_used,i) ) 
     .      write(iout,280) i,eq_codes(i), eq_dur(2,i),
     .              (eq_mar_post(j,i)/365.25d-6,j=1,6)
          end do
          write(iout,310)
 310      format(/,' POST-SEISMIC LOG Estimates',/,
     .         ' #  CODE    Tau              Static Log   ',
     .         '        Spatial Process (Depth/Dist)^2',/,
     .         '          (days)    North      East       Height  ',
     .         '      North     East       Height',/,
     .         '                                (m)           ',
     .         '             (m)')
          do i = 1, num_eq
            if( kbit(eq_used,i) ) 
     .      write(iout,280) i,eq_codes(i), eq_log_tau(i),
     .              (eq_log_sig(j,i),j=1,6)
          end do
      end if
 
****  Thats all
      return
      end
 
CTITLE REPORT_SMAR
 
      subroutine report_smar(iout)

      implicit none 
 
*     Routine to report the site markov process noise being used
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/globk_common.h'
 
* PASSED variables
 
*   iout        - Output unit
 
      integer*4 iout
 
* LOCAL variables
 
*   i,j     - Loop counters
*   date_start(5), date_end(5)  - Start and stop times of
*           - renames
 
      integer*4 i,j

      logical kbit
 
 
      write(iout,120)
 120  format(/,'SITE RANDOM WALK PROCESS NOISE VALUES      ',
     .       '                                         SMAR',/,
     .       '  SITE         N RW       E RW        U RW',
     .       '         NR RW       ER RW       UR RW    SMAR',/,
     .       '                         (m^2/yr)              ', 
     .       '             (m/yr)^2/yr             SMAR')

      do i = 1, gnum_sites
         if( kbit(guse_site,i) ) then
            write(iout,220) gsite_names(i), 
     .        (mar_neu(j,1,i)*1.d6,j=1,3),
     .        (mar_neu(j,2,i)*1.d6,j=1,3)
 220        format(1x,a8,1x,3(F8.2,'E-6 '),1x,3(F8.2,'E-6 '),
     .          ' SMAR')
         endif
      enddo

****  Thats all
      return
      end


 
 
