CTITLE WRITE_EPHEM
 
      subroutine write_ephem( unit )

      implicit none
 
*     This routine will add the ephemeris information from
*     this Q file to an epheris file connected to unit.
*
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/sln_def.h'
      include '../includes/glb_hdr_def.h'
 
* PASSED VARIABLES
 
*   unit        - unit number for output
 
      integer*4 unit
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   i,j     - Loop counters
*   trimlen - Length of string
*   date(5) - Date of the experiment. Lines output with this
*           - time so that GLOBK can match the ephemeris entry
*           - global file
*   eph_date(4) - ymdh for the ephemeris entry
 
      integer*4 ierr, i,j, trimlen, date(5), eph_date(4)
 
*   sectag  - Seconds tag for time
*   Dec_yrs - Deciminal years
 
      real*8 sectag, dec_yrs
 
*     Write some header information
      write(unit,100, iostat=ierr) hfile(1:trimlen(hfile))
 100  format('* EPHEMERIS INFORMATION FROM ',a)
      call report_error('IOSTAT',ierr,'head writ','empheris file',0,
     .    'write_ephem/htoglb')
 
*     Get the output time to use
 
      call jd_to_ymdhms( sepoch, date, sectag )
      do i = 1,4
          eph_date(i) = date(i)
      end do

      call jd_to_decyrs( sepoch, dec_yrs)
 
*     For the moment, write the station names and positions in globk
*     format.
 
      do i = 1, qnum_sites
          write(unit,200,iostat=ierr) qsite_names(i),
     .            (site_pos(j,i),j=1,3), 
     .            (site_vel(j,i),j=1,3), dec_yrs
 200      format('X ',a8, 3(1x,f15.5), 3(1x,f9.5),1x,f9.4)
      end do
 
*     Now write the satellite ephemeris
      do i = 1, qnum_svs
          write(unit,300, iostat=ierr) eph_date, qsvs_names(i),
     .            (svs_pos(j,i),j=1,max_svs_elem)
 300  format(i5,3i3,1x,a8,6(1x,f15.5),1x,17(1x,f8.5))
      end do
      call report_error('IOSTAT',ierr,'ephemeris writ',
     .        'Ephemeris file',0,'write_ephem/htoglb')
 
****  Thats all
      return
      end
 
 
 
