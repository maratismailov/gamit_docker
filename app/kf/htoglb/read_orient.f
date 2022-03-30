CTITLE READ_ORIENT
 
      subroutine read_orient ( unit )
 
      implicit none

*     This routine will read the apriori values for the
*     pole position, UT1 and nutation used in the GAMIT
*     run and save them in the appropriate locations
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
* unit  - Inoput file unit number
 
 
      integer*4 unit
 
* LOCAL VARIABLES
 
*      pmu_jd, pmu_secs     - JD and seconds tag for the
*                           - Julian date of the PMU epoch
*   dpmu_epoch      - Difference in days bewteen the midpoint
*                   - of the experiment and the PMU epoch
*   sectag          - Second tag when date from read
 
      real*8 pmu_jd, pmu_secs, dpmu_epoch, sectag
 
*   ierr            - IOSTAT error flag
*   date(5)         - Date when date form read.
 
      integer*4 ierr, date(5)
 
*   line            - line read from input file
 
      character*256 line
 
****  Get the next line containing the Polar motion epoch
      read(unit,'(a)', iostat=ierr) line
      if( line(1:11).eq.' E-rotation' ) then
          read(line(34:),*, iostat=ierr ) pmu_jd, pmu_secs
          pmu_jd = pmu_jd - 0.5d0 + pmu_secs/86400.d0
      else

*         New format with yy mm dd etc
          read(line,'(13x,i4,2x,i2,2x,i2,2x,i2,3x,i2,3x,f5.3)') 
     .                                    date,sectag 
          call ymdhms_to_jd( date, sectag, pmu_jd )
      end if

*     Get time difference between this and midpoint fo
*     data span
      dpmu_epoch = (qstart_epoch+qend_epoch)/2 - pmu_jd

      write(*,100) dpmu_epoch
 100  format(' Earth orientation data found, Epoch difference',
     .       ' from mid-point is ',F4.2,' days')
 
****  Get UT1
      read(unit,'(a)', iostat=ierr ) line
      read(line(34:),*,iostat=ierr) qut1_apr(1), qut1_apr(2)
 
*     Convert the units from secs and secs/day to mas and mas/day
*     and map values to the correct epoch.
      qut1_apr(1) = (qut1_apr(1)+dpmu_epoch*qut1_apr(2))*15.d03
      qut1_apr(2) = qut1_apr(2)*15.d03
 
****  Get the X-pole positions
      read(unit,'(a)', iostat=ierr ) line
      read(line(34:),*,iostat=ierr) qwob_apr(1,1), qwob_apr(1,2)
 
*     Convert the units from secs and secs/day to mas and mas/day
*     and map values to the correct epoch.
      qwob_apr(1,1) = (qwob_apr(1,1)+dpmu_epoch*qwob_apr(1,2))*1.d03
      qwob_apr(1,2) = qwob_apr(1,2)*1.d03
 
****  Get the Y-pole positions
      read(unit,'(a)', iostat=ierr ) line
      read(line(34:),*,iostat=ierr) qwob_apr(2,1), qwob_apr(2,2)
 
*     Convert the units from secs and secs/day to mas and mas/day
*     and map values to the correct epoch.
      qwob_apr(2,1) = (qwob_apr(2,1)+dpmu_epoch*qwob_apr(2,2))*1.d03
      qwob_apr(2,2) = qwob_apr(2,2)*1.d03
 
****  Now do the nutation angles.
 
      read(unit,'(a)', iostat=ierr ) line
      read(line(34:),*,iostat=ierr) qnut_ang_apr(1,1),
     .                              qnut_ang_apr(1,2)
 
*     Convert the units from secs and secs/day to mas and mas/day
*     and map values to the correct epoch.
      qnut_ang_apr(1,1) = (qnut_ang_apr(1,1) +
     .                     dpmu_epoch*qnut_ang_apr(1,2))*1.d03
      qnut_ang_apr(1,2) = qnut_ang_apr(1,2)*1.d03
 
****  Get the Nutatoin in obliquity positions
      read(unit,'(a)', iostat=ierr ) line
      read(line(34:),*,iostat=ierr) qnut_ang_apr(2,1),
     .                              qnut_ang_apr(2,2)
 
*     Convert the units from secs and secs/day to mas and mas/day
*     and map values to the correct epoch.
      qnut_ang_apr(2,1) = (qnut_ang_apr(2,1)+
     .                     dpmu_epoch*qnut_ang_apr(2,2))*1.d03
      qnut_ang_apr(2,2) = qnut_ang_apr(2,2)*1.d03
 
****  Thats all, get back to reading the rest of the file
      return
      end
 
 
 
