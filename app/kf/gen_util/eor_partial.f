
CTITLE EOR_PARTIAL
 
      subroutine eor_partial

      implicit none 
 
*     Routine to compute the diurnal and semi UT1 partials and
*     The retrograde diurnal PM, and the prograde and retrograde
*     semidiunal polar motion partials
*
*     The layout of the partials is:
*     EOR_UT1_PART:
*     1    - Diurnal UT1 cos(gst)
*     2    - Diurnal UT1 sin(gst)
*     3    - Semidiurnal UT1 cos(2*gst)
*     4    - Semidiurnal UT1 sin(2*gst)

*     EOR_XY_PART:
*     1    - Diurnal retrograde cos(gst)
*     2    - Diurnal retrograde sin(gst)
*     3    - Semidiurnal prograde cos(gst) [Same direction as nutation]
*     4    - Semidiurnal prograde sin(gst)  ""
*     5    - Semidiurnal retrograde cos(gst)
*     6    - Semidiurnal retrograde sin(gst)

      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
      include '../includes/const_param.h'
 
 
 
*          j          - Loop counter for delay and rate
 
      integer*4 j
 
*       arg           - GST of observation (rad)
*       arg_semi      - Semidiurnal argument (2*arg)
*       arg_diur      - Diurnal argument (arg)
 
 
      real*8 arg, arg_semi, arg_diur
 
***   calculate phase of mode (reference epoch JD 2451545.0)
      call gst_jd(epoch, arg)

      arg_diur = mod(arg,2*pi)
      arg_semi = mod(2*arg, 2*pi)

****  Do the UT1 partials
*                     ! Loop over delay and rate
      do j = 1, 2

          eor_ut1_part(1,j) = pmu_part(3,j) * cos(arg_diur)
          eor_ut1_part(2,j) = pmu_part(3,j) * sin(arg_diur)
          eor_ut1_part(3,j) = pmu_part(3,j) * cos(arg_semi)
          eor_ut1_part(4,j) = pmu_part(3,j) * sin(arg_semi)

      end do

****  Now do retrorgrade PM (opposite to nutation)
*                     ! Loop on delay and rate.
      do j = 1, 2
          eor_xy_part(1,j) = -pmu_part(1,j)* cos(arg_diur) +
     .                        pmu_part(2,j)* sin(arg_diur)
          eor_xy_part(2,j) =  pmu_part(2,j)* cos(arg_diur) +
     .                        pmu_part(1,j)* sin(arg_diur)

*         Renaming parials are semidiurnal.
*         Same direction as nutation first. Prograde.
          eor_xy_part(3,j) = -pmu_part(1,j)* cos(arg_semi) -
     .                        pmu_part(2,j)* sin(arg_semi)
          eor_xy_part(4,j) =  pmu_part(2,j)* cos(arg_semi) -
     .                        pmu_part(1,j)* sin(arg_semi)

          eor_xy_part(5,j) = -pmu_part(1,j)* cos(arg_semi) +
     .                        pmu_part(2,j)* sin(arg_semi)
          eor_xy_part(6,j) = +pmu_part(2,j)* cos(arg_semi) +
     .                        pmu_part(1,j)* sin(arg_semi)

      end do
 
****  Thats all
      return
      end
 
 
 
