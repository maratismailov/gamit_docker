CTITLE SCAN_PARAM
 
      subroutine scan_parm

      implicit none  
 
*     Routine to scan the estimated parameters for a solution and
*     save some information we need.  In particular, we check
*     to see if nutation angles, and PM/UT1 have been estimated
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
 
      include '../includes/glb_hdr_def.h'
 
*   i       - Loop counter
*   indx    - Index for parameter
*   type    - type of parameter
 
      integer*4 i, indx, type, subin
 
****  See if nutation angles or PM/UT1 have been estimated
 
      nut_ang_est = .false.
      pmu_est     =  0      
      tran_est    =  0
      rot_est     =  0
      scale_est   =  0
 
      do i = 1, cnum_parn
          call decode_code( gpar_codes(i), type, indx)
          subin       = mod(indx,256) 
 
*         Check nutation angles
          if( type.eq.15 ) nut_ang_est = .true.
*         Check nutation series coefficients
          if( type.eq.47 ) nut_ang_est = .true.
 
*         Check wobble
          if( type.eq.13 .and. indx.eq.1 ) call sbit(pmu_est,1,1)
          if( type.eq.13 .and. indx.eq.2 ) call sbit(pmu_est,2,1)
          if( type.eq.13 .and. indx.eq.3 ) call sbit(pmu_est,4,1)
          if( type.eq.13 .and. indx.eq.4 ) call sbit(pmu_est,5,1)
*         Check UT1
          if( type.eq.14 .and. indx.eq.1 ) call sbit(pmu_est,3,1)
          if( type.eq.14 .and. indx.eq.2 ) call sbit(pmu_est,6,1)

*         Check for multi-epoch polar motion
          if( type.eq.56 .and. subin.eq.1 ) call sbit(pmu_est,1,1)
          if( type.eq.57 .and. subin.eq.1 ) call sbit(pmu_est,2,1)
          if( type.eq.58 .and. subin.eq.1 ) call sbit(pmu_est,3,1)

          if( type.eq.56 .and. subin.eq.2 ) call sbit(pmu_est,4,1)
          if( type.eq.57 .and. subin.eq.2 ) call sbit(pmu_est,5,1)
          if( type.eq.58 .and. subin.eq.2 ) call sbit(pmu_est,6,1)

*         Check translation
          if( type.eq.52 .and. indx.eq.1 ) call sbit(tran_est,1,1)
          if( type.eq.52 .and. indx.eq.2 ) call sbit(tran_est,2,1)
          if( type.eq.52 .and. indx.eq.3 ) call sbit(tran_est,3,1)
*         Check translation rate of change
          if( type.eq.53 .and. indx.eq.1 ) call sbit(tran_est,4,1)
          if( type.eq.53 .and. indx.eq.2 ) call sbit(tran_est,5,1)
          if( type.eq.53 .and. indx.eq.3 ) call sbit(tran_est,6,1)

* MOD TAH 030615: Check to see if log are included in the solutions
*         Save as bits in the trsns_est values (invoked when combined
*         globals are read)
          if( type.eq.62 ) call sbit(tran_est,16,1)
          if( type.eq.63 ) call sbit(tran_est,17,1)
          if( type.eq.64 ) call sbit(tran_est,18,1)

*         Check rotation
          if( type.eq.59 .and. indx.eq.1 ) call sbit(rot_est,1,1)
          if( type.eq.59 .and. indx.eq.2 ) call sbit(rot_est,2,1)
          if( type.eq.59 .and. indx.eq.3 ) call sbit(rot_est,3,1)
*         Check rotation  rate of change
          if( type.eq.60 .and. indx.eq.1 ) call sbit(rot_est,4,1)
          if( type.eq.60 .and. indx.eq.2 ) call sbit(rot_est,5,1)
          if( type.eq.60 .and. indx.eq.3 ) call sbit(rot_est,6,1)

*         Check the scale
          if( type.eq.54 ) call sbit(scale_est,1,1)
          if( type.eq.55 ) call sbit(scale_est,2,1)
 
      end do
 
***** Thats all
      return
      end
 
