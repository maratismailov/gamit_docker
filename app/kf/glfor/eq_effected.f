CTITLE eq_effected
 
      logical function eq_effected(ns)

      implicit none 
 
*     Routine to test if site is effected by an earthquake
*     in some fashion.
* WARNING * This routine is meant only for controlling output
*     and should not be for actual control becuase the pre-
*     and post-seismic intevals are extended to get some points
*     not in the earthquake effected region.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* PASSED variables
 
*   ns      - Site number from global list

      integer*4 ns
 
* LOCAL VARIABLES
 
*   i       - Loop counter
 
      integer*4 i
 
*   dt_eq   - Time difference (days) between earthquake
*           - and current experiment mid-point
*   dist    - Distance from site to earthquake
*   coseismic_time(max_eq) - Keeps track of epoch at which
*             coesimic offset was output.  This is locally
*             saved
 
      real*8 dt_eq, dist, coseismic_time(max_eq)

*   First_call - Indicates that this is first call. USed to
*             clear the coseismic_time array
   
      logical first_call

      save coseismic_time, first_call

      data first_call / .true. /

****  If first call clear the cosemic time array
      if( first_call ) then
         do i = 1, num_eq
             coseismic_time(i) = 0.d0
         end do
         first_call = .false.
      end if
 
****  First scan the earthquakes to see if we need to anything

      eq_effected = .false.
 
      do i = 1, num_eq
 
*         Get time difference from earthquake
          dt_eq = gepoch_expt - eq_epoch(i)
 
*         See if we need to do anything with coseismic displacement
          if( sign(1.d0,deltat)*dt_eq.gt.0 .and.
     .       (coseismic_time(i).eq.0 .or. gepoch_expt.eq.
     .        coseismic_time(i) )       ) then
 
*             We just passed the Earthquake epoch so apply the
*             covariance increment to all sites effected by the
*             earthquake
              call eval_dist(eq_pos(1,i), apr_val_site(1,1,ns), dist)
              coseismic_time(i) = gepoch_expt
              if( dist.le.eq_rad(i) ) then
                  eq_effected = .true.
              end if 

*             Extend interval by two days 
          else if( dt_eq.lt.eq_dur(2,i)+2.0 .and.
     .             dt_eq.gt.0.d0 ) then
 
*             Post-seimic deformation process
              call eval_dist(eq_pos(1,i), apr_val_site(1,1,ns), dist)
              if( dist.le.eq_rad(i) .and. 
     .          (eq_mar_post(1,i)+eq_mar_post(2,i)+eq_mar_post(3,i)+
     .           eq_mar_post(4,i)+eq_mar_post(5,i)+eq_mar_post(6,i))
     .              .gt.0 ) then
                    eq_effected = .true.
              end if
 
*             Extend interval by two days 
          else if(  -dt_eq. lt. eq_dur(1,i)+2.0 .and.
     .               dt_eq.lt.0.d0  ) then
 
*             Pre-seismic deformation process.
              call eval_dist(eq_pos(1,i), apr_val_site(1,1,ns), dist)
              if( dist.le.eq_rad(i) .and. 
     .          (eq_mar_pre(1,i)+eq_mar_pre(2,i)+eq_mar_pre(3,i)+
     .           eq_mar_pre(4,i)+eq_mar_pre(5,i)+eq_mar_pre(6,i))
     .              .gt.0 ) then
                  eq_effected = .true.
              end if
          end if

* MOD TAH 980519: Set bit eq_used to show earthquake used.
          if( eq_effected ) then
              call sbit( eq_used, i, 1)
          end if
      end do
 
****  Thats all
      return
      end
 
