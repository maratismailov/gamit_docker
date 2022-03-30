      subroutine comp_rawomc( ep )

      implicit none

*     Routine to compute raw omc for L1/L2 phase (cycles) and P1/P2
*     range (meters).

      include '../includes/const_param.h'
      include '../../libraries/includes/const_freq.h'
      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep            ! Counter for number of epochs of data processed.

* LOCAL
      integer*4 i, j, k   ! counter
      integer*4 np        ! Counter for parameter number

* trn_time  -- Transmitt time as given by satellite clock (MJD)
* trn_sec   -- Transmit time as seconds from start of SP3 file
* site_pos(3) -- Approximate position of the site (XYZ, m) Earth fixed
* svs_ear(3)  -- Position of the satellite in earth fixed frame,
* range(2)    -- Range to the satellite at the L1 and L2 frequencies (m)
* phase(2)    -- Phase for the satellite at L1 and L2 freqs (cycles at each freq)
* elev        -- Elevation angle to satellite (deg)
* dry_map     -- Dry mapping function
* wet_map     -- Wet mapping function

      real*8 trn_time, trn_sec, site_pos(3), svs_ear(3), range(2), 
     .       phase(2), elev, az, dry_map, wet_map 

      do i = 1, num_site

*         Get the current coordinates of site i
          curr_site_xyz(:,i) = site_int(:,i) + 
     .        site_vel(:,i)*(RT_MJD_obs(sblk)-site_ep(i))/365.25d0

          do j = sblk, eblk
****         Get the approximate trn_sec (tranmission in seconds from
*            SP3 start MJD)
             if( rt_sitenum(j).eq.i ) then
*                OK correct station
                 trn_time = RT_MJD_obs(j)  - clk_error(i)/86400.d0

                 trn_sec = (RT_GPSWeek(j)*7+44244-sp3_refmjd)*86400 
     .                    + RT_GPSWeeks(j) - clk_error(i)
                         
*                Call theory with 'Short' option which just computes
*                the range
                 call theoryRT(trn_time, trn_sec, i, RT_satNum(j), 
     .                 curr_site_xyz(1,i), svs_ear,range, phase,  
     .                 elev, az, dry_map, wet_map,'S', ep)

*                Now do the full calculation with the correction
*                for the light-propagation time
                 
                 trn_time = RT_MJD_obs(j)  - sec_offset(i,1)/86400.d0
     .                    - range(1)/vel_light/86400.d0

                 trn_sec = (RT_GPSWeek(j)*7+44244-sp3_refmjd)*86400 
     .                    + RT_GPSWeeks(j) - sec_offset(i,1)
     .                    - range(1)/vel_light
                         
*                Call theory with 'Short' option which just computes
*                the range
                 call theoryRT(trn_time, trn_sec, i, RT_satNum(j), 
     .                 curr_site_xyz(1,i), svs_ear,range, phase,  
     .                 elev, az, dry_map, wet_map,'F', ep)

*                Now compute the raw omc values.  Note: We do not apply
*                station clock correction directly here.
                 RT_omcraw(1,j) = RT_L1(j) - phase(1)
     .                + svs_clk(RT_satNum(j))*fL1
                 RT_omcraw(2,j) = RT_L2(j) - phase(2)
     .                + svs_clk(RT_satNum(j))*fL2
                 RT_omcraw(3,j) = RT_P1(j) - range(1)
     .                + svs_clk(RT_satNum(j))*vel_light
                 RT_omcraw(4,j) = RT_P2(j) - range(1)
     .                + svs_clk(RT_satNum(j))*vel_light

*                Now check P1-P2 range and flag if too large
                 if( abs(RT_P1(j)-RT_P2(j)).gt.100.d0 ) then
                     call sbit(RT_errflag(j),2,1)
                 endif

*                Save AZ and Elevation
                 RT_azel(1,j) = az
                 RT_azel(2,j) = elev

*****            See if data below elevation limit
                 if( RT_azel(2,j).lt. elev_cutoff ) then
                     call sbit(RT_errflag(j),5,1)
C                     write(*,220) ep, site_names(i), RT_satNum(j), 
C     .                   RT_azel(2,j),RT_errflag(j) 
 220                 format('ElevCut Ep ',I6,1x,a4,' PRN ',I2.2,
     .                   F6.2,' ErrFlag ',o4) 
                 end if

****             See if PRN is in excluded list (-ve ss_exclude are
*                from not being in SP3 file).
                 do k = 1, num_exclude
                    if( RT_satNum(j).eq.abs(ss_exclude(k)) ) then
                       call sbit(RT_errflag(j),6,1)
                    end if
                 end do


****             Save the partial derivatives (some care need here with
*                differece btween max_obs and max_rtobs; partials 
*                dimensioned to max_obs
                 do k = 1,max_parm
                    ow_part(k,j) = 0.d0
                 end do
*                Now site partial
                 do k = 1,3
                    if( site_parn(k,i).ne.0 ) then
                        np = site_parn(k,i)
                        ow_part(np,j) = (curr_site_xyz(k,i)
     .                                  -svs_ear(k))/range(1)
                    end if
                 end do
                 np = atm_parn(i)
                 if( np.gt.0 ) then
                     ow_part(np,j) = wet_map
                 end if

             end if
          end do
      end do
                 
****  Thats all
      return
      end


