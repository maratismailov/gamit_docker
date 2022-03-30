      subroutine get_clkRT(ep) 

      implicit none

*     Routine to get the clock errors at each station at the current
*     epoch.

      include '../includes/const_param.h' 
      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED
      integer*4 ep      ! Counter with number of epochs processed.

* LOCAL
      integer*4 i, j, k   ! counter
      integer*4 iter   ! Iteration counter

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
      real*8 prev_range(max_rtobs) ! Estimate of L1 range on previous iteration

      real*8 av_clk, sd_clk    ! Average and standard deviation range residual (meters)
     .,      sum_clk, sqr_clk  ! Sum and sum**2 of clock errors
     .,      drng(max_rtobs)   ! Range error (by count)
     .,      tol, max_err      ! Tolerance RMS and max range error

      integer*4 num_clk  ! Number in sum clock
     .,      jmax        ! Index of largest residual
     .,      trimlen     ! Length of string

      logical done  ! Set true when clock converged
     .,       editing  ! Set true while editing range residuals

      character*80 line

***   Loop over the sites we are processing.
      call comp_svs_clkRT( RT_MJD_obs(sblk) )

      do j = sblk, eblk
         prev_range(j) = 20000.d3   ! Approximate range
      end do

      do i = 1, num_site

*         Get the current coordinates of site i
          curr_site_xyz(:,i) = site_int(:,i) + 
     .        site_vel(:,i)*(RT_MJD_obs(sblk)-site_ep(i))/365.25d0

*         Loop over data for this site
          done = .false.
          iter = 0
          clk_error(i) = 0.d0
*         Loop refining the clock error for max of 5 iterations. 
          do while ( .not.done .and. iter.lt.5 )
             sum_clk = 0.d0
             sqr_clk = 0.d0
             num_clk = 0
             iter = iter + 1
             do j = sblk, eblk
****            Get the approximate trn_sec (tranmission in seconds from
*               SP3 start MJD)
                if( rt_sitenum(j).eq.i ) then
*                   OK correct station
* MOD TAH 120615: Remove the site clock error from the receiver time tag.
*                 NOTE: Sign convention here for the clock error is opposite 
*                 to track.
                    trn_time = RT_MJD_obs(j)  - clk_error(i)/86400.d0
     .                       - prev_range(j)/vel_light/86400.d0

                    trn_sec = (RT_GPSWeek(j)*7+44244-sp3_refmjd)*86400 
     .                       + RT_GPSWeeks(j) - clk_error(i)
     .                       - prev_range(j)/vel_light
                            
*                   Call theory with 'Short' option which just computes
*                   the range
                    call theoryRT(trn_time, trn_sec, i, RT_satNum(j), 
     .                    curr_site_xyz(1,i), svs_ear,range, phase,  
     .                    elev, az, dry_map, wet_map,'F', ep)

* MOD TAH 120615: Changed sign of clk_error.  The clock error should 
*                 be subtracted from the theory range.                    
                    drng(j) = RT_P1(j) - range(1)
     .                       +svs_clk(RT_satNum(j))*vel_light
     .                       -clk_error(i)*vel_light
                    prev_range(j) = range(1)  
C                   print *,'OMC ',j,i,RT_satNum(j),drng(j), 
C    .                  RT_P1(j),range(1),
C    .                  svs_clk(RT_satNum(j))*vel_light, 
C    .                  clk_error(i)*vel_light 
                    if( RT_errflag(j).eq.0  ) then  ! add
                       sum_clk = sum_clk + drng(j)
                       sqr_clk = sqr_clk + drng(j)**2
                       num_clk = num_clk + 1
                    endif
                endif
             end do
*            OK Get estimate of clock error and RMS
             editing = .true.
             sd_clk = 999.99d0
             tol = 999.99d0
             do while ( editing )
                if( num_clk.gt.1 ) then
                    av_clk = sum_clk/num_clk
                    sd_clk = sqrt((sqr_clk-av_clk**2*num_clk)/
     .                            (num_clk-1))
                    tol = min(max(sd_clk,10.d0),100.d0)  ! Set 10<sd_tol<100 
*                   Now see if we need to edit data
                    max_err = 0
                    jmax = 0
                    do j = sblk, eblk
                       if( rt_sitenum(j).eq.i .and. 
     .                     RT_errflag(j).eq.0 ) then
                          if( abs((drng(j)-av_clk)/tol).gt.2.5 .and.
     .                        abs((drng(j)-av_clk)/tol).gt.max_err )then
                              max_err = abs((drng(j)-av_clk)/tol)
                              jmax = j
                          endif
                       end if
                    end do
*                   See if we edited
                    if( jmax.gt.0 ) then
                        if( ep.ge.debug(5) .and. ep.le.debug(6) ) 
     .                  write(*,220) jmax,RT_MJD_obs(jmax), 
     .                      site_names(i), RT_satNum(jmax),
     .                      tol, drng(jmax)
 220                    format('CLOCK EDIT Chan ',i3,' MJD ',F13.6, 1x,
     .                      A4,1x,' PRN ',I2.2,1x, ' Tol, error ',F6.2,
     .                      1x,E14.3,' m')
                        call sbit(RT_errflag(jmax),1,1)
*                       Now remove from sum
                        sum_clk = sum_clk - drng(jmax)
                        sqr_clk = sqr_clk - drng(jmax)**2
                        num_clk = num_clk - 1
                    else
                        editing = .false.
                    endif
                else
*                   Not enough data for clock estimate, remove
*                   all data at this time
                    write(line,240) RT_MJD_obs(sblk), site_names(i), 
     .                 tol
 240                format('No clock MJD ',F13.6,' Site ',a4,' Tol ',
     .                 F6.2,' m')
                    call report_stat('WARNING','TRACKRT','get_clk',
     .                 ' ',line,num_rtobs)
                    editing = .false.
                    done = .true.
                    do j = sblk, eblk
                       if( rt_sitenum(j).eq.i .and. 
     .                    RT_errflag(j).eq.0 ) then
                          call sbit(RT_errflag(j),1,1)
                       end if
                    end do
                    av_clk = 0.d0
                endif
             end do     !   Editing clock estimates
*            Now see if clock estimate needs updating
             if( abs(av_clk).gt.1.d-3 ) then 
*               Update clock
                clk_error(i) = clk_error(i) + av_clk/vel_light
             else
*               OK All done
                done = .true.
                clk_error(i) = clk_error(i) + av_clk/vel_light
                if( ep.ge.debug(3) .and. ep.le.debug(4) )
     .          write(*,320) RT_MJD_obs(sblk), site_names(i), num_clk,  
     .             av_clk, sd_clk, clk_error(i), iter
 320            format('CLOCK MJD ',F13.6,1x,a4,1x,I4,' AV/STD ',
     .             F13.2,1x,F13.2,' m, CLOCK ',E13.3,' sec',I2)

*               Save the clock error in the sec_offset array
                sec_offset(i,1) = clk_error(i)
 
*               Output residuals (if debug)
                if( ep.ge.debug(5) .and. ep.le.debug(6) ) then                
                   do j = sblk, eblk
                      if( rt_sitenum(j).eq.i ) then
                         if( RT_errflag(j).eq.0 ) then
                             write(*,340) RT_MJD_obs(j),site_names(i),
     .                           RT_satNum(j), drng(j),
     .                           RT_P1(j)-RT_P2(j) 
 340                         format('P1_RRES MJD ',F13.6,1x,a4,' PRN ',
     .                           I2.2,' Res,  P1-P2 ',2F13.3,' m')
                         else
                             write(*,360) RT_MJD_obs(j),site_names(i),
     .                           RT_satNum(j), drng(j),
     .                           RT_P1(j)- RT_P2(j) 
 360                         format('P1_BADR MJD ',F13.6,1x,a4,' PRN ',
     .                           I2.2,' Res, P1-P2 ',2G16.6,' m')
                         end if
                      endif
                   end do
                end if
             end if
          end do
       end do

*****  Thats all
       return
       end

CTITLE DATA_OK

      logical function data_ok(flag,mask)

      implicit none

*     Routine which returns true is none of the bits in mask are
*     set in flag.

* PASSED VARIABLES
* flag  -- data flag being checked
* mask  -- Masking bits to check only those for bad data

      integer*4 flag, mask

* cand  -- Common library version of the and function

      integer*4 cand

****  OK see if any bits match
      data_OK = .true.
      if( cand(flag,mask).gt.0 ) data_OK = .false.

****  Thats all
      return
      end


          
