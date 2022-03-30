ftn77 
      subroutine get_num_sites
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the number of sites in this experiment.
* 
*     T.Herring                    3:27 PM  TUE., 13  JAN., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
  
*     Scratch common
  
      common err, nd1, nd2, nd3 
  
***** See if # sites is available 
  
      call get_int(Lnum_sites  ,Fnum_sites  ,num_sites,1,1,1,err) 
  
***** See if too many 
  
      call check_max(Lnum_sites, num_sites, max_sites, icrt, log_lu)
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_num_sources
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the number of sources in this experiment.
* 
*     T.Herring                    3:27 PM  TUE., 13  JAN., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
  
*     Scratch common
  
      common err, nd1, nd2, nd3 
  
***** See if # sources is available 
  
      call get_int(Lnum_sources,Fnum_sources,num_sources,1,1,1,err) 
  
***** See if too many 
  
      call check_max(Lnum_sources,num_sources,max_sources, icrt, log_lu)
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_site_names 
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the names of the sources in this experiment. 
*     The names firstly have the blanks replaced with underscores 
*     and then the alias name of the sites is found.
* 
*     T.Herring                    3:27 PM  TUE., 13  JAN., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
     .,   work_array(4,max_sites) ! Array for temporary storage 
                          ! site names
*     Scratch common
  
      common err, nd1, nd2, nd3, work_array 
  
***** See if names of the sites are available 
  
      call get_asc(Lsite_names  ,Fsite_names  ,site_names,
     .             4,num_sites,1,work_array, err) 
  
***** Replace all blanks with under_score characters
  
      call underscore(site_names, num_sites)
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
  
      subroutine get_source_names 
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the names of the sites in this experiment. 
*     The names firstly have the blanks replaced with underscores 
*     and then the alias name of the sources is found.
* 
*     T.Herring                    3:27 PM  TUE., 13  JAN., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
     .,   work_array(4,max_sources) ! Array for temporary storage 
                          ! site names
  
*     Scratch common
  
      common err, nd1, nd2, nd3, work_array 
  
***** See if names of the sources are available 
  
      call get_asc(Lsource_names  ,Fsource_names  ,source_names,
     .             4,num_sources,1,work_array, err) 
  
***** Replace all blanks with under_score characters
  
      call underscore(source_names, num_sources)
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_CALC_ver 
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the version of CALC used to process this experiment
*     If the version is less than 6.0 then the atm_cont is subtracted 
*     from the theoretical delays 
* 
*     T.Herring                    3:27 PM  TUE., 13  JAN., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
  
      real*6
     .    work_6         ! Work array used to read database 
  
  
*     Scratch common
  
      common err, nd1, nd2, nd3, work_6 
  
***** See if CALC VER is available
  
      call get_6_to_4(LCALC_ver  ,FCALC_ver  ,1,1,1,err, work_6, 1, 
     .            CALC_ver, icrt, log_lu) 
  
***** Write out CALC ver
  
      call write_version(CALC_ver, num_obs, icrt, log_lu) 
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_ellip_hgt
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the ellipsoidal heighted of the sites for
*     this experiment.  (These values are used mainly for atm.
*     delay calculation.) 
* 
*     T.Herring                   08:52 PM TUE., 13 Jan., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
  
      real*6
     .    work_6(max_sites)   ! Used to retrieve data from database 
  
*     Scratch common
  
      common err, nd1, nd2, nd3, work_6 
  
***** See if SIT ELEV is available
  
C     call get_6_to_4(Lellip_hgt  ,Fellip_hgt  ,num_sites,1,1,err,
C    .             work_6, max_sites, ellip_hgt, icrt, log_lu ) 
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_site_pos 
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the site positions of the sites for
*     this experiment.
*     This routine will also compute the GEOD lat, long and height
*     of the site.
* 
*     T.Herring                   08:52 PM TUE., 13 Jan., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include const_param.ftni  , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
  
      integer*2 
     .    i               ! Loop counter
     .,   err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
  
      real*6
     .    work_6(3,max_sites)   ! Used to retrieve data from database 
  
*     Scratch common
  
      common err, nd1, nd2, nd3, work_6 
  
***** See if SITERECS is available
  
      call get_6_to_8(Lsite_pos  ,Fsite_pos  , 3,num_sites,1,err, 
     .             work_6, 3*max_sites, site_pos, icrt, log_lu )
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_geod_pos 
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to convert XYZ site positions into GEOD lat, long 
*     and height.  This routine must be called after the site positions 
*     have been sorted since the aplha_site routine will not sort 
*     the geodetic coordinates
* 
*     T.Herring                   08:52 PM TUE., 13 Jan., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include const_param.ftni  , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
  
      integer*2 
     .    i               ! Loop counter
  
      real*8
     .    geod_coords(3)  ! Geodetic coordinates of each site 
     .,   rot_mat(3,3)    ! Rotation matrix between XYZ and GEOD
  
*     Scratch common
  
      common geod_coords, rot_mat 
  
***** Now get the geodetic coordinates
  
      do i = 1, num_sites 
          call XYZ_to_GEOD( rot_mat, site_pos(1,i), geod_coords)
  
*         Save the values in the header 
  
          ellip_hgt(i)  = geod_coords(3)
          latitudes(i)  = pi/2.d0 - geod_coords(1)
          longitudes(i) = geod_coords(2)
  
      end do
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
      subroutine get_source_pos 
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the source positions of the sources for
*     this experiment.
* 
*     T.Herring                   08:52 PM TUE., 13 Jan., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
      include const_param.ftni  , nolist      ! Constants (rad_to_mas)
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   i               ! Loop counter
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
  
      real*6
     .    work_6(2,max_sources)   ! Used to retrieve data from database 
  
*     Scratch common
  
      common err, nd1, nd2, nd3, work_6 
  
***** See if STAR2000 is available
  
      call get_6_to_8(Lsource_pos  ,Fsource_pos  , 2,num_sources,1,err, 
     .             work_6, 2*max_sources, source_pos, icrt, log_lu )
  
***** Now convert these into Kalman filter units (mas for both RA and 
*     dec)
  
      do i = 1, num_sources 
          source_pos(1,i) = source_pos(1,i) * rad_to_mas
  
* MOD TAH 870921 Check for negative Right Ascension 
          if( source_pos(1,i).lt.0 ) then 
              source_pos(1,i) = source_pos(1,i) + 1296.0d6 ! 360 deg in mas 
          end if
          source_pos(2,i) = source_pos(2,i) * rad_to_mas
      end do
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_etd_apr
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the apriori values for the earth tide parameters     
*     this experiment.  
* 
*     T.Herring                   08:52 PM TUE., 13 Jan., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
      include const_param.ftni  , nolist      ! Constants (rad_to_mas)
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
  
      real*6
     .    work_6(3)       ! Used to retrieve data from database 
  
*     Scratch common
  
      common err, nd1, nd2, nd3, work_6 
  
***** See if EDT DATA is available
  
      call get_6_to_8(Lapr_etd  ,Fapr_etd  , 3, 1,1,err,
     .             work_6, 3, etd_apr, icrt, log_lu ) 
  
***** Convert the Lag angle into degs 
  
      etd_apr(1) = etd_apr(1) * rad_to_deg
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_gamma_apr
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the apriori value of gamma used in 
*     this experiment.
* 
*     T.Herring                   08:52 PM TUE., 13 Jan., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
  
      real*6
     .    work_6        ! Used to retrieve data from database 
  
*     Scratch common
  
      common err, nd1, nd2, nd3, work_6 
  
***** See if REL DATA is available
  
      call get_6_to_8(Lapr_gamma  ,Fapr_gamma  , 1,1,1,err, 
     .             work_6, 1, gamma_apr, icrt, log_lu ) 
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_delay_flag 
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the delay type flag for this experiment.  If 
*     the delay flag is one then the delays will be read as group 
*     delays.  If the delay is 2 then phase delays will be read.
* 
*     T.Herring                    3:27 PM  TUE., 13  JAN., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
  
  
*     Scratch common
  
      common err, nd1, nd2, nd3 
  
***** See if DELTFLAG is available
  
      call get_int(Ldelay_flag  ,Fdelay_flag  ,delay_flag,1,1,1,err)
  
***** See if it matches the expected values 
  
      call check_delay_flag( delay_flag, icrt, log_lu, 0) 
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_axis_offsets 
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the axis offset values and type for the sites  
*     in this experiment.   
*     The axis types are (from CALC AXOG) 
*     1   -- Equatorial mount 
*     2   -- XY mount fixed axis points North-South 
*     3   -- ALT-AZ mount 
*     4   -- XY mount fixed axis points East-West 
*     5   -- Special case for Richmond Fl.
* 
*     T.Herring                   08:52 PM TUE., 13 Jan., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
  
      real*6
     .    work_6(max_sites)   ! Used to retrieve data from database 
  
*     Scratch common
  
      common err, nd1, nd2, nd3, work_6 
  
  
***** See if AXISOFFS and AXISTYPS are available
  
      call get_6_to_8(Laxis_offsets  ,Faxis_offsets ,num_sites,1,1,err, 
     .             work_6, num_sites, axis_offsets, icrt, log_lu )
      call get_int(Laxis_types,Faxis_types,axis_types,num_sites,1,1,
     .             err )
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_UT1
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the UT1 value and rate of change at the
*     midepoch of the experiment.  To do this the UT1 tables
*     are read from the database and re-interpolated at the 
*     midpoint epoch using a 4 pt Lagrangian interpolation. 
*     The rate of change is computed using the derivative of
*     the Lagranian interpolator. 
*     REF: Abramowitz, M., and I.A. Stegun, Handbook of mathematical
*         functions, Dover Publ. Inc., New York, 1972. Pages 879 and 883. 
*         (pp1046)
* 
*     During this routine the sign the database entries is changed
*     so that we obtain UT1-AT (rather than the reverse convention) 
*     (Partials are also changed when the data is read) 
* 
* MOD JLD 870602 To change the sign on the entries in the UT1 table 
*     so that they also represent UT1-TAI, constistent with the 
*     other UT1 definitions 
* 
*     The PMU_EPOCH is set in this routine also.
* 
*     RESTRICTION:  This routine assumes that MID_EPOCH is available
*     Therefore the database must be read once before this routine is 
*     is called (Done in GET_NUM_OBS) 
* 
*     T.Herring                   10:18 AM WED., 14 Jan., 1987
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
      include const_param.ftni  , nolist
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
     .,   i               ! Loop counter
  
      real*6
     .    work_6(5)       ! Used to retrieve data from database 
                          ! (5 values for the data points)
      real*8
     .    tai_utc_array(3) ! Array of values for TAI-UTC.  First is JD, 
                           ! then AT-UT1 (sec) and then rate offset.
                           ! We only use the TAI-UTC value (saved as
                           ! integer) 
  
*     Scratch common
  
      common err, nd1, nd2, nd3, work_6 
  
  
***** See if FUT1 INF and FUT1 PTS are available
  
      call get_6_to_8(LFut1_inf  ,FFut1_inf , 4,1,1,err,
     .             work_6, 4, Fut1_inf, icrt, log_lu )
      call get_6_to_8(LFut1_pts  ,FFut1_pts , 5,1,1,err,
     .             work_6, 5, Fut1_pts, icrt, log_lu )
  
*     See if found final values, if not try preliminary 
  
      if( .not. FFut1_inf ) then  ! Try preliminary 
          call get_6_to_8(LPut1_inf  ,FPut1_inf , 4,1,1,err,
     .                 work_6, 4, Fut1_inf, icrt, log_lu )
          call get_6_to_8(LPut1_pts  ,FPut1_pts , 5,1,1,err,
     .                 work_6, 5, Fut1_pts, icrt, log_lu )
  
*         If we did not find preliminary try extrapolated 
  
          if( .not. FPut1_inf ) then  ! Try extrapolted 
              call get_6_to_8(LEut1_inf  ,FEut1_inf , 4,1,1,err,
     .                    work_6, 4, Fut1_inf, icrt, log_lu ) 
              call get_6_to_8(LEut1_pts  ,FEut1_pts , 5,1,1,err,
     .                    work_6, 5, Fut1_pts, icrt, log_lu ) 
  
              if( FEut1_inf ) call sbit(data_notes,5,1) 
          else
              call sbit(data_notes,4,1)   ! Say we found preliminary
          end if
      else
          call sbit(data_notes,3,1)       ! Say we found finals 
      end if
  
***** Now get the interpolated value and rate at the mid_epoch
*     of the experiment 
  
      call Lagrange_intp( Fut1_inf, Fut1_pts, mid_epoch, UT1_apr, 1)
  
***** Now convert the units to mas and set at UT1-TAI 
  
      UT1_apr(1) = -UT1_apr(1)*Fut1_inf(4) *15.d0*1000.d0 
      UT1_apr(2) = -UT1_apr(2)*Fut1_inf(4) *15.d0*1000.d0 
  
* MOD JLD 870602 Change the sign on the table entries to represent
*     UT1-TAI 
      do i = 1, fut1_inf(3) 
          fut1_pts(i) = - fut1_pts(i)   ! Units are still seconds 
      end do
  
***** Now get difference bewteen AT and UT1 
  
      call get_6_to_8(Ltai_utc    ,Ftai_utc , 3,1,1,err,
     .             work_6, 3, tai_utc_array, icrt, log_lu ) 
  
*     Now save the offset 
  
      tai_utc = anint(tai_utc_array(2)) 
  
***** Save the value of the epoch 
  
      pmu_epoch = mid_epoch 
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine Lagrange_intp( INF, PTS, epoch, value, dim)
  
      implicit none 
  
*-------------------------------------------------------------------- 
*     4 point Lagrange (equally spaced) interpolation routine.
*     Routine also computes the derivative of the function. 
*     The routine will simualtaneosly interpolate several 
*     difference type of data.  The number of types of data 
*     is given by variable DIM. 
* 
*     The derivative is returned in units of units of the 
*     tabluar points per unit of the spacing entries (usally days)
* 
*     REF: Abramowitz, M., and I.A. Stegun, Handbook of mathematical
*         functions, Dover Publ. Inc., New York, 1972. Pages 879 and 883. 
*         (pp1046)
* 
*     T.Herring                   11:02 AM WED., 14 Jan., 1987
*-------------------------------------------------------------------- 
  
      integer*2 
     .    dim         ! Number of types of data to interpolate. 
     .,   i,j,k       ! Loop counters 
     .,   st_index    ! the start index in the table to be used.  This
                      ! is the position in PTS of the first point 
                      ! imediately before EPOCH.
  
  
      real*8
     .    coeff(-1:2) ! The 4 coefficients used for the interpolation 
     .,   dt          ! Difference in time bewteen the first data point 
                      ! to be used in the interpolation and the epoch 
                      ! In units of the spacing of the table. 
     .,   epoch       ! Epoch at which we want intepolated value
  
     .,   INF(3)      ! Gives information about the spacing of
                      ! the data points.  The values are
                      ! 1 -- Julian date of first point.
                      ! 2 -- Spacing of the data in days
                      ! 3 -- Number of points in table
     .,   PTS(dim,1)  ! The tabluar values of the quantity to be
                      ! interolated.
  
     .,   value(dim,2)! Interpolated value and its derivative at
                      ! 'epoch'.
  
***** START, Determine the index of the first data point imediately 
*     before the epoch, and get the time difference 
  
      st_index = (epoch-INF(1))/INF(2) + 1
      if( st_index-1.lt.1      ) st_index = 2   ! Make sure we use valid entry
      if( st_index+2.gt.inf(3) ) st_index = inf(3) - 2
  
      dt       = (epoch - (INF(1)+(st_index-1)*INF(2)) )/INF(2) 
  
***** Now compute the coefficients for the intepolation 
  
      coeff(-1) = - dt*(dt-1)*(dt-2)/6
      coeff( 0) =  (dt**2-1)*(dt-2)/2 
      coeff( 1) = - dt*(dt+1)*(dt-2)/2
      coeff( 2) =  (dt**2-1)*dt/6 
  
*     Get interpolated value
  
      do i = 1,dim        ! Loop over different types 
  
          value(i,1) = 0.d0 
          do j = -1,2 
              value(i,1) = value(i,1) + coeff(j)*PTS(i,st_index+j)
          end do
      end do
  
***** Now do the derivative 
  
      coeff(-1) = -(3*dt**2-6*dt+2)/6 
      coeff( 0) =  (3*dt**2-4*dt-1)/2 
      coeff( 1) = -(3*dt**2-2*dt-2)/2 
      coeff( 2) =  (3*dt**2-1     )/6 
  
*     Compute the derivative
  
      do i = 1, dim 
  
          value(i,2) = 0.d0 
          do j = -1,2 
              value(i,2) = value(i,2) + coeff(j)*PTS(i,st_index+j)
          end do
  
*         Get into the correct units
          value(i,2) = value(i,2)/INF(2)   ! divide by data spacing 
  
      end do
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
      subroutine get_wobble 
  
      implicit none 
  
*---------------------------------------------------------------------- 
*     Routine to get the wobble valuse and their rates of change at the 
*     midepoch of the experiment.  To do this the wobble tables 
*     are read from the database and re-interpolated at the 
*     midpoint epoch using a 4 pt Lagrangian interpolation. 
*     The rate of change is computed using the derivative of
*     the Lagranian interpolator. 
*     REF: Abramowitz, M., and I.A. Stegun, Handbook of mathematical
*         functions, Dover Publ. Inc., New York, 1972. Pages 879 and 883. 
*         (pp1046)
* 
*     Units are changed to mas in this routine. 
* 
*     RESTRICTION:  This routine assumes that MID_EPOCH is available
*     Therefore the database must be read once before this routine is 
*     is called (Done in GET_NUM_OBS) 
* 
*     T.Herring                   12:21 PM WED., 14 Jan., 1987
* 
* MOD J. Davis 870225 To correct error made by multiplying PM values
*                     by 1000 to convert to mas.  PM values are already 
*                     in milliarcseconds in the database. 
*---------------------------------------------------------------------- 
  
      include kalman_param.ftni , nolist
      include readin_user.ftni  , nolist
      include obs_header.ftni   , nolist
      include const_param.ftni  , nolist
  
      integer*2 
     .    err             ! err returned from CHECH and GET_INT 
     .,   nd1,nd2,nd3     ! Dimensions returned from ASK command
  
      real*6
     .    work_6(10)      ! Used to retrieve data from database 
                          ! (5 values for the data points)
  
*     Scratch common
  
      common err, nd1, nd2, nd3, work_6 
  
  
***** See if FWOB INF and FWOB PTS are available
  
      call get_6_to_8(LFwob_inf  ,FFwob_inf , 3,1,1,err,
     .             work_6, 3, Fwob_inf, icrt, log_lu )
      call get_6_to_8(LFwob_pts  ,FFwob_pts , 2,5,1,err,
     .             work_6, 10, Fwob_pts, icrt, log_lu ) 
  
***** If they are not available then try preliminary values 
  
      if( .not.FFwob_inf ) then 
          call get_6_to_8(LPwob_inf, FPwob_inf, 3,1,1,err,
     .            work_6, 3, Fwob_inf, icrt, log_lu)
          call get_6_to_8(LPwob_pts, FPwob_pts, 2,5,1,err,
     .            work_6, 10, Fwob_pts, icrt, log_lu) 
  
*         See if preliminary values were available, if not get extroplated
*         values
          if( .not.FPwob_inf ) then 
              call get_6_to_8(LEwob_inf, FEwob_inf, 3,1,1,err,
     .                work_6, 3, Fwob_inf, icrt, log_lu)
              call get_6_to_8(LEwob_pts, FEwob_pts, 2,5,1,err,
     .                work_6, 10, Fwob_pts, icrt, log_lu) 
  
              if( FEwob_inf ) call sbit(data_notes,5,1) 
          else
              call sbit(data_notes,4,1)   ! Say we used PRELIMINARY values
          end if
      else
          call sbit(data_notes,3,1)       ! Say we used FINAL values
      end if
  
***** Now get the interpolated value and rate at the mid_epoch
*     of the experiment 
  
      call Lagrange_intp( Fwob_inf, Fwob_pts, mid_epoch, wob_apr, 2)
  
***** Save the epoch
  
      pmu_epoch = mid_epoch 
  
***** Thats all 
      return
      end 
  
*-----------------------------------------------------------------------
  
$title INIT_FCODES
      subroutine init_fcodes( fcodes, max_lcodes) 
  
      implicit none 
  
*     Routine to initalize the fcodes to .true. for all the 
*     the lcodes which will be searched.  If the Lcode is 
*     searched and not found in the data base, the lcode will 
*     be reported at the finish of the run.  (Many of the lcodes
*     set up will not be searched because they are only in old
*     data bases.  If these lcodes are not needed they will not be
*     reported as not being found.
*                                  3:40 PM  WED., 27  MAY , 1987
  
      integer*2 
     .    i           ! Loop counter
     .,   max_lcodes  ! Number of lcodes and fcode dimensioned in 
                      ! program.
  
      logical 
     .    fcodes(1)   ! Fcodes indicating that Lcode was not found
  
      do i = 1, max_lcodes
          fcodes(i) = .true.
      end do
  
****  Thats all 
      return
      end 
 
