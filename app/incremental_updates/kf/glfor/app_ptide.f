CTITLE APP_PTIDE

      subroutine app_ptide ( sol_obs, cstatus, cmpole )

      implicit none

*     Routine to apply the pole tide corrections based on
*     current status (true if appled) and current mean pole (cmpole).
*     Updates applied directly to the solution o-minus-c vector.

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* PASSED IN:
      logical cstatus   ! True if pole tide corrently applied
      character(*) cmpole  ! The Mean pole to which current model 
                       ! is applied.

*   sol_obs(cnum_parn)  - Solution vector.(Updated).

      real*8 sol_obs(cnum_parn) 

* LOCAL VARIABLES
* i,j  -- Loop counters
* dXYZ(3), dNEU(3) -- Pole tide contribution to XYZ and NEU
* type, indx  -- Type and index for parameters in solution
* dt -- Years since 2000.0

      integer*4 i, j,  type, indx, in 
      real*8 dXYZ(3), dNEU(3)

* pmx, pmy -- Pole positions.  Sign and algorithm checked by
*     comparing results with CALC pole-tide calculation.  (program
*     kalptd in kalupd directory of stdhp.)
      real*8 pmy, pmx

* mpxo, mpyo -- Mean pole postions for current orrection
      real*8 mpxo, mpyo
* mpxn, mpyn -- Mean pole postions selected mean pole (if the 
*               If the same no change needed
      real*8 mpxn, mpyn

* dxp, dyp -- Change in pole position needed to implement 
*             requested model change.
      real*8 dxp, dyp

      logical kbit     ! Test bit function

***** Start
*     First see if we need to do anything.  See if mean pole
*     models match and that correction is applied and we
*     are not requesting it be removed.
      if( cmpole(1:6).eq. mean_pole_def(1:6) .and.
     .    cstatus .and. .not. kbit(ptide_opt,3) ) RETURN 

****  Get the pole position we should use. Find values that
*     are non-zero.  (If only switching mean pole these 
*     values will not be needed.
      pmx =  cwob_apr(1)
      pmy =  cwob_apr(2)
      if( pmx.eq.0.d0 ) pmx =  gwob_apr(1)
      if( pmy.eq.0.d0 ) pmy =  gwob_apr(2)
      in = num_mul_pmu/2 + 1
      if( pmx.eq.0.d0 ) pmx =  apr_val_mul_pmu(1,1,in)
      if( pmy.eq.0.d0 ) pmy =  apr_val_mul_pmu(1,2,in)

****  Now get mean/secular pole for old and new models.
      call mean_pole( gepoch_expt,cmpole,        mpxo, mpyo)
      call mean_pole( gepoch_expt,mean_pole_def, mpxn, mpyn)

      dxp = 0 ; dyp = 0

*     Now compute the change in pole we need to implement model
      if( cstatus ) then
        if(  .not. kbit(ptide_opt,3) ) then
*          Currently applied and we are not requeting it be 
*          removced
           dxp = -mpxn + mpxo
           dyp = -mpyn + mpyo
        else   ! Removal of correction requested
           dxp = -( pmx - mpxo)
           dyp = -( pmy - mpyo)
        endif
      else
*       Currently not applied so if requested to be applied
*       apply to new mean pole
        if(  .not. kbit(ptide_opt,3) ) then
           dxp = pmx - mpxn
           dyp = pmy - mpyn
        end if
      endif

****  Now report and apply
      write(*,110) trim(cmpole), trim(mean_pole_def), dxp, dyp
      if( log_unit.ne.6 )
     .write(log_unit,110) trim(cmpole), trim(mean_pole_def), dxp, dyp
 110  format('Updating SE pole tide from ',a,' to ',a,
     .       ' dX/Y-pole ',2(F7.1,1x),' mas')


****  Update the solution vector
      do i = 1, cnum_parn
         call decode_code(gpar_codes(i ), type, indx )
         if( type.eq.7 ) then
              call comp_ptide(apr_val_site(1,1,ltog_sites(indx)),
     .                        dxp, dyp, dXYZ, dNEU)
              sol_obs(i) = sol_obs(i) - dXYZ(1)
              if( kbit(ptide_opt,16) )
     .        write(*,150) ptide_opt, gsite_names(ltog_sites(indx)), 
     .                   (dNEU(j)*1000,j=1,3), (dXYZ(j)*1000,j=1,3)
 150          format('Pole tide OPT ',o6,' Site ',a,' dNEU ',3f7.2,
     .               ' mm, dXYZ ',3F7.2,' mm')
         else if( type.eq.8 ) then
              sol_obs(i) = sol_obs(i) - dXYZ(2)
         else if( type.eq.9 ) then
              sol_obs(i) = sol_obs(i) - dXYZ(3)
         end if
      end do

****  Thats all
      return
      end

CTITLE APP_OPTD

      subroutine app_optd ( sol_obs, cstatus, cmpole )

      implicit none

*     Routine to add or remove the ocean pole tide correction from the
*     solution o-minus-c vector.  The variable opt determine the sign

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/const_param.h'
  
* PASSED IN:
      logical cstatus      ! True if pole tide corrently applied
      character(*) cmpole  ! The Mean pole to which current model 
                       ! is applied

*   sol_obs(cnum_parn)  - Solution vector.  

      real*8 sol_obs(cnum_parn) 

* LOCAL VARIABLES
      integer*4 i,j,k   ! Loop counters
      integer*4 in      ! Index in multi-day PMU is needed. Also used to insert
                        ! lines in ascending sort of reading coefficient file
      integer*4 ierr    ! IOSTAT error
      integer*4 lng_gridnum  ! Number of longitude grid lines (720 nominal)
      integer*4 nlng, nlat, nl  ! Index in long, lat and line number relative
                        ! to first data line (nl=1 is first line).
      integer*4 nc      ! Current line number in coefficent file
      integer*4 nread   ! Number of lines to read to get to next long/lat point
                        ! (maybe zero if sites in same cell)

      integer*4 rdline(max_glb_sites),  ! Line with long/lat for sites sorted in
                        ! ascending order).
     .          sitenm(max_glb_sites),  ! Global site number of corresponding site
     .          parnnm(max_glb_sites)   ! Parameter number of X-coordinate of
                        ! corresponding site

      integer*4 ns        ! Number of stations found in sol_obs decoding)
      integer*4 type, indx   ! Parameter type and index in sol_obs (type 7 is 
                          ! x-coordinate and indx in local site number)
      
      real*8 pmx,   pmy   ! Position (mas)
      real*8 mpxn,  mpyn  ! Mean pole position from IERS 2010 model (mas)
      real*8 mpxo,  mpyo  ! Mean pole position from IERS 1996 model (mas)
* dxp, dyp -- Change in pole position needed to implement 
*             requested model change.
      real*8 dxp, dyp      ! Converted to m1 and m2 at the end.
      real*8 m1, m2       ! dpole postion with -yp (rads)
      real*8 Kcoeff       ! coefficient to get displacement from m1/m2
                          ! due to ocean pole tide (m/rads)
      real*8 gamma_r, gamma_i  ! Real and Imag part of (1+k_2-h_2)

      real*8 slng, dlng, slat, dlat  ! Start and increment in long and latiude
                          ! (degs); lat assumed to start at south pole

      real*8 rot_mat(3,3)  ! Rotation from global to local coordinates
      real*8 loc_rad(3)    ! Co-lat, long (rads) and radius to site (spherical)
      real*8 glng, glat    ! Grid long and lat point (deg)
      real*8 ucoeff(6)     ! coeffs Radial (real/imag), North (real/imag), 
                           ! East (real/imag)

      real*8 dNEU(3), dXYZ(3)  ! Change in NEU and XYZ due to ocean pole tide

      character*256 home_dir   ! User home directory (gg link assumed)
      character*256 optd_file  ! Name of ocean coefficient file

      character*256 line       ! Line read from coefficient file

      logical kbit             ! Test bit function

***** Start
*     First see if we need to do anything.  See if mean pole
*     models match and that correction is applied and we
*     are not requesting it be removed.
      if( cmpole(1:6).eq. mean_pole_def(1:6) .and.
     .    cstatus .and. .not. kbit(ptide_opt,4) ) RETURN 


***** Get the pole position (try different methods to make sure OK)
*     (If only switching mean pole these values will not be needed)
      pmx =  cwob_apr(1)
      pmy =  cwob_apr(2)
      if( pmx.eq.0.d0 ) pmx =  gwob_apr(1)
      if( pmy.eq.0.d0 ) pmy =  gwob_apr(2)
      in = num_mul_pmu/2 + 1
      if( pmx.eq.0.d0 ) pmx =  apr_val_mul_pmu(1,1,in)
      if( pmy.eq.0.d0 ) pmy =  apr_val_mul_pmu(1,2,in)

****  Now get the pole positions we need,
      call mean_pole( gepoch_expt,cmpole,        mpxo, mpyo)
      call mean_pole( gepoch_expt,mean_pole_def, mpxn, mpyn)

      dxp = 0 
      dyp = 0

*     Now compute the change in pole we need to implement model
      if( cstatus ) then
        if(  .not. kbit(ptide_opt,4) ) then
*          Currently applied and we are not requeting it be 
*          removced
           dxp = -mpxn + mpxo
           dyp = -mpyn + mpyo
        else   ! Removal of correction requested
           dxp = -( pmx - mpxo)
           dyp = -( pmy - mpyo)
        endif
      else
*       Currently not applied so if requested to be applied
*       apply to new mean pole
        if(  .not. kbit(ptide_opt,4) ) then
           dxp = pmx - mpxn
           dyp = pmy - mpyn
        end if
      endif

      write(*,110) trim(cmpole), trim(mean_pole_def), dxp, dyp
      if( log_unit.ne.6 )
     .write(log_unit,110) trim(cmpole), trim(mean_pole_def), dxp, dyp
 110  format('Updating Ocean pole tide from ',a,' to ',a,
     .       ' dX/Y-pole ',2(F7.1,1x),' mas')

***** Now compute m1 = dxp and m2 = -dyp (rads
      m1 =  dxp/rad_to_mas
      m2 = -dyp/rad_to_mas

****  OPTD K factor
*     Kcoeff = (4*pi*G*Ae*rho_w*H_p)/(3*g_e)
*     H_p = (8*pi/15)^0.5*(Om^2*Ae^4)/(GM)
      Kcoeff = 5.3394043696d+03 ! meters/radian
*
*     Gamma = (1+k_2-h_2) 
      gamma_r = 0.6870d0   ! Real
      gamma_i = 0.0036d0   ! Imag

***** OK: We need to compute correction.  Work out which stations are
*     used and which line in the ~/gg/tables/opoleloadcoefcmcor.txt file
*     we need to read.  Lines are sorted in ascending order so that file
*     is only read once (per hfile being processed).
      call getenv('HOME',home_dir)
      optd_file = trim(home_dir) // '/gg/tables/opoleloadcoefcmcor.txt'
      open(110,file=optd_file,status='old',iostat=ierr)
      call report_error('IOSTAT',ierr,'open',optd_file,1,'APP_OPTD')
*     Read the header to get the start and step in longitude and latitude
      do i = 1, 5
          read(110,'(a)') line
      end do 
! First_longitude_degrees      =      0.25       -- Line  5
! Last_longitude_degrees       =    359.75
! Longitude_step_degrees       =      0.50       -- Line  7
! Number_latitude_grid_points  =       360
! First_latitude_degrees       =    -89.75       -- Line  9
! Last_latitude_degrees        =     89.75
! Latitude_step_degrees        =      0.50       -- Line 11
      read(line(32:),*,iostat=ierr) slng
      call report_error('IOSTAT',ierr,'read',line,1,
     .                  'APP_OPTD/First Long')
      read(110,'(a)') line
      read(110,'(a)') line
      read(line(32:),*,iostat=ierr) dlng
      call report_error('IOSTAT',ierr,'read',line,1,
     .                  'APP_OPTD/Step Long')
      read(110,'(a)') line
      read(110,'(a)') line
      read(line(32:),*,iostat=ierr) slat
      call report_error('IOSTAT',ierr,'read',line,1,
     .                  'APP_OPTD/First Lat')

      read(110,'(a)') line
      read(110,'(a)') line
      read(line(32:),*,iostat=ierr) dlat
      call report_error('IOSTAT',ierr,'read',line,1,
     .                  'APP_OPTD/Step Lat')
      write(*,120) trim(optd_file),slng, dlng, slat, dlat
 120  format('Ocean Pole Tide ',a,' Long ',2F10.2,' Lat ',2F10.2,' deg')
      write(*,140) pmx, pmy, trim(mean_pole_def), mpxn,  mpyn,
     .                       trim(cmpole), mpxo, mpxn, m1,m2, cstatus
 140  format('PM Values (mas) ',2F10.4,1x,a,' Mean ',2F10.4,
     .       1x,a,' Mean ',2F10.4,' m12 ',2E12.3,' rad, Applied ',L1)

****  Now read to start of data records
      do i = 1, 3 
         read(110,'(a)') line
      end do

***** OK now loop over sites and compute line number for nearest
*     node
      ns = 0
      lng_gridnum = nint(360/dlng)

      do i = 1, cnum_parn
         call decode_code(gpar_codes(i ), type, indx )
         if( type.eq.7 ) then  ! This is an X-ccordinate.  Used this
*                              ! to getr location and grid
            ns = ns + 1
            call XYZ_to_NEU( rot_mat, 
     .          apr_val_site(1,1,ltog_sites(indx)),loc_rad)

*           Get the longitude and latitude grid cell (since the
*           cells are centered with use int and not not nint to
*           the correct cell number.  The first line in the file
*           has nlng and nlat equal 0.
            if( loc_rad(2).lt.0 ) loc_rad(2) = loc_rad(2)+2*pi
            nlng = int((loc_rad(2)*180/pi)/dlng)
            nlat = int((180-loc_rad(1)*180/pi)/dlat)
            nl   = nlat*lng_gridnum + nlng + 1  ! nl = 1 for first line
            
*           OK, now put this line in the correct location
            in = ns   ! Add to end by default
            do j = 1, ns-1
               if( nl.le. rdline(j) ) then
*                  We need to insert the new line at this point
*                  Move all the higher entries up to make space
*                  Even when value is equal, we need to insert into
*                  list
                   in = j    ! Save index where we will insert
                   do k = ns, j+1, -1
                      rdline(k) = rdline(k-1)
                      sitenm(k) = sitenm(k-1)
                      parnnm(k) = parnnm(k-1)
                   end do
                   exit
               endif
            enddo
*           Now insert current values
            rdline(in) = nl
            sitenm(in) = ltog_sites(indx)
            parnnm(in) = i
         end if  
      end do     ! Looping over parameters.

****  Now read the file and for each site compute the ocean pole tide 
*     contribution.  Depending on option this may be added or substracted 
*     from the sol_obs vector
      nc = 0   ! Current line in table
      do i = 1, ns
         nread = rdline(i) - nc
         do j = 1, nread
            read(110,'(a)') line
         enddo
         nc =  rdline(i)
*        The last line read should be the one we want (if line numbers
*        are the same, we decode from the same line string
*        Ucoeffs Are Radial (real/imag), North (real/imag), East (real/imag)
         read(line,*) glng, glat, ucoeff

*****    Now compute dNEU for site
         dNEU(1) = Kcoeff*( (m1*gamma_r+m2*gamma_i)*ucoeff(3) +
     .                      (m2*gamma_r-m1*gamma_i)*ucoeff(4) )
         dNEU(2) = Kcoeff*( (m1*gamma_r+m2*gamma_i)*ucoeff(5) +
     .                      (m2*gamma_r-m1*gamma_i)*ucoeff(6) )
         dNEU(3) = Kcoeff*( (m1*gamma_r+m2*gamma_i)*ucoeff(1) +
     .                      (m2*gamma_r-m1*gamma_i)*ucoeff(2) )

****      Now convert dNEU to dXYZ
          call rotate_geod(dNEU, dXYZ, 'NEU', 'XYZ', 
     .                 apr_val_site(1,1,sitenm(i)),loc_rad,
     .                 rot_mat)

* NOTE: Geocentric latitude is used to compute the grid cell and
*        so differenceces from geodetic latitude at the cells can
*        be up to 0.4 degrees. 
         if( kbit(ptide_opt,16) )
     .   write(*,320) ptide_opt, gsite_names(sitenm(i)),parnnm(i),
     .       glng, glat, ucoeff, dNEU*1000, dXYZ*1000
 320     format('OPT ',o6,1x,a8,1x,I4,1x,2(F10.2),1x,6(F10.6,1x),
     .          ' dNEU ',3F10.2, ' dXYZ ',3F10.2,' mm') 

*****     Now apply correct (sign of m1 and m2 sets whether added
*         of subtrcated).  (Ths is O-C anf the model would be added
*         to C and thus we substract now.
          do j = 1,3
             sol_obs(parnnm(i)+j-1) = sol_obs(parnnm(i)+j-1) -
     .                                dXYZ(j)
          enddo

      end do

****  Thats all
      close(110)
      end
