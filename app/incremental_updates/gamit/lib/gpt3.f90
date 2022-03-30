Subroutine GPT3 ( lun,intdeg,mjd,lat,lon,h_ell,it,p,T,dT,Tm,e,ah,aw,la,undu,Gn_h,Ge_h,Gn_w,Ge_w )
            
! TU-WIEN subroutine from https://vmf.geo.tuwien.ac.at/codes/
! Mods 201216 by RWK; debugged 201218 by MAF:
!   - single station only (vectors become scalars, remove n_stat argument)
!   - change name and unit number for grid file
!   - added argument for unit number to allow opening in calling program
!   - added argument for grid spacing and made dimensions of grid arrays and
!     and limit tests made variable to allow combining gpt3_5 and gpt3_1 routines                               
!   - appended f77 subroutines sassthyd and asknewet, and change comment character to ! for f90
         
! MIT subroutine gpt3.f90 from a combination of TU-WIEN routines gpt3_5.f90
! and gpt3_1.f90, with saasthyd.f and asknewet.f appended.
!
! (c) Department of Geodesy and Geoinformation, Vienna University of
! Technology, 2018
!
!
! This subroutine determines pressure, temperature, temperature lapse rate, 
! mean temperature of the water vapor, water vapour pressure, hydrostatic 
! and wet mapping function coefficients ah and aw, water vapour decrease
! factor, geoid undulation and empirical tropospheric gradients for 
! specific sites near the earth's surface. All output values are valid for
! the specified ellipsoidal height h_ell. 
!
! GPT3_5 and GPT3_1 are based, respectively  on a 5°x5° or 1°x1° external grid
! file (e.g. gpt3_5.grd or gpt2_1.grd) with mean values as well as sine and
! cosine amplitudes for the annual and semiannual variation of the coefficients.
!
! In order to get the respective mapping function, the ah and aw coefficients 
! must be input to vmf3_ht.f90.
! In order to get the respective zenith hydrostatic delay, pressure p must be
! input to saasthyd.f.
! In order to get an approximate value for the respective zenith wet delay,
! water vapor pressure e, mean temperature Tm and water vapor decrease factor
! la must be input to asknewet.f.
!
!
! Reference:
! Landskron, D. and J. Böhm (2017), VMF3/GPT3: Refined Discrete and Empirical 
! Troposphere Mapping Functions, J Geod, doi: 10.1007/s00190-017-1066-2
!
!
! Input parameters:
!             
! lun:   logical unit number for grid file           
! intdeg: 1 or 5 to indicate the grid spacing in degrees
! mjd:   modified Julian date (scalar, only one epoch per call is possible)
! lat:   ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
! lon:   longitude in radians [-pi:pi] or [0:2pi] (vector)
! h_ell: ellipsoidal height in m (vector)
! it:    case 1: no time variation but static quantities
!        case 0: with time variation (annual and semiannual terms)
! 
! Output parameters:
!
! p:    pressure in hPa (vector) 
! T:    temperature in degrees Celsius (vector)
! dT:   temperature lapse rate in degrees per km (vector)
! Tm:   mean temperature weighted with the water vapor in degrees Kelvin 
!       (vector) 
! e:    water vapour pressure in hPa (vector)
! ah:   hydrostatic mapping function coefficient (VMF3) (vector)
! aw:   wet mapping function coefficient (VMF3) (vector)
! la:   water vapour decrease factor (vector)
! undu: geoid undulation in m (vector)
! Gn_h: hydrostatic north gradient in m (vector)
! Ge_h: hydrostatic east gradient in m (vector)
! Gn_w: wet north gradient in m (vector)
! Ge_w: wet east gradient in m (vector)
!
!
! File created by Daniel Landskron, 2018/02/22
!
! ==========================================================

 
! Input arguments
  
! MOD RWK 201211: No station dimension; add unit number and degree spacing
! integer, intent(in) :: n_stat
integer, intent(in) :: lun 
integer, intent(in) :: intdeg
integer, intent(in) :: it
double precision, intent(in) :: mjd
double precision, intent(in) :: lat
double precision, intent(in) :: lon
double precision, intent(in) :: h_ell


! Output arguments
                                    
! MOD RWK 201211: Remove n_stat dimension for all of these
double precision, intent(out) :: p
double precision, intent(out) :: T
double precision, intent(out) :: dT
double precision, intent(out) :: Tm
double precision, intent(out) :: e
double precision, intent(out) :: ah
double precision, intent(out) :: aw
double precision, intent(out) :: la
double precision, intent(out) :: undu
double precision, intent(out) :: Gn_h
double precision, intent(out) :: Ge_h
double precision, intent(out) :: Gn_w
double precision, intent(out) :: Ge_w


! Internal variables

! pi
double precision, parameter :: pi = 4 * atan(1.0d0)

! file open
! MOD RWK 201211: Unit number set in the calling program
!  integer :: file_unit_GPT3_grid = 1
!  integer :: open_status
character(len=1000) :: comment_line
integer :: i_grid, i_ind, n_grid
                         
! Set the grid dimensions; change to 2592 for smaller module with 5-deg grid     
                         
! variables for the gpt3_5.grd
double precision, dimension(64800) :: lat_grid
double precision, dimension(64800) :: lon_grid
double precision, dimension(64800,5) :: p_grid
double precision, dimension(64800,5) :: T_grid
double precision, dimension(64800,5) :: Q_grid
double precision, dimension(64800,5) :: dT_grid
double precision, dimension(64800) :: u_grid
double precision, dimension(64800) :: Hs_grid
double precision, dimension(64800,5) :: ah_grid
double precision, dimension(64800,5) :: aw_grid
double precision, dimension(64800,5) :: la_grid
double precision, dimension(64800,5) :: Tm_grid
double precision, dimension(64800,5) :: Gn_h_grid
double precision, dimension(64800,5) :: Gn_w_grid
double precision, dimension(64800,5) :: Ge_h_grid
double precision, dimension(64800,5) :: Ge_w_grid

! variable for the conversion from mjd to doy
double precision :: sec, jd
integer :: hour, minu, day, month, year, jd_int
double precision :: aa, bb, cc, dd, ee, mm
integer, dimension(12) :: days
integer :: leapYear
double precision :: doy

! constants
double precision :: gm
double precision :: dMtr
double precision :: Rg

! coefficients for amplitudes
double precision :: cosfy
double precision :: coshy
double precision :: sinfy
double precision :: sinhy

! variables for coordinate conversions
double precision :: plon, ppod, difflon, diffpod
integer :: ilon, ipod, ilon1, ipod1
double precision :: dnpod1, dnpod2, dnlon1, dnlon2

! variable to decide whether bilinear or nearest neighbor shall be applied
integer :: bilinear

! variable for the indices in the .grd file
integer, dimension(4) :: indx
integer :: ix

! further variables
double precision :: hgt, redh, c, Hs1
double precision :: Tv, T0, p0, Q, e0
double precision, dimension(4) :: undul, Ql, dTl, pl, Tl, &
                                  lal, ahl, awl, Tml, el, &
                                  Gn_hl, Ge_hl, Gn_wl, Ge_wl
double precision :: R1, R2

!------------------------------------
                                
! RWK MOD 201211: Set the grid length based on the spacing
if ( intdeg == 1 ) then
  n_grid = 64800
else 
  n_grid = 2592 
end if  
  
! MOD RWK 201211: Move this statement to calling program
! open and read the gpt3_5.grd file
! open(unit= file_unit_GPT3_grid, file= 'gpt3_5.grd' , action= 'read', status= 'old', iostat= open_status)
        
! Check if textfile can be opened (iostat == 0 --> no error)
! if (open_status /= 0) then
!     ! report error message
!     print '(a /)', 'Error: gpt3_5.grd can''t be opened...check if it is stored in the same directory as gpt3_5.f90!'
!     ! stop the program
!     stop
! end if
        
! The first line is the name of the grid file and is read in the calling program
! Read the second line, which is a comment giving the column headers
! MOD RWK 201211: change variable name of unit number
read(unit= lun, fmt=*) comment_line
        
! loop over all lines in the file up to the maximum index
! MOD RWK 122011: Make the number of grid points variable
! do i_grid = 1, 2592                                      
do i_grid = 1, n_grid
                 
    ! raise the counter of the index and read the data
    read(unit= lun, fmt=*) lat_grid(i_grid), &
        lon_grid(i_grid), &
        p_grid(i_grid,1), p_grid(i_grid,2), p_grid(i_grid,3),&
        p_grid(i_grid,4), p_grid(i_grid,5), &                  ! pressure in Pascal
        T_grid(i_grid,1), T_grid(i_grid,2), T_grid(i_grid,3),&
        T_grid(i_grid,4), T_grid(i_grid,5), &                  ! temperature in Kelvin
        Q_grid(i_grid,1), Q_grid(i_grid,2), Q_grid(i_grid,3),&
        Q_grid(i_grid,4), Q_grid(i_grid,5), &                  ! specific humidity in kg/kg
        dT_grid(i_grid,1), dT_grid(i_grid,2), &
        dT_grid(i_grid,3), dT_grid(i_grid,4), &
        dT_grid(i_grid,5), &                                   ! temperature lapse rate in Kelvin/m
        u_grid(i_grid),   &                                    ! geoid undulation in m
        Hs_grid(i_grid),   &                                   ! orthometric grid height in m
        ah_grid(i_grid,1), ah_grid(i_grid,2), &
        ah_grid(i_grid,3), ah_grid(i_grid,4), &
        ah_grid(i_grid,5), &                                   ! hydrostatic mapping function coefficient, dimensionless
        aw_grid(i_grid,1), aw_grid(i_grid,2), &
        aw_grid(i_grid,3), aw_grid(i_grid,4), &
        aw_grid(i_grid,5), &                                   ! wet mapping function coefficient, dimensionless
        la_grid(i_grid,1), la_grid(i_grid,2), &
        la_grid(i_grid,3), la_grid(i_grid,4), &
        la_grid(i_grid,5), &                                   ! water vapor decrease factor, dimensionless
        Tm_grid(i_grid,1), Tm_grid(i_grid,2), &
        Tm_grid(i_grid,3), Tm_grid(i_grid,4), &
        Tm_grid(i_grid,5), &                                   ! mean temperature in Kelvin
        Gn_h_grid(i_grid,1), Gn_h_grid(i_grid,2), &
        Gn_h_grid(i_grid,3), Gn_h_grid(i_grid,4), &
        Gn_h_grid(i_grid,5), &                                 ! hydrostatic north gradient in m
        Ge_h_grid(i_grid,1), Ge_h_grid(i_grid,2), &
        Ge_h_grid(i_grid,3), Ge_h_grid(i_grid,4), &
        Ge_h_grid(i_grid,5), &                                 ! hydrostatic east gradient in m
        Gn_w_grid(i_grid,1), Gn_w_grid(i_grid,2), &
        Gn_w_grid(i_grid,3), Gn_w_grid(i_grid,4), &
        Gn_w_grid(i_grid,5), &                                 ! wet north gradient in m
        Ge_w_grid(i_grid,1), Ge_w_grid(i_grid,2), &
        Ge_w_grid(i_grid,3), Ge_w_grid(i_grid,4), &
        Ge_w_grid(i_grid,5)                                    ! wet east gradient in m
                
end do
        
! Don't close the file here
! close(unit= file_unit_GPT3_grid)

! divide by constants in order to get the real values
Q_grid = Q_grid/1000
dT_grid = dt_grid/1000
ah_grid = ah_grid/1000
aw_grid = aw_grid/1000
Gn_h_grid = Gn_h_grid/100000
Ge_h_grid = Ge_h_grid/100000
Gn_w_grid = Gn_w_grid/100000
Ge_w_grid = Ge_w_grid/100000


! convert mjd to doy

hour = floor((mjd-floor(mjd))*24)   ! get hours
minu = floor((((mjd-floor(mjd))*24)-hour)*60)   ! get minutes
sec = (((((mjd-floor(mjd))*24)-hour)*60)-minu)*60   ! get seconds

! change secs, min hour whose sec==60 and days, whose hour==24
if (sec==60) then
    minu = minu +1
end if
if (minu==60) then
    hour = hour +1
end if

! if hr==24, correct jd and set hour==0
jd = mjd+2400000.5d0
if (hour==24) then
    jd = jd + 1
end if

! integer Julian date
jd_int = floor(jd+0.5d0)

aa = jd_int+32044
bb = floor((4*aa+3)/146097)
cc = aa-floor((bb*146097)/4)
dd = floor((4*cc+3)/1461)
ee = cc-floor((1461*dd)/4)
mm = floor((5*ee+2)/153)

day = ee-floor((153*mm+2)/5)+1
month = mm+3-12*floor(mm/10)
year = bb*100+dd-4800+floor(mm/10)

! first check if the specified year is leap year or not (logical output)
if ( (modulo(year,4) == 0   .AND.   modulo(year,100) /= 0)   .OR.   modulo(year,400) == 0 ) then
    leapYear = 1
else
    leapYear = 0
end if

days = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
doy = sum(days(1:month-1)) + day
if (leapYear == 1   .AND.   month > 2) then
    doy = doy + 1
end if
doy = doy + mjd-floor(mjd)   ! add decimal places


! mean gravity in m/s**2
gm = 9.80665d0
! molar mass of dry air in kg/mol
dMtr = 28.965d-3
! universal gas constant in J/K/mol
Rg = 8.3143d0


! factors for amplitudes
if (it==1) then   ! then  constant parameters
    cosfy = 0
    coshy = 0
    sinfy = 0
    sinhy = 0
else 
    cosfy = cos(doy/365.25d0*2*pi)   ! coefficient for A1
    coshy = cos(doy/365.25d0*4*pi)   ! coefficient for B1
    sinfy = sin(doy/365.25d0*2*pi)   ! coefficient for A2
    sinhy = sin(doy/365.25d0*4*pi)   ! coefficient for B2
end if

             
! RWK MOD 201211: Remove the station loop 
! loop over all stations
! do i_stat = 1,n_stat

    ! only positive longitude in degrees
    if (lon < 0) then
        plon = (lon + 2*pi)*180/pi
    else
        plon = lon*180/pi
    end if
    ! transform to polar distance in degrees
    ppod = (-lat + pi/2)*180/pi
    
! RWK MOD 201211: These calculations depend on the grid spacing

    if( intdeg == 5 ) then 

    ! find the index (line in the grid file) of the nearest point
    ipod = floor((ppod+5)/5)
    ilon = floor((plon+5)/5)
    
    ! normalized (to one) differences, can be positive or negative
    diffpod = (ppod - (ipod*5 - 2.5d0))/5
    difflon = (plon - (ilon*5 - 2.5d0))/5
    if (ipod == 37) then
        ipod = 36
    end if
    if (ilon == 73) then
		ilon = 1
    end if
    if (ilon == 0) then
        ilon = 72
    end if

    ! get the number of the corresponding line
    indx(1) = (ipod - 1)*72 + ilon
    
    ! near the poles: nearest neighbour interpolation, otherwise: bilinear
    bilinear = 0
    if (ppod > 2.5d0 .AND. ppod < 177.5d0 ) then
        bilinear = 1          
    end if

!  else assume 1-degree 
    else
                          
    ! find the index (line in the grid file) of the nearest point
    ipod = floor(ppod+1.d0)
    ilon = floor(plon+1.d0)
    
    ! normalized (to one) differences, can be positive or negative
    diffpod = ppod - (ipod - 0.5d0)
    difflon = plon - (ilon - 0.5d0)
    if (ipod == 181) then
        ipod = 180
    end if
    if (ilon == 361) then
		ilon = 1
    end if
    if (ilon == 0) then
        ilon = 360
    end if

    ! get the number of the corresponding line
    indx(1) = (ipod - 1)*360 + ilon
    
    ! near the poles: nearest neighbour interpolation, otherwise: bilinear
    bilinear = 0
    if (ppod > 0.5d0 .AND. ppod < 179.5d0 ) then
        bilinear = 1          
    end if
    
! Endif on degree spacing
    end if 

    ! case of nearest neighbourhood
    if (bilinear == 0) then

        ix = indx(1)
        
        ! transforming ellipsoidial height to orthometric height
        undu = u_grid(ix)
        hgt = h_ell-undu
            
        ! pressure, temperature at the heigtht of the grid
        T0 = T_grid(ix,1) &
             + T_grid(ix,2)*cosfy &
             + T_grid(ix,3)*sinfy &
             + T_grid(ix,4)*coshy &
             + T_grid(ix,5)*sinhy
        p0 = p_grid(ix,1) &
             + p_grid(ix,2)*cosfy &
             + p_grid(ix,3)*sinfy &
             + p_grid(ix,4)*coshy &
             + p_grid(ix,5)*sinhy
         
        ! specific humidity
        Q = Q_grid(ix,1) &
            + Q_grid(ix,2)*cosfy &
            + Q_grid(ix,3)*sinfy &
            + Q_grid(ix,4)*coshy &
            + Q_grid(ix,5)*sinhy
            
        ! lapse rate of the temperature
        dT = dT_grid(ix,1) &
             + dT_grid(ix,2)*cosfy &
             + dT_grid(ix,3)*sinfy &
             + dT_grid(ix,4)*coshy &
             + dT_grid(ix,5)*sinhy 

        ! station height - grid height
        redh = hgt - Hs_grid(ix)

        ! temperature at station height in Celsius
        T = T0 + dT*redh - 273.15d0
        
        ! temperature lapse rate in degrees / km
        dT = dT*1000

        ! virtual temperature in Kelvin
        Tv = T0*(1+0.6077d0*Q)
        
        c = gm*dMtr/(Rg*Tv)
        
        ! pressure in hPa
        p = (p0*exp(-c*redh))/100.d0
            
        ! hydrostatic and wet coefficients ah and aw 
        ah = ah_grid(ix,1) &
             + ah_grid(ix,2)*cosfy &
             + ah_grid(ix,3)*sinfy &
             + ah_grid(ix,4)*coshy &
             + ah_grid(ix,5)*sinhy
        aw = aw_grid(ix,1) &
             + aw_grid(ix,2)*cosfy &
             + aw_grid(ix,3)*sinfy &
             + aw_grid(ix,4)*coshy &
             + aw_grid(ix,5)*sinhy
		
		! water vapour decrease factor la
        la = la_grid(ix,1) &
             + la_grid(ix,2)*cosfy &
             + la_grid(ix,3)*sinfy &
             + la_grid(ix,4)*coshy &
             + la_grid(ix,5)*sinhy
		
		! mean temperature of the water vapor Tm
        Tm = Tm_grid(ix,1) &
             + Tm_grid(ix,2)*cosfy &
             + Tm_grid(ix,3)*sinfy &
             + Tm_grid(ix,4)*coshy &
             + Tm_grid(ix,5)*sinhy
            
        ! north and east gradients (total, hydrostatic and wet)
        Gn_h = Gn_h_grid(ix,1) &
               + Gn_h_grid(ix,2)*cosfy &
               + Gn_h_grid(ix,3)*sinfy &
               + Gn_h_grid(ix,4)*coshy &
               + Gn_h_grid(ix,5)*sinhy
        Ge_h = Ge_h_grid(ix,1) &
               + Ge_h_grid(ix,2)*cosfy &
               + Ge_h_grid(ix,3)*sinfy &
               + Ge_h_grid(ix,4)*coshy &
               + Ge_h_grid(ix,5)*sinhy
        Gn_w = Gn_w_grid(ix,1) &
               + Gn_w_grid(ix,2)*cosfy &
               + Gn_w_grid(ix,3)*sinfy &
               + Gn_w_grid(ix,4)*coshy &
               + Gn_w_grid(ix,5)*sinhy
        Ge_w = Ge_w_grid(ix,1) &
               + Ge_w_grid(ix,2)*cosfy &
               + Ge_w_grid(ix,3)*sinfy &
               + Ge_w_grid(ix,4)*coshy &
               + Ge_w_grid(ix,5)*sinhy
		
		! water vapor pressure in hPa
		e0 = Q*p0/(0.622d0+0.378d0*Q)/100.d0   ! on the grid
        e = e0*(100.d0*p/p0)**(la+1)   ! on the station height - (14) Askne and Nordius, 1987
		
    else   ! bilinear interpolation
        
        ipod1 = ipod + int(sign(1.d0,diffpod))
        ilon1 = ilon + int(sign(1.d0,difflon)) 
! RWK MOD 201211: These depend on the degree spacing
        if( intdeg == 5 ) then 

          if (ilon1 == 73) then
              ilon1 = 1
          end if
          if (ilon1 == 0) then
              ilon1 = 72
          end if                 
          ! get the number of the line
          indx(2) = (ipod1-1)*72 + ilon    ! along same longitude
          indx(3) = (ipod-1)*72 + ilon1   ! along same polar distance
          indx(4) = (ipod1-1)*72 + ilon1   ! diagonal

        else

          if (ilon1 == 361) then
              ilon1 = 1
          end if
          if (ilon1 == 0) then
              ilon1 = 360
          end if
          ! get the number of the line
          indx(2) = (ipod1-1)*360 + ilon    ! along same longitude
          indx(3) = (ipod-1)*360 + ilon1   ! along same polar distance
          indx(4) = (ipod1-1)*360 + ilon1   ! diagonal
        
        end if 

                
        do i_ind = 1,4
              
            ! transforming ellipsoidal height to orthometric height:
            ! Hortho = -N + h_ell
            undul(i_ind) = u_grid(indx(i_ind))
            hgt = h_ell-undul(i_ind)
        
            ! pressure, temperature at the height of the grid
            T0 = T_grid(indx(i_ind),1) &
                 + T_grid(indx(i_ind),2)*cosfy &
                 + T_grid(indx(i_ind),3)*sinfy &
                 + T_grid(indx(i_ind),4)*coshy &
                 + T_grid(indx(i_ind),5)*sinhy
            p0 = p_grid(indx(i_ind),1) &
                 + p_grid(indx(i_ind),2)*cosfy &
                 + p_grid(indx(i_ind),3)*sinfy &
                 + p_grid(indx(i_ind),4)*coshy &
                 + p_grid(indx(i_ind),5)*sinhy

            ! humidity 
            Ql(i_ind) = Q_grid(indx(i_ind),1) &
                        + Q_grid(indx(i_ind),2)*cosfy &
                        + Q_grid(indx(i_ind),3)*sinfy &
                        + Q_grid(indx(i_ind),4)*coshy &
                        + Q_grid(indx(i_ind),5)*sinhy
 
            ! reduction = stationheight - gridheight
            Hs1 = Hs_grid(indx(i_ind))
            redh = hgt - Hs1

            ! lapse rate of the temperature in degree / m
            dTl(i_ind) = dT_grid(indx(i_ind),1) &
                         + dT_grid(indx(i_ind),2)*cosfy &
                         + dT_grid(indx(i_ind),3)*sinfy &
                         + dT_grid(indx(i_ind),4)*coshy &
                         + dT_grid(indx(i_ind),5)*sinhy 

            ! temperature reduction to station height
            Tl(i_ind) = T0 + dTl(i_ind)*redh - 273.15d0

            ! virtual temperature
            Tv = T0*(1+0.6077d0*Ql(i_ind))  
            c = gm*dMtr/(Rg*Tv)
            
            ! pressure in hPa
            pl(i_ind) = (p0*exp(-c*redh))/100.d0
            
            ! hydrostatic coefficient ah
            ahl(i_ind) = ah_grid(indx(i_ind),1) &
                         + ah_grid(indx(i_ind),2)*cosfy &
                         + ah_grid(indx(i_ind),3)*sinfy &
                         + ah_grid(indx(i_ind),4)*coshy &
                         + ah_grid(indx(i_ind),5)*sinhy
            
            ! wet coefficient aw
            awl(i_ind) = aw_grid(indx(i_ind),1) &
                         + aw_grid(indx(i_ind),2)*cosfy &
                         + aw_grid(indx(i_ind),3)*sinfy &
                         + aw_grid(indx(i_ind),4)*coshy &
                         + aw_grid(indx(i_ind),5)*sinhy
					 
            ! water vapor decrease factor la - added by GP
            lal(i_ind) = la_grid(indx(i_ind),1) &
                         + la_grid(indx(i_ind),2)*cosfy &
                         + la_grid(indx(i_ind),3)*sinfy &
                         + la_grid(indx(i_ind),4)*coshy &
                         + la_grid(indx(i_ind),5)*sinhy
					 
            ! mean temperature of the water vapor Tm - added by GP
            Tml(i_ind) = Tm_grid(indx(i_ind),1) &
                         + Tm_grid(indx(i_ind),2)*cosfy &
                         + Tm_grid(indx(i_ind),3)*sinfy &
                         + Tm_grid(indx(i_ind),4)*coshy &
                         + Tm_grid(indx(i_ind),5)*sinhy
					 		 
            ! water vapor pressure in hPa - changed by GP
            e0 = Ql(i_ind)*p0/(0.622d0+0.378d0*Ql(i_ind))/100.d0      ! on the grid
            el(i_ind) = e0*(100.d0*pl(i_ind)/p0)**(lal(i_ind)+1.d0)   ! on the station height  (14) Askne and Nordius, 1987
            
            ! gradients
            Gn_hl(i_ind) = Gn_h_grid(indx(i_ind),1) &
                           + Gn_h_grid(indx(i_ind),2)*cosfy &
                           + Gn_h_grid(indx(i_ind),3)*sinfy &
                           + Gn_h_grid(indx(i_ind),4)*coshy &
                           + Gn_h_grid(indx(i_ind),5)*sinhy
            Ge_hl(i_ind) = Ge_h_grid(indx(i_ind),1) &
                           + Ge_h_grid(indx(i_ind),2)*cosfy &
                           + Ge_h_grid(indx(i_ind),3)*sinfy &
                           + Ge_h_grid(indx(i_ind),4)*coshy &
                           + Ge_h_grid(indx(i_ind),5)*sinhy
            Gn_wl(i_ind) = Gn_w_grid(indx(i_ind),1) &
                           + Gn_w_grid(indx(i_ind),2)*cosfy &
                           + Gn_w_grid(indx(i_ind),3)*sinfy &
                           + Gn_w_grid(indx(i_ind),4)*coshy &
                           + Gn_w_grid(indx(i_ind),5)*sinhy
            Ge_wl(i_ind) = Ge_w_grid(indx(i_ind),1) &
                           + Ge_w_grid(indx(i_ind),2)*cosfy &
                           + Ge_w_grid(indx(i_ind),3)*sinfy &
                           + Ge_w_grid(indx(i_ind),4)*coshy &
                           + Ge_w_grid(indx(i_ind),5)*sinhy
            
          end do
			
            
        dnpod1 = abs(diffpod)   ! distance nearer point
        dnpod2 = 1 - dnpod1     ! distance to distant point
        dnlon1 = abs(difflon)
        dnlon2 = 1 - dnlon1
        
        ! pressure
        R1 = dnpod2*pl(1)+dnpod1*pl(2)
        R2 = dnpod2*pl(3)+dnpod1*pl(4)
        p = dnlon2*R1+dnlon1*R2
            
        ! temperature
        R1 = dnpod2*Tl(1)+dnpod1*Tl(2)
        R2 = dnpod2*Tl(3)+dnpod1*Tl(4)
        T = dnlon2*R1+dnlon1*R2
        
        ! temperature in degree per km
        R1 = dnpod2*dTl(1)+dnpod1*dTl(2)
        R2 = dnpod2*dTl(3)+dnpod1*dTl(4)
        dT = (dnlon2*R1+dnlon1*R2)*1000.d0
            
        ! water vapor pressure in hPa
		R1 = dnpod2*el(1)+dnpod1*el(2)
        R2 = dnpod2*el(3)+dnpod1*el(4)
        e = dnlon2*R1+dnlon1*R2
            
        ! ah and aw
        R1 = dnpod2*ahl(1)+dnpod1*ahl(2)
        R2 = dnpod2*ahl(3)+dnpod1*ahl(4)
        ah = dnlon2*R1+dnlon1*R2
        R1 = dnpod2*awl(1)+dnpod1*awl(2)
        R2 = dnpod2*awl(3)+dnpod1*awl(4)
        aw = dnlon2*R1+dnlon1*R2
        
        ! undulation
        R1 = dnpod2*undul(1)+dnpod1*undul(2)
        R2 = dnpod2*undul(3)+dnpod1*undul(4)
        undu = dnlon2*R1+dnlon1*R2
		
		! water vapour decrease factor
        R1 = dnpod2*lal(1)+dnpod1*lal(2)
        R2 = dnpod2*lal(3)+dnpod1*lal(4)
        la = dnlon2*R1+dnlon1*R2
        
        ! gradients
        R1 = dnpod2*Gn_hl(1)+dnpod1*Gn_hl(2)
        R2 = dnpod2*Gn_hl(3)+dnpod1*Gn_hl(4)
        Gn_h = (dnlon2*R1 + dnlon1*R2)
        R1 = dnpod2*Ge_hl(1)+dnpod1*Ge_hl(2)
        R2 = dnpod2*Ge_hl(3)+dnpod1*Ge_hl(4)
        Ge_h = (dnlon2*R1 + dnlon1*R2)
        R1 = dnpod2*Gn_wl(1)+dnpod1*Gn_wl(2)
        R2 = dnpod2*Gn_wl(3)+dnpod1*Gn_wl(4)
        Gn_w = (dnlon2*R1 + dnlon1*R2)
        R1 = dnpod2*Ge_wl(1)+dnpod1*Ge_wl(2)
        R2 = dnpod2*Ge_wl(3)+dnpod1*Ge_wl(4)
        Ge_w = (dnlon2*R1 + dnlon1*R2)
		
		! mean temperature of the water vapor Tm
        R1 = dnpod2*Tml(1)+dnpod1*Tml(2)
        R2 = dnpod2*Tml(3)+dnpod1*Tml(4)
        Tm = dnlon2*R1+dnlon1*R2

    end if
           
! RWK MOD 201216: Remove i_stat loop
! end do

end subroutine gpt3

! ***************************************************************************

      subroutine saasthyd (p,dlat,hell,zhd)

! This subroutine determines the zenith hydrostatic delay based on the
! equation by Saastamoinen (1972) as refined by Davis et al. (1985)
!
! c Reference:
! Saastamoinen, J., Atmospheric correction for the troposphere and 
! stratosphere in radio ranging of satellites. The use of artificial 
! satellites for geodesy, Geophys. Monogr. Ser. 15, Amer. Geophys. Union, 
! pp. 274-251, 1972.
! Davis, J.L, T.A. Herring, I.I. Shapiro, A.E.E. Rogers, and G. Elgered, 
! Geodesy by Radio Interferometry: Effects of Atmospheric Modeling Errors 
! on Estimates of Baseline Length, Radio Science, Vol. 20, No. 6, 
! pp. 1593-1607, 1985.
!
! input parameters:
!
! p:     pressure in hPa
! dlat:  ellipsoidal latitude in radians 
! dlon:  longitude in radians 
! hell:  ellipsoidal height in m 
! 
! output parameters:
!
! zhd:  zenith hydrostatic delay in meter 
!
! Example 1 :
!
! p = 1000;
! dlat = 48d0*pi/180.d0
! hell = 200.d0
!
! output:
! zhd = 2.2695 m
!
! Johannes Boehm, 8 May 2013
! Johannes Boehm, 24 December 2014, converted to Fortran
! ---

      implicit double precision (a-h,o-z)

      !% calculate denominator f
      f = 1.d0-0.00266d0*dcos(2.d0*dlat) - 0.00000028d0*hell

      !% calculate the zenith hydrostatic delay
      zhd = 0.0022768d0*p/f

      end subroutine saasthyd

! *******************************************************************

      subroutine asknewet (e,Tm,lambda,zwd)

! This subroutine determines the zenith wet delay based on the
! equation 22 by Aske and Nordius (1987) 
!
! c Reference:
! Askne and Nordius, Estimation of tropospheric delay for microwaves from
! surface weather data, Radio Science, Vol 22(3): 379-386, 1987.
!
! input parameters:
!
! e:      water vapor pressure in hPa 
! Tm:     mean temperature in Kelvin
! lambda: water vapor lapse rate (see definition in Askne and Nordius 1987)
! 
! output parameters:
!
! zwd:  zenith wet delay in meter 
!
! Example 1 :
!
! e =  10.9621 hPa
! Tm = 273.8720
! lambda = 2.8071
!
! output:
! zwd = 0.1176 m
!
! Johannes Boehm, 3 August 2013
! Johannes Boehm, 24 December 2014, converted to Fortran
! ---

      implicit double precision (a-h,o-z)
      
      double precision k1,k2,k2p,k3,lambda

      !% coefficients
      k1  = 77.604d0                     !% K/hPa
      k2 = 64.79d0                       !% K/hPa
      k2p = k2 - k1*18.0152d0/28.9644d0  !% K/hPa
      k3  = 377600.d0                    !% KK/hPa

      !% mean gravity in m/s**2
      gm = 9.80665d0
      !% molar mass of dry air in kg/mol
      dMtr = 28.965d-3
      !% universal gas constant in J/K/mol
      R = 8.3143d0

      !% specific gas constant for dry consituents
      Rd = R/dMtr 

      zwd = 1.0d-6*(k2p + k3/Tm)*Rd/(lambda + 1.d0)/gm*e
      
      end subroutine asknewet




  

