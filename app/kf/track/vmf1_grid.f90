subroutine vmf1_grid( VMF_grid_file,mjd,lat,lon_old,h_ell,zd , mfh,mfw,zhd,zwd ) 
! MOD AZ 190305: use built Module to store VMF grids
    use vmf_type_build 
! MOD AZ 190305: tailored specifically for track implementation
! Original comments by the author kept for further elaboration
! Mapping functions plus zenith delays from the gridded VMF1 files,
! available from: http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/VMFG/ 
!
! MOD TAH 210408: location comments
! https://vmf.geo.tuwien.ac.at/trop_products/GRID/1x1/VMF3/VMF3_OP//<YYYY> 
!    where <YYYY> is the Year.  File names are of form VMF3_YYYYMMDD.Hhh
! The orography files are at
! https://vmf.geo.tuwien.ac.at/station_coord_files/
! Named https://vmf.geo.tuwien.ac.at/station_coord_files/orography_ell_1x1
! For 1x1 deg file
!
! On the temporal scale, the values from the two surrounding NWM epochs are 
! linearly interpolated to the respective mjd.
! In the horizontal, a bilinear interpolation is done for the mapping 
! function coefficients as well as for the zenith delays. In the vertical, 
! on the one hand the height correction by Niell (1996) is applied in order 
! to "lift" the hydrostatic mapping function from zero height to h_ell; on
! the other hand, specific formulae as suggested by Kouba (2008) are 
! applied in order to "lift" the zenith delays from the respective heights 
! of the grid points (orography_ell) to that of the desired location. 
! The respective (yearly subdivided) gridded VMF1 files and the 
! 'orography_ell' are assumed to be stored in a certain directory 'indir'.
!
! Reference for conversion of mapping functions:
! Niell, A.E. (1996), Global mapping functions for the atmosphere delay at 
! 310 radio wavelengths. J. Geophys. Res., 101, 3227-3246
! 
! Reference for conversion of zenith delays:
! Kouba, J. (2008), Implementation and testing of the gridded Vienna 
! Mapping Function 1 (VMF1). J. Geodesy, Vol. 82:193-205, 
! DOI: 10.1007/s00190-007-0170-0
!
! MOD AZ 190305: include Header for grid path
      include 'vmf_com.h'

! INPUT:
!         o indir_VMF1_grid ... input directory where the yearly subdivided VMF1 gridded files are stored
!         o indir_orography ... input directory where the orography_ell file is stored
!         o VMF1_grid_file: ... type containing filenames, VMF1 data and the orography, which is always passed with the function
!         o mjd ............... modified Julian date
!         o lat ............... ellipsoidal latitude (radians)
!         o lon ............... ellipsoidal longitude (radians)
!         o h_ell ............. ellipsoidal height (m)
!         o zd ................ zenith distance (radians)
!
! OUTPUT:
!         o mfh ............... hydrostatic mapping function, valid at h_ell
!         o mfw ............... wet mapping function, valid at h_ell
!         o zhd ............... zenith hydrostatic delay (m), valid at h_ell
!         o zwd ............... zenith wet delay (m), valid at h_ell
!         o VMF1_grid_file: ... cell containing filenames, VMF1 data and the orography, which is always passed with the function
!
!
! Make sure to set the following variables as actual arguments in the script from which you call vmf1_grid.f90:
!
! character(len=:), allocatable :: indir_VMF1_grid
! character(len=:), allocatable :: indir_orography
! double precision :: mjd, lat, lon, h_ell, zd
! type :: VMF1_grid_file_type
!     character(len=17), dimension(:), allocatable :: filename
!     double precision, dimension(2,13104,6) :: vmf_data_all
!     integer, dimension(13104) :: orography_ell
!     double precision :: lat, lon
! end type VMF1_grid_file_type
! type(VMF1_grid_file_type) :: VMF1_grid_file
! double precision :: mfh, mfw, zhd, zwd
!
! -------------------------------------------------------------------------
!
! written by Daniel Landskron (2018/01/30)
!
! =========================================================================



! Input arguments

!character(len=:), allocatable, intent(in) :: indir_VMF1_grid
!character(len=:), allocatable, intent(in) :: indir_orography

!type :: vmf_grid
!  character(len=17), dimension(:), allocatable :: filename
!  double precision, dimension(2,13104,6) :: vmf_data_all
!  integer, dimension(13104) :: orography_ell
!  double precision :: lat, lon
!end type vmf_grid

!type(vmf_grid), intent(inout) :: VMF1_grid_file
double precision, intent(in) :: mjd
double precision, intent(in) :: lat
double precision, intent(in) :: lon_old
double precision, intent(in) :: h_ell
double precision, intent(in) :: zd


! Output arguments

double precision, intent(out) :: mfh
double precision, intent(out) :: mfw
double precision, intent(out) :: zhd
double precision, intent(out) :: zwd


! Internal variables

! pi
double precision, parameter :: pi = 4 * atan(1.0d0)

! further variables for lat/lon
double precision :: lon, lat_deg, lon_deg
double precision, dimension(2) :: lat_int_deg, lon_int_deg, ind_lat_int_deg, ind_lon_int_deg
integer :: num_lon
double precision, dimension(4) :: lat_red

! the variables for the 1 or 3 mjd's
double precision, dimension(2) :: mjd_int
double precision, dimension(:), allocatable :: mjd_all, jd_all
integer, dimension(:), allocatable :: jd_all_int

! further time variables
double precision, dimension(:), allocatable :: hour, minu, sec, day, month, year, epoch
character(len=4) :: year_str
character(len=2) :: month_str, day_str, epoch_str
integer :: year_temp, month_temp, day_temp, epoch_temp

! indexing variables
integer :: i_mjd

! variables for the conversion of date to mjd
double precision, dimension(:), allocatable :: aa, bb, cc, dd, ee, mm

! filename
character(len=17), dimension(:), allocatable :: filename
    
! variable for storing the temporary orography data
integer, dimension(1365,10) :: orography_ell_temp

! variable for orography_ell and VMF1_data
integer :: load_new
integer :: file_unit_orography_ell = 1
integer :: file_unit_VMF1_data = 2
integer :: open_status, read_status
integer :: i, ind_data, i_line, i_col, i_file
integer, dimension(4) :: index_p
character :: first_char
double precision, dimension(2,4,6) :: VMF1_data
integer :: num_comment_lines

! VMF1 data after temporal interpolation
double precision, dimension(4,7) :: VMF1_data_int_h0
double precision, dimension(4,9) :: VMF1_data_int_h1

! mapping function variables
double precision :: phh, c11h, c10h
double precision, dimension(4) :: ah, aw, bh, bw, ch, cw, c0h
double precision :: a_ht, b_ht, c_ht
double precision :: ht_corr, ht_corr_coef
double precision :: zhd_lon1, zhd_lon2, zwd_lon1, zwd_lon2, mfh_lon1, mfh_lon2, mfw_lon1, mfw_lon2

double precision :: inter_mjd
! day of year
double precision :: doy

! elevation 
double precision :: el

!======================================================================


! save lat and lon also in degrees
lat_deg = lat*180/pi
lon_deg = lon_old*180/pi   ! one must be named "_old" because as a dummy argument it cannot be changed within the subroutine

! only positive longitude in degrees
if (lon_deg < 0) then
    lon = lon_old + 2*pi
    lon_deg = lon_deg + 360
end if



!----------------------------------------------------------------------
! (1) convert the mjd to year, month, day in order to find the correct files
!----------------------------------------------------------------------

! MOD AZ 190305: In some Linux system where mem control is strict,
!                it's neccessary to pre-allocate all these following
!                variables. Any "allocate" action in this file is probably
!                added by me, which is different from the original
!                source file.
! find the two surrounding epochs
if (modulo(mjd,0.25d0)==0) then
    allocate(mjd_all(1))
    allocate(jd_all(1))
    allocate(jd_all_int(1))
    allocate(hour(1))
    allocate(minu(1))
    allocate(sec(1))
    allocate(day(1))
    allocate(month(1))
    allocate(year(1))
    allocate(epoch(1))
    allocate(aa(1))
    allocate(bb(1))
    allocate(cc(1))
    allocate(dd(1))
    allocate(ee(1))
    allocate(mm(1))
    mjd_all = mjd
else
    allocate(mjd_all(3))
    allocate(jd_all(3))
    allocate(jd_all_int(3))
    allocate(hour(3))
    allocate(minu(3))
    allocate(sec(3))
    allocate(day(3))
    allocate(month(3))
    allocate(year(3))
    allocate(epoch(3))
    allocate(aa(3))
    allocate(bb(3))
    allocate(cc(3))
    allocate(dd(3))
    allocate(ee(3))
    allocate(mm(3))
    mjd_int = (/ floor(mjd*4.0d0)/4.0d0, ceiling(mjd*4.0d0)/4.0d0 /)
    mjd_all = (/ mjd, mjd_int /)
end if


hour = floor((mjd_all-floor(mjd_all))*24)   ! get hours
minu = floor((((mjd_all-floor(mjd_all))*24)-hour)*60)   ! get minutes
sec = (((((mjd_all-floor(mjd_all))*24)-hour)*60)-minu)*60   ! get seconds

! change secs, min hour whose sec==60 and days, whose hour==24
do i_mjd = 1,size(hour)
    if (sec(i_mjd)==60) then
        minu(i_mjd) = minu(i_mjd) +1
    end if
    if (minu(i_mjd)==60) then
        hour(i_mjd) = hour(i_mjd) +1
    end if
    if (hour(i_mjd)==24) then
        mjd_all(i_mjd) = mjd_all(i_mjd) +1
    end if
end do

! calc jd (yet wrong for hour==24)
jd_all = mjd_all+2400000.5d0

! integer Julian date
jd_all_int = floor(jd_all+0.5d0)

aa = jd_all_int+32044
bb = floor((4*aa+3)/146097)
cc = aa-floor((bb*146097)/4)
dd = floor((4*cc+3)/1461)
ee = cc-floor((1461*dd)/4)
mm = floor((5*ee+2)/153)

day = ee-floor((153*mm+2)/5)+1
month = mm+3-12*floor(mm/10)
year = bb*100+dd-4800+floor(mm/10)

epoch = (mjd_all-floor(mjd_all))*24

! derive related VMFG filename(s)
if (size(mjd_all)==1) then   ! if the observation epoch coincides with an NWM epoch
    year_temp = year(1)
    month_temp = month(1)
    day_temp = day(1)
    epoch_temp = epoch(1)
    write(year_str,'(i4.4)') year_temp
    write(month_str, '(i2.2)') month_temp
    write(day_str, '(i2.2)') day_temp
    write(epoch_str, '(i2.2)') epoch_temp
    allocate(filename(1))
    filename(1) = 'VMFG_' // year_str // month_str // day_str // '.H' // epoch_str
else
    allocate(filename(2))
    do i_mjd = 2,size(mjd_all)
        year_temp = year(i_mjd)
        month_temp = month(i_mjd)
        day_temp = day(i_mjd)
        epoch_temp = epoch(i_mjd)
        write(year_str,'(i4.4)') year_temp
        write(month_str, '(i2.2)') month_temp
        write(day_str, '(i2.2)') day_temp
        write(epoch_str, '(i2.2)') epoch_temp
        filename(i_mjd-1) = 'VMFG_' // year_str // month_str // day_str // '.H' // epoch_str
    end do
end if



!----------------------------------------------------------------------
! (2) check if new files have to be loaded or if the overtaken ones are sufficient
!----------------------------------------------------------------------


if (.not. allocated(VMF1_grid_file % filename)) then   ! in the first run, 'VMF1_file' is always empty and the orography_ell file has to loaded
    
    load_new = 1
    allocate(VMF1_grid_file % filename(2))
    VMF1_grid_file % filename = filename   ! replace the empty cell by the current filenames
    
    ! load the orography_ell file
    ! Open the file
    open(unit= file_unit_orography_ell,file= indir_orography // &
        '/orography_ell', action= 'read', &
         status= 'old', iostat= open_status)
    ! Check if textfile can be opened (iostat == 0 --> no error)
    if (open_status /= 0) then
        ! report error message
        print '(a /)', 'Error: Problem with opening the orography_ell file! Program stopped!'
        ! stop the program
        stop
    end if
            
    ! first line is comment
    read(unit= file_unit_orography_ell, fmt= '(a)', iostat= read_status) first_char
        
    ! loop over all lines in the file
    do i= 1, 1365
        
        if (modulo(i,15)==0) then
            read(unit= file_unit_orography_ell, fmt=*) orography_ell_temp(i,1), &
                                                       orography_ell_temp(i,2), &
                                                       orography_ell_temp(i,3), &
                                                       orography_ell_temp(i,4), &
                                                       orography_ell_temp(i,5)
            orography_ell_temp(i,5) = 99999   ! overwrite this number as it is a repitition
            orography_ell_temp(i,6) = 99999   ! fill the remaining indicex with dummy numbers so that one later knows, which indices were empty
            orography_ell_temp(i,7) = 99999
            orography_ell_temp(i,8) = 99999
            orography_ell_temp(i,9) = 99999
            orography_ell_temp(i,10) = 99999
        else
            read(unit= file_unit_orography_ell, fmt=*) orography_ell_temp(i,1), &
                                                       orography_ell_temp(i,2), &
                                                       orography_ell_temp(i,3), &
                                                       orography_ell_temp(i,4), &
                                                       orography_ell_temp(i,5), &
                                                       orography_ell_temp(i,6), &
                                                       orography_ell_temp(i,7), &
                                                       orography_ell_temp(i,8), &
                                                       orography_ell_temp(i,9), &
                                                       orography_ell_temp(i,10) 
        end if
                
    end do
    
    ! convert to the real column format
    ind_data = 1
    do i_line = 1,1365
        do i_col = 1,10
            if (orography_ell_temp(i_line,i_col) /= 99999) then
                VMF1_grid_file % orography_ell(ind_data) = orography_ell_temp(i_line,i_col)
                ind_data = ind_data+1
            end if
        end do
    end do
    
    ! close the file
    close(unit= file_unit_orography_ell)
    
    VMF1_grid_file % lat = lat
    VMF1_grid_file % lon = lon
    
elseif ( size(filename)==2 ) then

    if  (  VMF1_grid_file % filename(1) == filename(1)   .AND.   &
          VMF1_grid_file % filename(2) == filename(2)   .AND.   &
          ( lat > VMF1_grid_file % lat  .OR.  &
          (lat == VMF1_grid_file % lat .AND. &
           lon <= VMF1_grid_file % lon) )  ) then

! if the current filenames are the same as in the forwarded files
    
    load_new = 0
    else 
    
    load_new = 1
    VMF1_grid_file % filename = filename
    VMF1_grid_file % lat = lat
    VMF1_grid_file % lon = lon

    end if
    
elseif ( size(filename)==1 ) then

    if  (  VMF1_grid_file % filename(1) == filename(1)   .AND.   &
          ( lat > VMF1_grid_file % lat  .OR.  &
          (lat == VMF1_grid_file % lat .AND. &
           lon <= VMF1_grid_file % lon) )  ) then

! if the current filenames are the same as in the forwarded files
    
    load_new = 0
    
    else 
    
    load_new = 1
    VMF1_grid_file % filename(1) = filename(1)
    VMF1_grid_file % lat = lat
    VMF1_grid_file % lon = lon

    end if
    
else   ! if new VMF1 files are required, then they must be loaded anew
    
    load_new = 1
    VMF1_grid_file % filename = filename
    VMF1_grid_file % lat = lat
    VMF1_grid_file % lon = lon
    
end if



!----------------------------------------------------------------------
! (3) find the indices of the 4 surrounding grid points
!----------------------------------------------------------------------


! find the coordinates (lat,lon) of the surrounding grid points
lat_int_deg = (/ floor(lat_deg/2)*2, ceiling(lat_deg/2)*2 /)
lon_int_deg = (/ floor(lon_deg/2.5d0)*2.5d0, ceiling(lon_deg/2.5d0)*2.5d0 /)


! find the indices of these surrounding grid points

num_lon = 144
ind_lat_int_deg = -(lat_int_deg-90)/2   ! auxiliary variable
ind_lon_int_deg = lon_int_deg/2.5d0+1   ! auxiliary variable

index_p(1) = ind_lat_int_deg(1)*num_lon+ind_lon_int_deg(1)
index_p(2) = ind_lat_int_deg(1)*num_lon+ind_lon_int_deg(2)
index_p(3) = ind_lat_int_deg(2)*num_lon+ind_lon_int_deg(1)
index_p(4) = ind_lat_int_deg(2)*num_lon+ind_lon_int_deg(2)



!----------------------------------------------------------------------
! (4) read the correct data and perform a linear time interpolation from the surrounding two epochs
!----------------------------------------------------------------------


if (load_new == 1) then
    
    do i_file = 1,size(filename)
        
        ! read the files and collect the data
        if (size(mjd_all)==1) then   ! if the observation epoch coincides with an NWM epoch
            year_temp = year(1)
            write(year_str,'(i4.4)') year_temp
        else
            year_temp = year(i_file+1)
            write(year_str,'(i4.4)') year_temp
        end if
            open(unit= file_unit_VMF1_data, file= indir_VMF1_grid // '/'  &
             // year_str // '/' // filename(i_file) , action= 'read', &
             status= 'old', iostat= open_status)
        
        ! Check if textfile can be opened (iostat == 0 --> no error)
        if (open_status /= 0) then
            ! report error message
            print '(a /)', 'Error: Problem with opening the gridded VMF1 data! Program stopped!'
            ! stop the program
            stop
        end if
        
        ind_data = 0 
        num_comment_lines = 0   ! there are comment lines in the VMF1 file
        
        ! determine the number of comment lines
        do
            ! read next line
            read(unit= file_unit_VMF1_data, fmt= '(a)', iostat= read_status) first_char
            
            ! test if end of file is reached and exit loop in this case
            if (is_iostat_end(read_status)) exit
            
            ! test if the first character of the line is a comment sign ("%" or "!")
            if (first_char == '%' .or. first_char == '!') then
                num_comment_lines = num_comment_lines + 1
            else
                exit
            end if
            
        end do
        
        ! rewind the file
        rewind(unit = file_unit_VMF1_data)
        
        ! loop over all lines in the file up to the maximum index
        do i= 1, maxval(index_p)+num_comment_lines 
            
            ! skip all comment lines
            if (i<=num_comment_lines) then
                ! read next line
                read(unit= file_unit_VMF1_data, fmt= '(a)', iostat= read_status) first_char
                cycle
            end if
                    
            ! raise the counter of the index and read the data
            ind_data= ind_data + 1
            read(unit= file_unit_VMF1_data, fmt=*) VMF1_grid_file % vmf_data_all(i_file,ind_data,:)
                
        end do
        
        ! close the file
        close(unit= file_unit_VMF1_data)
    
        ! reduce data to the respective indices
        VMF1_data(i_file,:,:) = VMF1_grid_file % vmf_data_all(i_file,index_p,:)
        
    end do

else
   
    ! reduce data to the respective indices
    do i_file = 1,size(filename)
        VMF1_data(i_file,:,:) = VMF1_grid_file % vmf_data_all(i_file,index_p,:)
    end do
    
end if


! perform the linear temporal interpolation
if (size(mjd_all)==1) then   ! if the observation epoch coincides with an NWM epoch
    VMF1_data_int_h0(:,1:6) = VMF1_data(1,:,1:6)
else   ! else perform the linear interpolation
    VMF1_data_int_h0(:,1:6) = VMF1_data(1,:,:) + (VMF1_data(2,:,:)-VMF1_data(1,:,:))*(mjd-mjd_int(1))/(mjd_int(2)-mjd_int(1))   ! the appendix 'h0' means that the values are valid at zero height
end if

! the first four columns are equal
VMF1_data_int_h1(:,1:4) = VMF1_data_int_h0(:,1:4)



!----------------------------------------------------------------------
! (5) bring mfh, mfw, zhd and zwd of the surrounding grid points to the respective height of the location
!----------------------------------------------------------------------


! (a) zhd
! to be exact, the latitudes of the respective grid points would have to be used instead of the latitude of the station (lat). However, the loss of accuracy is only in the sub-micrometer range.
VMF1_data_int_h0(:,7) = (VMF1_data_int_h0(:,5)/0.0022768d0) * &
                        (1-0.00266d0*cos(2*lat)-0.00000028d0*VMF1_grid_file % orography_ell(index_p))
! (1) convert the hydrostatic zenith delay at grid height to the respective pressure value
VMF1_data_int_h1(:,7) = VMF1_data_int_h0(:,7)*(1-0.0000226d0* &
                        (h_ell-VMF1_grid_file % orography_ell(index_p)))**5.225d0
! (2) lift the pressure each from grid height to site height
VMF1_data_int_h1(:,5) = 0.0022768d0*VMF1_data_int_h1(:,7) / (1-0.00266d0*cos(2*lat)-0.00000028d0*h_ell)
! (3) convert the lifted pressure to zhd again (as proposed by Kouba, 2008)

! (b) zwd
! simple exponential decay approximation function
VMF1_data_int_h1(1:4,6) = VMF1_data_int_h0(1:4,6) * exp(-(h_ell-VMF1_grid_file % orography_ell(index_p))/2000)

! convert zenith distance to elevation
el = pi/2 - zd

! (c) ah => mfh, aw => mfw

! reference day is 28 January (this is taken from Niell (1996) to be consistent)
doy = mjd  - 44239.d0 + 1 - 28
      
bh = 0.0029d0
c0h = 0.062d0
if (lat .lt. 0.d0) then   ! southern hemisphere
    phh  = pi
    c11h = 0.007d0
    c10h = 0.002d0
else                     ! northern hemisphere
    phh  = 0.d0
    c11h = 0.005d0
    c10h = 0.001d0
end if

! read the latitudes of the neighboring grid points
lat_red = VMF1_data_int_h0(1:4,1)*pi/180

ah = VMF1_data_int_h0(:,3)
aw = VMF1_data_int_h0(:,4)
bh = spread(0.0029d0,1,size(ah))
c0h = spread(0.062d0,1,size(ah))
ch = c0h + ((cos(doy/365.25d0*2*pi + phh)+1)*c11h/2 + c10h)*(1-cos(lat_red))
bw = spread(0.00146d0,1,size(ah))
cw = spread(0.04391d0,1,size(ah))

! calculating the hydrostatic and wet mapping function
VMF1_data_int_h1(:,8) = (1+(ah/(1+bh/(1+ch))))   /   (sin(el)+(ah/(sin(el)+bh/(sin(el)+ch))))
VMF1_data_int_h1(:,9) = (1+(aw/(1+bw/(1+cw))))   /   (sin(el)+(aw/(sin(el)+bw/(sin(el)+cw))))

! height correction for the hydrostatic part [Niell, 1996]
a_ht = 2.53d-5
b_ht = 5.49d-3
c_ht = 1.14d-3
ht_corr_coef = 1/sin(el)   -   (1+(a_ht/(1+b_ht/(1+c_ht))))  /  (sin(el)+(a_ht/(sin(el)+b_ht/(sin(el)+c_ht))))
ht_corr      = ht_corr_coef * h_ell/1000
VMF1_data_int_h1(:,8) = VMF1_data_int_h1(:,8) + ht_corr



!----------------------------------------------------------------------
! (6) perform the bilinear interpolation
!----------------------------------------------------------------------

if ( index_p(1)==index_p(2) .AND. index_p(2)==index_p(3) .AND. index_p(3)==index_p(4) ) then   ! if the point is directly on a grid point
    
    zhd = VMF1_data_int_h1(1,5)
    zwd = VMF1_data_int_h1(1,6)
    mfh = VMF1_data_int_h1(1,8)
    mfw = VMF1_data_int_h1(1,9)
    
else

    ! bilinear interpolation (interpreted as two 1D linear interpolations for lat and lon, but programmed without subfunctions)
    
    ! (a) linear interpolation for longitude
    if (VMF1_data_int_h1(1,2) /= VMF1_data_int_h1(2,2)) then   ! if longitude must be interpolated (that is, the point does not have a longitude on the interval [0:2.5:357.5])
        zhd_lon1 = VMF1_data_int_h1(1,5) + (VMF1_data_int_h1(2,5)-VMF1_data_int_h1(1,5))* &
                  (lon_deg-VMF1_data_int_h1(1,2))/(VMF1_data_int_h1(2,2)-VMF1_data_int_h1(1,2))
        zhd_lon2 = VMF1_data_int_h1(3,5) + (VMF1_data_int_h1(4,5)-VMF1_data_int_h1(3,5))* &
                  (lon_deg-VMF1_data_int_h1(3,2))/(VMF1_data_int_h1(4,2)-VMF1_data_int_h1(3,2))
        zwd_lon1 = VMF1_data_int_h1(1,6) + (VMF1_data_int_h1(2,6)-VMF1_data_int_h1(1,6))* &
                  (lon_deg-VMF1_data_int_h1(1,2))/(VMF1_data_int_h1(2,2)-VMF1_data_int_h1(1,2))
        zwd_lon2 = VMF1_data_int_h1(3,6) + (VMF1_data_int_h1(4,6)-VMF1_data_int_h1(3,6))* &
                  (lon_deg-VMF1_data_int_h1(3,2))/(VMF1_data_int_h1(4,2)-VMF1_data_int_h1(3,2))
        mfh_lon1 = VMF1_data_int_h1(1,8) + (VMF1_data_int_h1(2,8)-VMF1_data_int_h1(1,8))* &
                  (lon_deg-VMF1_data_int_h1(1,2))/(VMF1_data_int_h1(2,2)-VMF1_data_int_h1(1,2))
        mfh_lon2 = VMF1_data_int_h1(3,8) + (VMF1_data_int_h1(4,8)-VMF1_data_int_h1(3,8))* &
                  (lon_deg-VMF1_data_int_h1(3,2))/(VMF1_data_int_h1(4,2)-VMF1_data_int_h1(3,2))
        mfw_lon1 = VMF1_data_int_h1(1,9) + (VMF1_data_int_h1(2,9)-VMF1_data_int_h1(1,9))* &
                  (lon_deg-VMF1_data_int_h1(1,2))/(VMF1_data_int_h1(2,2)-VMF1_data_int_h1(1,2))
        mfw_lon2 = VMF1_data_int_h1(3,9) + (VMF1_data_int_h1(4,9)-VMF1_data_int_h1(3,9))* &
                  (lon_deg-VMF1_data_int_h1(3,2))/(VMF1_data_int_h1(4,2)-VMF1_data_int_h1(3,2))
    else   ! if the station coincides with the longitude of the grid
        zhd_lon1 = VMF1_data_int_h1(1,5)
        zhd_lon2 = VMF1_data_int_h1(3,5)
        zwd_lon1 = VMF1_data_int_h1(1,6)
        zwd_lon2 = VMF1_data_int_h1(3,6)
        mfh_lon1 = VMF1_data_int_h1(1,8)
        mfh_lon2 = VMF1_data_int_h1(3,8)
        mfw_lon1 = VMF1_data_int_h1(1,9)
        mfw_lon2 = VMF1_data_int_h1(3,9)
    end if
    
    ! linear interpolation for latitude
    if (VMF1_data_int_h1(1,1) /= VMF1_data_int_h1(3,1)) then
        zhd = zhd_lon1 + (zhd_lon2-zhd_lon1)*(lat_deg-VMF1_data_int_h1(1,1))/ &
              (VMF1_data_int_h1(3,1)-VMF1_data_int_h1(1,1))
        zwd = zwd_lon1 + (zwd_lon2-zwd_lon1)*(lat_deg-VMF1_data_int_h1(1,1))/ &
              (VMF1_data_int_h1(3,1)-VMF1_data_int_h1(1,1))
        mfh = mfh_lon1 + (mfh_lon2-mfh_lon1)*(lat_deg-VMF1_data_int_h1(1,1))/ &
              (VMF1_data_int_h1(3,1)-VMF1_data_int_h1(1,1))
        mfw = mfw_lon1 + (mfw_lon2-mfw_lon1)*(lat_deg-VMF1_data_int_h1(1,1))/ &
              (VMF1_data_int_h1(3,1)-VMF1_data_int_h1(1,1))
    else   ! if the station coincides with the latitude of the grid
        zhd = zhd_lon1
        zwd = zwd_lon1
        mfh = mfh_lon1
        mfw = mfw_lon1
    end if
      
end if
! MOD AZ 190305: Remember to deallocate to avoid mem leak
   deallocate(mjd_all)
   deallocate(jd_all) 
   deallocate(jd_all_int)
   deallocate(filename)
   deallocate(hour)
   deallocate(minu)
   deallocate(sec)
   deallocate(day)
   deallocate(month)
   deallocate(year)
   deallocate(epoch)    
   deallocate(aa)
   deallocate(bb)
   deallocate(cc)
   deallocate(dd)
   deallocate(ee)
   deallocate(mm)    

end subroutine
