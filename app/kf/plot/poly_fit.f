CTITLE POLY_FIT
 
      subroutine poly_fit( x_data, y_data, pt_data )

      implicit none 
 
 
*     Routine to fit a polynomial to the data, and compute the RMS
*     scatter about the polynomial.
 
      include 'plot_param.h'
      include 'plot_com.h'
 
 
*   dates(5)        - Date for calender reference
*   ierr            - IOSTAT err
*   iuse(0:max_poly)    - Pivot index used by SMNV8
*   i,j,k           - Loop counters
*   jel             - Function to return index in lower triangular
*                   - matrix
*   len_poly        - Length of label for poly_labels
*   pt_data(2,1)    - Point information (Second entry UNW flag is
*                   - used if it is given in P_FIELD)
*   trimlen         - HP function for length of string
*   fit_point       - Type of point to used in fitting polynomial
*                    (0 will use all points)
 
 
      integer*4 dates(5), ierr, iuse(0:max_poly), i,j,k, jel, len_poly,
     .    pt_data(2,1), trimlen, indx, fit_point
 
*   use             - Indicates we should used data point in solution
 
      logical use
 
*   x_data(2,1)     - The x data information (only the value is used)
*   y_data(2,1)     - The y data information (value and sigma are used
*                   - if sigma has been given in P_field)
*   chisq           - Chi**2/f for the fit of the data
*   minform         - Minimum value for format
*   maxform         - Maximum value for format
*   result          - Result for polynomimial coefficient (scaled by
*                   - user factor)
*   sigma           - Sigma for cooefficient, scaled by NRMS and
*                   - user factor)
*   wrms            - WRMS scatter
*   nrms            - NRMS scatter (used to scale sigmas)
 
      real*4 x_data(2,1), y_data(2,1), chisq, minform, maxform, result,
     .    sigma, wrms, nrms

      real*4 resy_min, resy_max  ! Min and Max residuals
 
*   deriv(0:max_poly)   - Partial for each data point
*   estimate        - Polynomial evaluated at each data point
*   pvrow(0:max_poly), pvcol(0:max_poly), pvrwb - Pivot data used
*                   - by SMNV8.
*   scaling         - Scaling factor for conversion of polynomial
*                   - coefficients
*   value           - Temporary real*8 value
*   wgh             - Weight to be used for each data
 
      real*8 deriv(0:max_poly), estimate,
     .    pvrow(0:max_poly), pvcol(0:max_poly), pvrwb, scaling, value,
     .    wgh, sectag
 
*   labform         - Format for output values
*   detrend         - D if data are be detrended

      character*10 labform, detrend
      character*1  cd

* MOD TAH 131021: Save refval so that total is printed for offset
      real*8 refsave_y   ! Save reference value of y to output totals
 
 
 
***** Firstly get the polynomial order
      indx = 9
      call read_line(buffer,indx,'I4',ierr, poly_order, cd)
*                                 ! Use an offset and slope
      if( ierr.ne.0 ) then
          poly_order = 1
      end if
 
      if( poly_order.gt. max_poly ) then
          call report(' Polynomial order too large')
          poly_order = max_poly
      end if

*     See if point type for fit is passed
      call read_line(buffer,indx,'I4',ierr, fit_point,cd)
      if( ierr.ne.0 ) fit_point = 0
      
*     See if residuals are to be detrended
      call read_line(buffer,indx,'CH',ierr, fit_point,detrend)
      if( ierr.ne.0 ) then
           detrend = 'N'
      else
           call casefold(detrend)
      end if

***** Now clear the poly solution and normal equations
 
      do i = 0, poly_order
          poly_vec(i) = 0.d0
      end do
 
      do i = 1, (poly_order+1)*(poly_order+2)/2
          poly_mat(i) = 0.d0
      end do
 
*                         ! Statistics sums
      sum_poly(1) = 0.d0
      sum_poly(2) = 0.d0
 
      num_poly_data = 0
 
***** Now loop over the data
 
      do i = 1, num_data
 
*         See if in window first
 
          if( x_data(1,i).ge. poly_window(1) .and.
     .        x_data(1,i).le. poly_window(2) .and.
     .        y_data(1,i).ge. poly_window(3) .and.
*                                                         ! Inside window
     .        y_data(1,i).le. poly_window(4)       ) then
 
*             See if we should use
              USE = .true.
*                                          ! UNW field given
              if( p_field(2).ne.0 ) then
                  if( pt_data(2,i).ne.0 ) then
                      use = .false.
                  end if
              end if

*             See if fails point type
              if( p_field(1).ne.0 .and.
     .            pt_data(1,i).ne.fit_point .and.
     .            fit_point.ne.0                 ) then
                     use = .false.
              end if
 
*             Now get the weight
              wgh = 1.d0
*                                         ! Errbar given, so use
              if( y_field(3).ne.0 ) then
                  wgh = 1.d0/y_data(2,i)**2
              end if
 
*             If we are using this point then incerment normal equations
              if( use ) then
                  num_poly_data = num_poly_data + 1
c                 write(*,*) num_poly_data, x_data(1,i), y_data(1,i),
c    .                       y_data(2,i)
                  do j = 0, poly_order
                      if( x_data(1,i).eq.0 .and. j.eq.0 ) then
*                                         ! Trap 0**0 error
                          deriv(j) = 1.d0
                      else
                          deriv(j) = x_data(1,i)**j
                      end if
                  end do
 
 
******            Increment normal equations
                  do j = 0, poly_order
                      poly_vec(j) = poly_vec(j) + y_data(1,i)*deriv(j)*
     .                                            wgh
                      do k = 0, j
                          poly_mat(jel(j+1,k+1)) =
     .                         poly_mat(jel(j+1,k+1)) + deriv(j)*
     .                                                  deriv(k)*wgh
                      end do
                  end do
 
*                         ! We should use
              end if
*                         ! Inside the window
          end if
*                         ! Looping over the data
      end do
 
***** Now finish the solution
 
      if( num_poly_data.ge.poly_order+1 ) then
          call SMNV8( poly_mat, poly_vec, poly_order+1, 1, pvrow, pvcol,
     .                iuse, pvrwb )
      else
          write( termlu,100)
 100      format(/' POLY_FIT Error: No usable data in window')
          RETURN
      end if
 
***** Now get the statistics
      resy_min =  1.d20
      resy_max = -1.d20

      do i = 1, num_data
 
*         See if in window first
 
          if( x_data(1,i).ge. poly_window(1) .and.
     .        x_data(1,i).le. poly_window(2) .and.
     .        y_data(1,i).ge. poly_window(3) .and.
*                                                         ! Inside window
     .        y_data(1,i).le. poly_window(4)       ) then
 
*             See if we should use
              USE = .true.
*                                          ! UNW field given
              if( p_field(2).ne.0 ) then
                  if( pt_data(2,i).ne.0 ) then
                      use = .false.
                  end if
              end if

*             See if fails point type
              if( p_field(1).ne.0 .and.
     .            pt_data(1,i).ne.fit_point .and.
     .            fit_point.ne.0                 ) then
                     use = .false.
              end if

 
*             Now get the weight
              wgh = 1.d0
*                                         ! Errbar given, so use
              if( y_field(3).ne.0 ) then
                  wgh = 1.d0/y_data(2,i)**2
              end if
 
*             If we are using this point then incerment normal equations
              do j = 0, poly_order
                  if( x_data(1,i).eq.0 .and. j.eq.0 ) then
*                                     ! Trap 0**0 error
                      deriv(j) = 1.d0
                  else
                      deriv(j) = x_data(1,i)**j
                  end if
              end do

******        Increment statistics
              estimate = 0.d0
              do j = 0, poly_order
                  estimate = estimate + poly_vec(j)*deriv(j)
              end do
 
              if( use ) then
                  sum_poly(1) = sum_poly(1) + (y_data(1,i)-estimate)**2
     .                                       *wgh
                  sum_poly(2) = sum_poly(2) + wgh
*                         ! We should use
              end if
*                         ! Inside the window
*             See if we should detrend just data in window
              if( detrend(1:1).eq.'W' ) then
                  y_data(1,i) = y_data(1,i)-estimate
              end if
          end if
*                         ! Looping over the data

*         See if we should detrend just data in window
          if( detrend(1:1).eq.'A' ) then
*             If we are using this point then incerment normal equations
              do j = 0, poly_order
                  if( x_data(1,i).eq.0 .and. j.eq.0 ) then
*                                     ! Trap 0**0 error
                      deriv(j) = 1.d0
                  else
                      deriv(j) = x_data(1,i)**j
                  end if
              end do

******        Increment statistics
              estimate = 0.d0
              do j = 0, poly_order
                  estimate = estimate + poly_vec(j)*deriv(j)
              end do
              y_data(1,i) = y_data(1,i)-estimate
          end if

* MID TAH 131021: See if data is in X-scale range
          use = .false.
          if( x_data(1,i).ge. scale_size(1) .and.
     .        x_data(1,i).le. scale_size(2)  ) then
 
*             See if we should use
              use = .true.
*                                          ! UNW field given
              if( p_field(2).ne.0 ) then
                  if( pt_data(2,i).ne.0 ) then
                      use = .false.
                  end if
              end if

*             See if fails point type
              if( p_field(1).ne.0 .and.
     .            pt_data(1,i).ne.fit_point .and.
     .            fit_point.ne.0                 ) then
                     use = .false.
              end if
             
              if( use .and. y_data(1,i).lt. resy_min ) 
     .                           resy_min = y_data(1,i)
              if( use .and. y_data(1,i).gt. resy_max ) 
     .                           resy_max = y_data(1,i)
           end if

      end do
      
****  Data has been detrended then set the reference value to
*     zero
      refsave_y = ref_valy 
      if( detrend(1:1).eq.'A' .or. detrend(1:1).eq.'W' ) then
         ref_valy = 0.0d0
      end if
 
***** Now build the polynomial labels
      if( num_poly_data-poly_order-1.gt.0 ) then 
          chisq = sum_poly(1)/(num_poly_data-poly_order-1)
      else
          chisq = 1.0
      end if
      wrms  = sqrt( (num_poly_data/sum_poly(2) )*chisq ) * conv_poly(2)
*                                 ! Used for scaling
      nrms  = sqrt( chisq )
      if( nrms.eq.0 ) nrms = 1.d0
 
      poly_labels(1) = 'Reference '
      value = ref_valx
 
*                                 ! Write calanender data
      if( x_field(1).eq.0 ) then
          call mjd_to_ymdhms(value, dates, sectag)
c         call epoc_8(dates(2), dates(3), dates(1),
c    .                dates(4), dates(5), value )
 
          if( x_field(3).lt.6 ) then
              write(poly_labels(1)(12:),200)
     .            (dates(j),j=1, 5)
          else
              write(poly_labels(1)(12:),200)
     .            (dates(j),j=1, 5), sectag
 200          format(i4,"/",i2,"/",i2,1x,i2,":",i2,1x,F6.2)
          endif
      else
 
*         Generate format for this value
          call format_label(scale_size(1), scale_size(2), value,
     .                      labform )
          write(poly_labels(1)(12:), labform, iostat=ierr) value
      end if
 
      len_poly = trimlen( poly_labels(1) )
 
      write( poly_labels(1)(len_poly+1:),210, iostat=ierr) num_poly_data
 210  format(' Num ',i8)
 
***** Now do rms line, first generate format
      call format_label( 0.0, wrms, 0.d0, labform)
      poly_labels(2) = 'WRMS '
      write(poly_labels(2)(6:), labform, iostat=ierr ) wrms
 
      len_poly = trimlen(poly_labels(2))
      poly_labels(2)(len_poly+2:) =
     .    poly_units(2)(1:max(1,trimlen(poly_units(2)))) //
     .    ' Chi**2/f'
 
*     Get format for Chi**2/f value
      minform = chisq
      maxform = chisq*1.001
      call format_label(minform,maxform, 0.d0, labform)
      len_poly = trimlen(poly_labels(2))
      write(poly_labels(2)(len_poly+2:),labform, iostat=ierr)
     .      chisq
 
***** Do the coefficents of the polynomial
      do j = 0, poly_order
          if( j.gt.1 ) then
              write(poly_labels(j+3), 300, iostat=ierr) j
  300         format('Coeff',i2)
          else
              if( j.eq.0 ) poly_labels(j+3) = 'Intercept'
              if( j.eq.1 ) poly_labels(j+3) = 'Slope '
          end if
 
*         Get format for the value
*                                                ! Convert to user units
          scaling = conv_poly(2)*conv_poly(1)**j
          result  = poly_vec(j)*scaling
          sigma   = sqrt(poly_mat(jel(j+1,j+1)))*nrms*scaling
 
          minform = poly_vec(j)* scaling
          maxform = minform  + sigma
 
          if( j.eq.0 ) then
              value = refsave_y*scaling
          else
              value = 0.d0
          end if
 
          call format_label(minform, maxform, value, labform)
 
*         Write out the value
          write(poly_labels(j+3)(12:), labform, iostat=ierr)
     .        result+value
          len_poly = trimlen(poly_labels(j+3))
          poly_labels(j+3)(len_poly+2:) = '+-'
 
*         Now format for the sigma
          call format_label( minform, maxform, 0.d0, labform)
          write( poly_labels(j+3)(len_poly+4:), labform, iostat=ierr)
     .           sigma
 
*         Now add labels (For offset do not add x units (not needed))
          len_poly = trimlen(poly_labels(j+3))
          if( j.ne.0 ) then
              poly_labels(3+j)(len_poly+2:) =
     .           poly_units(2)(1:max(1,trimlen(poly_units(2)))) // '/'
     .        // poly_units(1)(1:max(1,trimlen(poly_units(1))))
*                                ! Only add power for greater than slope
              if( j.gt.1 ) then
                  len_poly = trimlen(poly_labels(j+3))
                  write(poly_labels(j+3)(len_poly+1:),330,iostat=ierr) j
  330             format('**',i1)
              end if
          else
              poly_labels(3+j)(len_poly+2:) =
     .            poly_units(2)(1:max(1,trimlen(poly_units(2))))
          end if
      end do
 
***** Now write out the labels
      do i = 1, poly_order+3
          len_poly = max(1,trimlen(poly_labels(i)))
          write(termlu,400, iostat=ierr) i, poly_labels(i)(1:len_poly)
 400      format(1x,i2,2x,a)
      end do

* MOD TAH 131016: See if we should reset scales.  (Don't reset reset_scales
*     is false.
* MOD TAH 131029: Only reset scale if we are removing trend.
      if ( reset_scales .and. 
     .    (detrend(1:1).eq.'A' .or. detrend(1:1).eq.'W') ) then
         scale_set = .false.
         use_def_scale = .false.
         scale_size(3) = resy_min - (resy_max-resy_min)*0.05
         scale_size(4) = resy_max + (resy_max-resy_min)*0.05
         call set_scale 
      end if   

 
***** Thats all
      return
      end
 
