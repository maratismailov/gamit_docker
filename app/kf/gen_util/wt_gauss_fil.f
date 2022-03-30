CTITLE 'WT_GAUSS_FILTER'
 
      SUBROUTINE WT_GAUSS_FILTER ( x_data,  y_data,  sigma , num_data,
     .                       x_value, y_value, sig_value, FWHM      )

      implicit none 

*    .,      29MAY87                <871210.1317>
 
*     Modifier: J. Davis           2:20 PM  FRI., 29  MAY , 1987
*
*     Based on the routine GAUSS_FILTER by T. Herring.  Modified to use
*     sigmas of Y values.
*
*     Author: T.Herring           10:05 AM  WED.,  9  JULY, 1986
*
*$S "WT_GAUSS_FILTER"
*----------------------------------------------------------------------
*     MOD JLD 870630 To return the sigma of the y value and to check
*     that the number of calling arguments is correct
*
*     Gaussian Filter utility subroutine to smooth data.  This routine
*     takes two sets of data (x_data and y_data) and will return the
*     gaussian filter smoothed estimate of the y_value at some x_value.
*     The Full width at half maximum is supplied by the user (FWHM)
*     The formulation of this filter is slightly different to the usual
*     formulation in that this filter will bias the results slightly
*     towards the mean value of y_data.  Consquently, when filtered
*     estimates are obtained for x_values well away from the span
*     of data (>10 FWHM's), the mean value of the y_data will returned
*     rather than the first or last value of y_data (the usual return
*     from Gaussian filters).
*
*     The filter is formulated as:
*
*                 sum[ y_data*weight]
*     y_value  =  -------------------
*                     sum[weight]
*
*                 sqrt{sum[weight**2 * sigma**2]}
*     sig_value = -------------------------------
*                          sum[weight]
*
*     where
*     tau     = FWHM/2.35482     (tau is the equivalent of the standard
*                                 deviation of Gaussian curve)
*     weight  =   sum[ exp(-(x_value-x_data)**2/2*tau**2) ] / sigma**2
*
*     When the value of the exponent is less than 1d-5 we take 1d-5 as
*     the value of the exponent.  This procedure ensures that the filter
*     calculations never result in 0 divided by 0 which occurs if all
*     the data is further than 13.2 tau away from the x_value at which
*     the smoothed value is to be obtained.
*
*     Each call to the filter returns a single estimate of y at the
*     value of x passed to the routine.
*
*     CALLING SEQUENCE:
*     =================
*     CALL wt_Gauss_filter ( x_data,  y_data, sigma, num_data,
*    .                       x_value, y_value, sig_value, FWHM    )
*
*     WHERE:
*     x_data      is an array of x_values (e.g., times)
*                 (REAL*8 array INPUT)
*     y_data      is an array of y_values to be smoothed
*                 (REAL*8 array INPUT)
*     sigma       is an array of sigmas corresponding to Y_DATA
*                 (REAL*8 array INPUT)
*     num_data    is the number of data samples in the x_data and y_data
*                 arrays.
*                 (I*2 INPUT)
*     x_value     is the value of x at which the smoothed estimate of
*                 y is to be obtained.
*                 (REAL*8 INPUT)
*     y_value     is the smoothed estimate of y at the x_value
*                 (REAL*8 OUTPUT)
*     sig_value   is the sigma of the returned Y_VALUE for assumed
*                 uncorrelated data
*                 (REAL*8 OUTPUT)
*     FWHM        is the full width half maximum of the filter (This
*                 value is 2.35482 times the standard deviation of
*                 Gaussian curve)
*                 (REAL*8 INPUT)
*
*----------------------------------------------------------------------
*$E
 
 
*          i    - Loop counter
*   num_data    - Number of data in x_data and y_data arrays
 
      integer*4 i, num_data
 
*   exp_term    - The exponential term
*   FWHM    - Full width at half maximum of the filter (same units
*           - as x_data
*   tau     - the sigma of the filter (FWHM/2.35482)
 
*   sig_value   - Sigma of Y_VALUE
*   sum_data    - sum of the data by the gaussian weights
*   sum_weight  - sum of the weights
 
*   weight  - value of the exp of the filter
 
 
*   sigma(1)    - The array of sigmas
*   x_data(1)   - the array of x values
*   y_data(1)   - the array of y values to be smoothed
 
*   x_value     - the value of x at the smoothed estimate is to be
*               - obtained
*   y_value     - the smoothed estimate of y at x_value
 
 
 
      real*8 exp_term, FWHM, tau, sig_value, sum_data, sum_weight,
     .    weight, sigma(num_data), x_data(num_data), y_data(num_data), 
     .    x_value, y_value
 
***** START, compute tau from the FWHM
 
      tau = FWHM / 2.35482d0
 
***** Clear summation variables and do the sum
 
      sum_data   = 0.d0
      sum_weight = 0.d0
      sig_value  = 0
 
      do i = 1, num_data
 
          exp_term = exp( -((x_data(i)-x_value)/tau)**2 / 2.d0 )
          if( exp_term.lt.1.d-5 ) then
              exp_term = 1.d-5
          end if
 
          weight = exp_term / sigma(i) ** 2
 
          sum_data   = sum_data   + y_data(i)*weight
          sum_weight = sum_weight + weight
          sig_value  = sig_value  + (weight * sigma(i)) ** 2
 
      end do
 
***** Now compute smoothed value
 
      if( num_data.gt.0 ) then
          y_value   = sum_data/sum_weight
          sig_value = sqrt(sig_value) / sum_weight
      else
*                             ! no data, as good a guess as any
          y_value   = 0
          sig_value = 0
      end if
 
      end
 
