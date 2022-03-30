CTITLE 'GAUSS_WRMS'
 
      subroutine gauss_wrms(x_data,y_data,sig,num_data,x_bar,y_bar,
     .    fwhm,wrms)

      implicit none 
 
*     J.L. Davis                   6:09 PM  TUE., 30  JUNE, 1987
*
*     This routine returns the WRMS value of a series of points about
*     th value Y_BAR occuring at X_BAR.  The weights include not only the
*     sigmas, but a Gaussian filter shape (see below).
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
*                      |sum[weight*delta_y**2]|
*     WRMS     =  sqrt|----------------------- |
*                      |      sum[weight]     |
*
*     where
*
*     tau     = FWHM/2.35482     (tau is the equivalent of the standard
*                                 deviation of Gaussian curve)
*
*     weight  =   sum[ exp(-(x_data-x_bar)**2/2*tau**2) ] / sigma**2
*
*     When the value of the exponent is less than 1d-5 we take 1d-5 as
*     the value of the exponent.  This procedure ensures that the filter
*     calculations never result in 0 divided by 0 which occurs if all
*     the data is further than 13.2 tau away from the x_value at which
*     the smoothed value is to be obtained.
 
 
*         i         - Loop counter
*   ,   num_data    - Number of data in x_data and y_data arrays
 
      integer*4 i, num_data
 
*       exp_term    - The exponential term
*   ,   FWHM        - Full width at half maximum of the filter
*                   -   (same units as x_data
*   ,   tau         - The sigma of the filter (FWHM/2.35482)
*   ,   sum_square  - Sum of the weighted residuals
*   ,   sum_weight  - Sum of the weights
*   ,   weight      - Value of the exp of the filter
*   ,   sig(1)      - The array of sigmas
*   ,   x_data(1)   - The array of x values
*   ,   y_data(1)   - The array of y values to be smoothed
*   ,   x_bar       - The value of x corresponding to Y_BAR
*   ,   y_bar       - The value around which the WRMS is to be taken
*   ,   wrms        - Weighted root-mean-square
 
      real*8 exp_term, FWHM, tau, sum_square, sum_weight, weight,
     .    sig(1), x_data(1), y_data(1), x_bar, y_bar, wrms
 
***** Compute tau from the FWHM
      tau = FWHM / 2.35482d0
 
***** Clear summation variables and do the sum
      sum_weight = 0
      sum_square = 0
 
      do i = 1, num_data
 
          exp_term = exp(-((x_data(i) - x_bar) / tau) ** 2 / 2.d0)
          if (exp_term .lt .1.d-5) then
              exp_term = 1.d-5
          end if
 
          weight = exp_term / sig(i) ** 2
 
          sum_weight = sum_weight + weight
          sum_square = sum_square + weight * (y_data(i) - y_bar) ** 2
 
      end do
 
***** Now compute WRMS
      if (num_data .gt. 0) then
          wrms = sqrt(sum_square / sum_weight)
      else
          wrms = 0
      end if
 
      end
 
