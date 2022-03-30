CTITLE 'WT_FILTER_FWHM'
 
      SUBROUTINE WT_FILTER_FWHM ( x_data, sigma, num_data, FWHM )

      implicit none 
*    .,      29MAY87                <871210.1317>
 
*     Modifier: J. Davis           2:36 PM  FRI., 29  MAY , 1987
*
*     Based on routine FILTER_FWHM by T. Herring.  Modified to
*     include sigmas in calculation of FWHM.
*
*     Author: T.Herring           11:02 AM  WED.,  9  JULY, 1986
*
*$S "WT_FILTER_FWHM"
*------------------------------------------------------------------------
*     Routine to compute an appropriate Full Width at Half Maximum
*     value of a set of data sampled at values given by x_data.
*     The value computed corresponds to a sigma in the filter of
*     0.5 times the rms of the x_data spacings.
*
*     The routine assumes that the x_data are arranged in a monotonic
*     order (either increasing or decreasing).
*
*     CALLING SEQUENCE:
*     =================
*     CALL wt_Filter_FWHM ( x_data, sigma, num_data, FWHM )
*
*     WHERE:
*     x_data      is an array of x values at which samples have
*                 been taken (e.g., times).  The x_data values
*                 are assumed to be monotonic.
*                 (REAL*8 array INPUT)
*     sigma       is an array of sigmas for the corresponding Y
*                 values
*                 (REAL*8 array INPUT)
*     num_data    is the number of data in x_data.
*                 (I*2 INPUT)
*     FWHM        is the appropriate FWHM of the filter to be used with
*                 data so that it is not excessively smoothed i.e.,
*                 the filter can be used as an interpolator rather than
*                 as a filter.
*                 (REAL*8 OUTPUT)
*
*---------------------------------------------------------------------
*$E
 
 
*          i    - loop counter
*   num_data    - number of data in x_data
 
      integer*4 i, num_data
 
*       FWHM        - Full width half maximum to be used.
*   rms_diff    - rms scatter of the differences in x_data
*   sum_diffsq  - sum of the squares of the differences
*   sum_weights - sum of the weights
*   weight      - Weight for this square difference
 
*   sigma(1)    - the array of sigmas
*   x_data(1)   - the array of x_values.
 
      real*8 FWHM, rms_diff, sum_diffsq, sum_weights, weight, 
     .    sigma(num_data), x_data(num_data)
 
***** START, compute the rms scatter of the differences in x_data
 
      sum_diffsq  = 0
      sum_weights = 0
 
      do i = 1, num_data - 1
 
          weight      = 1.0D0 / (sigma(i) ** 2 + sigma(i+1) ** 2)
          sum_diffsq  = sum_diffsq + weight*(x_data(i)-x_data(i+1))**2
          sum_weights = sum_weights + weight
 
      end do
 
*                                   ! compute rms
      if ( sum_weights.gt.0 ) then
          rms_diff = sqrt( sum_diffsq/sum_weights)
*                                ! set an arbitrary value
      else
          rms_diff = 1.d0
      end if
 
***** Now compute FWHM
 
      FWHM = 0.5d0* rms_diff* 2.35482d0
 
      end
 
