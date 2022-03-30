*     Include file seasonal.f  This file containes the common block
*     to be used with this program.
 
* PARAMETERS
 
*   max_rsites  - Maximum number of sites
*   max_robs        - Maximum number of data per site.
 
      integer*4 max_rsites, max_robs
 
      parameter ( max_rsites = 14 )
*                                     ! Enough for 4 years of data.
      parameter ( max_robs = 3000 )
 
* DECLARATIONS
 
*   epochs(max_robs, max_rsites)    - Times of the observations
 
      real*8 epochs(max_robs, max_rsites)
 
*   num_robs(max_rsites)    - Number of observations per site.
*   num_used(max_rsites) - Actual number of observations used
*                       - in fit
*   edit_flag(max_robs, max_rsites) - Edit flag on data
 
 
      integer*4 num_robs(max_rsites), num_used(max_rsites),
     .    edit_flag(max_robs, max_rsites)
 
*   site_lat(max_rsites)        - Latitudes of the sites
*   site_ht(max_rsites)     - Heights of the sites.
*   temps(max_robs, max_rsites) - Projected surface
*                               - temperature
*   tbias(max_robs, max_rsites) - Biases in the surface
*                           - temperatures
*   seas_temp(3, max_rsites)    - Mean, cos(month), sin(month)
*                           - fit to temperature.
*   seas_bias(3, max_rsites)    - Mean, cos(month), sin(month)
*                           - fit to temperature bias.
*   rms_temp(max_rsites)        - RMS fit to temperature
*   rms_bias(max_rsites)        - RMS fit to temperature bias
*   temp_mean(3)                - mean, sin(latitude), ht (km)
*                           - fit to temperature
*   temp_cosm(3)                - mean, sin(latitude), ht (km)
*                           - fit to cos(month) term
*   temp_sinm(3)                - mean, sin(latitude), ht (km)
*                           - fit to sin(month) term
*   tbias_mean(3)           - mean, sin(latitude), ht (km)
*                           - fit to tbias
*   tbias_cosm(3)           - mean, sin(latitude), ht (km)
*                           - fit to cos(month) tbias
*   tbias_sinm(3)           - mean, sin(latitude), ht (km)
*                           - fit to sin(month) tbias
 
*   rms_global(6)           - RMS fit of all of the terms in
*                           - temp and bias seasonal model
 
      real*4 site_lat(max_rsites), site_ht(max_rsites),
     .    temps(max_robs, max_rsites), tbias(max_robs, max_rsites),
     .    seas_temp(3, max_rsites), seas_bias(3, max_rsites),
     .    rms_temp(max_rsites), rms_bias(max_rsites), temp_mean(3),
     .    temp_cosm(3), temp_sinm(3), tbias_mean(3), tbias_cosm(3),
     .    tbias_sinm(3), rms_global(6)
 
*   site_codes(max_rsites)  - These are the WMO number (eg
*                           - 72518)
*   site_names(max_rsites)  - Names of the sites
 
 
      character*8 site_codes(max_rsites), site_names(max_rsites)
 
*---------------------------------------------------------------------
*     COMMON DECLARATION
 
      common / seas_com / epochs, num_robs, num_used, edit_flag,
     .    site_lat, site_ht, temps, tbias, seas_temp, seas_bias,
     .    rms_temp, rms_bias, temp_mean, temp_cosm, temp_sinm,
     .    tbias_mean, tbias_cosm, tbias_sinm, rms_global, site_codes,
     .    site_names
 
*----------------------------------------------------------------------
 
 
