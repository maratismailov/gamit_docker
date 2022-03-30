
 
CTITLE SD_GET_TIDES
 
      subroutine sd_get_tides( site_num, fund_arg, tide_vals )

      implicit none 
 
*     This routine scans the list of tide corrections to be
*     applied and computes the contributions at site 'site_num'
*     and saves the values in tide_vals.  (The diurnal tides
*     are in the first six elements and the semidiurnals in the
*     last six elements.
*     Here we call sd_arg with only 5 arguments because GST is
*     already built into the tidal partials.
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include '../includes/sd_common.h'
 
* PASSED VARIABLES
 
*   site_num        - Number of current site
 
      integer*4 site_num
 
*   fund_arg(6)     - Values of the fundamental arguments
*   tide_vals(12)   - Values of the diurnal and semidiurnal
*                   - tide Spherical harmonic coefficients
*                   - for radial, north and East (all m)
 
      real*8 fund_arg(6), tide_vals(12)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   offset  - Change in index between diurnal and semidiurnal
 
      integer*4 i,j, offset
 
*   arg     - Argument of the tide (computed ignoring GST)
 
      real*8 arg
 
*   reported(max_sd_tides)  - Set true when an error in the
*           - GST argumented has been reported.
 
      logical reported(max_sd_tides)
 
      data reported / max_sd_tides*.false. /

* MOD TAH 910130: Quick fix to not having enough sites for tide
*     values

      if( site_num.gt.max_sites ) RETURN
 
 
***** Clear all of the values first
 
      do i = 1,12
          tides_val(i) = 0.0d0
      end do
 
***** Scan the list of tides and evaluate those which apply to this
*     site.
 
      do i = 1, sd_tides_num
          if( site_num.eq.sd_tides_site(i) ) then
              call sd_arg( sd_tides_arg(1,i), fund_arg, arg, 5)
 
*             Now sum up the values
*                                                     ! Diurnal
              if( sd_tides_arg(6,i).eq.-1 ) then
                  offset = 0
*                                                     ! Semidiurnal
              else if( sd_tides_arg(6,i).eq.-1 ) then
                  offset = 6
*                                                 ! Invalid
              else
*                                                 ! argument
                  if( .not.reported(i) ) then
                      write(*,100) i, (sd_tides_arg(j,i),j=1,6)
 100                  format('SD_GET_TIDES ERROR: Invalid GST ',
     .                    ' argument for tide # ',i3,
     .                    '. Arguments are ',6I3,/,
      .                20x,'Remaining tide correctoions ignored')
                      reported(i) = .true.
                  end if
 
                  RETURN
              end if
 
*                             ! Loop in pairs for cos and sin
              do j = 1,6, 2
                  tides_val(j+offset) = tides_val(j+offset) +
     .                sd_tides_val(j,i)* cos(arg)
                  tides_val(j+offset+1) = tides_val(j+offset+1) +
     .                sd_tides_val(j,i)* sin(arg)
              end do
*                         ! Site matches
          end if
*                         ! Looping over all tides input.
      end do
 
****  Thats all
      return
      end
 
 
* additional code for ADD_SD
 
*   tides_val(12)   - Values for the spherical harmonic
*                   - coefficients for the 12 extended
*                   - tidal displacements.
 
      real*8 tides_val(12)
 
***** Compute the tidal contribtions.  Add contributions for both sites.
*                     ! Loop over the two sites
      do i = 1,2
 
*         Get the values of the tides for this site
          call sd_get_tides(site(i), fund_arg, tides_val)
 
*         Now compute contribution
*                         ! Loop over 12 components
          do j = 1,12
*                         ! Loop over delay and rate
              do k = 1,2
                  dt(k) = dt(k) + tides_val(j)*etd_ext_part(j,i,k)
              end do
          end do
      end do
 
****  Code to be added to sd_get_mid
 
****  Get the tidal values at each of the sites
      do i = 1, num_sites
          call sd_get_tides( i, fund_arg, sd_tides_mid(1,i))
      end do
 
 
 
 
