CTITLE ADD_ION
 
      subroutine add_ion
 

      implicit none 
 
*     Routine to add the ionospheric delay to the theoretical
*     delays and rates.  For this routine the data sigmas are
*     also updated when the correction is applied.  If the quality
*     of the ion correction does not pass the ion_mask bits
*     then the data_flag will be set to show ion bad.
*
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
      include '../includes/obs_data.h'
 
*   in      - Index to be used in ion_corr to set group
*           - or phase correcion (in=1 for group; in=2
*           - for phase)
 
      integer*4 in, cand
 
*   kbit    - Function to check if a bit is set.
 
 
      logical kbit
 
***** See if we want to apply ion correction 
 
      in = 0
 
      IF( kbit(cont_baseline,6) .or. kbit(cont_baseline,7) ) THEN
*                                             apply ion correction
*         See if it is available
          IF( (kbit(cont_baseline,6).and.kbit(avail_baseline,6)).or.
     .        (kbit(cont_baseline,7).and.kbit(avail_baseline,7)) ) THEN
*                                      desired ion correction available
 
*             Now see if ion correction passed the ion_mask condition
              if( kbit(cont_baseline,6) ) in = 1
              if( kbit(cont_baseline,7) ) in = 2
 
*             Now apply correction
              theoretical(1) = theoretical(1) + ion_corr(in)
              theoretical(2) = theoretical(2) - ion_corr(in)
              theoretical(3) = theoretical(3) + ion_corr(in)
*                                                            ! Rate
              theoretical(4) = theoretical(4) + ion_corr(3)
 
*             Now update data sigmas for noise in ion correction
* MOD TAH 900129: Added check to make sure ion sigma is not too large
*             Otherwize a floating pointing overflow occurrs
              if( ion_sigma(in).gt.1.d6 ) then
                  ion_sigma(in) = 1.d6
              end if
              sigma(1) = sqrt(sigma(1)**2 + ion_sigma(in)**2)
              sigma(2) = sqrt(sigma(2)**2 + ion_sigma(in)**2)
              sigma(3) = sqrt(sigma(3)**2 + ion_sigma(in)**2)

              if( ion_sigma(3).gt.1.d6 ) then
                  ion_sigma(3) = 1.d6
              end if

              sigma(4) = sqrt(sigma(4)**2 + ion_sigma(3)**2)
 
*             Now see if we even use this data i.e., iand of
*             ion_mask and ion_flag must be zero for us to use it.
 
*                                 ! Group delay
              if( in.eq.1 ) then
                  if( cand(ion_mask,ion_flag).ne.0 ) then
                      call sbit(data_flag,2,1)
                  else
*                                                ! Reset if OK
                      call sbit(data_flag,2,0)
                  end if
*                                 ! Phase delay
              else
                  if( cand(ion_mask*64,ion_flag).ne.0 ) then
                      call sbit(data_flag,2,1)
                  else
*                                                ! Reset if OK
                      call sbit(data_flag,2,0)
                  end if
              end if
 
*                         ! Ion not available, reset contribution so
          ELSE
*                         ! that it will not be tried again
*                                           ! Reset both, only one
              call sbit(cont_baseline,6,0)
*                                           ! should be on anyway
              call sbit(cont_baseline,7,0)
          END IF
*                         ! if we were to apply
      END IF
 
***** Thats all
      return
      end
 
