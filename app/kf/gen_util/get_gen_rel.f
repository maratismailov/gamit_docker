CTITLE GET_GEN_REL_CONT
 
      subroutine get_gen_rel_cont ( temp_cont )

      implicit none 
 
 
*     Routine to get the general relativity contribution.  This routine is
*     needed to fix the sign of the potential term in CALC Ver 6.0.

* MOD TAH 900710: Changed routine so that Calc Ver 6.0 delay, with relativity
*     coorection applied, will yield a delay the same as that fromnm Calc 7.0
*     (Nominally the SHAP DEL value)
 
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
 
      include '../includes/obs_values.h'
      include '../includes/obs_apr.h'
      include '../includes/obs_data.h'
 
*   temp_cont(2)    - Return values for the contribution.  Set equal
*                   - to gen_rel_cont for non-Calc Ver 6.0 KalObs
*                   - files.
 
      real*4 temp_cont(2)
 
*   C1              - Same as Calc Ver 6.0 C1 term
*   rsun            - Distance to sun from earth (m)
*   sdotb           - Unit vector to source dotted with site 1-site 2
*                   - baseline (same definition as CALC)
*   sdbdt           - Rate of change of sdotb (same as Calc)
 
      real*8 C1, rsun, sdotb, sdbdt
 
      if( Calc_ver.ne. 6.0 ) then
          temp_cont(1) = gen_rel_cont(1)
          temp_cont(2) = gen_rel_cont(2)
*                         ! Correct the sign on the potential term
      else
 
          rsun = sqrt( ephem_sun(1,1,2)**2 + ephem_sun(2,1,2)**2 +
*                                            ! Distance to sun a Mid_epoch
     .                 ephem_sun(3,1,2)**2 )
 
          C1 = 2.d0*Gm_sun/( vel_light**3*rsun)
 
          sdotb = -db_theoretical(1)*vel_light
          sdbdt = -db_theoretical(4)*vel_light
 
*         Now change the sign of the C1 term in the contribution
 
          temp_cont(1) = gen_rel_cont(1) - 2*( -C1*sdotb )
          temp_cont(2) = gen_rel_cont(2) - 2*( -C1*sdbdt )

* MOD TAH 900719: Add on the additional bending term.  Effectively
*         undoes one of C1*sdotb terms above.
          temp_cont(1) = temp_cont(1) - C1*sdotb
          temp_cont(2) = temp_cont(2) - C1*sdbdt
 
      end if
 
****  That all
      return
      end
 
 
