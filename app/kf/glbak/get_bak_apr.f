CTITLE GET_BAK_APR
 
      subroutine get_bak_apr( iel, code, sol_parm, obs_corr )

      implicit none  
 
*     Routine to decode the codes in code and to copy the element
*     from sol_parm to obs_corr.  This routine will allow us later
*     to generate an observed parameter value from a local solution
*     with the translations and orientations removed.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   code    - The current code we are decoding
 
*   iel     - Element in obs_corr into which the results should
*           - be saved.
*   indx    - Indx from current code
*   jel     - element number from sol_parm to use
*   type    - type of parameter
 
      integer*4 code, iel, indx, jel, type
 
*   obs_corr(1) - correction which will be applied to the observation
*   parm_to_obs - Real*8 function to copy a parameter adjustement
*               - estimate to the observation correction vecotor.
*   sol_parm(1) - Current solution vector
 
      real*8 obs_corr(1), parm_to_obs, sol_parm(1)
 
 
***** Get the parameter type and then go to correct code
 
      call decode_code ( code, type, indx )

* MOD TAH 040703: See if atm is estimated, then set the type to 28.
      if( type.eq.61 ) type = 28
 
*                                             ! These values should not
      if( type.lt.7 .or. type.gt.28 ) return
*                                             ! appear
 
***** Use a GOTO
 
      GOTO (  100,  200,  300,  400,  500,  600,  700,  800,  900, 1000,
     .       1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
     .       2100, 2200 ) type - 6
 
***** X-site position
  100 continue
          jel = parn_site(1,1,ltog_sites(indx))
          obs_corr(iel) = parm_to_obs( jel, sol_parm)
          RETURN
 
***** Y-site position
  200 continue
          jel = parn_site(2,1,ltog_sites(indx))
          obs_corr(iel) = parm_to_obs(jel, sol_parm)
          RETURN
 
***** Z-site position
  300 continue
          jel = parn_site(3,1,ltog_sites(indx))
          obs_corr(iel) = parm_to_obs(jel, sol_parm)
          RETURN
 
***** Axis offset
  400 continue
          jel = parn_axo(1,ltog_sites(indx))
          obs_corr(iel) = parm_to_obs(jel, sol_parm)
          RETURN
 
***** RA of source
  500 continue
          jel = parn_source(1,1,ltog_sources(indx))
          obs_corr(iel) = parm_to_obs(jel, sol_parm)
          RETURN
 
***** Dec of source
  600 continue
          jel = parn_source(2,1,ltog_sources(indx))
          obs_corr(iel) = parm_to_obs(jel, sol_parm)
          RETURN
 
***** Wobble components
  700 continue
          obs_corr(iel) = parm_to_obs(parn_wob(indx), sol_parm)
          RETURN
 
***** UT1 values
  800 continue
          obs_corr(iel) = parm_to_obs(parn_ut1(indx), sol_parm)
          RETURN
 
***** Nutation angles
  900 continue
          obs_corr(iel) = parm_to_obs(parn_nut_ang(indx), sol_parm)
          RETURN
 
***** Tide l
 1000 continue
          jel = parn_tid(1,ltog_sites(indx))
          obs_corr(iel) = parm_to_obs(jel, sol_parm)
          RETURN
 
***** Tide h
 1100 continue
          jel = parn_tid(2,ltog_sites(indx))
          obs_corr(iel) = parm_to_obs(jel, sol_parm)
          RETURN
 
***** Lag angle
 1200 continue
          jel = parn_tid(3,ltog_sites(indx))
          obs_corr(iel) = parm_to_obs(jel, sol_parm)
          RETURN
 
***** Extendented earth tides 1
 1300 continue
          obs_corr(iel) = parm_to_obs(0, sol_parm)
          RETURN
 
***** Extendented earth tides 2
 1400 continue
          obs_corr(iel) = parm_to_obs(0, sol_parm)
          RETURN
 
***** Extendented earth tides 3
 1500 continue
          obs_corr(iel) = parm_to_obs(0, sol_parm)
          RETURN
 
***** Extendented earth tides 4
 1600 continue
          obs_corr(iel) = parm_to_obs(0, sol_parm)
          RETURN
 
***** Extendented earth tides 5
 1700 continue
          obs_corr(iel) = parm_to_obs(0, sol_parm)
          RETURN
 
***** Extendented earth tides 6
 1800 continue
          obs_corr(iel) = parm_to_obs(0, sol_parm)
          RETURN
 
***** Extendented earth tides 7
 1900 continue
          obs_corr(iel) = parm_to_obs(0, sol_parm)
          RETURN
 
***** Extendented earth tides 8
 2000 continue
          obs_corr(iel) = parm_to_obs(0, sol_parm)
          RETURN
 
***** Gamma
 2100 continue
          obs_corr(iel) = parm_to_obs(parn_gamma, sol_parm)
          RETURN

****  Atm delay
 2200 continue
          jel = parn_atm(ltog_sites(indx))
          obs_corr(iel) = parm_to_obs(jel, sol_parm)
          RETURN
         
 
***** Thats all
      end
 
