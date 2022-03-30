CTITLE GLB_SAVE_APR
 
      subroutine glb_save_apr( type, indx, apr )

      implicit none  
 
*     Routine to save the apriori/site source or axis offset parameter
*     in the list of apriri values.
 
      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/obs_header.h'
 
      include '../includes/globk_common.h'
 
*   comp        - Component for nutation series terms
*   indx        - Index in array of values to be saved (eg. site
*               - number or source number)
*   phase       - 1 if inphase and 2 if out of phase
*   type        - type of parameter to be saved
 
      integer*4 comp, indx, phase, type

*   orb_el      - Orbital element decode from index
*   sv_num      - Number of satellite.

      integer*4 sv_num, orb_el
 
*   apr         - Place to save the result
*   dt          - Time difference between reference and solution
*               - epoch
 
      real*8 apr, dt
 
 
***** Use a series of if statements.  We only want to do site
*     positions and rates, axis offsets and rates, and source
*     positions and rates
 
*     X site position
      if( type.eq.7 ) then
          dt = (gepoch_out-site_epoch(indx))/365.25d0
          apr = apr_val_site(1,1,indx) + apr_val_site(1,2,indx)*dt
      end if
 
*     Y site position
      if( type.eq.8 ) then
          dt = (gepoch_out-site_epoch(indx))/365.25d0
          apr = apr_val_site(2,1,indx) + apr_val_site(2,2,indx)*dt
      end if
 
*     Z site position
      if( type.eq.9 ) then
          dt = (gepoch_out-site_epoch(indx))/365.25d0
          apr = apr_val_site(3,1,indx) + apr_val_site(3,2,indx)*dt
      end if

* MOD TAH 040703: Atmospheric delay
      if( type.eq.61 ) then
          apr = apr_val_atm(indx)
      endif
 
*     Axis offset
      if( type.eq.10 ) then
          dt = (gepoch_out-axo_epoch(indx))/365.25d0
          apr = apr_val_axo(1,indx) + apr_val_axo(2,indx)*dt
      end if
 
*     RA of source
      if( type.eq.11 ) then
          dt = (gepoch_out-source_epoch(indx))/365.25d0
          apr = apr_val_source(1,1,indx) + apr_val_source(1,2,indx)*dt
      end if
 
*     Dec of source
      if( type.eq.12 ) then
          dt = (gepoch_out-source_epoch(indx))/365.25d0
          apr = apr_val_source(2,1,indx) + apr_val_source(2,2,indx)*dt
      end if
 
*     X-position rate
      if( type.eq.42 ) then
          apr= apr_val_site(1,2,indx)
      end if
 
*     Y-position rate
      if( type.eq.43 ) then
          apr= apr_val_site(2,2,indx)
      end if
 
*     Z-position rate
      if( type.eq.44 ) then
          apr= apr_val_site(3,2,indx)
      end if
 
*     Axis offset rate
C     if( type.eq.35 ) then
C         apr= apr_val_axo(2,indx)
C     end if
 
*     RA rate
      if( type.eq.45 ) then
          apr = apr_val_source(1,2,indx)
      end if
 
*     Declination rate
      if( type.eq.46 ) then
          apr = apr_val_source(2,2,indx)
      end if
 
*     Nutation series coefficients
      if( type.eq.47 ) then
          comp = mod(indx,128)
          phase = (indx-comp)/256 + 1
          apr = apr_val_nut_coeff(phase,comp)
      end if

*     Orbital elements
      if( type.eq.51 ) then
          call decode_code(indx, orb_el, sv_num)
          apr = apr_val_svs(orb_el, sv_num)
      end if

 
***** Thats all
      return
      end
 
