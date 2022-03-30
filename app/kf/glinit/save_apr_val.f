CTITLE save_apr_val
 
      subroutine save_apr_val( type, indx, pos, apr, check, rad_ck)

      implicit none 
 
 
*     Routine to save the value from the apr array in the GLOBK
*     common, if the type matches those which we save. (See
*     GLB_HEADER_DEF.FTNI for definitions of the types.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'
 
*   check       - if value is 1 then new value will be checked
*               - against any non-zero previous error.  An error
*               - will be generated if they do not match and GLINIT
*               - will stop.
*   indx        - site or source number for the type.  For satellites, the
*                 index contains the both the sv number and the orbital
*                 element
*   pos         - Position in the apr array that the value
*               - stored.
*   type        - Parameter or apririo type. eg. type 7 is X
*               - site coordinate.
*   rad_ck      - If the set to 1 then the radiation parameters will be
*                 checked to see which have been estimated.
 
      integer*4 check, indx, pos, type, rad_ck
 
*   apr(1j)     - List of apriori values to be extracted and saved
*               - GLOBK common.
      real*8 apr(*)

      real*8 old_x  ! Old  X coordinates to see if we should
                 !  update epoch.

*   sv_num      - Satellite number decoded from index
*   orb_el      - Orbital element decode from index

      integer*4 sv_num, orb_el
 
***** X-site coordinate
      if( type.eq.7 ) then
          old_x = apr_val_site(1,1,ltog_sites(indx))
          call check_apr( apr_val_site(1,1,ltog_sites(indx)),apr(pos),
     .                    check, type)
          apr_val_site(1,1,ltog_sites(indx)) = apr(pos)
          if( apr_val_site(1,1,ltog_sites(indx)).ne.old_x ) then
              site_epoch(ltog_sites(indx)) = cepoch_expt
          endif 
          RETURN
      end if
 
***** Y-site coordinate
      if( type.eq.8 ) then
          call check_apr( apr_val_site(2,1,ltog_sites(indx)),apr(pos),
     .                    check, type)
          apr_val_site(2,1,ltog_sites(indx)) = apr(pos)
          RETURN
      end if
 
***** Z-site coordinate
      if( type.eq.9 ) then
          call check_apr( apr_val_site(3,1,ltog_sites(indx)),apr(pos),
     .                    check, type )
          apr_val_site(3,1,ltog_sites(indx)) = apr(pos)
          RETURN

      end if
 
***** Axis offset
      if( type.eq.10) then
          call check_apr( apr_val_axo(1,ltog_sites(indx)),apr(pos),
     .                    check, type)
          apr_val_axo(1,ltog_sites(indx)) = apr(pos)
          RETURN
      end if
 
***** RA of source
      if( type.eq.11) then
          call check_apr(apr_val_source(1,1,ltog_sources(indx)),
     .                   apr(pos),check, type )
          apr_val_source(1,1,ltog_sources(indx)) = apr(pos)
          RETURN
      end if
 
***** Dec of source
      if( type.eq.12) then
          call check_apr( apr_val_source(2,1,ltog_sources(indx)),
     .                    apr(pos), check, type )
          apr_val_source(2,1,ltog_sources(indx)) = apr(pos)
          RETURN
      end if

***** See if satellite element
      if( type.eq.51 ) then
          call decode_code( indx, orb_el, sv_num )
          if( rad_ck.eq.1 ) then
             call sbit(type_svs_elem, orb_el, 1)
          end if

*         For radiation parameters, set to a small value
*         if the apriori is zero.  This way they will not
*         be updated with observed values later.
          if( orb_el.gt.6 .and. abs(apr(pos)).eq.0.d0 ) then
              apr(pos) = 1.d-12
          end if
          call check_apr( apr_val_svs(orb_el,ltog_svs(sv_num)),
     .                    apr(pos), check, type )
*         Only save the the last value of it is not axis offests
          if( orb_el.le. max_svs_elem - 3 ) then
             apr_val_svs(orb_el,ltog_svs(sv_num)) = apr(pos)
          end if
*         Save the axis offsets if this is the apriori pass (rad_ck = 0)
          if( orb_el.gt. max_svs_elem - 3 .and. rad_ck.eq. 0 ) then
             apr_val_svs(orb_el,ltog_svs(sv_num)) = apr(pos)
          endif
          RETURN
      end if

* MOD TAH 040703: Check for atmospheric delay apriori
      if( type.eq.61) then
          call check_apr( apr_val_atm(ltog_sites(indx)),
     .                    apr(pos), check, type )
         apr_val_atm(ltog_sites(indx)) = apr(pos)
          RETURN
      end if
 
***** For the monent that is all.  Later we can add saving extended
*     tide model, and nutation series coefficients. It is not critical
*     that we save these values since they will be small and we can
*     start with apriori values of zero.  (NOTE: Total estimated values
*     are saved in the global files, so we donot need apriori values to
*     interprett these.)
 
      end



